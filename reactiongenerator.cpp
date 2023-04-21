#include "reactiongenerator.h"

#include <GraphMol/GraphMol.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/MolDraw2D/MolDraw2DCairo.h>

#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/ChemReactions/ReactionParser.h>

#include <GraphMol/new_canon.h>
#include <unordered_map>

#include <GraphMol/ChemTransforms/MolFragmenter.h>
#include <GraphMol/MolStandardize/Fragment.h>

#define RXN_SPTR boost::shared_ptr<RDKit::ChemicalReaction>

namespace {
    RDKit::UINT_VECT unique_atoms(RDKit::ROMOL_SPTR mol){
        RDKit::UINT_VECT rank;
        RDKit::Canon::rankMolAtoms(*mol, rank, false);

        std::unordered_map<RDKit::UINT, int> uniqueIds;
        for(unsigned int i = 0; i < rank.size(); i ++){
            uniqueIds.insert(std::make_pair(rank[i], i));
            // std::cout << i  << " " << mol->getAtomWithIdx(i)->getAtomicNum() << " " << result[i] << "\n";
        }

        RDKit::UINT_VECT result;
        for(auto &kv: uniqueIds){
            result.push_back(kv.second);
        }

        return result;
    }

    RDKit::MOL_SPTR_VECT run_reaction(RXN_SPTR rxn, RDKit::ROMOL_SPTR mol){
        RDKit::ROMOL_SPTR mol1(new RDKit::ROMol(*RDKit::MolOps::removeAllHs(*mol)));
        RDKit::ROMOL_SPTR mol2(new RDKit::ROMol(*RDKit::MolOps::removeAllHs(*mol)));
        RDKit::MOL_SPTR_VECT rVect = {mol1, mol2};
        // std::cout << "I am here\n";
        std::vector<RDKit::MOL_SPTR_VECT> products = rxn->runReactants(rVect);
        // std::cout << "Finished reaction\n";
        std::vector<std::string> smiles;
        for(auto p: products){
            smiles.push_back(RDKit::MolToSmiles(*p[0]));
        }

        RDKit::MOL_SPTR_VECT res;
        for(auto s: smiles){
            res.push_back(RDKit::ROMOL_SPTR(RDKit::SmilesToMol(s)));
        }
        return res;
    }

    std::unordered_map<std::string, RDKit::ROMOL_SPTR> run_reaction_with_symm(RXN_SPTR rxn, RDKit::ROMOL_SPTR mol, const int expected_molecules = 128){
        auto uniqueIds = unique_atoms(mol);
        std::unordered_map<std::string, RDKit::ROMOL_SPTR> product; // The final reaction products
        product.reserve(expected_molecules);

        for(auto uId: uniqueIds){
            // Protect all atoms
            for(auto atom: mol->atoms()){
                atom->setProp("_protected", "1");
            }
            mol->getAtomWithIdx(uId)->clearProp("_protected");

            for(auto &new_mol: run_reaction(rxn, mol)){
                product.insert(std::make_pair(RDKit::MolToSmiles(*new_mol), new_mol)); // Removes the duplicates if there are any
            }
        }
        return product;
    }

    RXN_SPTR mol_to_reaction(const RDKit::ROMOL_SPTR mol, const unsigned int idx1, const unsigned int idx2){

        RDKit::ROMOL_SPTR tmp(RDKit::MolOps::addHs(*mol)), r; // benzene
        RDKit::MolStandardize::LargestFragmentChooser neZnam;
        r.reset(neZnam.choose(*RDKit::MolFragmenter::fragmentOnBonds(*tmp, RDKit::UINT_VECT({idx1,idx2}))));

        //Create reaction SMARTS
        std::string smiles = RDKit::MolToSmiles(*RDKit::MolOps::removeAllHs(*r));
        std::string idx1_name = "[" + std::to_string(idx1) + "*]";
        std::string idx2_name = "[" + std::to_string(idx2) + "*]";
        // std::cout << smiles << " mol fragment smiles\n";

        std::string smarts = smiles;
        smarts.replace(smarts.find(idx1_name), idx1_name.size(), "[cH1:1]");
        smarts.replace(smarts.find(idx2_name), idx2_name.size(), "[cH1:2]");
        smarts = "([cH1:1]).([cH1:2])>>" + smarts;

        RXN_SPTR react(RDKit::RxnSmartsToChemicalReaction(smarts));
        react->initReactantMatchers();

        return react;
    }

    int get_bond_idx(RDKit::ROMOL_SPTR mol, RDKit::UINT Cidx){
        RDKit::ROMOL_SPTR new_mol(RDKit::MolOps::addHs(*mol));
        RDKit::Atom *carbon = new_mol->getAtomWithIdx(Cidx);
        auto neighbours = new_mol->getAtomNeighbors(carbon);
        int Hidx = -1;
        for(const auto &nbri: make_iterator_range(new_mol->getAtomBonds(carbon))){
            const RDKit::Bond *bond =  (*new_mol)[nbri];
            if(bond->getOtherAtom(carbon)->getAtomicNum() == 1){
                if(Hidx == -1){
                    Hidx = bond->getIdx();
                }
                else{
                    return -1;
                }
            }
        }
        return Hidx;
    }

    std::string get_reaction_key(RXN_SPTR react){
        RDKit::ROMOL_SPTR benzene(RDKit::SmilesToMol("C1=CC=CC=C1"));
        auto product = run_reaction_with_symm(react, benzene);
        std::string key = RDKit::MolToSmiles(*(product.begin()->second)); // There should be only one product
        return key;
    }

    std::unordered_map<std::string,RXN_SPTR> generate_bridges(RDKit::ROMOL_SPTR mol){
        std::unordered_map<std::string, RXN_SPTR> uniqueReactions;
        RDKit::UINT_VECT rank;
        RDKit::Canon::rankMolAtoms(*mol, rank, false);

        std::unordered_map<int,RDKit::UINT_VECT> symmAtoms;
        for(unsigned int i = 0; i < rank.size(); i ++){
            if( ! symmAtoms.contains(rank[i])){
                symmAtoms[rank[i]] = {};
            }
            symmAtoms[rank[i]].push_back(i);
        }

        // What to do if the molecule is not symmetric
        if(symmAtoms.size() == rank.size()){
            return uniqueReactions; // Empty map
        }

        RDKit::ROMOL_SPTR new_mol(RDKit::MolOps::addHs(*mol));
        for(auto &atomSet : symmAtoms){
            for(unsigned int i = 0; i < atomSet.second.size(); i ++){
                for(unsigned int j = i+1; j < atomSet.second.size(); j++){
                    int idx1 = get_bond_idx(new_mol, atomSet.second[i]);
                    int idx2 = get_bond_idx(new_mol, atomSet.second[j]);
                    if(idx1 == -1 || idx2 == -1){
                        break;
                    }

                    RXN_SPTR r = mol_to_reaction(mol, idx1, idx2);
                    std::string key = get_reaction_key(r);

                    uniqueReactions.insert(std::make_pair(key, r));
                }
            }
        }
        return uniqueReactions;
    }
}
reactionGenerator::reactionGenerator(){

}

std::unordered_map<std::string,ChemicalReactionWidget*> reactionGenerator::generate_reactions(RDKit::ROMOL_SPTR mol, QWidget* parent){
    std::unordered_map<std::string,RXN_SPTR> reactions = generate_bridges(mol);
//    int counter = 1;

    std::unordered_map<std::string,ChemicalReactionWidget*> ans;
    for(const auto &p: reactions){
        ChemicalReactionWidget *widget = new ChemicalReactionWidget(parent);
        widget->set_reaction(p.second);
//        widget->set_smarts("bridge" + std::to_string(counter));
//        counter ++;

        ans.insert(std::make_pair(p.first, widget));
    }

    return ans;

}

std::unordered_map<std::string,ChemicalReactionWidget*> reactionGenerator::generate_reactions(std::vector<RDKit::ROMOL_SPTR> mol_vect, QWidget* parent){
    std::unordered_map<std::string,ChemicalReactionWidget*> result;
    for(const auto &m: mol_vect){
        for(const auto &p: generate_reactions(m, parent)){
            result.insert(p);
        }
    }
    return result;
}

