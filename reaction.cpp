#include "reaction.h"

#include <QPainter>

#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/new_canon.h>

#define RXN_SPTR boost::shared_ptr<RDKit::ChemicalReaction>

RDKit::UINT_VECT unique_atoms(RDKit::ROMOL_SPTR mol){
    RDKit::UINT_VECT rank;
    RDKit::Canon::rankMolAtoms(*mol, rank, false);

    std::unordered_map<RDKit::UINT, int> uniqueIds;
    for(int i = 0; i < rank.size(); i ++){
        uniqueIds.insert(std::make_pair(rank[i], i));
        // std::cout << i  << " " << mol->getAtomWithIdx(i)->getAtomicNum() << " " << result[i] << "\n";
    }

    RDKit::UINT_VECT result;
    for(auto &kv: uniqueIds){
        result.push_back(kv.second);
    }

    return result;
}

RDKit::MOL_SPTR_VECT generate_all_products(RXN_SPTR rxn, RDKit::ROMOL_SPTR mol){
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

std::unordered_map<std::string, RDKit::ROMOL_SPTR> run_reaction(RXN_SPTR rxn, RDKit::ROMOL_SPTR mol, const int expected_molecules = 128){
    auto uniqueIds = unique_atoms(mol);
    std::unordered_map<std::string, RDKit::ROMOL_SPTR> product; // The final reaction products
    product.reserve(expected_molecules);


    // for(auto &id: uniqueIds){ std::cout << id << " ";} std::cout << std::endl;
    for(auto uId: uniqueIds){
        // Protect all atoms
        for(auto atom: mol->atoms()){
            atom->setProp("_protected", "1");
        }
        mol->getAtomWithIdx(uId)->clearProp("_protected");
        // std::cout << uId << " index\n";
        // for(auto atom: mol->atoms()){
        //     std::cout << atom->getIdx() << " " << ((atom->hasProp("_protected"))? "":"not ") << "protected\n";
        // }
        // std::cout << "****\n";

        for(auto &new_mol: generate_all_products(rxn, mol)){

         auto flag = product.insert(std::make_pair(RDKit::MolToSmiles(*new_mol), new_mol)); // Removes the duplicates if there are any
            // if(!flag.second){
            //     std::cout << "Failed to insert " << RDKit::MolToSmiles(*new_mol)  << " on id " << uId << "\n";
            // }
        }

    }

    return product;
}

ChemicalReactionWidget::ChemicalReactionWidget(QWidget  *parent, Qt::WindowFlags f)
    : QWidget{parent, f}
{
    setMinimumSize(70, 50);
}

void ChemicalReactionWidget::mousePressEvent(QMouseEvent *event){
    return;
}

void ChemicalReactionWidget::paintEvent(QPaintEvent *event){
    QPainter qp;

    qp.begin(this);

    qp.setRenderHint(QPainter::Antialiasing, true);
    qp.setRenderHint(QPainter::TextAntialiasing, true);

    if(!m_reaction){ return; }

    draw_reaction(qp);
}

void ChemicalReactionWidget::draw_reaction(QPainter &qp){

    int h = rect().height();
    if(!m_smarts.empty()){
        h = int(float(rect().height()) * 0.85);
    }

    m_react_drawer.reset(new RDKit::MolDraw2DQt(rect().width(), h, &qp));

    m_react_drawer->drawReaction(*m_reaction);
//    add_reaction_title(qp, m_smarts, h);

}

void ChemicalReactionWidget::add_reaction_title(QPainter &qp,
                                  const std::string &mol_name,
                                  int label_box_height) {
  qp.setPen("Black");
  if (!mol_name.empty()) {
    int box_height = height() - label_box_height;
    qp.fillRect(0, label_box_height, width(), box_height, qp.background());
    while (true) {
      QRect br = qp.boundingRect(0, label_box_height, width(), box_height,
                                 Qt::AlignHCenter | Qt::AlignVCenter,
                                 mol_name.c_str());
      if (br.height() > box_height) {
        float scale = float(box_height) / float(br.height());
        QFont ft = qp.font();
        float new_fs = scale * ft.pointSizeF();
        ft.setPointSizeF(new_fs);
        qp.setFont(ft);
      } else {
        break;
      }
    }
    qp.drawText(0, label_box_height, width(), height() - label_box_height,
                Qt::AlignHCenter | Qt::AlignVCenter, mol_name.c_str());
  }
}

void ChemicalReactionWidget::set_smarts(const std::string &smarts){
    m_smarts = smarts;
}

std::string ChemicalReactionWidget::smarts() const{
    return m_smarts;
}

void ChemicalReactionWidget::set_reaction(boost::shared_ptr<RDKit::ChemicalReaction> new_react){
    if (!new_react) {
      m_reaction.reset();
    }
    else {
        m_reaction.reset();
        m_reaction = new_react;
        m_reaction->initReactantMatchers();
    }

    update();
}

boost::shared_ptr<RDKit::ChemicalReaction> ChemicalReactionWidget::reaction() const{
    return m_reaction;
}


std::vector<Molecule*> ChemicalReactionWidget::run(RDKit::ROMOL_SPTR mol,
                                                   QWidget *parent, Qt::WindowFlags f){
    auto products = run_reaction(reaction(), mol);
    std::vector<Molecule*> result;
    result.reserve(products.size());
    for(auto &p: products){
        Molecule *new_mol = new Molecule(parent, f);
        new_mol->set_display_mol(p.second);
        new_mol->set_title(p.first);
        result.push_back(new_mol);
    }

    return result;

}
