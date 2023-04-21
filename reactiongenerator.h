#pragma once

#include <unordered_map>
#include <vector>

#include "reaction.h"

class reactionGenerator
{
public:
    reactionGenerator();
    std::unordered_map<std::string,ChemicalReactionWidget*> generate_reactions(RDKit::ROMOL_SPTR mol, QWidget *parent = nullptr);
    std::unordered_map<std::string,ChemicalReactionWidget*> generate_reactions(std::vector<RDKit::ROMOL_SPTR> mol, QWidget *parent = nullptr);
};

