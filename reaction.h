#pragma once

#include "molecule.h"

#include <QWidget>

#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolDraw2D/Qt/MolDraw2DQt.h>

class ChemicalReactionWidget : public QWidget
{
    Q_OBJECT
public:
    explicit ChemicalReactionWidget(QWidget *parent = nullptr, Qt::WindowFlags flags = Qt::WindowFlags(0));

    void set_reaction(boost::shared_ptr<RDKit::ChemicalReaction> new_react);
    boost::shared_ptr<RDKit::ChemicalReaction> reaction() const;

    void set_smarts(const std::string &smarts);
    std::string smarts() const;

    std::vector<Molecule*> run(RDKit::ROMOL_SPTR mol, QWidget *parent = nullptr, Qt::WindowFlags f = Qt::WindowFlags(0));

protected:
    void mousePressEvent(QMouseEvent *event) override;
    void paintEvent(QPaintEvent *event) override;

    void draw_reaction( QPainter &qp );
    void add_reaction_title( QPainter &qp, const std::string &title, int label_box_height);
private:
    boost::shared_ptr<RDKit::ChemicalReaction> m_reaction;
    boost::shared_ptr<RDKit::MolDraw2D> m_react_drawer;
    std::string m_smarts;

signals:

};
