#pragma once

#include <QWidget>

#include <GraphMol/GraphMol.h>
#include <GraphMol/MolDraw2D/Qt/MolDraw2DQt.h>

#include <memory>
#include <vector>

class Molecule : public QWidget
{
    Q_OBJECT

public:
    explicit Molecule(QWidget *parent = nullptr, Qt::WindowFlags flags = Qt::WindowFlags(0));

    void set_display_mol(boost::shared_ptr<RDKit::ROMol> new_mol );
    boost::shared_ptr<RDKit::ROMol> display_mol() { return m_mol; }

    void set_title(const std::string &mol_name);
    std::string title() const;

    boost::shared_ptr<RDKit::MolDraw2D> mol_drawer() const { return m_mol_drawer; }

    int pick_circle_rad() const;

    QSize minimumSize() const;
    QSize sizeHint() const override;

    void add_molecule_title(QPainter &qp, const std::string &mol_name, int label_box_height);

protected:
    void mousePressEvent(QMouseEvent *event) override;
    void paintEvent(QPaintEvent *event) override;

    void drawMolecule( QPainter &qp );
    void identify_selected_atoms( QPainter &qp );

private:
    int find_nearest_atom(int x_screen_pos, int y_screen_pos) const;

    mutable int m_pick_circle_rad;

    boost::shared_ptr<RDKit::MolDraw2D> m_mol_drawer;
    boost::shared_ptr<RDKit::ROMol> m_mol;
    std::string m_title;
    std::vector<int> m_selected_atoms;
//    std::vector<std::string> mols = {"C1=CC=CC=C1", "C1=CC=C2C(=C1)C(=O)C3=CC=CC=C3C2=O", "C1=CC=C2C=C3C=CC=CC3=CC2=C1"};
};

