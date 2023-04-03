#include "molecule.h"

#include <GraphMol/GraphMol.h>

#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolFileStereochem.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/MolDraw2D/Qt/MolDraw2DQt.h>
#include <GraphMol/Depictor/RDDepictor.h>

#include <QPainter>
#include <QSize>

#include <iostream>
#include <random>

Molecule::Molecule(QWidget *parent, Qt::WindowFlags f)
    : QWidget{parent, f},
      m_pick_circle_rad(-1)
{
    setMinimumSize(70, 50);
}

void Molecule::mousePressEvent(QMouseEvent *event){

}

int Molecule::find_nearest_atom(int x_screen_pos,
                                int y_screen_pos) const {
    int nearest_at = -1, nearest_dist = std::numeric_limits<int>::max();

    for(int i = 0, is = m_mol->getNumAtoms(); i < is; i++){
        Point2D screen_cds = m_mol_drawer->getDrawCoords(i);
        int x = screen_cds.x - x_screen_pos;
        int y = screen_cds.y - y_screen_pos;
        int dist = x*x + y*y;

        if(dist < nearest_dist){
            nearest_dist = dist;
            nearest_at = i;
        }
    }

    if(nearest_dist < m_pick_circle_rad * m_pick_circle_rad){
        return nearest_at;
    }
    else {
        return -1;
    }
}

void Molecule::paintEvent(QPaintEvent *event){
    QPainter qp;

    qp.begin(this);

    qp.setRenderHint(QPainter::Antialiasing, true);
    qp.setRenderHint(QPainter::TextAntialiasing, true);

    if(!m_mol){ return; }

    drawMolecule(qp);
}

void Molecule::drawMolecule(QPainter &qp){

    int h = rect().height();
    if(!m_title.empty()){
        h = int(float(rect().height()) * 0.85);
    }

    m_mol_drawer.reset(new RDKit::MolDraw2DQt(rect().width(), h, &qp));

    m_mol_drawer->drawMolecule(*m_mol, &m_selected_atoms);
//    add_molecule_title(qp, m_title, h);

    identify_selected_atoms(qp);
}

void Molecule::identify_selected_atoms(QPainter &qp) {
  static QPen sel_pen(QColor("Orange"));
  qp.setPen(sel_pen);

  // put an orange square round selected atoms.

  for (auto sa: m_selected_atoms) {
    Point2D at_cds = mol_drawer()->getDrawCoords(sa);
    qp.drawRect(at_cds.x - m_pick_circle_rad,
                at_cds.y - m_pick_circle_rad,
                2 * m_pick_circle_rad,
                2 * m_pick_circle_rad);
  }
}

QSize Molecule::sizeHint() const { return QSize(400, 400); }

QSize Molecule::minimumSize() const {return QSize(100, 100);}

int Molecule::pick_circle_rad() const {
  if (-1 == m_pick_circle_rad) {
    m_pick_circle_rad = std::min(width() / 100, height() / 100);
  }

  return m_pick_circle_rad;
}

void Molecule::set_display_mol(boost::shared_ptr<RDKit::ROMol> new_mol){
    if (!new_mol) {
      m_mol.reset();
    } else {
      boost::shared_ptr<RDKit::RWMol> mol(new RDKit::RWMol(*new_mol));
      RDKit::MolOps::Kekulize(*mol);
      RDDepict::compute2DCoords(*mol);
      if (!mol->hasProp("_drawingBondsWedged")) {
        RDKit::Conformer conf = mol->getConformer();
        RDKit::WedgeMolBonds(*mol, &conf);
      }

      m_mol = boost::shared_ptr<RDKit::ROMol>(new RDKit::ROMol(*mol));
      mol.reset();
    }
    update();
}

void Molecule::add_molecule_title(QPainter &qp,
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

void Molecule::set_title(const std::string &mol_name){
    m_title = mol_name;
}

std::string Molecule::title() const{
    return m_title;
}
