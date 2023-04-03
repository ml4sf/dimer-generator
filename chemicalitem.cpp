#include "chemicalitem.h"
#include "ui_chemicalitem.h"

ChemicalItem::ChemicalItem(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::ChemicalItem)
{
    ui->setupUi(this);
    ui->header->setAlignment(Qt::AlignCenter);
}

ChemicalItem::~ChemicalItem()
{
    delete ui;
}

void ChemicalItem::set_text(const std::string &new_text){
    m_text = new_text;
    ui->header->setText(QString::fromStdString(new_text));
}

void ChemicalItem::set_widget(QWidget *widget){
    m_widget = widget;
    ui->item_layout->addWidget(widget);
    update();
}

QWidget* ChemicalItem::widget(){
    return m_widget;
}

std::string ChemicalItem::text() const{
    return m_text;
}
