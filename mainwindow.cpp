#include "mainwindow.h"
#include "molecule.h"
#include "reaction.h"
#include "chemicalitem.h"

#include "./ui_mainwindow.h"

#include <QTextStream>
#include <QDebug>
#include <QDir>
#include <QHBoxLayout>

#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/FileParserUtils.h>
#include <GraphMol/ChemReactions/ReactionParser.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/MolOps.h>


void add_item(QTableWidget *table, QWidget* item, const std::string &header){
    int rows = table->rowCount();
    table->setRowCount(rows+1);

    ChemicalItem *new_item = new ChemicalItem(table);
    new_item->set_text(header);
    new_item->set_widget(item);

    table->setCellWidget(rows, 0, new_item);
}


void add_molecule(QTableWidget *table, const std::string &mol_path){

    boost::shared_ptr<RDKit::RWMol> new_mol(RDKit::MolFileToMol(mol_path));
    Molecule *mol = new Molecule(table);
    mol->set_display_mol(new_mol);
    mol->set_title(mol_path);
    add_item(table, mol, mol_path);
}

void add_reaction(QTableWidget *table, const std::string &smarts){
    boost::shared_ptr<RDKit::ChemicalReaction> new_react(RDKit::RxnSmartsToChemicalReaction(smarts));

    ChemicalReactionWidget *react = new ChemicalReactionWidget(table);
    react->set_smarts(smarts);
    react->set_reaction(new_react);

    add_item(table, react, smarts);
}

void save_all(QTableWidget *table, int row, int column){

}



MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
    , fileDialog(nullptr)
    , messageBox(nullptr)
{
    ui->setupUi(this);

    fileDialog = new QFileDialog(this);
    fileDialog->setNameFilter("MOL file (*.mol; *.mdl)");
    messageBox = new QMessageBox(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}


void MainWindow::on_actionExit_triggered()
{
    close();
}

void MainWindow::on_actionOpen_triggered()
{
    QString path = fileDialog->getOpenFileName(this, "Select a molecule", "Mol files (*.mdl *.mol)");
    if (path.isEmpty()) return;

    QString title;
    QTextStream(&title) << windowTitle() << "(" << path << ")";
    setWindowTitle(title);

    if (path.endsWith(".mol") || path.endsWith(".mdl")){
        try{
            add_molecule(ui->input_table, path.toStdString());
        }
        catch(std::exception){
            messageBox->setText("Cannot read file correctly");
            messageBox->setWindowTitle("Invalid file");
            messageBox->exec();
            return;
        }
    }
    //readFileToTable(ui->separatorBox->currentIndex());
//    readFilesFrom(path);
}

void MainWindow::on_actionSave_triggered()
{

    QString path = fileDialog->getSaveFileName(this, "Save");
    if(path.isEmpty()){
        return;
    }

    QFile file(path);

    auto mols = ui->output_table->children();
    std::string text = "";
    for(auto m: mols){
        auto h = ((ChemicalItem*)m)->text();
        text += h;
        text += " ";
    }
    messageBox->setText(QString::fromStdString(text));
    messageBox->show();
//    Molecule* mol = (Molecule*)ui->output_table->cellWidget(0, 0);
//    if(!mol){
//        messageBox->setText(QString::fromStdString("There is no molecule"));
//        messageBox->show();
//        return;
//    }
//    if(!mol->display_mol()){
//        messageBox->setText(QString::fromStdString("There is no molecule in the molecule"));
//        messageBox->show();
//        return;
//    }
//    RDKit::MolToMolFile(*(mol->display_mol()), file.fileName().toStdString());

}

void MainWindow::on_actionRun_triggered()
{
//    QLayoutItem* item;
//    while( (item = ui->output_table->layout()->takeAt(0)) != nullptr){
//        delete item->widget();
//        delete item;
//    }

    int numberOfReactions = ui->react_table->rowCount();
    int numberOfMolecules = ui->input_table->rowCount();
    int counter = 0;
    for(int i = 0; i < numberOfReactions; i ++){
        ui->progressBar->setValue((i+1) * 100.f / numberOfReactions);
        ChemicalItem *item = (ChemicalItem*)ui->react_table->cellWidget(i, 0);
        ChemicalReactionWidget *react = ((ChemicalReactionWidget*)item->widget());
        for(int j = 0; j < numberOfMolecules; j ++){
            ChemicalItem *m = (ChemicalItem*)ui->input_table->cellWidget(j, 0);
            RDKit::ROMOL_SPTR curr_mol = ((Molecule*)m->widget())->display_mol();
            // Ne pitai, ne znam
//            boost::shared_ptr<RDKit::ROMol> mol1(new RDKit::ROMol(*((Molecule*)m->widget())->display_mol()));
//            boost::shared_ptr<RDKit::ROMol> mol2(new RDKit::ROMol(*((Molecule*)m->widget())->display_mol()));
//            std::vector<boost::shared_ptr<RDKit::ROMol>> vec{mol1, mol2};

            auto res = react->run(curr_mol, ui->output_table, this->windowFlags());
            for(Molecule* new_mol: res){
                add_item(ui->output_table, new_mol, new_mol->title());
            }
        }
    }
    return;
}

void MainWindow::on_actionAdd_Reaction_triggered()
{
    add_reaction(ui->react_table, "([RH1:1]).([RH1:2])>>[RH1:1]/C=C/[RH1:2]");
    return;
}

