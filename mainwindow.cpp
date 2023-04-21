#include "mainwindow.h"
#include "molecule.h"
#include "reaction.h"
#include "reactiongenerator.h"
#include "chemicalitem.h"

#include "./ui_mainwindow.h"

#include <QTextStream>
#include <QDebug>
#include <QDir>
#include <QDirIterator>
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

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
    , fileDialog(nullptr)
    , messageBox(nullptr)
    , reactionDialog(nullptr)

{
    ui->setupUi(this);
    saveFileName = "";
    fileDialog = new QFileDialog(this);
    messageBox = new QMessageBox(this);
    reactionDialog = new ReactionDialog(this);
    connect(reactionDialog, SIGNAL(accepted()), this, SLOT(reactionDialogAccepted()));
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
    QString path = fileDialog->getExistingDirectory(this, "Select folder", "");
    if (path.isEmpty()) return;

    QStringList filters;
    filters << "*.mol" << "*.mdl";
    QDirIterator it(path, filters, QDir::Files);
    QString badFiles = "";
    while (it.hasNext()){
        QString filePath = it.next();
        try{
            add_molecule(ui->input_table, filePath.toStdString());
        }
        catch(RDKit::BadFileException()){
            badFiles += filePath + "\n";
        }
    }
    if(!badFiles.isEmpty()){
        messageBox->setText("Cannot read file correctly\n" + badFiles);
        messageBox->setWindowTitle("Invalid file");
        messageBox->exec();
    }
}

void MainWindow::on_actionSave_triggered()
{

    QString path = fileDialog->getExistingDirectory(this, "Select save directory", "");
    if(path.isEmpty()){
        return;
    }
    reactionDialog->show();

    QEventLoop loop;
    connect(this, &MainWindow::name_accepted, &loop, &QEventLoop::quit);
    loop.exec(); //exec will delay execution until the signal has arrived
    if(saveFileName.isEmpty()){
        return;
    }
    qDebug() << saveFileName;

    int numberOfMolecules = ui->output_table->rowCount();
    for(int i = 0; i < numberOfMolecules; i ++){
        QFile file(path + "/" + saveFileName + "_" + QString::fromStdString(std::to_string(i)) + ".mol");
        file.open(QFile::WriteOnly);

        ChemicalItem* item = (ChemicalItem*)ui->output_table->cellWidget(i, 0);
        Molecule* mol = (Molecule*)item->widget();
        mol->display_mol()->setProp("_Name", item->text());
        file.write(RDKit::MolToMolBlock(*mol->display_mol()).data());

        file.close();
    }
    messageBox->setWindowTitle("Saving");
    messageBox->setText("Saving done");
    messageBox->show();
}

void MainWindow::on_actionRun_triggered()
{
    ui->output_table->setRowCount(0);

    int numberOfReactions = ui->react_table->rowCount();
    int numberOfMolecules = ui->input_table->rowCount();
    if(!numberOfMolecules || !numberOfReactions){
        return;
    }

    for(int i = 0; i < numberOfReactions; i ++){
        ui->progressBar->setValue((i+1) * 100.f / numberOfReactions);
        ChemicalItem *item = (ChemicalItem*)ui->react_table->cellWidget(i, 0);
        ChemicalReactionWidget *react = ((ChemicalReactionWidget*)item->widget());
        for(int j = 0; j < numberOfMolecules; j ++){
            ChemicalItem *m = (ChemicalItem*)ui->input_table->cellWidget(j, 0);
            RDKit::ROMOL_SPTR curr_mol = ((Molecule*)m->widget())->display_mol();

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
    std::vector<RDKit::ROMOL_SPTR> mols;

    //add_reaction(ui->react_table, "([RH1:1]).([RH1:2])>>[RH1:1]/C=C/[RH1:2]");
    QString path = fileDialog->getExistingDirectory(this, "Select folder", "");
    if (path.isEmpty()) return;

    QStringList filters;
    filters << "*.mol" << "*.mdl";
    QDirIterator it(path, filters, QDir::Files);
    QString badFiles;
    while (it.hasNext()){
        QString filePath = it.next();
        try{
            mols.push_back(RDKit::ROMOL_SPTR(RDKit::MolFileToMol(filePath.toStdString())));
        }
        catch(RDKit::BadFileException()){
            badFiles += filePath + "\n";
        }
    }
    if(!badFiles.isEmpty()){
        messageBox->setText("Cannot read file correctly\n" + badFiles);
        messageBox->setWindowTitle("Invalid file");
        messageBox->exec();
    }

    reactionGenerator gen;
    auto res = gen.generate_reactions(mols, this);

    for(auto &p: res){
        add_item(ui->react_table, p.second, p.second->smarts());
    }

    return;
}

void MainWindow::reactionDialogAccepted()
{
    saveFileName = reactionDialog->getReaction();
    emit name_accepted();
}

