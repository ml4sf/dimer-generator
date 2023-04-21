#include "reactiondialog.h"
#include "ui_reactiondialog.h"

ReactionDialog::ReactionDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::ReactionDialog)
{
    ui->setupUi(this);
}

ReactionDialog::~ReactionDialog()
{
    delete ui;
}

QString ReactionDialog::getReaction()
{
    QString str = ui->reactionLineEdit->text();
    return str.simplified();
}
