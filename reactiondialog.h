#ifndef REACTIONDIALOG_H
#define REACTIONDIALOG_H

#include <QDialog>

namespace Ui {
class ReactionDialog;
}

class ReactionDialog : public QDialog
{
    Q_OBJECT

public:
    explicit ReactionDialog(QWidget *parent = nullptr);
    ~ReactionDialog();
    QString getReaction();

private:
    Ui::ReactionDialog *ui;
};

#endif // REACTIONDIALOG_H
