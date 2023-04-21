#pragma once

#include <QMainWindow>
#include <QFileDialog>
#include <QMessageBox>

#include "reactiondialog.h"

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

signals:
    void proccessStarted();
    QString name_accepted();

private slots:
//    void on_proccessStared();

    void on_actionExit_triggered();

    void on_actionOpen_triggered();

    void on_actionSave_triggered();

    void on_actionRun_triggered();

    void on_actionAdd_Reaction_triggered();

    void reactionDialogAccepted();

private:
    Ui::MainWindow *ui;
    QFileDialog *fileDialog;
    QString filePath, saveFileName;
    QMessageBox *messageBox;
    ReactionDialog *reactionDialog;

    void handleResults();
};
