#ifndef CHEMICALITEM_H
#define CHEMICALITEM_H

#include <QWidget>

namespace Ui {
class ChemicalItem;
}

class ChemicalItem : public QWidget
{
    Q_OBJECT

public:
    explicit ChemicalItem(QWidget *parent = nullptr);
    ~ChemicalItem();

    void set_text(const std::string &new_text);
    void set_widget(QWidget *widget);

    std::string text() const;
    QWidget* widget();

private:
    Ui::ChemicalItem *ui;
    QWidget* m_widget;
    std::string m_text;
};

#endif // CHEMICALITEM_H
