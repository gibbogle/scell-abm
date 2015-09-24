#ifndef PLOTWIN_H
#define PLOTWIN_H

#include <QMainWindow>

namespace Ui {
class PlotWin;
}

class PlotWin : public QMainWindow
{
    Q_OBJECT
    
public:
    explicit PlotWin(QWidget *parent = 0);
    ~PlotWin();
    
private:
    Ui::PlotWin *ui;
};

#endif // PLOTWIN_H
