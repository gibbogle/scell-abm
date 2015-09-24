#include "plotwin.h"
#include "ui_plotwin.h"

PlotWin::PlotWin(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::PlotWin)
{
    ui->setupUi(this);
}

PlotWin::~PlotWin()
{
    delete ui;
}
