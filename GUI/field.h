#ifndef FIELD_H
#define FIELD_H

#include <QDialog>
#include <QMessageBox>
#include <QtGui>
#include <QMouseEvent>
#include "global.h"
#include "plot.h"
#include "myqgraphicsview.h"

#define CANVAS_WIDTH 696

struct field_data {
    int NX, NY, NZ;
    int NCONST;
    double DX;
    double *Cave;   // Caverage(NX,NY,NZ,NCONST) where NCONST = MAX_CHEMO
    // NOT conc[MAX_CONC+N_EXTRA+1] with added CFSE, dVdt, volume, O2byVol
    int ncells;
    CELL_DATA *cell_data;
};
typedef field_data FIELD_DATA;

#define X_AXIS 1
#define Y_AXIS 2
#define Z_AXIS 3

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

extern "C" {
    void get_fielddata(int *, double *, FIELD_DATA *, int *);
}

//class Field : public QMainWindow
class Field : public QWidget
{
public:
    Field(QWidget *);
    ~Field();
    void displayField(int, int *);
    void displayField1();
    void setSliceChanged();
    void chooseFieldColor(double c, double cmin, double cmax, bool use_log, int rgbcol[]);
    void chooseRateColor(double fr, int rgbcol[]);
    void getTitle(int iconst, QString *title);
    bool isConcPlot();
    void setConcPlot(bool);
    bool isVolPlot();
    void setVolPlot(bool);
    bool isOxyPlot();
    void setOxyPlot(bool);
    void selectCellConstituent();
    void setExecuting(bool);
    void setSaveImages(bool);
    void setUseLogScale(bool);
    void setCellConstituentButtons(QGroupBox *gbox, QButtonGroup *bg, QVBoxLayout **vbox, QList<QRadioButton *> *rb_list, QString tag);
    void setFieldConstituentButtons(QGroupBox *gbox, QButtonGroup *bg, QVBoxLayout **vbox, QList<QRadioButton *> *rb_list, QString tag);

    QWidget *field_page;
    bool save_images;
    bool use_log;
    MyQGraphicsView* view;
    int axis;
    double fraction;
    int hour;
    int ifield;
    QString const_name[16];     // not used
    int cell_constituent;
    int field_constituent;
    bool slice_changed;
    bool useConcPlot;
    bool useVolPlot;
    bool useOxyPlot;
    FIELD_DATA fdata;
    Plot *pGconc, *pGvol, *pGoxy;
    bool executing;
    char msg[1024];

    QButtonGroup *buttonGroup_cell_constituent;
    QButtonGroup *buttonGroup_field_constituent;
    QVBoxLayout *vbox_cell_constituent;
    QVBoxLayout *vbox_field_constituent;
//    QRadioButton **cell_constituent_rb_list;
//    QRadioButton **field_constituent_rb_list;

    QList<QRadioButton *> cell_constituent_rb_list;
    QList<QRadioButton *> field_constituent_rb_list;
    QList<QLineEdit *> line_maxConc_list;
    QVBoxLayout *vbox_cell_max_concentration;

    void setCellConstituent(QAbstractButton* button);
    void setFieldConstituent(QAbstractButton* button);
    void setPlane(QAbstractButton* button);
	void setFraction(QString text);
    void setMaxConcentrations(QGroupBox *gbox);
};

#endif // FIELD_H
