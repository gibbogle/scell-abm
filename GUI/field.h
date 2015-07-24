#ifndef FIELD_H
#define FIELD_H

#include <QDialog>
#include <QMessageBox>
#include <QtGui>
#include <QMouseEvent>
#include "plot.h"
#include "myqgraphicsview.h"

#define CANVAS_WIDTH 696
#define MAX_CONC 9  // must = MAX_CHEMO in DLL
#define NEXTRA 3    // must = N_EXTRA in DLL

struct field_data {
    int site[3];
    int state;
    double volume;
    double conc[MAX_CONC+NEXTRA+1];    // added CFSE, dVdt, volume, O2byVol
};

typedef field_data FIELD_DATA;

#define X_AXIS 1
#define Y_AXIS 2
#define Z_AXIS 3

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

extern "C" {
    void get_fieldinfo(int *, int *, double *, int *, int *, int *, int *);
    void get_fielddata(int *, double *, int *, int *, FIELD_DATA *, int *);
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
//    void makeConcPlot(QMdiArea *);
//    void updateConcPlot();
    bool isVolPlot();
    void setVolPlot(bool);
//    void makeVolPlot(QMdiArea *);
//    void updateVolPlot();
    bool isOxyPlot();
    void setOxyPlot(bool);
//    void makeOxyPlot(QMdiArea *);
//    void updateOxyPlot();
    void selectConstituent();
    void setExecuting(bool);
    void setSaveImages(bool);
    void setUseLogScale(bool);
//    void setConstUsage(int, int *);
    void setConstituentButtons(QGroupBox *gbox, QButtonGroup *bg, QVBoxLayout **vbox, QRadioButton ***rb_list, QString tag);

    QWidget *field_page;
    bool save_images;
    bool use_log;
    MyQGraphicsView* view;
    int NX;
    int axis;
    double fraction;
    int hour;
    int ifield;
    int nsites, nconst, const_used[MAX_CONC+NEXTRA+1];
    int nvars_used;
    int cvar_index[32];
    QString const_name[16];
    QString constituentText;
    int constituent;
    bool slice_changed;
	bool constituent_changed;
    bool useConcPlot;
    bool useVolPlot;
    bool useOxyPlot;
    FIELD_DATA *data;
    Plot *pGconc, *pGvol, *pGoxy;
    bool executing;
    char msg[1024];

    QButtonGroup *buttonGroup_constituent;
    QVBoxLayout *vbox_constituent;
    QRadioButton **constituent_rb_list;

    void setConstituent(QAbstractButton* button);
    void setPlane(QAbstractButton* button);
	void setFraction(QString text);
};

#endif // FIELD_H
