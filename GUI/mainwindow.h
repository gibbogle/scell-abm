/****************************************************************************

****************************************************************************/

#ifndef MAINWINDOW_H
#define MAINWINDOW_H


#include <string>
#include <fstream>
#include <QTcpServer>
#include <QTcpSocket>
#include <QWaitCondition>

using namespace std;

//#include "ui_scell_GUI.h"
#include "ui_spheroid_GUI.h"
#include <qwt_plot_curve.h>
#include <qwt_plot_marker.h>
#include <qwt_symbol.h>
#include "params.h"
#include "qmylabel.h"
#include "qmycheckbox.h"
#include "misc.h"
#include "plot.h"
#include "myvtk.h"
#include "field.h"
#include "result_set.h"
#include "log.h"
#include "SimpleView3DUI.h"
#include "SimpleView2DUI.h"
#include "qvideooutput.h"

#include <qwt_plot.h>
#include <qwt_plot_grid.h>
#include <qwt_interval_data.h>
#include "histogram_item.h"

QT_BEGIN_NAMESPACE
class QAction;
class QMenu;
class QPlainTextEdit;
class QMdiArea;
class QTcpServer;
class QTcpSocket;
QT_END_NAMESPACE

class SliderPlus 
{
	QString name;
	int pindex;
	int windex;
	double vmin;
	double vmax;
	double dv;
	int n;

public:

	SliderPlus(QString, double, double, int,int, int);
	~SliderPlus();
	int val_to_int(double);
	double int_to_val(int);
	QString val_to_str(double);
	double str_to_val(QString);
	int pIndex();
	int wIndex();
	int nTicks();
};


class MainWindow : public QMainWindow, private Ui::MainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = 0);

	char msg[2048];

protected:
    void closeEvent(QCloseEvent *event);

private slots:
    void on_action_show_gradient2D_triggered();
    void on_action_show_gradient3D_triggered();
    void on_action_FACS_triggered();
    void on_checkBox_FACS_PLOT_toggled(bool checked);
    void newFile();
    void open();
    void about();
    void documentWasModified();

    bool save();
    bool saveAs();
	void readInputFile();
	void loadResultFile();
    void goToInputs();
    void goToOutputs();
    void goToVTK();
    void goToField();
    void goToFACS();
    void runServer();
    void pauseServer();
    void stopServer();
	void changeParam();
	void redrawDistPlot();
	void showMore(QString);
	void updateSliderBox();

	double getMaximum(RESULT_SET *, double *);
	void addGraph();
	void removeGraph();
	void removeAllGraphs();
	void playVTK();
	void setVTKSpeed();
	void saveSnapshot();
    void saveProfileData();
    void showGradient3D();
    void showGradient2D();
    void setSavePosStart();

//    void on_radioButton_oxygen_clicked();
//    void on_radioButton_glucose_clicked(bool checked);

    void onSelectConstituent();

    void on_verticalSliderTransparency_sliderMoved(int position);
    void on_checkBox_CELLDISPLAY_1_toggled(bool display);
    void on_checkBox_CELLDISPLAY_2_toggled(bool display);
    void on_comboBox_CELLCOLOUR_1_currentIndexChanged(int index);
    void on_comboBox_CELLCOLOUR_2_currentIndexChanged(int index);

    void on_cbox_SAVE_PROFILE_DATA_toggled(bool checked);

    void LoadProtocol(QString);
    void SaveProtocol(QString);
    void ProtocolChanged(int, int);

public slots:
	void preConnection();
	void outputData(QString);
	void postConnection();
	void timer_update();
	void errorPopup(QString);
	void displayScene();
    void reportBadDLL(QString);
    void showSummary(int);
    void showFACS();
    void showHisto();
    void startRecorderVTK();
    void stopRecorderVTK();
    void startRecorderFACS();
    void stopRecorderFACS();
    bool getVideoFileInfo(int *nframes, QString *itemFormat, QString *itemCodec, QString *videoFileName);

    void buttonClick_constituent(QAbstractButton* button);
    void buttonClick_plane(QAbstractButton* button);
    void buttonClick_canvas(QAbstractButton* button);
    void textChanged_fraction(QString text);
	void textEdited_fraction(QString text);
    void setupConstituents();

    void on_cbox_USE_TPZ_DRUG_toggled(bool checked);
//    void on_cbox_TPZ_DRUG_SIMULATE_METABOLITE_toggled(bool checked);
    void on_cbox_USE_DNB_DRUG_toggled(bool checked);
//    void on_cbox_DNB_DRUG_SIMULATE_METABOLITE_toggled(bool checked);
//    void on_cbox_USE_RADIATION_toggled(bool checked);
    void on_line_CELLPERCENT_1_textEdited(QString pc1_str);
    void on_line_CELLPERCENT_2_textEdited(QString pc2_str);
    void radioButtonChanged(QAbstractButton *b);
    void killModelChanged();
    void on_buttonGroup_celltype_buttonClicked(QAbstractButton* button);
    void on_buttonGroup_histotype_buttonClicked(QAbstractButton* button);
    void on_checkBox_histo_logscale_toggled();

    void on_comb_TPZ_currentIndexChanged(int);
    void on_comb_DNB_currentIndexChanged(int);
// For Kd computed in the GUI
//    void on_pushButton_SN30K_Kd_1_clicked();
//    void on_pushButton_SN30K_Kd_2_clicked();

    void processGroupBoxClick(QString);
signals:
    void facs_update();
    void histo_update();

private:
    void createActions();
	void createLists();
	void drawDistPlots();
    void initFACSPlot();
    void initHistoPlot();
    void setupParamList();
	void loadParams();
	void reloadParams();

    void trackError();

    void enableUseOxygen();
    void disableUseOxygen();
    void enableUseGlucose();
    void disableUseGlucose();
    void enableUseTracer();
    void disableUseTracer();
//    void enableUseSN30K();
//    void disableUseSN30K();
//    void enableUseTPZ();
//    void disableUseTPZ();
//    void enableUseDNB();
//    void disableUseDNB();
//    void enableUseTreatmentFile();
//    void disableUseTreatmentFile();
//    void setTreatmentFileUsage();

	void writeout();
	void execute_para();
	void init_VTK();
    void read_cell_positions();
	void close_sockets();
	void compareOutputs();
	void clearAllGraphs();
	void initializeGraphs(RESULT_SET *);
	void drawGraphs();
	QString selectResultSet();
	int selectGraphCase();
    void setupCellColours();
    void setupGraphSelector();
    void setGraphsActive();
    void initDrugComboBoxes();
    void test_histo();
    void makeHistoPlot(int numValues, double xmin, double width, QVector<double> values);
    void showBool(QString, bool);

	double erf(double z);
    double pnorm(double x1, double x2, double mu, double sig);
    double plognorm(double x1, double x2, double mu, double sig);
    void create_lognorm_dist(double p1, double p2,int n, double *x, double *prob);
	int dist_limit(double *p, int n);
	QString parse_rbutton(QString wtag, int *rbutton_case);
	void setBdryRadioButton(QRadioButton *w_rb, int val);
	void setLineEditVisibility(QString wname, int val);

	PARAM_SET get_param(int);

    void createMenus();
    void createToolBars();
    void createStatusBar();
    void readSettings();
    void writeSettings();
    bool maybeSave();
    void loadFile(const QString &fileName);
//    bool saveFile(const QString &fileName);
    void setCurrentFile(const QString &fileName);
    QString strippedName(const QString &fullFileName);

    void SetupProtocol();
    void setField(QTableWidget *, int, int, QString);
    int getField(QTableWidget *, int, int, QString *);

    QPlainTextEdit *textEdit;
    QString curFile;
	QList<QLineEdit *> lineEdit_list;
	QList<QSpinBox *> spin_list;
	QList<QComboBox *> combo_list;
	QList<QCheckBox *> checkbox_list;
	QList<QRadioButton *> radiobutton_list;
	QList<QSlider *> slider_list;
	QList<QLabel *> label_list;
	QList<SliderPlus *> sliderplus_list;
	QList<QWidget *> sliderParam;
	QList<RESULT_SET *> result_list;

	QwtPlot *distplot_list[5];
	QwtPlotCurve *curve_list[5];

	QList<QWidget *> widget_list;

    QwtPlot *qpFACS;
    QwtPlotCurve *curveFACS;
    QwtPlot *qpHistoBar, *qpHistoLine;
    HistogramItem *histogram;

	int nDistPts;
	int nTicks;
	int nParams;
    int ndistplots;
    int nSliders;
	int nWidgets;
	int nLabels;
    int nCheckBoxes;
	int *param_to_sliderIndex;
	bool paramSaved;
	bool paused;
	bool posdata;
	bool DCmotion;
	bool done;
	bool first;
    bool firstDisplay;
    bool started;
	bool firstVTK;
	bool playingVTK;
	int tickVTK;
	int currentDescription;
	QString defaultInputFile;
	QString inputFile;
	QString cellfile;
	QString vtkfile;
	QTextBrowser *box_outputData;
	SocketHandler *sthread0;
	SocketHandler *sthread1;
	QTimer *timer;

	int step;
	int ntimes;
	int savepos_start;
	int ncpu;
	double hours;
	double hour;
	int progress;
	int nGraphs;		// act, ntot_LN, ncog_PER, ...
	int nGraphCases;
    QColor comboColour[30];

	RESULT_SET *newR;

//	Plot *graph_act;
//	Plot *graph_ntot_LN;
//	Plot *graph_ncog_PER;
//	Plot *graph_ncog_LN;
//	Plot *graph_ncog;
//	Plot *graph_ncogseed;
//	Plot *graph_nDC;
//	Plot *graph_teffgen;
//	Plot *graph_nbnd;
	Plot *graph_dummy;	// placeholder

    Plot *pGraph[32];
    QMyCheckBox *checkBox_conc;
    QMyCheckBox *checkBox_vol;
    QMyCheckBox *checkBox_oxy;
    QMyCheckBox **cbox_ts;
//    QCheckBox **cbox_ts;  // GRAPHS_CHANGE

    QButtonGroup *buttonGroup_histo;
    QVBoxLayout *vbox_histo;
    QRadioButton **histo_rb_list;
    QButtonGroup *buttonGroup_FACS_x_vars;
    QVBoxLayout *vbox_FACS_x_vars;
    QRadioButton **FACS_x_vars_rb_list;
    QButtonGroup *buttonGroup_FACS_y_vars;
    QVBoxLayout *vbox_FACS_y_vars;
    QRadioButton **FACS_y_vars_rb_list;

	QString graphCaseName[Plot::ncmax];
	RESULT_SET *graphResultSet[Plot::ncmax];
	static const bool show_outputdata = false;
	static const bool use_CPORT1 = false;

	static const int CPORT0 = 5000;
	static const int CPORT1 = 5001;
	static const bool USE_RANGES = false;

	MyVTK *vtk;
    Field *field;
	ExecThread *exthread;

    QVideoOutput   *videoVTK;
    QVideoOutput   *videoFACS;

signals:
    void pause_requested();

};

class MyDoubleValidator : public QDoubleValidator
{
public:
	MyDoubleValidator( double bottom, double top, int decimals, QObject* parent = 0)
		: QDoubleValidator( bottom, top, decimals, parent)
	{}

	QValidator::State validate ( QString &input, int &pos ) const
	{
		if ( input.isEmpty() || input == "." ) {
			return Intermediate;
		}
		bool ok;
		double entered = input.toDouble(&ok);
		if (!ok) return Invalid;
		if (entered < bottom())
			return Intermediate;
		if ( QDoubleValidator::validate( input, pos ) != Acceptable ) {
			return Invalid;
		}
		return Acceptable;
	}
};

#endif
