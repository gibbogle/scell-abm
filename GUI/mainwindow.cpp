/****************************************************************************
 scell_GUI
****************************************************************************/

#include <QtGui>

#include "mainwindow.h"
#include "log.h"
#include "params.h"
#include "graphs.h"
#include "misc.h"
#include "plot.h"
#include "myvtk.h"
#include "field.h"
#include "transfer.h"
#include <QDebug>

#include "dialog.h"
#include "drug.h"
#include "global.h"

#include "../src/version.h"

#ifdef linux
#include <QTcpServer>
#else
#include <QTcpServer.h>
#endif

LOG_USE();

Params *parm;	// I don't believe this is the right way, but it works
Graphs *grph;

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
MainWindow::MainWindow(QWidget *parent)
   : QMainWindow(parent)
{
    LOG_MSG("Executable build versions:");
    Global::GUI_build_version = GUI_BUILD_VERSION;
    Global::DLL_build_version = DLL_BUILD_VERSION;
    LOG_QMSG("GUI build version: " + Global::GUI_build_version);
    LOG_QMSG("DLL build version: " + Global::DLL_build_version);
    LOG_MSG("");

    LOG_MSG("Started MainWindow");
	setupUi(this);
    LOG_MSG("did setupUi");
    showMaximized();
    setWindowTitle("Tumour Spheroid ABM - off-lattice");

    QString currPath = QDir::currentPath();
    LOG_QMSG("starting path: " + currPath);
    QString newPath = currPath + "/execution";
    LOG_QMSG("newPath: " + newPath);
    QDir::setCurrent(newPath);
    currPath = QDir::currentPath();
    LOG_QMSG("current path: " + currPath);

    makeDrugFileLists();
    initDrugComboBoxes();

    // Some initializations
    nDistPts = 200;
	nTicks = 1000;
	tickVTK = 100;	// timer tick for VTK in milliseconds
    ndistplots = 2;
    paramSaved = false;
	paused = false;
	posdata = false;
    DCmotion = false;
    done = false;
    first = true;
    firstDisplay = true;
    started = false;
    firstVTK = true;
    exthread = NULL;
    histogram = NULL;
    Global::recordingVTK = false;
    Global::showingVTK = false;
    Global::recordingFACS = false;
    Global::showingFACS = false;
    Global::casename = "";
    nGraphCases = 0;
	for (int i=0; i<Plot::ncmax; i++) {
		graphResultSet[i] = 0;
	}
	vtkfile = "basecase.pos";
	savepos_start = 0;
	ntimes = 0;
    step = 0;
	hour = 0;

	param_to_sliderIndex = NULL;
	defaultInputFile = "basecase.inp";
    inputFile = defaultInputFile;

	parm = new Params();
	nParams = parm->nParams;
    field = new Field(page_2D);
    grph = new Graphs();
//    pal = new Colours();
//    histo_rb_list = NULL;
    histo_rb_list.clear();
    vbox_histo = NULL;
    buttonGroup_histo = new QButtonGroup;
//    FACS_x_vars_rb_list = NULL;
    FACS_x_vars_rb_list.clear();
    vbox_FACS_x_vars = NULL;
    buttonGroup_FACS_x_vars = new QButtonGroup;
//    FACS_y_vars_rb_list = NULL;
    FACS_y_vars_rb_list.clear();
    vbox_FACS_y_vars = NULL;
    buttonGroup_FACS_y_vars = new QButtonGroup;

    setupGraphSelector();
    setGraphsActive();

    for(int i=0; i<16; i++)
        pGraph[i] = NULL;
    LOG_QMSG("did Graphs");

    createLists();
    LOG_QMSG("did createLists");
    createActions();
    LOG_QMSG("did createActions");
    drawDistPlots();
    LOG_QMSG("did drawDistPlots");
    initFACSPlot();
    LOG_QMSG("did initFACSPlot");
    initHistoPlot();
    LOG_QMSG("did initHistoPlot");
//    initDrugComboBoxes();
    loadParams();
    LOG_QMSG("Did loadparams");

    setFields();

    SetupProtocol();

    writeout();

    timer = new QTimer(this);
    vtk = new MyVTK(mdiArea_VTK, widget_key);
    vtk->init();

    LOG_QMSG("do setupCellColours");
    setupCellColours();
    LOG_QMSG("did setupCellColours");
    QRect rect;
    rect.setX(50);
    rect.setY(30);
#ifdef __DISPLAY768
    rect.setHeight(642);
    rect.setWidth(642);
#else
    rect.setHeight(786);
    rect.setWidth(786);
#endif
    mdiArea_VTK->setGeometry(rect);

    videoVTK = new QVideoOutput(this, VTK_SOURCE, vtk->renWin, NULL);
    videoFACS = new QVideoOutput(this, QWT_SOURCE, NULL, qpFACS);

    tabs->setCurrentIndex(6);
    setupPopup();
    goToInputs();
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
//void MainWindow::initDrugComboBoxes()
//{
//    comb_TPZ->addItem("SN30000");
//    comb_TPZ->addItem("TPZ");
//    comb_TPZ->addItem("Drug3");
//    comb_DNB->addItem("PR104A");
//    comb_DNB->addItem("DNB2");
//    comb_DNB->addItem("DNB3");
//    comb_DNB->addItem("DNB4");
//}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::createActions()
{
	action_stop->setEnabled(false);
    action_pause->setEnabled(false);
    action_inputs->setEnabled(false);
    action_outputs->setEnabled(false);
    action_save_3D_snapshot->setEnabled(false);
    action_save_profile_data->setEnabled(false);
    action_show_gradient3D->setEnabled(false);
    action_show_gradient2D->setEnabled(false);
    action_field->setEnabled(false);
//    action_start_recording->setEnabled(true);
//    action_stop_recording->setEnabled(false);
    text_more->setEnabled(false);
    connect(action_open_input, SIGNAL(triggered()), this, SLOT(readInputFile()));
    connect(action_load_results, SIGNAL(triggered()), this, SLOT(loadResultFile()));
    connect(action_saveAs, SIGNAL(triggered()), this, SLOT(saveAs()));
    connect(action_save, SIGNAL(triggered()), this, SLOT(save()));
    connect(action_inputs, SIGNAL(triggered()), SLOT(goToInputs()));
    connect(action_outputs, SIGNAL(triggered()), SLOT(goToOutputs()));
	connect(action_VTK, SIGNAL(triggered()), SLOT(goToVTK()));
    connect(action_FACS, SIGNAL(triggered()), SLOT(goToFACS()));
    connect(action_field, SIGNAL(triggered()), SLOT(goToField()));
    connect(action_run, SIGNAL(triggered()), SLOT(runServer()));
    connect(action_pause, SIGNAL(triggered()), SLOT(pauseServer()));
    connect(action_stop, SIGNAL(triggered()), SLOT(stopServer()));
    connect(action_play_VTK, SIGNAL(triggered()), SLOT(playVTK()));
    connect(action_set_speed, SIGNAL(triggered()), SLOT(setVTKSpeed()));
    connect(buttonGroup_FACS_x_vars, SIGNAL(buttonClicked(QAbstractButton*)), this, SIGNAL(facs_update()));
    connect(buttonGroup_FACS_y_vars, SIGNAL(buttonClicked(QAbstractButton*)), this, SIGNAL(facs_update()));
    connect(checkBox_FACS_log_x, SIGNAL(stateChanged(int)), this, SIGNAL(facs_update()));
    connect(checkBox_FACS_log_y, SIGNAL(stateChanged(int)), this, SIGNAL(facs_update()));
    connect(buttonGroup_histo, SIGNAL(buttonClicked(QAbstractButton*)), this, SIGNAL(histo_update()));

//    connect(line_TPZ_KILL_MODEL_CELL1,SIGNAL(textChanged(QString)),this,SLOT(killModelChanged()));
//    connect(line_TPZ_KILL_MODEL_CELL2,SIGNAL(textChanged(QString)),this,SLOT(killModelChanged()));
//    connect(line_DNB_KILL_MODEL_CELL1_MET1,SIGNAL(textChanged(QString)),this,SLOT(killModelChanged()));
//    connect(line_DNB_KILL_MODEL_CELL1_MET2,SIGNAL(textChanged(QString)),this,SLOT(killModelChanged()));
//    connect(line_DNB_KILL_MODEL_CELL2_MET1,SIGNAL(textChanged(QString)),this,SLOT(killModelChanged()));
//    connect(line_DNB_KILL_MODEL_CELL2_MET2,SIGNAL(textChanged(QString)),this,SLOT(killModelChanged()));

    connect(line_NXB,SIGNAL(textChanged(QString)),this,SLOT(setFields()));
    connect(line_DELTA_X,SIGNAL(textChanged(QString)),this,SLOT(setFields()));

    connect(this,SIGNAL(pause_requested()),SLOT(pauseServer()));

    for (int i=0; i<parm->nInfolabel; i++) {
        QString tag;
        parm->get_labeltag(i, &tag);
        QString objName = "infolabel_" + tag;
        QLabel *label = findChild<QLabel *>(objName);
        connect((QObject *)label, SIGNAL(labelClicked(QString)), this, SLOT(showMore(QString)));
    }

	for (int i=0; i<nLabels; i++) {
		QLabel *label = label_list[i];
		QString label_str = label->objectName();
//		LOG_QMSG(label_str);
		if (label_str.startsWith("label_")) {
			connect((QObject *)label, SIGNAL(labelClicked(QString)), this, SLOT(showMore(QString)));
//			LOG_QMSG(label_str);
		}
	}

    for (int i=0; i<nCheckBoxes; i++) {
        QCheckBox *cbox = checkbox_list[i];
        QString cbox_str = cbox->objectName();
//		LOG_QMSG(label_str);
        if (cbox_str.startsWith("cbox_")) {
            connect((QObject *)cbox, SIGNAL(checkBoxClicked(QString)), this, SLOT(showMore(QString)));
//			LOG_QMSG(label_str);
        }
    }

    // Graph menu
//    connect(action_add_graph, SIGNAL(triggered()), this, SLOT(addGraph()));
//    connect(action_remove_graph, SIGNAL(triggered()), this, SLOT(removeGraph()));
//    connect(action_remove_all, SIGNAL(triggered()), this, SLOT(removeAllGraphs()));
    connect(action_save_3D_snapshot, SIGNAL(triggered()), this, SLOT(saveSnapshot()));
    connect(action_save_profile_data, SIGNAL(triggered()), this, SLOT(saveProfileData()));
    connect(actionStart_recording_VTK, SIGNAL(triggered()), this, SLOT(startRecorderVTK()));
    connect(actionStop_recording_VTK, SIGNAL(triggered()), this, SLOT(stopRecorderVTK()));
    connect(actionStart_recording_FACS, SIGNAL(triggered()), this, SLOT(startRecorderFACS()));
    connect(actionStop_recording_FACS, SIGNAL(triggered()), this, SLOT(stopRecorderFACS()));

//    connect(action_show_gradient3D, SIGNAL(triggered()), this, SLOT(showGradient3D()));
//    connect(action_show_gradient2D, SIGNAL(triggered()), this, SLOT(showGradient2D()));
    connect(field->buttonGroup_cell_constituent, SIGNAL(buttonClicked(QAbstractButton*)), this, SLOT(buttonClick_cell_constituent(QAbstractButton*)));
    connect(field->buttonGroup_field_constituent, SIGNAL(buttonClicked(QAbstractButton*)), this, SLOT(buttonClick_field_constituent(QAbstractButton*)));
    connect(buttonGroup_plane, SIGNAL(buttonClicked(QAbstractButton*)), this, SLOT(buttonClick_plane(QAbstractButton*)));
//    connect(buttonGroup_constituent, SIGNAL(buttonClicked(QAbstractButton*)), field, SLOT(setConstituent(QAbstractButton*)));
//	  connect(lineEdit_fraction, SIGNAL(textChanged(QString)), this, SLOT(textChanged_fraction(QString)));
	connect(lineEdit_fraction, SIGNAL(textEdited(QString)), this, SLOT(textEdited_fraction(QString)));
    connect(action_select_constituent, SIGNAL(triggered()), SLOT(onSelectConstituent()));
    connect(line_CELLPERCENT_1, SIGNAL(textEdited(QString)), this, SLOT(on_line_CELLPERCENT_1_textEdited(QString)));
    connect(line_CELLPERCENT_2, SIGNAL(textEdited(QString)), this, SLOT(on_line_CELLPERCENT_2_textEdited(QString)));

// For Kd computed in the GUI
//    connect(buttonGroup_SN30K_killmodel_1, SIGNAL(buttonClicked(QAbstractButton*)), this, SLOT(radioButtonChanged(QAbstractButton*)));
//    connect(buttonGroup_SN30K_killmodel_2, SIGNAL(buttonClicked(QAbstractButton*)), this, SLOT(radioButtonChanged(QAbstractButton*)));
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::createLists()
{        
	lineEdit_list = findChildren<QLineEdit *>();
	spin_list = findChildren<QSpinBox *>();
	combo_list = findChildren<QComboBox *>();
	checkbox_list = findChildren<QCheckBox *>();
	radiobutton_list = findChildren<QRadioButton *>();
	slider_list = findChildren<QSlider *>();
	label_list = findChildren<QLabel *>();

	for (int i=0; i<lineEdit_list.length(); i++) {
		widget_list.append(lineEdit_list[i]);
    }
	for (int i=0; i<spin_list.length(); i++) {
		widget_list.append(spin_list[i]);
	}
	for (int i=0; i<combo_list.length(); i++) {
		widget_list.append(combo_list[i]);
	}
	for (int i=0; i<checkbox_list.length(); i++) {
		widget_list.append(checkbox_list[i]);
	}
	for (int i=0; i<radiobutton_list.length(); i++) {
		widget_list.append(radiobutton_list[i]);
	}

	nWidgets = widget_list.length();
	nSliders = slider_list.length();
    nLabels = label_list.length();
    nCheckBoxes = checkbox_list.length();

	for (int i=0; i<nWidgets; i++) {
		QWidget *w = widget_list[i];
        QString wname = w->objectName();
		if (wname.startsWith("line_")) {
			connect(w, SIGNAL(textChanged(QString)), this, SLOT(changeParam()));
			connect(w, SIGNAL(textChanged(QString)), this, SLOT(redrawDistPlot()));
		}
		if (wname.startsWith("text_")) {
			connect(w, SIGNAL(textChanged(QString)), this, SLOT(changeParam()));
		}
		if (wname.startsWith("spin_")) {
			connect(w, SIGNAL(valueChanged(int)), this, SLOT(changeParam()));
		}
		if (wname.startsWith("comb_")) {
			connect(w, SIGNAL(activated(QString)), this, SLOT(changeParam()));
		}
		if (wname.startsWith("cbox_")) {
			connect(w, SIGNAL(toggled(bool)), this, SLOT(changeParam()));
		}
		if (wname.startsWith("rbut_")) {
			connect(w, SIGNAL(toggled(bool)), this, SLOT(changeParam()));
		}
	}

	QwtPlot *qp;

//    qp = (QwtPlot *)qFindChild<QObject *>(this, "qwtPlot_DIVIDE_TIME");
//    distplot_list[0] = qp;
    qp = (QwtPlot *)qFindChild<QObject *>(this, "qwtPlot_DIVIDE_TIME_1");
    distplot_list[0] = qp;
    qp = (QwtPlot *)qFindChild<QObject *>(this, "qwtPlot_DIVIDE_TIME_2");
    distplot_list[1] = qp;
}


//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
bool MainWindow::getVideoFileInfo(int *nframes, QString *itemFormat, QString *itemCodec, QString *videoFileName)
{
    bool ok;

    int i = QInputDialog::getInteger(this, tr("Set nframes"),tr("Number of frames to capture: "), *nframes, 0, 10000, 1, &ok);
    if (ok) {
        *nframes = i;
    }
    if (!ok || nframes == 0) {
        return false;
    }

    QStringList formatItems;
    formatItems << tr("avi") << tr("mov") << tr("mpg");
    *itemFormat = QInputDialog::getItem(this, tr("QInputDialog::getItem()"),
                                         tr("Video file format:"), formatItems, 0, false, &ok);
    QStringList codecItems;
    codecItems << tr("h264") << tr("mpeg4") << tr("mpeg");
    *itemCodec = QInputDialog::getItem(this, tr("QInputDialog::getItem()"),
                                              tr("Codec:"), codecItems, 0, false, &ok);

    const char *prompt;
    if (itemFormat->contains("avi")) {
        prompt = "Videos (*.avi)";
    } else if (itemFormat->contains("mov")) {
        prompt = "Videos (*.mov)";
    } else if (itemFormat->contains("mpg")) {
        prompt = "Videos (*.mpg)";
    }
    *videoFileName = QFileDialog::getSaveFileName(this,
                                                    tr("Save File"),
                                                    QString(),
                                                    tr(prompt));
    return true;
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow:: startRecorderVTK()
{
    bool ok;
    int nframes=0;
    QString itemFormat, itemCodec, videoFileName;

    ok = getVideoFileInfo(&nframes, &itemFormat, &itemCodec, &videoFileName);
    if (!ok) return;
    goToVTK();
    videoVTK->startRecorder(videoFileName,itemFormat,itemCodec,nframes);
    actionStart_recording_VTK->setEnabled(false);
    actionStop_recording_VTK->setEnabled(true);
    Global::recordingVTK = true;
    started = true;
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow:: stopRecorderVTK()
{
    videoVTK->stopRecorder();
    actionStart_recording_VTK->setEnabled(true);
    actionStop_recording_VTK->setEnabled(false);
    Global::recordingVTK = false;
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow:: startRecorderFACS()
{
    bool ok;
    int nframes=0;
    QString itemFormat, itemCodec, videoFileName;

    ok = getVideoFileInfo(&nframes, &itemFormat, &itemCodec, &videoFileName);
    if (!ok) return;
    goToFACS();
    videoFACS->startRecorder(videoFileName,itemFormat,itemCodec,nframes);
    actionStart_recording_FACS->setEnabled(false);
    actionStop_recording_FACS->setEnabled(true);
    Global::recordingFACS = true;
    LOG_QMSG("startRecorderFACS");
    LOG_QMSG(videoFileName);
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow:: stopRecorderFACS()
{
    videoFACS->stopRecorder();
    actionStart_recording_FACS->setEnabled(true);
    actionStop_recording_FACS->setEnabled(false);
    Global::recordingFACS = false;
    LOG_QMSG("stopRecorderFACS");
}

/*
//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow:: drawDistPlots()
{
	double *x, *prob;
	x = new double[nDistPts];
	prob = new double[nDistPts];
	QwtPlot *qp;
	string name_str;
	QString median_qstr, shape_qstr;
	double median, shape;

    for (int j=0; j<ndistplots; j++) {
		qp = distplot_list[j];
        QString name = qp->objectName();
        if (j == 0) {
            qp->setTitle("Tumour cell division time (hrs)");
            median_qstr = line_DIVIDE_TIME_MEDIAN->text();
            shape_qstr = line_DIVIDE_TIME_SHAPE->text();
        }
		median = median_qstr.toDouble();
		shape = shape_qstr.toDouble();
        create_lognorm_dist(median,shape,nDistPts,x,prob);

        int n = dist_limit(prob,nDistPts);
        double xmax = x[n];
		sprintf(msg,"%f %f %d",median,shape,n);
		for (int i=0;i<40;i++) {
			sprintf(msg,"%d %f %f",i,x[i],prob[i]);
		}
        qp->setAxisScale(QwtPlot::xBottom, 0.0, xmax, 0.0);
        QwtPlotCurve *curve = new QwtPlotCurve("title");
        curve->attach(qp);
        curve->setData(x, prob, n);
        curve_list[j] = curve;
		qp->replot();
	}
	delete [] x;
    x = NULL;
	delete [] prob;
	prob = NULL;
}
*/

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow:: drawDistPlots()
{
    double *x, *prob;
    x = new double[nDistPts];
    prob = new double[nDistPts];
    QwtPlot *qp;
    string name_str;
    QString median_qstr, shape_qstr;
    double median, shape;

    for (int j=0; j<ndistplots; j++) {
        qp = distplot_list[j];
        QString name = qp->objectName();
        if (j == 0) {
            qp->setTitle("Type 1 division time (hrs)");
            median_qstr = line_DIVIDE_TIME_1_MEDIAN->text();
            shape_qstr = line_DIVIDE_TIME_1_SHAPE->text();
        } else if (j == 1) {
            qp->setTitle("Type 2 division time (hrs)");
            median_qstr = line_DIVIDE_TIME_2_MEDIAN->text();
            shape_qstr = line_DIVIDE_TIME_2_SHAPE->text();
        }
        median = median_qstr.toDouble();
        shape = shape_qstr.toDouble();
        create_lognorm_dist(median,shape,nDistPts,x,prob);

        int n = dist_limit(prob,nDistPts);
        double xmax = x[n];
        sprintf(msg,"%f %f %d",median,shape,n);
        for (int i=0;i<40;i++) {
            sprintf(msg,"%d %f %f",i,x[i],prob[i]);
        }
        qp->setAxisScale(QwtPlot::xBottom, 0.0, xmax, 0.0);
        QwtPlotCurve *curve = new QwtPlotCurve("title");
        curve->attach(qp);
        curve->setData(x, prob, n);
        curve_list[j] = curve;
        qp->replot();
    }
    delete [] x;
    x = NULL;
    delete [] prob;
    prob = NULL;
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow:: initFACSPlot()
{
    qpFACS = (QwtPlot *)qFindChild<QObject *>(this, "qwtPlot_FACS");
    qpFACS->clear();
    qpFACS->setTitle("FACS");
//    QwtSymbol symbol = QwtSymbol( QwtSymbol::Diamond, Qt::blue, Qt::NoPen, QSize( 3,3 ) );
    qpFACS->replot();
    connect((QObject *)groupBox_FACS,SIGNAL(groupBoxClicked(QString)),this,SLOT(processGroupBoxClick(QString)));
}

//--------------------------------------------------------------------------------------------------------
// Possible variables to plot are Global::vars_used[]
// Use this to trigger FACS plot.
//--------------------------------------------------------------------------------------------------------
void MainWindow::showFACS()
{
    double xmin, xmax, ymin, ymax, x, y, xscale, yscale;
    int i, ivar, kvar_x, kvar_y, ichemo;
    bool x_logscale, y_logscale;
    QString xlabel, ylabel;
    QRadioButton *rb;

//    LOG_MSG("showFACS");
//    if (Global::showingFACS) LOG_MSG("showingFACS");
//    if (Global::recordingFACS) LOG_MSG("recordingFACS");
    qpFACS = (QwtPlot *)qFindChild<QObject *>(this, "qwtPlot_FACS");
    qpFACS->size();
    qpFACS->clear();
    qpFACS->setTitle("FACS");
//    QwtSymbol symbol = QwtSymbol( QwtSymbol::Diamond, Qt::blue, Qt::NoPen, QSize( 3,3 ) );
    QwtSymbol symbol = QwtSymbol( QwtSymbol::Rect, Qt::blue, Qt::NoPen, QSize( 2,2 ) );

    // Determine which button is checked:
    for (ivar=0; ivar<Global::nvars_used; ivar++) {
        rb = FACS_x_vars_rb_list[ivar];
        if (rb->isChecked()) {
            ichemo = Global::GUI_to_DLL_index[ivar];
            break;
        }
    }
    kvar_x = ivar;
    switch(ichemo) {
    case CFSE:
        xscale = 1000;
        xlabel = "CFSE";
        xmin = 0.1;
        xmax = 1500;
        break;
    case OXYGEN:
        xscale = 1;
        xlabel = "Oxygen";
        xmin = 1.0e-4;
        xmax = 1.0;
        break;
    case GLUCOSE:
        xscale = 1;
        xlabel = "Glucose";
        xmin = 1.0e-3;
        xmax = 10.0;
        break;
    case TRACER:
        break;
    case DRUG_A_PARENT:
        break;
    case DRUG_A_METAB_1:
        break;
    case DRUG_A_METAB_2:
        break;
    case DRUG_B_PARENT:
        break;
    case DRUG_B_METAB_1:
        break;
    case DRUG_B_METAB_2:
        break;
    case GROWTH_RATE:
        xscale = 1;
        xlabel = "Growth rate";
        xmin = 1.0e-7;
        xmax = 1.0e-5;
        break;
    case CELL_VOLUME:
        xscale = 1;
        xlabel = "Cell volume";
        xmin = 0;
        xmax = 2;
        break;
    case O2_BY_VOL:
        xscale = 1;
        xlabel = "O2 x volume";
        xmin = 0;
        xmax = 2;
        break;
    }
    for (ivar=0; ivar<Global::nvars_used; ivar++) {
        rb = FACS_y_vars_rb_list[ivar];
        if (rb->isChecked()) {
            ichemo = Global::GUI_to_DLL_index[ivar];
            break;
        }
    }
    kvar_y = ivar;
    switch(ichemo) {
    case CFSE:
        yscale = 1000;
        ylabel = "CFSE";
        ymin = 0.1;
        ymax = 1500;
        break;
    case OXYGEN:
        yscale = 1;
        ylabel = "Oxygen";
        ymin = 1.0e-4;
        ymax = 1.0;
        break;
    case GLUCOSE:
        yscale = 1;
        ylabel = "Glucose";
        ymin = 1.0e-3;
        ymax = 10.0;
        break;
    case TRACER:
        break;
    case DRUG_A_PARENT:
        break;
    case DRUG_A_METAB_1:
        break;
    case DRUG_A_METAB_2:
        break;
    case DRUG_B_PARENT:
        break;
    case DRUG_B_METAB_1:
        break;
    case DRUG_B_METAB_2:
        break;
    case GROWTH_RATE:
        yscale = 1;
        ylabel = "Growth rate";
        ymin = 1.0e-7;
        ymax = 1.0e-5;
        break;
    case CELL_VOLUME:
        yscale = 1;
        ylabel = "Cell volume";
        ymin = 0;
        ymax = 2;
        break;
    case O2_BY_VOL:
        yscale = 1;
        ylabel = "O2 x volume";
        ymin = 0;
        ymax = 2;
        break;
    }

    x_logscale = checkBox_FACS_log_x->isChecked();
    y_logscale = checkBox_FACS_log_y->isChecked();
//    xmin = 0.1;
//    xmax = 1500;
    double cfse_min = 1.0e20;
    for (i=0; i<Global::nFACS_cells; i++) {
        x = Global::FACS_data[Global::nvars_used*i+kvar_x];
        cfse_min = MIN(x,cfse_min);
        y = Global::FACS_data[Global::nvars_used*i+kvar_y];
        x = xscale*x;
        y = yscale*y;
        y = max(y,1.01*ymin);
//        ymax = max(y,ymax);
        if (x >= xmin) {
            QwtPlotMarker* m = new QwtPlotMarker();
            m->setSymbol( symbol );
            m->setValue( QPointF( x,y ) );
            m->attach( qpFACS );
        }
    }
    qpFACS->setAxisScale(QwtPlot::yLeft, ymin, ymax, 0);
    qpFACS->setAxisTitle(QwtPlot::yLeft, ylabel);
    if (y_logscale) {
        qpFACS->setAxisScaleEngine(QwtPlot::yLeft, new QwtLog10ScaleEngine);
    } else {
        qpFACS->setAxisScaleEngine(QwtPlot::yLeft, new QwtLinearScaleEngine);
    }
    qpFACS->setAxisMaxMinor(QwtPlot::yLeft, 10);
    qpFACS->setAxisMaxMajor(QwtPlot::yLeft, 5);

    qpFACS->setAxisScale(QwtPlot::xBottom, xmin, xmax, 0);
    qpFACS->setAxisTitle(QwtPlot::xBottom, xlabel);
    if (x_logscale) {
        qpFACS->setAxisScaleEngine(QwtPlot::xBottom, new QwtLog10ScaleEngine);
    } else {
        qpFACS->setAxisScaleEngine(QwtPlot::xBottom, new QwtLinearScaleEngine);
    }
    qpFACS->setAxisMaxMinor(QwtPlot::xBottom, 10);
    qpFACS->setAxisMaxMajor(QwtPlot::xBottom, 5);

    qpFACS->replot();
//    sprintf(msg,"x range: %f %f  y range: %f %f",xmin,xmax,ymin,ymax);
//    LOG_MSG(msg);
//    QSize size = qpFACS->size();
//    sprintf(msg,"qpFACS size: %d %d",size.width(),size.height());
//    LOG_MSG(msg);
    if (videoFACS->record) {
        videoFACS->recorder();
    } else if (actionStop_recording_FACS->isEnabled()) {
        actionStart_recording_FACS->setEnabled(true);
        actionStop_recording_FACS->setEnabled(false);
    }
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::test_histo()
{
    int numValues = 20;
    double width = 10, xmin = 0;
    QwtArray<double> values(numValues);
    for (int i=0; i<numValues; i++) {
        values[i] = rand() %100;
    }
    makeHistoPlot(numValues,xmin,width,values);
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow:: initHistoPlot()
{
    qpHistoBar = (QwtPlot *)qFindChild<QObject *>(this, "qwtPlot_Histo");
    qpHistoBar->setTitle("Histogram");
//    QwtSymbol symbol = QwtSymbol( QwtSymbol::Diamond, Qt::blue, Qt::NoPen, QSize( 3,3 ) );
    qpHistoBar->replot();

    qpHistoLine = (QwtPlot *)qFindChild<QObject *>(this, "qwtPlot_HistoLine");
    qpHistoLine->hide();

    connect((QObject *)groupBox_Histo,SIGNAL(groupBoxClicked(QString)),this,SLOT(processGroupBoxClick(QString)));
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::makeHistoPlot(int numValues, double xmin, double width,  QwtArray<double> values)
{
    QwtPlot *plot;
    double pos;
//    LOG_MSG("makeHistoPlot");
    bool use_HistoBar = radioButton_histotype_1->isChecked();
    if (use_HistoBar) {
        plot = qpHistoBar;
        qpHistoLine->hide();
    } else {
        plot = qpHistoLine;
        qpHistoBar->hide();
    }
    plot->clear();
    plot->setCanvasBackground(QColor(Qt::white));
    plot->setTitle("Histogram");

    QwtPlotGrid *grid = new QwtPlotGrid;
    grid->enableXMin(true);
    grid->enableYMin(true);
    grid->setMajPen(QPen(Qt::black, 0, Qt::DotLine));
    grid->setMinPen(QPen(Qt::gray, 0 , Qt::DotLine));
    grid->attach(plot);

    if (use_HistoBar) {
        if (histogram) {
            histogram->detach();
        } else {
            histogram = new HistogramItem();
        }
        histogram->setColor(Qt::darkCyan);

        QwtArray<QwtDoubleInterval> intervals(numValues);

        pos = xmin;
        for ( int i = 0; i < numValues; i++ )
        {
            intervals[i] = QwtDoubleInterval(pos, pos + width);
            pos += width;
        }

        histogram->setData(QwtIntervalData(intervals, values));
        histogram->attach(plot);
    } else {
        double x[100], y[100];
        for ( int i = 0; i < numValues; i++ ) {
            x[i] = xmin + (i + 0.5)*width;
            y[i] = values[i];
        }
        pos = x[numValues-1] + width/2;
        QwtPlotCurve *curve = new QwtPlotCurve("");
        QPen *pen = new QPen();
        pen->setColor(Qt::black);
        curve->attach(plot);
        curve->setPen(*pen);
        curve->setData(x, y, numValues);
    }

    plot->setAxisScale(QwtPlot::yLeft, 0.0, 100.0);
    plot->setAxisScale(QwtPlot::xBottom, xmin, pos);
    plot->replot();
//    plot->resize(600,400);
    plot->show();

}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::on_buttonGroup_celltype_buttonClicked(QAbstractButton* button)
{
    LOG_MSG("on_buttonGroup_celltype_buttonClicked");
    if (button->text() == "Cell type 1") {
        Global::histo_celltype = 1;
    } else if (button->text() == "Cell type 2") {
        Global::histo_celltype = 2;
    } else {
        Global::histo_celltype = 0; // both cell types
    }
    showHisto();
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::on_checkBox_histo_logscale_toggled()
{
    showHisto();
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::on_buttonGroup_histotype_buttonClicked(QAbstractButton* button)
{
    showHisto();
}
//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow:: showHisto()
{
    int ivar, k, k0, numValues;
    QRadioButton *rb;
    QString xlabel;
    double width, xmin;
    bool log_scale;

    log_scale = checkBox_histo_logscale->isChecked();
    numValues = Global::nhisto_bins;
    QwtArray<double> values(numValues);

    // Determine which button is checked:
    for (ivar=0; ivar<Global::nvars_used; ivar++) {
        rb = histo_rb_list[ivar];
        if (rb->isChecked()) {
            break;
        }
    }
    xlabel = Global::var_string[ivar];
    k0 = Global::histo_celltype*numValues*Global::nvars_used;
//    sprintf(msg,"histo_celltype: %d numValues: %d nvars_used: %d k0: %d",Global::histo_celltype,numValues,Global::nvars_used,k0);
//    LOG_MSG(msg);
    if (!Global::histo_data) {
        LOG_MSG("No histo_data");
        return;
    }
    for (int i=0; i<numValues; i++) {
        k = k0 + ivar*numValues + i;
        if (log_scale)
            values[i] = Global::histo_data_log[k];
        else
            values[i] = Global::histo_data[k];
    }
    if (log_scale) {
        xmin = Global::histo_vmin_log[ivar];
        width = (Global::histo_vmax_log[ivar] - Global::histo_vmin_log[ivar])/numValues;
    } else {
        xmin = Global::histo_vmin[ivar];
        width = (Global::histo_vmax[ivar] - Global::histo_vmin[ivar])/numValues;
    }
    makeHistoPlot(numValues,xmin,width,values);
}

//-------------------------------------------------------------
// Loops through the workingParameterList and fills in the GUI.
//-------------------------------------------------------------
void MainWindow::loadParams()
{
	sprintf(msg,"nParams: %d nSliders: %d",nParams,nSliders);
	for (int k=0; k<nSliders; k++) {		// there must be a neater way to do this
		SliderPlus *splus = 0;
		sliderplus_list.append(splus);
		QWidget *w = 0;
		sliderParam.append(w);		
	}
	if (param_to_sliderIndex == NULL) {
		param_to_sliderIndex = new int[nParams];
		for (int i=0; i<nParams; i++)
			param_to_sliderIndex[i] = -1;
	}
	for (int i=0; i<nWidgets; i++) {
		QWidget *w = widget_list[i];							// w = widget_list[i] is the ith widget in the UI
        QString qsname = w->objectName();
		if (qsname.startsWith("line_") || qsname.startsWith("spin_")
			|| qsname.startsWith("comb_") || qsname.startsWith("cbox_")
			|| qsname.startsWith("rbut_") || qsname.startsWith("text_")) {
			QString wtag = qsname.mid(5);
			int rbutton_case = 0;
			if (qsname.startsWith("rbut_")) {
//				parse_rbutton(wtag,&rbutton_case);
                wtag = parse_rbutton(qsname,&rbutton_case);
            }
            // Find corresponding data in workingParameterList
            bool found = false;
			for (int k=0; k<nParams; k++) {
				PARAM_SET p = parm->get_param(k);
				QString ptag = p.tag;		// ptag is the tag of the kth parameter in the list
				if (wtag.compare(ptag) == 0) {
					double vmax = p.maxvalue;
					double vmin = p.minvalue;
                    // Update the widget (line_, spin_ or comb_) with data from the parameter list
                    // ---LineEdits
					if (qsname.startsWith("line_")) {
                        double val = p.value;
						QString val_str = QString::number(val);
						QLineEdit *w_l = (QLineEdit *)w;
                        w_l->setText(val_str);
						if (USE_RANGES) {
							// Set max and min values. If min=max=0, there're no restrictions.
							if (!(vmin == 0 && vmax == 0)) {
//	                            QValidator *aValidator = new QDoubleValidator(vmin, vmax, 10, w_l);
								QValidator *aValidator = new MyDoubleValidator(vmin, vmax, 8, w_l);
								w_l->setValidator(aValidator);
							}
						}						
					} else if (qsname.startsWith("spin_")) {
						double val = p.value;
						QSpinBox *w_s = (QSpinBox *)w;
                        w_s->setValue(val);
						if (!(vmin == 0 && vmax == 0)) {
                            w_s->setMinimum(vmin);
                            w_s->setMaximum(vmax);
						}
						if (qsname.contains("NCPU")) {
							ncpu = p.value;
						}
                    } else if (qsname.startsWith("comb_")) {
                        int val = p.value - 1;	//0-based indexing
						QComboBox *w_c = (QComboBox *)w;
                        w_c->setCurrentIndex(val);
					} else if (qsname.startsWith("cbox_")) {
						QCheckBox *w_cb = (QCheckBox *)w;

                        bool use_OXYGEN = qsname.contains("USE_OXYGEN");
                        if (p.value == 1) {
                            w_cb->setChecked(true);
                            if (use_OXYGEN)
                                enableUseOxygen();
                        } else {
                            w_cb->setChecked(false);
                            if (use_OXYGEN)
                                disableUseOxygen();
                        }
                        bool use_GLUCOSE = qsname.contains("USE_GLUCOSE");
                        if (p.value == 1) {
                            w_cb->setChecked(true);
                            if (use_GLUCOSE)
                                enableUseGlucose();
                        } else {
                            w_cb->setChecked(false);
                            if (use_GLUCOSE)
                                disableUseGlucose();
                        }
                        bool use_TRACER = qsname.contains("USE_TRACER");
                        if (p.value == 1) {
                            w_cb->setChecked(true);
                            if (use_TRACER)
                                enableUseTracer();
                        } else {
                            w_cb->setChecked(false);
                            if (use_TRACER)
                                disableUseTracer();
                        }
//                        bool use_TPZ = qsname.contains("USE_TPZ_DRUG");
                        if (p.value == 1) {
                            w_cb->setChecked(true);
                        } else {
                            w_cb->setChecked(false);
                        }
//                        bool use_DNB = qsname.contains("USE_DNB_DRUG");
                        if (p.value == 1) {
                            w_cb->setChecked(true);
                        } else {
                            w_cb->setChecked(false);
                        }
//                        bool use_TREATMENT_FILE = qsname.contains("USE_TREATMENT_FILE");
//                        if (p.value == 1) {
//                            w_cb->setChecked(true);
//                            if (use_TREATMENT_FILE)
//                                enableUseTreatmentFile();
//                        } else {
//                            w_cb->setChecked(false);
//                            if (use_TREATMENT_FILE)
//                                disableUseTreatmentFile();
//                        }
					} else if (qsname.startsWith("rbut_")) {
						QRadioButton *w_rb = (QRadioButton *)w;
                        if (int(p.value) == rbutton_case) {
							w_rb->setChecked(true);
							LOG_QMSG("loadParams: setChecked true")
						} else {
							w_rb->setChecked(false);
							LOG_QMSG("loadParams: setChecked false")
						}
						setLineEditVisibility(qsname,int(p.value));
					} else if (qsname.startsWith("text_")) {
						QLineEdit *w_l = (QLineEdit *)w;
						w_l->setText(p.label);
					}
					
					// Update Label text (except for "text_" variables)
                    // Get the corresponding label from the label list
                    QString labelString = "label_" + wtag;
					QLabel *label = NULL;
					bool foundLabel = false;
					for (int j=0; j<nLabels; j++) {
						label = label_list[j];
						if (!qsname.startsWith("text_") && labelString.compare(label->objectName()) == 0) {
							foundLabel = true;
							break;
						}
					}										// label is the pointer to the UI label for wtag and ptag
                    QString labelText = p.label;
                    
                    // Hardcode the distribution label names for now
                    if (wtag.compare("DIVIDE_TIME_MEDIAN") == 0)
                        labelText = "Median";
                    else if (wtag.compare("DIVIDE_TIME_SHAPE") == 0)
                        labelText = "Shape";

					bool is_slider = false;
					int j;
					QSlider *s=0;
					QString sliderString;
					for (j=0; j<nSliders; j++) {
						sliderString = "slider_" + wtag;
						s = slider_list[j];
						if (sliderString.compare(s->objectName()) == 0) {
							is_slider = true;					// the jth slider in the list corresponds to wtag and ptag
							break;
						}
					}

					// Try this change to eliminate sliders except for distributions
					if (labelText.compare("Shape") != 0 && labelText.compare("Median") != 0) {
						is_slider = false;
					}

					if (is_slider) {
                        // If there is a slider corresponding to wtag, then just use the label.
						if (foundLabel)
	                        label->setText(labelText);
					} else {
						if (!(vmin == 0 && vmax == 0)) {
                            // If there is no slider, then add min and max values to the label text.
							QString min_str = QString::number(vmin);
							QString max_str = QString::number(vmax);
							if (foundLabel)
		                        label->setText(labelText + "  [ " + min_str + "-" + max_str + " ]");
						} else {
							if (foundLabel)
		                        label->setText(labelText);
						}
					}
						
                    // If there is a corresponding slider for this parameter, then apply settings.
					if (is_slider) {						
                        SliderPlus *splus = new SliderPlus(wtag,vmin,vmax,nTicks,k,i);
						sliderplus_list[j] = splus;
                        int ival = splus->val_to_int(p.value);
                        s->setMinimum(0);
                        s->setMaximum(splus->nTicks());
						s->setSliderPosition(ival);
						sliderParam[j] = w;
                        connect(s, SIGNAL(valueChanged(int)), this, SLOT(updateSliderBox())); //sliderReleased               
                        param_to_sliderIndex[k] = j;
					}                  
                    found = true;
                    break;

					if (!found) {
						sprintf(msg,"%s was not found in the parameter list",(wtag.toStdString()).data());
						LOG_MSG(msg);
					}
				}
			}
		}
	}
//    setTreatmentFileUsage();
    text_GUI_VERSION_NAME->setText(Global::GUI_build_version);
    text_DLL_VERSION_NAME->setText(Global::DLL_build_version);
}

//--------------------------------------------------------------------------------------------------------
// This is to disable unused fields (because spheroid_GUI.ui is shared with spheroid_GUI).
// The MEDIUM_VOLUME field value is computed from the specified DELTA_X and NXB.
//--------------------------------------------------------------------------------------------------------
void MainWindow::setFields()
{
    line_NT_CONC->setEnabled(false);
    line_NMM3->setEnabled(false);
    line_FLUID_FRACTION->setEnabled(false);
    line_UNSTIRRED_LAYER->setEnabled(false);
    cbox_USE_RELAX->setEnabled(false);
    cbox_USE_PAR_RELAX->setEnabled(false);
    groupBox_drop->setEnabled(false);
    int nxb = line_NXB->text().toInt();
    double dx = line_DELTA_X->text().toDouble();
    double vol_cm3 = pow(4*nxb*dx,3)*1.0e-12;
    QString str = QString::number(vol_cm3,'g',3);
    line_MEDIUM_VOLUME->setText(str);
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
//QString MainWindow::parse_rbutton(QString wtag, int *rbutton_case)
//{
//	// parse wtag into part before '_' and part after '_'
//	int j = wtag.indexOf('_');
//	QString suffix = wtag.mid(j+1);
//	// the prefix becomes wtag, the suffix becomes rbutton_case, an integer 0,1,2,...
//	wtag = wtag.mid(0,j);
//	bool ok;
//	*rbutton_case = suffix.toInt(&ok);
//	return wtag;
//}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
QString MainWindow::parse_rbutton(QString qsname, int *rbutton_case)
{
    // parse wtag into part before '_' and part after '_'
    QString wtag = qsname.mid(5);   // strips off "rbut_"
    int j = wtag.lastIndexOf('_');  // position of last '_'
    QString suffix = wtag.mid(j+1);
    // the prefix becomes wtag0, the suffix becomes rbutton_case, an integer 0,1,2,...
    QString wtag0 = wtag.mid(0,j);
    bool ok;
    *rbutton_case = suffix.toInt(&ok);
    return wtag0;
}

//--------------------------------------------------------------------------------------------------------
// Only change the LineEdit visibility is one is already enabled
//--------------------------------------------------------------------------------------------------------
void MainWindow::setLineEditVisibility(QString wname, int val)
{
    QString lineRateName, lineConcName;
	// Need to locate the LineEdit widgets from their names
	QWidget *w_rate=0, *w_conc=0;	//, *w_cb, *w_cb2=0;
	for (int i=0; i<nWidgets; i++) {
		QWidget *w = widget_list[i];							// w = widget_list[i] is the ith widget in the UI
		QString qsname = w->objectName();
		if (qsname.contains(lineRateName)) {
			w_rate = w;
		}
		if (qsname.contains(lineConcName)) {
			w_conc = w;
		}
	}
	if (w_rate->isEnabled() || w_conc->isEnabled()) {
		if (val == 0) {
			w_rate->setEnabled(false);
			w_conc->setEnabled(true);
		} else {
			w_rate->setEnabled(true);
			w_conc->setEnabled(false);
		}
	}
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::reloadParams()
{
	for (int i=0; i<nWidgets; i++) {
		QWidget *w = widget_list[i];							// w = widget_list[i] is the ith widget in the UI
        QString qsname = w->objectName();
		if (qsname.startsWith("line_") || qsname.startsWith("spin_") 
			|| qsname.startsWith("comb_") || qsname.startsWith("cbox_")
			|| qsname.startsWith("rbut_") || qsname.startsWith("text_")) {
			QString wtag = qsname.mid(5);
			int rbutton_case = 0;
			if (qsname.startsWith("rbut_")) {
//				parse_rbutton(wtag,&rbutton_case);
                wtag = parse_rbutton(qsname,&rbutton_case);
            }
            // Find corresponding data in workingParameterList
            bool found = false;
			for (int k=0; k<nParams; k++) {
				PARAM_SET p = parm->get_param(k);
				QString ptag = p.tag;		// ptag is the tag of the kth parameter in the list
				if (wtag.compare(ptag) == 0) {
					found = true;
                    // Update the widget (line_, spin_ or comb_) with data from the parameter list
					if (qsname.startsWith("line_")) {
                        double val = p.value;
						QString val_str = QString::number(val);
						QLineEdit *w_l = (QLineEdit *)w;
                        w_l->setText(val_str);
					} else if (qsname.startsWith("text_")) {
						QLineEdit *w_l = (QLineEdit *)w;
						w_l->setText(p.label);
					} else if (qsname.startsWith("spin_")) {
						double val = p.value;
						QSpinBox *w_s = (QSpinBox *)w;
                        w_s->setValue(val);
						if (qsname.contains("NCPU")) {
							ncpu = p.value;
						}
					} else if (qsname.startsWith("comb_")) {
                        int val = p.value - 1;	//0-based indexing
						QComboBox *w_c = (QComboBox *)w;
                        w_c->setCurrentIndex(val);
					} else if (qsname.startsWith("cbox_")) {
						QCheckBox *w_cb = (QCheckBox *)w;
                        LOG_QMSG(qsname);

                        bool use_OXYGEN = qsname.contains("USE_OXYGEN");
                        if (p.value == 1) {
                            w_cb->setChecked(true);
                            if (use_OXYGEN)
                                enableUseOxygen();
                        } else {
                            w_cb->setChecked(false);
                            if (use_OXYGEN)
                                disableUseOxygen();
                        }
                        bool use_GLUCOSE = qsname.contains("USE_GLUCOSE");
                        if (p.value == 1) {
                            w_cb->setChecked(true);
                            if (use_GLUCOSE)
                                enableUseGlucose();
                        } else {
                            w_cb->setChecked(false);
                            if (use_GLUCOSE)
                                disableUseGlucose();
                        }
                        bool use_TRACER = qsname.contains("USE_TRACER");
                        if (p.value == 1) {
                            w_cb->setChecked(true);
                            if (use_TRACER)
                                enableUseTracer();
                        } else {
                            w_cb->setChecked(false);
                            if (use_TRACER)
                                disableUseTracer();
                        }
//                        bool use_TPZ = qsname.contains("USE_TPZ_DRUG");
                        if (p.value == 1) {
                            w_cb->setChecked(true);
//                            if (use_TPZ)
//                                enableUseTPZ();
                        } else {
                            w_cb->setChecked(false);
//                            if (use_TPZ)
//                                disableUseTPZ();
                        }
//                        bool use_DNB = qsname.contains("USE_DNB_DRUG");
                        if (p.value == 1) {
                            w_cb->setChecked(true);
//                            if (use_DNB)
//                                enableUseDNB();
                        } else {
                            w_cb->setChecked(false);
//                            if (use_DNB)
//                                disableUseDNB();
                        }
//                        bool use_TREATMENT_FILE = qsname.contains("USE_TREATMENT_FILE");
//                        if (p.value == 1) {
//                            w_cb->setChecked(true);
//                            if (use_TREATMENT_FILE)
//                                enableUseTreatmentFile();
//                        } else {
//                            w_cb->setChecked(false);
//                            if (use_TREATMENT_FILE)
//                                disableUseTreatmentFile();
//                        }
					} else if (qsname.startsWith("rbut_")) {
                        QRadioButton *w_rb = (QRadioButton *)w;
                        if (p.value == rbutton_case) {
                            w_rb->setChecked(true);
                        } else {
                            w_rb->setChecked(false);
                        }
					}
				}
			}
//			if (!found) {
//				LOG_MSG("Widget tag not found:");
//				LOG_QMSG(qsname);
//				LOG_QMSG(wtag);
//			}
		}
    }
    text_GUI_VERSION_NAME->setText(Global::GUI_build_version);
    text_DLL_VERSION_NAME->setText(Global::DLL_build_version);
    setFields();
}

//--------------------------------------------------------------------------------------------------------
// When we want to uncheck a RadioButton, it is necessary to check some other RB.
// In this case there is only one other RB.
//--------------------------------------------------------------------------------------------------------
void MainWindow::setBdryRadioButton(QRadioButton *w_rb, int val)
{
	if (val == 1) {
		w_rb->setChecked(true);
	} else {
		QButtonGroup *bg = (QButtonGroup *)w_rb->group();
		QRadioButton *rb1 = (QRadioButton *)bg->buttons().last();
		rb1->setChecked(true);
	}
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::showMore(QString moreText)
{
//    LOG_QMSG("showMore " + moreText);
	
    text_more->setEnabled(true);
    text_more->setText(moreText); // text_description
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::writeout()
{
    int ndrugs;
    QString line;
    QFile file(inputFile);
    if (!file.open(QFile::WriteOnly | QFile::Text)) {
        QMessageBox::warning(this, tr("Application"),
                             tr("Cannot write file %1:\n%2.")
                             .arg(inputFile)
                             .arg(file.errorString()));
        LOG_MSG("File open failed");
        return;
    }
    QTextStream out(&file);
    for (int k=0; k<parm->nParams; k++) {
        PARAM_SET p = parm->get_param(k);
        double val = p.value;
        bool is_text = false;
        if (p.tag.contains("_NAME")) {
            is_text = true;
            line = p.label;
        } else if (val == int(val)) 	// whole number, write as integer
            line = QString::number(int(val));
        else
            line = QString::number(val);
        int nch = line.length();
        for (int i=0; i<max(12-nch,1); i++)
            line += " ";
        line += p.tag;
        nch = line.length();
        for (int i=0; i<max(50-nch,1); i++)
            line += " ";
        if (is_text)
            line += p.text;
        else
            line += p.label;
        line += "\n";
        out << line;
        if (p.tag.contains("SAVE_PROFILE_DATA_NUMBER")) {   // insert the drug data here, before plot data
            ndrugs = 0;
            if (ProtocolUsesDrug()) {
                if (cbox_USE_DRUG_A->isChecked()) ndrugs++;
                if (cbox_USE_DRUG_B->isChecked()) ndrugs++;
            }
            line = QString::number(ndrugs);
            int nch = line.length();
            for (int k=0; k<max(16-nch,1); k++)
                line += " ";
            line += "NDRUGS_USED\n";
            out << line;
            if (ndrugs > 0) {
                if (cbox_USE_DRUG_A->isChecked()) {
                    writeDrugParams(&out,DRUG_A);
                }
                if (cbox_USE_DRUG_B->isChecked()) {
                    writeDrugParams(&out,DRUG_B);
                }
            }
        }
    }
    SaveProtocol(&out,ndrugs);
    file.close();
    paramSaved = true;
    LOG_MSG("Input data saved");
}

/*
//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::writeout()
{
	QString line;
    QFile file(inputFile);
    if (!file.open(QFile::WriteOnly | QFile::Text)) {
        QMessageBox::warning(this, tr("Application"),
                             tr("Cannot write file %1:\n%2.")
                             .arg(inputFile)
                             .arg(file.errorString()));
		LOG_MSG("File open failed");
        return;
    }
    QTextStream out(&file);
	for (int k=0; k<parm->nParams; k++) {
        PARAM_SET p = parm->get_param(k);
		double val = p.value;
        bool is_text = false;
        if (p.tag.contains("_NAME")) {
            is_text = true;
            line = p.label;
        } else if (val == int(val)) 	// whole number, write as integer
			line = QString::number(int(val));
		else
			line = QString::number(val);
		int nch = line.length();
		for (int i=0; i<max(12-nch,1); i++)
			line += " ";
        line += p.tag;
        nch = line.length();
        for (int i=0; i<max(50-nch,1); i++)
            line += " ";
        if (is_text)
            line += p.text;
        else
            line += p.label;
        line += "\n";
		out << line;
        if (p.tag.contains("SAVE_PROFILE_DATA_NUMBER")) {   // insert the drug data here, before plot data
            int ndrugs = 0;
            if (cbox_USE_DRUG_A->isChecked()) ndrugs++;
            if (cbox_USE_DRUG_B->isChecked()) ndrugs++;
            line = QString::number(ndrugs);
            int nch = line.length();
            for (int k=0; k<max(16-nch,1); k++)
                line += " ";
            line += "NDRUGS_USED\n";
            out << line;
            if (cbox_USE_DRUG_A->isChecked()) {
                writeDrugParams(&out,DRUG_A);
            }
            if (cbox_USE_DRUG_B->isChecked()) {
                writeDrugParams(&out,DRUG_B);
            }
        }
	}
    file.close();
    SaveProtocol(inputFile);
    paramSaved = true;
	LOG_MSG("Input data saved");
}
*/
//--------------------------------------------------------------------------------------------------------
// Note: To be treated as text, a parameter tag must contain "_NAME"
//--------------------------------------------------------------------------------------------------------
void MainWindow::readInputFile()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open ..."), ".", tr("Input Files (*.inp)"));
	if (fileName.compare("") == 0)
		return;
    paramSaved = false;
    QFile file(fileName);
    if (!file.open(QFile::ReadOnly | QFile::Text)) {
        QMessageBox::warning(this, tr("Application"),
                             tr("Cannot read file %1:\n%2.")
                             .arg(fileName)
                             .arg(file.errorString()));
        return;
    }

    QTextStream in(&file);
	QString line;
	for (int k=0; k<parm->nParams; k++) {
		line = in.readLine();
		QStringList data = line.split(" ",QString::SkipEmptyParts);
		PARAM_SET p = parm->get_param(k);
		QString ptag = p.tag;
        if (ptag.contains("GUI_VERSION")) {
            if (data[0] != Global::GUI_build_version) {
                QMessageBox::warning(this, tr("Application"),
                                     tr("Incorrect GUI version in input file: %1\nthis program is version: %2")
                                     .arg(data[0])
                                     .arg(Global::GUI_build_version));
                return;
            }
        }
        if (ptag.contains("DLL_VERSION")) {
            if (data[0] != Global::DLL_build_version) {
                QMessageBox::warning(this, tr("Application"),
                                     tr("Incorrect DLL version in input file: %1\nthis program was built with DLL version: %2")
                                     .arg(data[0])
                                     .arg(Global::DLL_build_version));
                return;
            }
        }
        if (ptag.contains("_NAME")) {
            parm->set_label(k,data[0]);
        } else {
			parm->set_value(k,data[0].toDouble());
		}
        if (p.tag.contains("SAVE_PROFILE_DATA_NUMBER")) {   // drug data follows, before plot data
            readDrugData(&in);
        }
    }

    LoadProtocol(fileName);

    reloadParams();
    paramSaved = true;
	inputFile = fileName;
    alabel_casename->setText(inputFile);
}

//--------------------------------------------------------------------------------------------------------
// NOT MAINTAINED
//--------------------------------------------------------------------------------------------------------
void MainWindow::loadResultFile()
{
	RESULT_SET *R;

	R = new RESULT_SET;
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open ..."), ".", tr("Result Files (*.res)"));
	if (fileName.compare("") == 0)
		return;
	R->casename = QFileInfo(fileName).baseName();
	for(int i=0; i<result_list.size(); i++) {
		if (R->casename.compare(result_list[i]->casename) == 0) {
			QMessageBox::warning(this, tr("Open results"),
                             tr("This result file is already loaded"));
			return;
		}
	}
    QFile file(fileName);
    if (!file.open(QFile::ReadOnly | QFile::Text)) {
        QMessageBox::warning(this, tr("Application"),
                             tr("Cannot read file %1:\n%2.")
                             .arg(fileName)
                             .arg(file.errorString()));
        return;
    }

    QTextStream in(&file);
	R->nsteps = 0;
	bool indata = false;
	QString line;
	do {
		line = in.readLine();
		if (line.length() > 0) {
			QStringList datalist = line.split(" ",QString::SkipEmptyParts);
			if (indata) {
				R->nsteps++;
			}
			if (datalist[0].contains("========")) {
				indata = true;
			}
		}
	} while (!line.isNull());
	in.seek(0);

	R->tnow = new double[R->nsteps];
    for (int i=0; i<nGraphs; i++) {
        if (!grph->isTimeseries(i)) continue;
        if (!grph->isActive(i)) continue;
        R->pData[i] = new double[R->nsteps];
    }
    indata = false;
	do {
		line = in.readLine();
		if (line.length() > 0) {
			QStringList dataList = line.split(" ",QString::SkipEmptyParts);
			if (indata) {
				int ndata = dataList.length();
				double *data = new double[ndata];
				for (int k=0; k<ndata; k++)
					data[k] = dataList[k].toDouble();
				step++;
				if (step >= R->nsteps) {
					LOG_MSG("ERROR: loadResultFile: step >= nsteps_p");
					return;
				}
				R->tnow[step] = step;		//data[1];step

				for (int i=0; i<nGraphs; i++) {
                    if (!grph->isTimeseries(i)) continue;
                    if (!grph->isActive(i)) continue;
					int k = grph->get_dataIndex(i);
					R->pData[i][step] = data[k]*grph->get_scaling(i);
				}
			}
			if (dataList[0].contains("========")) {
				indata = true;
			}
		}
	} while (!line.isNull());

	// Compute the maxima
    for (int i=0; i<nGraphs; i++) {
        if (!grph->isTimeseries(i)) continue;
        if (!grph->isActive(i)) continue;
        double maxval = getMaximum(R,R->pData[i]);
        grph->set_maxValue(i,maxval);
        R->maxValue[i] = maxval;
    }
	// Now add the result set to the list
	result_list.append(R);
	if (nGraphCases == 0) {
		if (show_outputdata)
			box_outputData = 0;
		initializeGraphs(R);
		drawGraphs();
		goToOutputs();
	}
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
bool MainWindow::save()
{
	writeout();
	return true;
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
bool MainWindow::saveAs()
{
    // show the file dialog
    const QString fileName = QFileDialog::getSaveFileName(this, tr("Select Input File"), ".", tr("Input Files (*.inp)"));
	if (fileName.compare("") != 0) {
		LOG_MSG("Selected file:");
		LOG_QMSG(fileName);
		inputFile = fileName;
        alabel_casename->setText(inputFile);
        writeout();
//        QDir dir = QDir(fileName);
//        QString currPath = dir.absoluteFilePath(fileName);
        QString currPath = QFileInfo(fileName).absolutePath();
        QDir::setCurrent(currPath);
        LOG_QMSG("currPath: " + currPath);
	}
    // Otherwise if user chooses cancel ...
	return true;
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
double MainWindow::getMaximum(RESULT_SET *R, double *x)
{
	double maxx = 0;
	for (int i=0; i<R->nsteps; i++)
		maxx = max(maxx,x[i]);
	return maxx;
}

//-------------------------------------------------------------
// Switches to the input screen
// For when the inputs are being displayed
//-------------------------------------------------------------
void MainWindow::goToInputs()
{
    stackedWidget->setCurrentIndex(0);
    Global::showingVTK = false;
    Global::showingFACS = false;
    action_inputs->setEnabled(false);
    action_outputs->setEnabled(true);
    action_VTK->setEnabled(true);
    action_FACS->setEnabled(true);
    action_field->setEnabled(true); // was commented
}

//-------------------------------------------------------------
// Switches to the output screen
//-------------------------------------------------------------
void MainWindow::goToOutputs()
{
    stackedWidget->setCurrentIndex(1);    
    Global::showingVTK = false;
    Global::showingFACS = false;
    action_outputs->setEnabled(false);
    action_inputs->setEnabled(true);
    action_VTK->setEnabled(true);
    action_FACS->setEnabled(true);
    action_field->setEnabled(true);
}

//-------------------------------------------------------------
// Switches to the VTK screen
//-------------------------------------------------------------
void MainWindow::goToVTK()
{
    stackedWidget->setCurrentIndex(2);
    action_outputs->setEnabled(true);
    action_inputs->setEnabled(true);
    action_field->setEnabled(true);
    action_VTK->setEnabled(false);
    action_FACS->setEnabled(true);
    Global::showingVTK = true;
    Global::showingFACS = false;
}

//-------------------------------------------------------------
// Switches to the FACS screen
//-------------------------------------------------------------
void MainWindow::goToFACS()
{
    stackedWidget->setCurrentIndex(4);
    action_outputs->setEnabled(true);
    action_inputs->setEnabled(true);
    action_field->setEnabled(true);
    action_VTK->setEnabled(true);
    action_FACS->setEnabled(false);
    Global::showingVTK = false;
    Global::showingFACS = true;
}

//-------------------------------------------------------------
// Switches to the field screen
//-------------------------------------------------------------
void MainWindow::goToField()
{
    stackedWidget->setCurrentIndex(3);
    action_outputs->setEnabled(true);
    action_inputs->setEnabled(true);
    action_field->setEnabled(false);
    action_VTK->setEnabled(true);
    action_FACS->setEnabled(true);
    Global::showingVTK = false;
    Global::showingFACS = false;
    LOG_MSG("goToField");
    int res;
    field->displayField(-1,&res);
    if (res != 0) {
        LOG_QMSG("displayField error")
//        stopServer();
    }
}

//-------------------------------------------------------------
// Load and play stored cell position data
//-------------------------------------------------------------
void MainWindow::playVTK()
{
	LOG_MSG("playVTK");
	// Select a file to play
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open ..."), ".", tr("Cell Path Files (*.pos)"));
	if (fileName.compare("") == 0)
		return;
    QMessageBox::StandardButton reply;
    reply = QMessageBox::question(this, tr("Animation player"),
									tr("Save animation frames to image files (.jpg)?"),
                                    QMessageBox::Yes | QMessageBox::No | QMessageBox::Cancel);
	bool save_image;
	if (reply == QMessageBox::Yes)
		save_image = true;
	else
		save_image = false;
	started = true;
	goToVTK();
	if (!vtk->startPlayer(QFileInfo(fileName).absoluteFilePath(), timer, save_image)) {
		LOG_MSG("startPlayer failed");
		errorPopup("Open failure on this file");
		return;
//		exit(1);
	}
    connect(timer, SIGNAL(timeout()), this, SLOT(timer_update()));
    timer->start(tickVTK);
	action_stop->setEnabled(true);
    action_pause->setEnabled(true);
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::timer_update()
{
	ntimes ++;
//	if (!vtk->nextFrame()) {
//		LOG_MSG("Player completed");
//		action_stop->setEnabled(false);
//		action_pause->setEnabled(false);
//	}
}


//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::setVTKSpeed()
{
    bool ok;
    int i = QInputDialog::getInt(this, tr("Set speed"),
		tr("Player timer tick (ms): "), tickVTK, 10, 10000, 1, &ok);
	if (ok) {
		tickVTK = i;
	}
	if (started) {
		timer->stop();
		timer->start(tickVTK);
	}
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::setSavePosStart()
{
    bool ok;
    int i = QInputDialog::getInt(this, tr("Set savepos start"),
		tr("Start recording cell positions at (hours): "), savepos_start, 0, 1000, 1, &ok);
	if (ok) {
		savepos_start = i;
		sprintf(msg,"savepos_start: %d",savepos_start);
		LOG_MSG(msg);
	}
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::saveSnapshot()
{
	LOG_MSG("saveSnapshot");
    QString imageFileName = QFileDialog::getSaveFileName(this, tr("Select image file"), ".",
		tr("Image files (*.png *.jpg *.tif *.bmp)"));    
    if (imageFileName.compare("") == 0) {
		goToVTK();
		return;
	}
    QFileInfo fi(imageFileName);
	QString imgType = fi.suffix();
	goToVTK();

    // Optionally save the snapshot as cell locations in a text file
    QString locationFileName;
    QMessageBox::StandardButton reply;
    reply = QMessageBox::question(this, "Cell positions", "Save cell positions?", QMessageBox::Yes|QMessageBox::No);
    if (reply == QMessageBox::Yes) {
        qDebug() << "Yes was clicked";
        locationFileName = QFileDialog::getSaveFileName(this, tr("Select cell location file"), ".",
            tr("Text files (*.out *.txt)"));
    } else {
        qDebug() << "Yes was *not* clicked";
        locationFileName = "";
    }

    vtk->saveSnapshot(imageFileName,imgType,locationFileName);
}

//--------------------------------------------------------------------------------------------------------
// Note that constituent data is currently hard-wired!!
// The names are NOT in field->cell_const_name[]
//--------------------------------------------------------------------------------------------------------
void MainWindow::saveProfileData()
{
    int i, ichemo;
    double x, c;
    QString line;

    LOG_MSG("saveProfileData");
    QString dataFile = QFileDialog::getSaveFileName(this, tr("Select profile data file"), ".", tr("Data files (*.txt *.dat)"));
    if (dataFile.compare("") == 0) {
        return;
    }
    QFile file(dataFile);
    if (!file.open(QFile::WriteOnly | QFile::Text)) {
        QMessageBox::warning(this, tr("saveProfileData"),
                             tr("Cannot write data file %1:\n%2.")
                             .arg(dataFile)
                             .arg(file.errorString()));
        LOG_MSG("Profile data file open failed");
        return;
    }
    QTextStream out(&file);
//    for (ichemo=0; ichemo<Global::MAX_CHEMO+1; ichemo++) {
//        if (!field->const_used[ichemo]) continue;
//        out << field->cell_const_name[ichemo];
//        out << "\n";
//    }
//    out << "\n";
    for (i=0; i<Global::conc_nc; i++) {
        x = i*Global::conc_dx*1.0e4;
        line = QString::number(x,'g',4);
        line += " ";
        for (ichemo=0; ichemo<Global::MAX_CHEMO+1; ichemo++) {
//            if (!field->const_used[ichemo]) continue;
            c = Global::concData[i*(Global::MAX_CHEMO+1)+ichemo];
            line += QString::number(c,'g',4);
            line += " ";
        }
        line += "\n";
        out << line;
    }
    file.close();
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::showGradient2D()
{
    LOG_MSG("showGradient2D");
    SimpleView2D *mySimpleView2D = new SimpleView2D();
//    QSize size = mySimpleView2D->size();
//    sprintf(msg,"mySimpleView2D size: %d %d",size.height(),size.width());
//    LOG_MSG(msg);
    mySimpleView2D->chooseParameters();
    mySimpleView2D->create();
    mySimpleView2D->show();
    mySimpleView2D->aimCamera();
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::showGradient3D()
{
    LOG_MSG("showGradient3D");
    SimpleView3D *mySimpleView3D = new SimpleView3D();
//    QSize size = mySimpleView3D->size();
//    sprintf(msg,"mySimpleView3D size: %d %d",size.height(),size.width());
//    LOG_MSG(msg);
    mySimpleView3D->chooseParameters();
    mySimpleView3D->create();
    mySimpleView3D->show();
    mySimpleView3D->aimCamera();
//    mySimpleView3D->GetRenderWindow()->SetSize(768,768);
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::runServer()
{
	if (paused) {
		if (vtk->playing) {
			vtk->playon();
		} else {
			exthread->unpause();
		}
        action_run->setEnabled(false);
        action_pause->setEnabled(true);
        action_stop->setEnabled(true);
        action_save_3D_snapshot->setEnabled(false);
        action_save_profile_data->setEnabled(false);
        action_show_gradient3D->setEnabled(false);
        action_show_gradient2D->setEnabled(false);
        paused = false;
        return;
    } else {
 //       field = new Field(page_2D, checkBox_record2D->isChecked());
        field->setSaveImages(checkBox_record2D->isChecked());
        field->setUseLogScale(checkBox_O2logscale->isChecked());
        if (radioButton_hypoxia_1->isChecked())
            Global::i_hypoxia_cutoff = 1;
        else if (radioButton_hypoxia_2->isChecked())
            Global::i_hypoxia_cutoff = 2;
        else if (radioButton_hypoxia_3->isChecked())
            Global::i_hypoxia_cutoff = 3;
        if (radioButton_growthfraction_1->isChecked())
            Global::i_growth_cutoff = 1;
        else if (radioButton_growthfraction_2->isChecked())
            Global::i_growth_cutoff = 2;
        else if (radioButton_growthfraction_3->isChecked())
            Global::i_growth_cutoff = 3;
    }

	if (!paramSaved) {
		int response = QMessageBox::critical(this, tr("ABM Model GUI"), \
					tr("The document has been modified.\nPlease save changes before continuing."), \
					QMessageBox::Save | QMessageBox::Cancel); // | Qt.QMessageBox.Discard
		if (response == QMessageBox::Save) {
            save();
		} else if (response == QMessageBox::Cancel) {
            return;
		}
	}
	
	if (!first) {
		int response = QMessageBox::question(this, tr("ABM Model GUI"),
						tr("Would you like to clear the graphs from the previous run?"),
						QMessageBox::Yes | QMessageBox::No);
		if (response == QMessageBox::Yes)
            mdiArea->closeAllSubWindows();
		else if (response == QMessageBox::Cancel)
            return;
	}
	
    // Display the outputs screen
    if (Global::showingVTK) {
        goToVTK();
    } else if(Global::showingFACS) {
        goToFACS();
    } else {
        goToOutputs();
    }
    // Disable parts of the GUI
    action_run->setEnabled(false);
    action_pause->setEnabled(true);
    action_stop->setEnabled(true);
    action_inputs->setEnabled(true);
    action_VTK->setEnabled(true);
    action_FACS->setEnabled(true);
    action_save_3D_snapshot->setEnabled(false);
    action_save_profile_data->setEnabled(false);
    action_show_gradient3D->setEnabled(false);
    action_show_gradient2D->setEnabled(false);
    action_field->setEnabled(true);
    tab_tumour->setEnabled(false);
//    tab_DC->setEnabled(false);
    tab_chemo->setEnabled(false);
    tab_run->setEnabled(false);

	if (show_outputdata)
	    box_outputData = new QTextBrowser();
	else
		box_outputData = 0;

	if (use_CPORT1) {

		// Port 5001
		sthread1 = new SocketHandler(CPORT1);
		connect(sthread1, SIGNAL(sh_output(QString)), this, SLOT(outputData(QString)));
		sthread1->start();
	}

	// Port 5000
	sthread0 = new SocketHandler(CPORT0);
	connect(sthread0, SIGNAL(sh_output(QString)), box_outputLog, SLOT(append(QString))); //self.outputLog)
	connect(sthread0, SIGNAL(sh_connected()), this, SLOT(preConnection()));
	connect(sthread0, SIGNAL(sh_disconnected()), this, SLOT(postConnection()));
	sthread0->start();
	vtk->cleanup();
    sleep(100);

	hours = 0;
    Global::nt_vtk = 0;
	for (int k=0; k<parm->nParams; k++) {
		PARAM_SET p = parm->get_param(k);
		if (p.tag.compare("NDAYS") == 0) {
			hours = p.value*24;
		}
		if (p.tag.compare("NT_ANIMATION") == 0) {
            Global::nt_vtk = p.value;
		}
	}
	started = true;
	exthread = new ExecThread(inputFile);
	connect(exthread, SIGNAL(display()), this, SLOT(displayScene()));
    connect(exthread, SIGNAL(summary(int)), this, SLOT(showSummary(int)));
    connect(exthread, SIGNAL(facs_update()), this, SLOT(showFACS()));
    connect(this, SIGNAL(facs_update()), this, SLOT(showFACS()));
    connect(exthread, SIGNAL(histo_update()), this, SLOT(showHisto()));
    connect(this, SIGNAL(histo_update()), this, SLOT(showHisto()));
    connect(exthread, SIGNAL(setupC()), this, SLOT(setupConstituents()));
    connect(exthread, SIGNAL(badDLL(QString)), this, SLOT(reportBadDLL(QString)));
    exthread->ncpu = ncpu;
    exthread->nsteps = int(hours*60/Global::DELTA_T);
	exthread->paused = false;
	exthread->stopped = false;
    firstDisplay = true;
    LOG_MSG("exthread->start");
    exthread->start();
}


//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::preConnection()
{
	LOG_MSG("preConnection");

    double hours = 0;
	for (int k=0; k<parm->nParams; k++) {
		PARAM_SET p = parm->get_param(k);
		if (p.tag.compare("NDAYS") == 0) {
			hours = p.value*24;
			break;
		}
	}
	// We assume that the model output is at hourly intervals
	newR = new RESULT_SET;
	QString casename = QFileInfo(inputFile).baseName();
	vtkfile = casename + ".pos";
	newR->casename = casename;
    LOG_QMSG(newR->casename);
	int nsteps = int(hours+1.5);
	newR->nsteps = nsteps;
	newR->tnow = new double[nsteps];

    for (int i=0; i<grph->nGraphs; i++) {
 //       if (!grph->isActive(i)) continue;
		newR->pData[i] = new double[nsteps];
		newR->pData[i][0] = 0;
	}
	LOG_MSG("preconnection: Allocated result set arrays");

	newR->tnow[0] = 0;	// These are not the right initial values
    step = -1;

	// Initialize graphs
    initializeGraphs(newR);
    LOG_MSG("did initializeGraphs");
    Global::nhisto_bins = lineEdit_nhistobins->text().toInt();
    posdata = false;
	LOG_MSG("preconnection: done");
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::reportBadDLL(QString dll_version)
{
    QMessageBox::warning(this, tr("Application"),
                         tr("The version of the DLL linked: %1\nis different from the version the GUI was built with: %2")
                         .arg(dll_version)
                         .arg(Global::DLL_build_version));
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::errorPopup(QString errmsg)
{
	LOG_QMSG(errmsg);
	QMessageBox msgBox;
	msgBox.setText(errmsg);
	msgBox.exec();
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::initializeGraphs(RESULT_SET *R)
{
    LOG_MSG("initializeGraphs");
	mdiArea->closeAllSubWindows();
	mdiArea->show();
    setGraphsActive();
    int non_ts = 0;
//    if (field->isConcPlot()) non_ts++;
//    if (field->isVolPlot()) non_ts++;
//    if (field->isOxyPlot()) non_ts++;
    grph->makeGraphList(non_ts);
    nGraphs = grph->nGraphs;
    if (nGraphCases > 0) {
        clearAllGraphs();
    }
    QString tag;
    QString title;
    QString yAxisTitle;
    for (int i=0; i<nGraphs; i++) {
        if (!grph->isTimeseries(i) && !grph->isProfile(i) && !grph->isDistribution(i)) continue;   // ???
        tag = grph->get_tag(i);
        title = grph->get_title(i);
        yAxisTitle = grph->get_yAxisTitle(i);
        if (pGraph[i] != NULL) {
            pGraph[i]->deleteLater();
            pGraph[i] = NULL;
        }
        if (pGraph[i] == NULL) {
            pGraph[i] = new Plot(tag,R->casename);
            pGraph[i]->setTitle(title);
            pGraph[i]->setAxisTitle(QwtPlot::yLeft, yAxisTitle);
//            LOG_QMSG(title);
//            LOG_QMSG(tag);
        }
    }

	nGraphCases = 1;
    graphResultSet[0] = R;

    for (int i=0; i<nGraphs; i++) {
        if (!grph->isTimeseries(i) && !grph->isProfile(i) && !grph->isDistribution(i)) continue;
        mdiArea->addSubWindow(pGraph[i]);
		pGraph[i]->show();
    }

	if (show_outputdata) {
		mdiArea->addSubWindow(box_outputData);	// Need another way of creating this window - should be floating
		box_outputData->show();
	}
/*
    if (field->isConcPlot())
        field->makeConcPlot(mdiArea);
    if (field->isVolPlot())
        field->makeVolPlot(mdiArea);
    if (field->isOxyPlot())
        field->makeOxyPlot(mdiArea);
*/
    mdiArea->tileSubWindows();

	for (int i=0; i<nGraphs; i++) {
        if (!grph->isTimeseries(i)) continue;
        pGraph[i]->setAxisScale(QwtPlot::xBottom, 0, R->nsteps, 0);
	}
    Global::dist_nv = 20;
}


//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::drawGraphs()
{
	RESULT_SET *R;
	for (int kres=0; kres<Plot::ncmax; kres++) {
		R = graphResultSet[kres];
		if (R != 0) {
			for (int i=0; i<nGraphs; i++) {
                if (!grph->isTimeseries(i)) continue;
                if (!grph->isActive(i)) continue;
				int k = grph->get_dataIndex(i);
                QString tag = grph->get_tag(i);
                double yscale = grph->get_yscale(i);
                pGraph[i]->redraw(R->tnow, R->pData[i], R->nsteps, R->casename, tag, yscale, false);
//                pGraph[i]->redraw(R->tnow, R->pData[i], R->nsteps, R->casename, tag);
				if (k == 0) {
					grph->set_maxValue(i,R->maxValue[i]);
				} else {
					double maxval = grph->get_maxValue(i);
					double newmax = R->maxValue[i];
					if (newmax > maxval) {
						grph->set_maxValue(i,newmax);
					}
				}
			}
		}
	}
	for (int i=0; i<nGraphs; i++) {
        if (!grph->isTimeseries(i)) continue;
        if (!grph->isActive(i)) continue;
		double maxval = grph->get_maxValue(i);
		pGraph[i]->setYScale(maxval);
		pGraph[i]->replot();
	}
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::displayScene()
{
//    LOG_MSG("displayScene");
//	bool redo = false;	// need to understand this
	started = true;
    if (firstDisplay) {
        vtk->setCameraPosition();
        firstDisplay = false;
    }
    vtk->get_cell_positions();
    vtk->renderCells();
    if (videoVTK->record) {
        videoVTK->recorder();
    } else if (actionStop_recording_VTK->isEnabled()) {
        actionStart_recording_VTK->setEnabled(true);
        actionStop_recording_VTK->setEnabled(false);
    }
}

//--------------------------------------------------------------------------------------------------------
// Currently summaryData[] holds istep,ntot,nborn.  Hourly intervals, i.e. every 240 timesteps
//--------------------------------------------------------------------------------------------------------
void MainWindow::showSummary(int hr)
{
    double val;
    int n, res;
    QString tag;

    sprintf(msg,"showSummary: step: %d",step);
    LOG_MSG(msg);
    step++;
    if (step >= newR->nsteps) {
		LOG_MSG("ERROR: step >= nsteps");
        stopServer();
		return;
	}
    hour = hr;
    exthread->mutex1.lock();

//    hour = summaryData[0]*DELTA_T/(60*60);
//    hour = summaryData[1]*DELTA_T/60;

    progress = int(100.*hour/hours);
	progressBar->setValue(progress);
	QString hourstr = QString::number(int(hour));
	hour_display->setText(hourstr);

    Global::casename = newR->casename;
    newR->tnow[step] = step;

    // TS plots
	for (int i=0; i<nGraphs; i++) {
        if (!grph->isTimeseries(i)) continue;
        if (!grph->isActive(i)) continue;
		int k = grph->get_dataIndex(i);
        val = Global::summaryData[k];
        newR->pData[i][step] = val*grph->get_scaling(i);
        tag = grph->get_tag(i);
//        pGraph[i]->redraw(newR->tnow, newR->pData[i], step+1, casename, tag);
        double yscale = grph->get_yscale(i);
        pGraph[i]->redraw(newR->tnow, newR->pData[i], step+1, Global::casename, tag, yscale, false);
    }
    LOG_QMSG("did ts graphs");

    // Profile plots
    updateProfilePlots();
    LOG_QMSG("did profile graphs");

    field->setSliceChanged();
    if (step > 0 && !action_field->isEnabled()) {
        field->displayField(hour,&res);
        if (res != 0) {
            sprintf(msg,"displayField returned res: %d",res);
            LOG_MSG(msg);
            stopServer();
        }
    }
    exthread->mutex1.unlock();
    exthread->summary_done.wakeOne();
}

//--------------------------------------------------------------------------------------------------------
// nvars = 1 + Global::MAX_CHEMO + Global::N_EXTRA;
//--------------------------------------------------------------------------------------------------------
void MainWindow::updateProfilePlots()
{
    if (Global::casename == "") return;
    int ivar=0;
    for (int i=0; i<nGraphs; i++) {
        if (!grph->isActive(i)) continue;
        if (Global::conc_nc > 0 && grph->isProfile(i)) {
            double x[100], y[100];
            double xscale, yscale;
            QString tag = grph->get_tag(i);
            int k = grph->get_dataIndex(i);
            if (k == MULTI) {
                ivar = field->cell_constituent;
                QString title;
                field->getTitle(ivar,&title);
                pGraph[i]->setTitle(title);
                k = Global::GUI_to_DLL_index[ivar];
            }
            int n = Global::conc_nc;
            int offset = k*n;
            for (int j=0; j<n; j++) {
                x[j] = j*Global::conc_dx*1.0e4;
                y[j] = Global::concData[offset+j];
            }
            xscale = grph->get_xscale(x[n-1]);
            double maxval = 0;
            for (int j=0; j<n; j++) {
                if (y[j] > maxval) maxval = y[j];
            }
            yscale = pGraph[i]->calc_yscale(maxval);
            pGraph[i]->setAxisScale(QwtPlot::xBottom, 0, xscale, 0);
            pGraph[i]->setAxisScale(QwtPlot::yLeft, 0, yscale, 0);
//            if (k == CFSE){
//                pGraph[i]->setAxisScale(QwtPlot::xBottom, -20.0, 1.0, 0);
//            }
            pGraph[i]->setAxisTitle(QwtPlot::xBottom, "Distance (microns)");
            pGraph[i]->setAxisTitle(QwtPlot::yLeft, grph->get_yAxisTitle(i));
            pGraph[i]->redraw(x, y, n, Global::casename, tag, yscale, true);
        }
    }
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::outputData(QString qdata)
{
//	if (qdata.compare("VTK") == 0) {
	if (qdata.startsWith("VTK")) {
		qdata.replace(0,3,"");
//        bool savepos = false;
//		bool savepos = cbox_savepos->isChecked();
//		if (savepos) {
//			if (step < savepos_start) {
//				savepos = false;
//			}
//		}
//		vtk->read_cell_positions(cellfile, vtkfile, savepos);
		started = true;
        if (Global::showingVTK || firstVTK) {
			firstVTK = false;
//			bool redo = false;
//            if (showingVTK ) {
//				redo = true;
//			}
            vtk->renderCells();
		} 
	    posdata = true;
		if (qdata.length() == 0)
			return;
	}
//	if (qdata.contains("__EXIT__",Qt::CaseSensitive) || qdata.contains("Fortran") ) {
	if (quitMessage(qdata) || qdata.contains("Fortran") ) {
		return;
	}
	if (show_outputdata)
	    box_outputData->append(qdata);

    QStringList dataList = qdata.split(" ",QString::SkipEmptyParts);
	double data[11];
	for (int k=0; k<11; k++)
		data[k] = dataList[k].toDouble();
	step++;
	if (step >= newR->nsteps) {
		LOG_MSG("ERROR: step >= nsteps");
		return;
	}
	QString casename = newR->casename;
    newR->tnow[step] = step;		//data[1];
	for (int i=0; i<nGraphs; i++) {
        if (!grph->isTimeseries(i)) continue;
        if (!grph->isActive(i)) continue;
		int k = grph->get_dataIndex(i);
		newR->pData[i][step] = data[k]*grph->get_scaling(i);
	}

	for (int i=0; i<nGraphs; i++) {
        if (!grph->isTimeseries(i)) continue;
        if (!grph->isActive(i)) continue;
        QString tag = grph->get_tag(i);
        double yscale = grph->get_yscale(i);
        pGraph[i]->redraw(newR->tnow, newR->pData[i], step+1, casename, tag, yscale, true);
//        pGraph[i]->redraw(newR->tnow, newR->pData[i], step+1, casename, grph->get_tag(i));
    }
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::postConnection()
{
	LOG_MSG("postConnection");
	if (use_CPORT1) {
		sthread1->socket->close();
		sthread1->tcpServer->close();
		sthread1->quit();
		sthread1->wait(100);
		if (sthread1->isRunning()) {
			LOG_MSG("sthread1 did not terminate");
		}
	}

    action_run->setEnabled(true);
    action_pause->setEnabled(false);
    action_stop->setEnabled(false);
    action_save_3D_snapshot->setEnabled(true);
    action_save_profile_data->setEnabled(true);
    action_show_gradient3D->setEnabled(true);
    action_show_gradient2D->setEnabled(true);
    action_field->setEnabled(true);
//    tab_B->setEnabled(true);
//    tab_DC->setEnabled(true);
    tab_tumour->setEnabled(true);
    tab_chemo->setEnabled(true);
    tab_run->setEnabled(true);

	// Check if a result set of this name is already in the list, if so remove it
	for (int i=0; i<result_list.size(); i++) {
		if (newR->casename.compare(result_list[i]->casename) == 0) {
			result_list.removeAt(i);
		}
	}
	// Compute the maxima
	for (int i=0; i<nGraphs; i++) {
        if (!grph->isTimeseries(i)) continue;
        if (!grph->isActive(i)) continue;
		double maxval = getMaximum(newR,newR->pData[i]);
		newR->maxValue[i] = maxval;
	}

	// Add the new result set to the list
//	result_list.append(newR);
//	vtk->renderCells(true,true);		// for the case that the VTK page is viewed only after the execution is complete
    if (actionStop_recording_VTK->isEnabled()) {
        stopRecorderVTK();
    }
    if (actionStop_recording_FACS->isEnabled()) {
        stopRecorderFACS();
    }
    posdata = false;
	LOG_MSG("completed postConnection");
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::close_sockets()
{
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::pauseServer()
{
//	if (vtk->playing) {
		vtk->pause();
		LOG_MSG("Paused the player");
//	} else {
		exthread->pause();
		LOG_MSG("Paused the ABM program.");
//	}
	paused = true;
	action_run->setEnabled(true); 
	action_pause->setEnabled(false);
	action_stop->setEnabled(true);
    action_save_3D_snapshot->setEnabled(true);
    action_save_profile_data->setEnabled(true);
    action_show_gradient3D->setEnabled(true);
    action_show_gradient2D->setEnabled(true);
    action_field->setEnabled(true);
}

//--------------------------------------------------------------------------------------------------------
// Need to check that DLL has completed.
//--------------------------------------------------------------------------------------------------------
void MainWindow::stopServer()
{
//    goflag = true;
    exthread->summary_done.wakeOne();
    if (vtk->playing) {
		vtk->stop();
		LOG_MSG("Stopped the player");
	} else {
		LOG_MSG("MainWindow::stopServer: stop requested");
		if (paused) {
			LOG_MSG("was paused, runServer before stopping");
			runServer();
		}
		exthread->snapshot();
		exthread->stop();
		sleep(1);		// delay for Fortran to wrap up (does this help?)
//		if (use_CPORT1) {
//			sthread1->quit();
//			sthread1->terminate();
//		}
//		sthread0->stop();
		newR->nsteps = step+1;
	}
    action_run->setEnabled(true); 
    action_pause->setEnabled(false);
    action_stop->setEnabled(false);
    action_save_3D_snapshot->setEnabled(true);
    action_save_profile_data->setEnabled(true);
    action_show_gradient3D->setEnabled(true);
    action_show_gradient2D->setEnabled(true);
    action_field->setEnabled(true);
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::clearAllGraphs()
{
    if (nGraphCases > 0) {
//        if (nGraphs > 0) {
//            for (int i=0; i<nGraphs; i++) {
//                if (!grph->isTimeseries(i)) continue;
//                if (!grph->isActive(i)) continue;
//            }
//        }
        nGraphCases = 0;
	}
    for (int i=0; i<Plot::ncmax; i++) {
		graphCaseName[i] = "";
		graphResultSet[i] = 0;
	}
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
QString MainWindow::selectResultSet()
{
    QStringList items;
	for (int k=0; k<result_list.size(); k++) {
		RESULT_SET *R = result_list[k];
		if (R == 0) continue;
		bool inlist = false;
		for (int i=0; i<Plot::ncmax; i++) {
			if (graphResultSet[i] == 0) continue;
			if (R->casename.compare(graphResultSet[i]->casename) == 0) {
				inlist = true;
				break;
			}
		}
		if (!inlist)
			items << R->casename;
	}
	if (items.size() == 0) {
		QMessageBox::warning(this, tr("Select result case"),
			tr("No result sets available - use 'File > Load results'"));
		return QString("");
	}

    bool ok;
    QString item = QInputDialog::getItem(this, tr("Select result case"),
		tr("Case:"), items, 0, false, &ok);
    if (ok && !item.isEmpty())
		return item;
	else
		return QString("");
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::addGraph()
{
	// Need to select a result set from result_list then add the corresponding curves
        RESULT_SET *R = NULL;

	if (nGraphCases == Plot::ncmax) {
		QString mess = QString("The maximum number of cases is %1").arg(Plot::ncmax);
		QMessageBox::warning(this, tr("Add graph"),	tr((mess.toStdString()).data())); 
		return;
	}
	QString casename = selectResultSet();
	if (casename.compare("") == 0)
		return;

	for (int k=0; k<result_list.size(); k++) {
		if (casename.compare(result_list[k]->casename) == 0) {
			R = result_list[k];		// OK after doing a run or a load, followed by another load
			break;
		}
	}

	graphResultSet[nGraphCases] = R;
	nGraphCases++;
	// First add the curves
	for (int i=0; i<nGraphs; i++) {
        if (!grph->isTimeseries(i)) continue;
        if (!grph->isActive(i)) continue;
		pGraph[i]->addCurve(R->casename);
		pGraph[i]->setAxisAutoScale(QwtPlot::xBottom);
	}
	// Now redraw with the data
	drawGraphs();
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
int MainWindow::selectGraphCase()
{
    QStringList items;
	for (int i=0; i<Plot::ncmax; i++) {
		if (graphResultSet[i] == 0) continue;
		items << graphResultSet[i]->casename;
	}
	if (items.size() == 0) {
		QMessageBox::warning(this, tr("Select graph case"),
			tr("No graph cases to remove"));
		return -1;
	}

    bool ok;
    QString item = QInputDialog::getItem(this, tr("Select graph case"),
		tr("Case:"), items, 0, false, &ok);
	if (ok && !item.isEmpty()) {
		for (int i=0; i<Plot::ncmax; i++) {
			if (graphResultSet[i] == 0) continue;
			if (item.compare(graphResultSet[i]->casename) == 0) {
				return i;
			}
		}
	}
	return -1;
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::removeGraph()
{
	int i = selectGraphCase();
	if (i == -1) return;
	RESULT_SET *R = graphResultSet[i];
	// First remove the curves
	for (int i=0; i<nGraphs; i++) {
        if (!grph->isTimeseries(i)) continue;
        if (!grph->isActive(i)) continue;
		pGraph[i]->removeCurve(R->casename);
	}
	// Then remove the graph case
	graphResultSet[i] = 0;
	nGraphCases--;
	drawGraphs();
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::removeAllGraphs()
{
    LOG_MSG("removeAllGraphs");
    clearAllGraphs();
}

//---------------------------------------------------------------------
// Updates input parameter (QLineEdit) widgets according to the slider
//---------------------------------------------------------------------
void MainWindow::updateSliderBox()
{
    paramSaved = false;          // Keeps track of the fact that a param has been changed but not saved.
    // ---Get value from slider
    QString slider_str = sender()->objectName();
    QString stag = slider_str.mid(7);
    int ival = ((QSlider *)sender())->value();

    // --- Get param index from workingParameterList
	int k;
	for (k=0; k<nParams; k++) {
		PARAM_SET p = parm->get_param(k);
		if (stag.compare(p.tag) == 0)
			break;
	}
    int j = param_to_sliderIndex[k];
    SliderPlus *sp = sliderplus_list[j];
    double v = sp->int_to_val(ival);
    QString vstr = sp->val_to_str(v);
    ((QLineEdit *)sliderParam[j])->setText(vstr);
}

//------------------------------------------------------------------------------------------------------
// changeParam() is invoked in response to a signal sent when a value in a QLineEdit etc. widget
// is changed.  Note that when a QRadioButton widget is changed, signals are sent both the radiobuttons
// that change, but only one signal is used to change the parameter value.
//------------------------------------------------------------------------------------------------------
void MainWindow::changeParam()
{
	paramSaved = false;
    QObject *w = sender(); // Gets the pointer to the object that invoked the changeParam slot.
	if (w->isWidgetType()) {
		QString wname = w->objectName();
//        LOG_QMSG("changeParam:" + wname);
        if (wname.contains("_PARENT_") || wname.contains("_METAB1_") || wname.contains("_METAB2_")) {
            changeDrugParam(w);
            return;
        }
		if (wname.contains("line_")) {
			QString wtag = wname.mid(5);
			QLineEdit *lineEdit = (QLineEdit *)w;
            QString text = lineEdit->displayText();
			// Determine if there is a slider associated with the sender widget
			for (int k=0; k<parm->nParams; k++) {
				PARAM_SET p = parm->get_param(k);
				if (wtag.compare(p.tag) == 0) {
					int j = param_to_sliderIndex[k];
					if (j >= 0) {
                        QSlider *slider = slider_list[j];
                        SliderPlus *sp = sliderplus_list[j];
                        double v = sp->str_to_val(text);
                        int ival = sp->val_to_int(v);
                        int ival_old = sp->val_to_int(p.value);
						if (ival != ival_old) {
                            slider->setSliderPosition(ival);
						}
					}
					parm->set_value(k,text.toDouble());
					break;
				}
			}
		} else if (wname.contains("text_")) {
			QString wtag = wname.mid(5);
			QLineEdit *lineEdit = (QLineEdit *)w;
			QString text = lineEdit->displayText();
			for (int k=0; k<parm->nParams; k++) {
				PARAM_SET p = parm->get_param(k);
				if (wtag.compare(p.tag) == 0) {
					parm->set_label(k,text);
					break;
				}
			}
		} else if (wname.contains("spin_")) {
			QSpinBox *spinBox = (QSpinBox *)w;
            int v = spinBox->value();
			QString wtag = wname.mid(5);
			for (int k=0; k<parm->nParams; k++) {
				PARAM_SET p = parm->get_param(k);
				if (wtag.compare(p.tag) == 0) {
					parm->set_value(k,v);
                    break;
				}
				if (wname.contains("NCPU")) {
					ncpu = v;
				}
			}
		} else if (wname.contains("cbox_")) {
			QCheckBox *checkBox = (QCheckBox *)w;
			int v;

            bool use_OXYGEN = wname.contains("USE_OXYGEN");
            if (checkBox->isChecked()) {
                v = 1;
                if (use_OXYGEN)
                    enableUseOxygen();
            } else {
                v = 0;
                if (use_OXYGEN)
                    disableUseOxygen();
            }

            bool use_GLUCOSE = wname.contains("USE_GLUCOSE");
            if (checkBox->isChecked()) {
                v = 1;
                if (use_GLUCOSE)
                    enableUseGlucose();
            } else {
                v = 0;
                if (use_GLUCOSE)
                    disableUseGlucose();
            }

            bool use_TRACER = wname.contains("USE_TRACER");
            if (checkBox->isChecked()) {
                v = 1;
                if (use_TRACER)
                    enableUseTracer();
            } else {
                v = 0;
                if (use_TRACER)
                    disableUseTracer();
            }

//            bool use_TPZ = wname.contains("USE_TPZ_DRUG");
            if (checkBox->isChecked()) {
                v = 1;
            } else {
                v = 0;
            }

//            bool use_DNB = wname.contains("USE_DNB_DRUG");
            if (checkBox->isChecked()) {
                v = 1;
            } else {
                v = 0;
            }

//            bool use_TREATMENT_FILE = wname.contains("USE_TREATMENT_FILE");
//            if (checkBox->isChecked()) {
//                v = 1;
//                if (use_TREATMENT_FILE)
//                    enableUseTreatmentFile();
//            } else {
//                v = 0;
//                if (use_TREATMENT_FILE)
//                    disableUseTreatmentFile();
//            }

			QString wtag = wname.mid(5);
			for (int k=0; k<parm->nParams; k++) {
				PARAM_SET p = parm->get_param(k);
				if (wtag.compare(p.tag) == 0) {
					parm->set_value(k,v);
                    break;
				}
			}
			if (wname.contains("savepos")) {
				if (checkBox->isChecked()) {
					setSavePosStart();
				}
			}
		} else if (wname.contains("comb_")) {
			QComboBox *comboBox = (QComboBox *)w;
            int v = comboBox->currentIndex();
			QString wtag = wname.mid(5);
//			sprintf(msg,"combo: %s  currentIndex: %d",wtag,v);
//			LOG_MSG(msg);
			for (int k=0; k<parm->nParams; k++) {
				PARAM_SET p = parm->get_param(k);
				if (wtag.compare(p.tag) == 0) {
					parm->set_value(k,v+1);
                    break;
				}
			}
		} else if (wname.contains("rbut_")) {
			QRadioButton *radioButton = (QRadioButton *)w;
            QString wtag = wname.mid(5);
            for (int k=0; k<parm->nParams; k++) {
                PARAM_SET p = parm->get_param(k);
                if (wtag.compare(p.tag) == 0) {
                    if (radioButton->isChecked()) {
                        parm->set_value(k,1);
                    } else {
                        parm->set_value(k,0);
                    }
                    break;
                }
            }
		}
	}
}

/*
//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::killModelChanged()
{
    int model;
    QString name;
    QLineEdit *line;

    line = (QLineEdit *)sender();
    model = line->text().toInt();
    name = line->objectName();
    if (name.contains("TPZ")) {
        if (model > 3) {
            errorPopup("For TPZ-type drug only kill models 1, 2 and 3 are possible");
            line->setText("1");
        }
     } else if (name.contains("DNB")) {
        if (model < 4) {
            errorPopup("For DNB-type drug only kill models 4 and 5 are possible");
            line->setText("4");

        }
    }
}
*/

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::enableUseOxygen()
{
    for (int i=0; i<lineEdit_list.length(); i++) {
        QLineEdit *w = lineEdit_list[i];
        QString wname = w->objectName();
        if (wname.contains("line_OXYGEN")) {
            w->setEnabled(true);
        }
    }
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::disableUseOxygen()
{
    for (int i=0; i<lineEdit_list.length(); i++) {
        QLineEdit *w = lineEdit_list[i];
        QString wname = w->objectName();
        if (wname.contains("line_OXYGEN")) {
            w->setEnabled(false);
        }
    }
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::enableUseGlucose()
{
    for (int i=0; i<lineEdit_list.length(); i++) {
        QLineEdit *w = lineEdit_list[i];
        QString wname = w->objectName();
        if (wname.contains("line_GLUCOSE")) {
            w->setEnabled(true);
        }
    }
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::disableUseGlucose()
{
    for (int i=0; i<lineEdit_list.length(); i++) {
        QLineEdit *w = lineEdit_list[i];
        QString wname = w->objectName();
        if (wname.contains("line_GLUCOSE")) {
            w->setEnabled(false);
        }
    }
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::enableUseTracer()
{
    for (int i=0; i<lineEdit_list.length(); i++) {
        QLineEdit *w = lineEdit_list[i];
        QString wname = w->objectName();
        if (wname.contains("line_TRACER")) {
            w->setEnabled(true);
        }
    }
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::disableUseTracer()
{
    for (int i=0; i<lineEdit_list.length(); i++) {
        QLineEdit *w = lineEdit_list[i];
        QString wname = w->objectName();
        if (wname.contains("line_TRACER")) {
            w->setEnabled(false);
        }
    }
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
//void MainWindow::enableUseTreatmentFile()
//{
//    for (int i=0; i<lineEdit_list.length(); i++) {
//        QLineEdit *w = lineEdit_list[i];
//        QString wname = w->objectName();
//        if (wname.contains("text_TREATMENT_FILE")) {
//            w->setEnabled(true);
//        }
//    }
//    cbox_USE_TREATMENT_FILE->setChecked(true);
//}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
//void MainWindow::disableUseTreatmentFile()
//{
//    for (int i=0; i<lineEdit_list.length(); i++) {
//        QLineEdit *w = lineEdit_list[i];
//        QString wname = w->objectName();
//        if (wname.contains("text_TREATMENT_FILE")) {
//            w->setEnabled(false);
//        }
//    }
//    cbox_USE_TREATMENT_FILE->setChecked(false);
//}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow:: on_cbox_SAVE_PROFILE_DATA_toggled(bool checked)
{
    text_SAVE_PROFILE_DATA_FILE_NAME->setEnabled(checked);
    line_SAVE_PROFILE_DATA_INTERVAL->setEnabled(checked);
    line_SAVE_PROFILE_DATA_NUMBER->setEnabled(checked);
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::redrawDistPlot()
{
    int i_m = 0, i_s = 0;
    QString sname = sender()->objectName();
    for (int k=0; k<ndistplots; k++) {
		QwtPlot *qp = distplot_list[k];
        QString tag = qp->objectName().mid(8);
        QString tag_m = tag + "_MEDIAN";
        QString tag_s = tag + "_SHAPE";
		if (sname.endsWith(tag_m) || sname.endsWith(tag_s)) {   
			for (int i=0; i<nWidgets; i++) {
				QString wname = widget_list[i]->objectName();
                if (wname.endsWith(tag_m))
                    i_m = i;
                else if (wname.endsWith(tag_s))
                    i_s = i;
			}

            QString median_str = ((QLineEdit *)widget_list[i_m])->text();
            QString shape_str = ((QLineEdit *)widget_list[i_s])->text() ;
			if (median_str.compare("") == 0) return;
			if (shape_str.compare("") == 0) return;
            double median = median_str.toDouble();
            median = max(0.001, median);
            double shape = shape_str.toDouble();
            shape = max(1.001, shape);
			
			double *x = new double[nDistPts];
			double *prob = new double[nDistPts];
            create_lognorm_dist(median,shape,nDistPts,x,prob);
            int n = dist_limit(prob,nDistPts);
            double xmax = x[n];
            curve_list[k]->setData(x, prob, n);
            qp->setAxisScale(QwtPlot::xBottom, 0.0, xmax, 0.0);
            qp->replot();
			delete [] x;
			delete [] prob;
			x = 0;
			prob = 0;
		}
	}
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
int MainWindow::dist_limit(double *p, int n)
{
	int i, imax;
    double pmax = 0;

	imax = n-1;
	for (i=0; i<n; i++) {
		if (p[i] > pmax) {
			pmax = p[i];
			imax = i;
		}
	}
    double plim = 0.01*pmax;
	for (i=n-1; i>0; i--) {
		if (p[i] > plim) {
            return min(n-1,max(i,2*imax));
		}
	}
	return 1;
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
double MainWindow::erf(double z)
{
   double t = 1.0 / (1.0 + 0.5 * fabs(z));
   // use Horner's method
   double ans = 1 - t * exp( -z*z -  1.26551223 +
         t * ( 1.00002368 +
         t * ( 0.37409196 +
         t * ( 0.09678418 +
         t * (-0.18628806 +
         t * ( 0.27886807 +
         t * (-1.13520398 +
         t * ( 1.48851587 +
         t * (-0.82215223 +
         t * ( 0.17087277))))))))));
   if (z >= 0.0)
       return ans;
   else
       return -ans;
}

//-----------------------------------------------------------------------------------------
// When X is distributed N(mu,sig), this gives Prob{x1 < X <= x2}
//-----------------------------------------------------------------------------------------
double MainWindow::pnorm(double x1, double x2, double mu, double sig)
{
    double z1, z2, e1, e2;
		
	z1 = (x1-mu)/sig;
    e1 = erf(z1/sqrt(2.0))/2;
    z2 = (x2-mu)/sig;
    e2 = erf(z2/sqrt(2.0))/2;
    return e2 - e1;
}

//-----------------------------------------------------------------------------------------
// When log(X) is distributed N(mu,sig), this gives Prob{x1 < X <= x2}
//-----------------------------------------------------------------------------------------
double MainWindow::plognorm(double x1, double x2, double mu, double sig)
{
    double z1, z2, e1, e2;

    z1 = 0;
    z2 = 0;
    if (x1 == 0)
        e1 = -0.5;
	else {
        z1 = (log(x1)-mu)/sig;
        e1 = erf(z1/sqrt(2.0))/2;
	}
    if (x2 == 0)
        e2 = -0.5;
	else {
        z2 = (log(x2)-mu)/sig;
        e2 = erf(z2/sqrt(2.0))/2;
	}
    return e2 - e1;
}

//-----------------------------------------------------------------------------------------
// Create the lognormal distribution with median = p1, shape = p2
// at n points stored in x[], probability values stored in prob[].
// Note that x[0] = 0.
// The range of x is currently just less than 4*median.  This should be
// OK for values of shape < 2.
// Convert probability into probability density
//-----------------------------------------------------------------------------------------
void MainWindow::create_lognorm_dist(double p1, double p2,int n, double *x, double *prob)
{
	double xmax, dx, mu_l, sig_l, x1, x2;

    if (p1 >= 0.5)
        xmax = p1*4;
    else
        xmax = p1*8;
        
    dx = xmax/n;
    mu_l = log(p1);
    sig_l = log(p2);
	for (int ix=0; ix<n; ix++) {
        x1 = (ix - 0.5)*dx;
        x2 = x1 + dx;
		x1 = max(x1,0.0);
        x[ix] = (x1+x2)/2;
        prob[ix] = plognorm(x1,x2,mu_l,sig_l)/(x2-x1);
	}
}


//======================================================================================================
//------------------------------------------------------------------------------------------------------
// SliderPlus class member definitions
//------------------------------------------------------------------------------------------------------
SliderPlus::SliderPlus(QString aname, double valmin, double valmax, int nval, int iparam, int kwidget)
{
	int i;
    name = aname;
    pindex = iparam;
    windex = kwidget;
    dv = (valmax - valmin)/nval;
	for (i=10; i>-10; i--) {
		if (pow(10.0,i) < dv) {
            dv = pow(10.0,i);
            break;
		}
	}
    i = int(valmin/dv);
    vmin = dv*(i+1);
    int n1 = (int)((valmax - vmin)/dv + 0.5);	//round
	if (n1 > 5*nval) {
        dv = 5*dv;
        n1 = n1/5;
	}
	else if (n1 > 2*nval) {
        dv = 2*dv;
        n1 = n1/2;
	}
    i = int(valmin/dv);
    vmin = dv*i;
    n = n1;
    vmax = vmin + n*dv;
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
SliderPlus::~SliderPlus()
{}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
int SliderPlus::val_to_int(double v) {
    int i = (int)((v-vmin)/dv + 0.5);	//round
    return i;
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
double SliderPlus::int_to_val(int i) {
    double v = vmin + i*dv;
    return v;
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
QString SliderPlus::val_to_str(double v) {
	int i = SliderPlus::val_to_int(v);
	QString vstr = QString::number(int_to_val(i));
    return vstr;
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
double SliderPlus::str_to_val(QString vstr) {
    double v = vstr.toDouble();
    return v;
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
int SliderPlus::pIndex() {
    return pindex;
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
int SliderPlus::wIndex() {
    return windex;
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
int SliderPlus::nTicks() {
    return n;
}

void MainWindow::on_action_show_gradient3D_triggered()
{
    showGradient3D();
}

void MainWindow::on_action_show_gradient2D_triggered()
{
    showGradient2D();
}

//-----------------------------------------------------------------------------------------
// Shows how to fetch active constituent names and use them to initialise a radiobutton
// group, and also to record their DLL index values (ichemo):
// 0 CFSE
// 1 Oxygen
// 2 Glucose
// 3 Tracer
// 4 TPZ drug
// 5 TPZ drug metabolite 1
// 6 TPZ drug metabolite 2
// 7 DNB drug
// 8 DNB drug metabolite 1
// 9 DNB drug metabolite 2
// ...
//
// GUI_to_DLL_index[ivar], ivar=0,..,nvars_used-1 = base DLL index (0,MAX_CHEMO+N_EXTRA)
// DLL_to_GUI_index[ichemo],ichemo=0,..,MAX_CHEMO+N_EXTRA = index in list of variables in use (0,nvars_used-1)
//
// Note that while the cell constituents have DLL index values ranging from 0 to MAX_CHEMO+N_EXTRA,
// the field constituent DLL indexes are in the range 1 to MAX_CHEMO, since the extra
// variables not applicable in the extracellular compartment.
//-----------------------------------------------------------------------------------------
void MainWindow::setupConstituents()
{
    int nvarlen, narraylen;
    char *name_array;
    char name[25];
    QString str, tag;
    int ivar, ichemo;

    narraylen = 1000;
    name_array = (char *)malloc(narraylen*sizeof(char));
    get_constituents(&Global::nvars_used, Global::GUI_to_DLL_index, &nvarlen, name_array, &narraylen);
    for (ichemo=0; ichemo<32; ichemo++)
        Global::DLL_to_GUI_index[ichemo] = -1;
    for (ivar=0; ivar<Global::nvars_used; ivar++)
        Global::DLL_to_GUI_index[Global::GUI_to_DLL_index[ivar]] = ivar;
    int k = 0;
    for (ivar=0; ivar<Global::nvars_used; ivar++) {
        for (int i=0; i<nvarlen; i++) {
            name[i] = name_array[k];
            k++;
        }
        name[nvarlen] = NULL;
        str = name;
        Global::var_string[ivar] = str.trimmed();
        LOG_QMSG(name);
    }
    free(name_array);
    LOG_MSG("set up Global");
    tag = "cell";
    field->setCellConstituentButtons(groupBox_cell_constituent, field->buttonGroup_cell_constituent, &field->vbox_cell_constituent, &field->cell_constituent_rb_list, tag);
    LOG_MSG("did setCellCellConstituentButtons: field");
    tag = "field";
    field->setFieldConstituentButtons(groupBox_field_constituent, field->buttonGroup_field_constituent, &field->vbox_field_constituent, &field->field_constituent_rb_list, tag);
    LOG_MSG("did setCellCellConstituentButtons: field");
    tag = "histo";
    field->setCellConstituentButtons(groupBox_Histo_y_vars, buttonGroup_histo, &vbox_histo, &histo_rb_list, tag);
    LOG_MSG("did setCellConstituentButtons: histo");
    tag = "FACS_x";
    field->setCellConstituentButtons(groupBox_FACS_x_vars, buttonGroup_FACS_x_vars, &vbox_FACS_x_vars, &FACS_x_vars_rb_list, tag);
    LOG_MSG("did setCellConstituentButtons: FACS_x");
//    FACS_x_vars_rb_list[0]->setChecked(true);
    tag = "FACS_y";
    field->setCellConstituentButtons(groupBox_FACS_y_vars, buttonGroup_FACS_y_vars, &vbox_FACS_y_vars, &FACS_y_vars_rb_list, tag);
    LOG_MSG("did setCellConstituentButtons: FACS_y");
    field->setMaxConcentrations(groupBox_maxconc);
    LOG_MSG("did setMaxConcentrations");
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::setupCellColours()
{
    QComboBox *combo;
    int i, k;

    for (i=0; i<2; i++) {
        if (i == 0)
            combo = comboBox_CELLCOLOUR_1;
        else
            combo = comboBox_CELLCOLOUR_2;
        k = 0;
        combo->addItem("red");
        comboColour[k] = QColor(Qt::red);
        k++;
        combo->addItem("orange");
        comboColour[k] = QColor(255,130,0);
        k++;
        combo->addItem("yellow");
        comboColour[k] = QColor(Qt::yellow);
        k++;
        combo->addItem("green");
        comboColour[k] = QColor(Qt::green);
        k++;
        combo->addItem("blue");
        comboColour[k] = QColor(Qt::blue);
        k++;
        combo->addItem("magenta");
        comboColour[k] = QColor(Qt::magenta);
        k++;
        combo->addItem("cyan");
        comboColour[k] = QColor(Qt::cyan);
        k++;

//        combo->addItem("red");
//        combo->addItem("orange");
//        combo->addItem("yellow");
//        combo->addItem("green");
//        combo->addItem("blue");
//        combo->addItem("purple");
//        combo->addItem("brown");
    }
//    QString colstr = "red";
//    comboBox_CELLCOLOUR_1->addItem(colstr);
//    comboBox_CELLCOLOUR_1->addItem("orange");
//    comboBox_CELLCOLOUR_1->addItem("yellow");
//    comboBox_CELLCOLOUR_1->addItem("green");
//    comboBox_CELLCOLOUR_1->addItem("blue");
//    comboBox_CELLCOLOUR_1->addItem("purple");
//    comboBox_CELLCOLOUR_1->addItem("brown");
//    comboBox_CELLCOLOUR_2->addItem("red");
//    comboBox_CELLCOLOUR_2->addItem("orange");
//    comboBox_CELLCOLOUR_2->addItem("yellow");
//    comboBox_CELLCOLOUR_2->addItem("green");
//    comboBox_CELLCOLOUR_2->addItem("blue");
//    comboBox_CELLCOLOUR_2->addItem("purple");
//    comboBox_CELLCOLOUR_2->addItem("brown");
    int k1 = 0;
    comboBox_CELLCOLOUR_1->setCurrentIndex(k1);
    vtk->celltype_colour[1] = comboColour[k1];
    int k2 = 1;
    comboBox_CELLCOLOUR_2->setCurrentIndex(k2);
    vtk->celltype_colour[2] = comboColour[k2];
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::setupGraphSelector()
{
    QGridLayout *grid = new QGridLayout;
    int row[3];
    row[0] = row[1] = row[2] = -1;

/* // superceded by profile plots
    checkBox_conc = new QMyCheckBox();
    checkBox_conc->setText("Concentration Profile");
    checkBox_conc->setChecked(field->isConcPlot());
    checkBox_conc->setObjectName("checkBox_conc");
    checkBox_conc->description = "Concentration along a line through the centre of the blob";
    row[0]++;
    grid->addWidget(checkBox_conc,row[0],0);
    connect((QObject *)checkBox_conc, SIGNAL(checkBoxClicked(QString)), this, SLOT(showMore(QString)));

    checkBox_vol = new QMyCheckBox();
    checkBox_vol->setText("Cell Volume Distribution");
    checkBox_vol->setChecked(field->isVolPlot());
    checkBox_vol->setObjectName("checkBox_vol");
    checkBox_vol->description = "Probability distribution (histogram) of cell volume";
    row[0]++;
    grid->addWidget(checkBox_vol,row[0],0);
    connect((QObject *)checkBox_vol, SIGNAL(checkBoxClicked(QString)), this, SLOT(showMore(QString)));

    checkBox_oxy = new QMyCheckBox();
    checkBox_oxy->setText("Cell Oxygen Distribution");
    checkBox_oxy->setChecked(field->isOxyPlot());
    checkBox_oxy->setObjectName("checkBox_oxy");
    checkBox_oxy->description = "Probability distribution (histogram) of cell oxygen concentration";
    row[0]++;
    grid->addWidget(checkBox_oxy,row[0],0);
    connect((QObject *)checkBox_oxy, SIGNAL(checkBoxClicked(QString)), this, SLOT(showMore(QString)));
*/
    cbox_ts = new QMyCheckBox*[grph->n_tsGraphs];
    for (int i=0; i<grph->n_tsGraphs; i++) {
        int col = grph->tsGraphs[i].type;
        row[col]++;
        QString text = grph->tsGraphs[i].title;
        cbox_ts[i] = new QMyCheckBox;
        cbox_ts[i]->setText(text);
//        cbox_ts[i]->setChecked(grph->tsGraphs[i].active);
        cbox_ts[i]->setObjectName("cbox_"+grph->tsGraphs[i].tag);
        grid->addWidget(cbox_ts[i],row[col],col);
        connect((QObject *)cbox_ts[i], SIGNAL(checkBoxClicked(QString)), this, SLOT(showMore(QString)));
//        if (col == 2) {
//            LOG_QMSG(grph->tsGraphs[i].tag);
//            sprintf(msg,"graph: %d %d",i,grph->tsGraphs[i].active);
//            LOG_MSG(msg);
//        }
    }
    groupBox_graphselect->setLayout(grid);

    QRect rect = groupBox_graphselect->geometry();
#ifdef __DISPLAY768
    rect.setHeight(460);
#else
    rect.setHeight(500);
#endif
    groupBox_graphselect->setGeometry(rect);
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::processGroupBoxClick(QString text)
{
    LOG_QMSG("processGroupBoxClick: " + text);
    QwtPlot *plot;

    if (text.compare("Histo") == 0) {
        LOG_MSG("save Histo plot");
        bool use_HistoBar = radioButton_histotype_1->isChecked();
        if (use_HistoBar) {
            plot = qpHistoBar;
            qpHistoLine->hide();
        } else {
            plot = qpHistoLine;
            qpHistoBar->hide();
        }
    } else if (text.compare("FACS") == 0) {
        LOG_MSG("save FACS plot");
        plot = qpFACS;
    } else {
        return;
    }

    int w = plot->width();
    int h = plot->height();
    QPixmap pixmap(w, h);
    pixmap.fill(Qt::white); // Qt::transparent ?

    QwtPlotPrintFilter filter;
    int options = QwtPlotPrintFilter::PrintAll;
    options &= ~QwtPlotPrintFilter::PrintBackground;
    options |= QwtPlotPrintFilter::PrintFrameWithScales;
    filter.setOptions(options);

    plot->print(pixmap, filter);

//		QString fileName = getImageFile();
    QString fileName = QFileDialog::getSaveFileName(0,"Select image file", ".",
        "Image files (*.png *.jpg *.tif *.bmp)");
    if (fileName.isEmpty()) {
        return;
    }
    pixmap.save(fileName,0,-1);
}

//--------------------------------------------------------------------------------------------------------
// Note that the initial selection of active graphs is now set in params.cpp
//--------------------------------------------------------------------------------------------------------
void MainWindow::setGraphsActive()
{
    for (int i=0; i<grph->n_tsGraphs; i++) {
        grph->tsGraphs[i].active = cbox_ts[i]->isChecked();
    }
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::showBool(QString qstr, bool result)
{
    if (result) {
        LOG_QMSG(qstr + "= true");
    } else {
        LOG_QMSG(qstr + "= false");
    }
}

/*
//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::on_comb_TPZ_currentIndexChanged(int index)
{
    text_TPZ_DRUG_NAME->setText(comb_TPZ->currentText());
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::on_comb_DNB_currentIndexChanged(int index)
{
    text_DNB_DRUG_NAME->setText(comb_DNB->currentText());
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::on_cbox_USE_TPZ_DRUG_toggled(bool checked)
{
    LOG_MSG("cbox_use_TPZ_drug toggled");
    QLineEdit *leb = findChild<QLineEdit *>("line_TPZ_DRUG_BDRY_CONC");
//    QCheckBox *cbd = findChild<QCheckBox *>("cbox_TPZ_DRUG_DECAY");
    QCheckBox *cbm = findChild<QCheckBox *>("cbox_TPZ_DRUG_SIMULATE_METABOLITE");
    leb->setEnabled(checked);
    cbm->setEnabled(checked);
    setTreatmentFileUsage();
    comb_TPZ->setEnabled(checked);
    text_TPZ_DRUG_NAME->setEnabled(checked);
    text_TPZ_DRUG_NAME->setText(comb_TPZ->currentText());
//    int indexTPZ = comb_TPZ->currentIndex();
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::on_cbox_USE_DNB_DRUG_toggled(bool checked)
{
//    LOG_MSG("cbox_use_DNB_drug toggled");
    QLineEdit *leb = findChild<QLineEdit *>("line_DNB_DRUG_BDRY_CONC");
//    QCheckBox *cbd = findChild<QCheckBox *>("cbox_DNB_DRUG_DECAY");
    QCheckBox *cbm = findChild<QCheckBox *>("cbox_DNB_DRUG_SIMULATE_METABOLITE");
    leb->setEnabled(checked);
    cbm->setEnabled(checked);
    setTreatmentFileUsage();
    comb_DNB->setEnabled(checked);
    text_DNB_DRUG_NAME->setEnabled(checked);
    text_DNB_DRUG_NAME->setText(comb_DNB->currentText());
//    int indexDNB = comb_DNB->currentIndex();
}
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void MainWindow::on_checkBox_CELLDISPLAY_1_toggled(bool display)
{
    vtk->display_celltype[1] = display;
//    vtk->cleanup();
    vtk->renderCells(false,false);
//    LOG_QMSG("toggled display_celltype[1]");
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void MainWindow::on_checkBox_CELLDISPLAY_2_toggled(bool display)
{
    vtk->display_celltype[2] = display;
    vtk->renderCells(false,false);
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void MainWindow::on_comboBox_CELLCOLOUR_1_currentIndexChanged(int index)
{
    QColor qcolor;
//    vtk->celltype_colour[1] = comboBox_CELLCOLOUR_1->currentText();
    qcolor = comboColour[index];
    vtk->celltype_colour[1] = qcolor;
    vtk->renderCells(false,false);
    sprintf(msg,"changed celltype_colour[1]: index: %d r,g,b: %d %d %d",index,qcolor.red(),qcolor.green(),qcolor.blue());
    LOG_MSG(msg);
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void MainWindow::on_comboBox_CELLCOLOUR_2_currentIndexChanged(int index)
{
//    vtk->celltype_colour[2] = comboBox_CELLCOLOUR_2->currentText();
    vtk->celltype_colour[2] = comboColour[index];
    vtk->renderCells(false,false);
}

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
void MainWindow::on_checkBox_FACS_PLOT_toggled(bool checked)
{
    line_FACS_INTERVAL->setEnabled(checked);
    if (!checked) {
        line_FACS_INTERVAL->setText("0");
    }
}

//-------------------------------------------------------------
// Switches to the FACS screen
//-------------------------------------------------------------
void MainWindow::on_action_FACS_triggered()
{
    stackedWidget->setCurrentIndex(4);
    action_outputs->setEnabled(true);
    action_inputs->setEnabled(true);
    action_VTK->setEnabled(true);
    action_FACS->setEnabled(false);
    Global::showingFACS = true;
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::on_cbox_USE_RADIATION_toggled(bool checked)
{
    setTreatmentFileUsage();
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::setTreatmentFileUsage()
{
    if (cbox_USE_TPZ_DRUG->isChecked() || cbox_USE_DNB_DRUG->isChecked() || cbox_USE_RADIATION->isChecked()) {
        enableUseTreatmentFile();
    } else {
        disableUseTreatmentFile();
    }
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::on_line_CELLPERCENT_1_textEdited(QString pc1_str)
{
    double pc1 = pc1_str.toDouble();
    double pc2 = 100 - pc1;
    QString pc2_str = QString::number(pc2);
    line_CELLPERCENT_2->setText(pc2_str);
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::on_line_CELLPERCENT_2_textEdited(QString pc2_str)
{
    double pc2 = pc2_str.toDouble();
    double pc1 = 100 - pc2;
    QString pc1_str = QString::number(pc1);
    line_CELLPERCENT_1->setText(pc1_str);
}


//------------------------------------------------------------------------------------------------------
// This should be used for any radioButtonGroups for model input parameters
//------------------------------------------------------------------------------------------------------
void MainWindow::radioButtonChanged(QAbstractButton *b)
{
    QString wtag = b->objectName();
    int rbutton_case;
    if (b->isChecked()) {
        QString ptag = parse_rbutton(wtag,&rbutton_case);
        // Now need to reflect the change in the workingParameterList
        // Need to locate ptag
        for (int k=0; k<nParams; k++) {
            PARAM_SET p = parm->get_param(k);
            if (ptag.compare(p.tag) == 0) {
                parm->set_value(k,double(rbutton_case));
                break;
            }
        }
    }
}

void MainWindow::buttonClick_constituent(QAbstractButton* button)
{
    LOG_MSG("buttonClick_constituent");
    field->setConstituent(button);
}

void MainWindow::buttonClick_plane(QAbstractButton* button)
{
    LOG_MSG("buttonClick_plane");
    field->setPlane(button);
}

void MainWindow::buttonClick_canvas(QAbstractButton* button)
{
    LOG_MSG("buttonClick_canvas");
}

void MainWindow::textChanged_fraction(QString text)
{
    LOG_MSG("textChanged_fraction");
    field->setFraction(text);
}

void MainWindow::textEdited_fraction(QString text)
{
    LOG_MSG("textEdited_fraction");
    field->setFraction(text);
}

void MainWindow::onSelectConstituent()
{
    if (exthread != NULL)
        field->selectConstituent();
}

void MainWindow::on_verticalSliderTransparency_sliderMoved(int position)
{
    vtk->setOpacity(position);
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::on_cbox_TPZ_DRUG_SIMULATE_METABOLITE_toggled(bool checked)
{
    LOG_MSG("cbox_use_drugA_metabolite toggled");
    QRadioButton *rbm = findChild<QRadioButton*>("radioButton_drugA_metabolite");
    QCheckBox *cbm = findChild<QCheckBox *>("cbox_DRUG_A_METABOLITE_DECAY");
    rbm->setEnabled(checked);
    cbm->setEnabled(checked);
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::on_cbox_DNB_DRUG_SIMULATE_METABOLITE_toggled(bool checked)
{
    LOG_MSG("cbox_use_drugB_metabolite toggled");
    QRadioButton *rbm = findChild<QRadioButton*>("radioButton_drugB_metabolite");
    QCheckBox *cbm = findChild<QCheckBox *>("cbox_DRUG_B_METABOLITE_DECAY");
    if (checked) {
        rbm->setEnabled(true);
        cbm->setEnabled(true);
    } else {
        rbm->setEnabled(false);
        cbm->setEnabled(false);
    }
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::initializeGraphs(RESULT_SET *R)
{
    mdiArea->closeAllSubWindows();
    mdiArea->show();
    setGraphsActive();
    grph->makeGraphList();
    LOG_QMSG("did makeGraphList");
    nGraphs = grph->nGraphs;
    if (nGraphCases > 0) {
        clearAllGraphs();
        LOG_QMSG("did clearAllGraphs");
    }

    QString tag;
    QString title;
    QString yAxisTitle;
    int k = 0;
    for (int i=0; i<nGraphs; i++) {
        if (grph->isActive(i)) {
            tag = grph->get_tag(i);
            title = grph->get_title(i);
            yAxisTitle = grph->get_yAxisTitle(i);
            k++;
        } else {
            tag = "";
            title = "";
            yAxisTitle = "";
        }
        if (k > maxGraphs) break;
            if (pGraph[i] != NULL) {
                pGraph[i]->deleteLater();
                pGraph[i] = NULL;
            }
            pGraph[i] = new Plot(tag,R->casename);
            pGraph[i]->setTitle(title);
            pGraph[i]->setAxisTitle(QwtPlot::yLeft, yAxisTitle);
    }
    LOG_QMSG("did setTitles");

    for (int i=0; i<nGraphs; i++) {
//        if (grph->isTimeseries(i)) {
            LOG_QMSG("addSubWindow: " + grph->get_tag(i));
            mdiArea->addSubWindow(pGraph[i]);
            pGraph[i]->show();
//        }
    }

    if (show_outputdata) {
        mdiArea->addSubWindow(box_outputData);	// Need another way of creating this window - should be floating
        box_outputData->show();
    }

    mdiArea->tileSubWindows();

    for (int i=0; i<nGraphs; i++) {
        if (!grph->isActive(i)) continue;
//        sprintf(msg,"i: %d isActive",i);
//        LOG_MSG(msg);
        if (grph->isTimeseries(i)) {
            pGraph[i]->setAxisScale(QwtPlot::xBottom, 0, R->nsteps, 0);
            pGraph[i]->setAxisTitle(QwtPlot::xBottom, "Time (hours)");
        }
    }
//    showmdiAreaSize();
}
*/

/*
//------------------------------------------------------------------------------------------------------
// Need to convert Kmet0 from /min to /sec, and kill_duration from mins to secs.
//------------------------------------------------------------------------------------------------------
void MainWindow::on_pushButton_SN30K_Kd_1_clicked()
{
//    kmet = (SN30K%C1(i) + SN30K%C2(i)*SN30K%KO2(i)/(SN30K%KO2(i) + SN30K%kill_O2(i)))*SN30K%Kmet0(i)
//	if (SN30K%kill_model(i) == 1) then
//		SN30K%Kd(i) = -log(1-SN30K%kill_fraction(i))/(SN30K%kill_duration(i)*kmet*SN30K%kill_drug(i))
//	elseif (SN30K%kill_model(i) == 2) then
//		SN30K%Kd(i) = -log(1-SN30K%kill_fraction(i))/(SN30K%kill_duration(i)*kmet*SN30K%kill_drug(i)**2)
//	elseif (SN30K%kill_model(i) == 3) then
//		SN30K%Kd(i) = -log(1-SN30K%kill_fraction(i))/(SN30K%kill_duration(i)*(kmet*SN30K%kill_drug(i))**2)
//	endif
    double Kd;
    double C1 = line_SN30K_C1_1->text().toDouble();
    double C2 = line_SN30K_C2_1->text().toDouble();
    double KO2 = line_SN30K_KO2_1->text().toDouble();
    double Kmet0 = line_SN30K_KMET0_1->text().toDouble();
    double kill_O2 = line_SN30K_KILL_O2_CONC_1->text().toDouble();
    double kill_drug = line_SN30K_KILL_DRUG_CONC_1->text().toDouble();
    double kill_duration = line_SN30K_KILL_DURATION_1->text().toDouble();
    double kill_fraction = line_SN30K_KILL_FRACTION_1->text().toDouble();
    Kmet0 /= 60;            // /min -> /sec
    kill_duration *= 60;    // min -> sec
    KO2 *= 1.0e-3;          // um -> mM
    double kmet = (C1 + C2*KO2/(KO2 + kill_O2))*Kmet0;
    if (rbut_SN30K_KILL_MODEL_1_0->isChecked()) {
        Kd = -log(1-kill_fraction)/(kill_duration*kmet*kill_drug);
    } else if (rbut_SN30K_KILL_MODEL_1_1->isChecked()) {
        Kd = -log(1-kill_fraction)/(kill_duration*kmet*qPow(kill_drug,2));
    } else if (rbut_SN30K_KILL_MODEL_1_2->isChecked()) {
        Kd = -log(1-kill_fraction)/(kill_duration*qPow(kmet*kill_drug,2));
    }
    line_SN30K_KD_1->setText(QString::number(Kd,'g',4));
}

//------------------------------------------------------------------------------------------------------
// Need to convert Kmet0 from /min to /sec, and kill_duration from mins to secs.
//------------------------------------------------------------------------------------------------------
void MainWindow::on_pushButton_SN30K_Kd_2_clicked()
{
    double Kd;
    double C1 = line_SN30K_C1_2->text().toDouble();
    double C2 = line_SN30K_C2_2->text().toDouble();
    double KO2 = line_SN30K_KO2_2->text().toDouble();
    double Kmet0 = line_SN30K_KMET0_2->text().toDouble();
    double kill_O2 = line_SN30K_KILL_O2_CONC_2->text().toDouble();
    double kill_drug = line_SN30K_KILL_DRUG_CONC_2->text().toDouble();
    double kill_duration = line_SN30K_KILL_DURATION_2->text().toDouble();
    double kill_fraction = line_SN30K_KILL_FRACTION_2->text().toDouble();
    Kmet0 /= 60;            // /min -> /sec
    kill_duration *= 60;    // min -> sec
    KO2 *= 1.0e-3;          // um -> mM
    double kmet = (C1 + C2*KO2/(KO2 + kill_O2))*Kmet0;
    if (rbut_SN30K_KILL_MODEL_2_0->isChecked()) {
        Kd = -log(1-kill_fraction)/(kill_duration*kmet*kill_drug);
    } else if (rbut_SN30K_KILL_MODEL_2_1->isChecked()) {
        Kd = -log(1-kill_fraction)/(kill_duration*kmet*qPow(kill_drug,2));
    } else if (rbut_SN30K_KILL_MODEL_2_2->isChecked()) {
        Kd = -log(1-kill_fraction)/(kill_duration*qPow(kmet*kill_drug,2));
    }
    line_SN30K_KD_2->setText(QString::number(Kd,'g',4));
}
*/

/*
//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::on_cbox_drugA_decay_toggled(bool checked)
{
    LOG_MSG("cbox_drugA_decay toggled");
    QLineEdit *leh = findChild<QLineEdit *>("line_DRUG_A_HALFLIFE");
    if (checked) {
        leh->setEnabled(true);
    } else {
        leh->setEnabled(false);
    }
}

/*
//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::on_cbox_drugA_metabolite_decay_toggled(bool checked)
{
    LOG_MSG("cbox_drugA_metabolite_decay toggled");
    QLineEdit *leh = findChild<QLineEdit *>("line_DRUG_A_METABOLITE_HALFLIFE");
    if (checked) {
        leh->setEnabled(true);
    } else {
        leh->setEnabled(false);
    }
}
*/

/*
//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::on_cbox_drugB_decay_toggled(bool checked)
{
    LOG_MSG("cbox_drugB_decay toggled");
    QLineEdit *leh = findChild<QLineEdit *>("line_DRUG_B_HALFLIFE");
    if (checked) {
        leh->setEnabled(true);
    } else {
        leh->setEnabled(false);
    }
}
*/

/*
//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::on_cbox_drugB_metabolite_decay_toggled(bool checked)
{
    LOG_MSG("cbox_drugB_metabolite_decay toggled");
    QLineEdit *leh = findChild<QLineEdit *>("line_DRUG_B_METABOLITE_HALFLIFE");
    if (checked) {
        leh->setEnabled(true);
    } else {
        leh->setEnabled(false);
    }
}
*/

//==================================================================================================================
// Code below here is not used
//----------------------------
/*
void MainWindow::on_cbox_SHOW_NONCOGNATE_toggled(bool checked)
{
    QLineEdit *le = findChild<QLineEdit*>("line_DISPLAY_FRACTION");
    if (checked) {
        le->setEnabled(true);
    } else {
        le->setEnabled(false);
    }
}
*/
void MainWindow::closeEvent(QCloseEvent *event)
{
    LOG_MSG("Close event");
    event->accept();

//    if (maybeSave()) {
//        writeSettings();
//        event->accept();
//    } else {
//        event->ignore();
//    }
}

void MainWindow::newFile()
{
//    if (maybeSave()) {
//        textEdit->clear();
//        setCurrentFile("");
//    }
}

void MainWindow::open()
{
//    if (maybeSave()) {
//        QString fileName = QFileDialog::getOpenFileName(this);
//        if (!fileName.isEmpty())
//            loadFile(fileName);
//    }
}

void MainWindow::about()
{
//   QMessageBox::about(this, tr("About Application"),
//            tr("The <b>Application</b> example demonstrates how to "
//               "write modern GUI applications using Qt, with a menu bar, "
//               "toolbars, and a status bar."));
}

void MainWindow::documentWasModified()
{
//    setWindowModified(textEdit->document()->isModified());
}

/*
void MainWindow::createMenus()
{
    fileMenu = menuBar()->addMenu(tr("&File"));
    fileMenu->addAction(newAct);
    fileMenu->addAction(openAct);
    fileMenu->addAction(saveAct);
    fileMenu->addAction(saveAsAct);
    fileMenu->addSeparator();
    fileMenu->addAction(exitAct);

    editMenu = menuBar()->addMenu(tr("&Edit"));
    editMenu->addAction(cutAct);
    editMenu->addAction(copyAct);
    editMenu->addAction(pasteAct);

    menuBar()->addSeparator();

    helpMenu = menuBar()->addMenu(tr("&Help"));
    helpMenu->addAction(aboutAct);
    helpMenu->addAction(aboutQtAct);
}

void MainWindow::createToolBars()
{
    fileToolBar = addToolBar(tr("File"));
    fileToolBar->addAction(newAct);
    fileToolBar->addAction(openAct);
    fileToolBar->addAction(saveAct);

    editToolBar = addToolBar(tr("Edit"));
    editToolBar->addAction(cutAct);
    editToolBar->addAction(copyAct);
    editToolBar->addAction(pasteAct);
}

void MainWindow::createStatusBar()
{
    statusBar()->showMessage(tr("Ready"));
}

void MainWindow::readSettings()
{
    QSettings settings("Trolltech", "Application Example");
    QPoint pos = settings.value("pos", QPoint(200, 200)).toPoint();
    QSize size = settings.value("size", QSize(400, 400)).toSize();
    resize(size);
    move(pos);
}

void MainWindow::writeSettings()
{
    QSettings settings("Trolltech", "Application Example");
    settings.setValue("pos", pos());
    settings.setValue("size", size());
}

bool MainWindow::maybeSave()
{
    if (textEdit->document()->isModified()) {
        QMessageBox::StandardButton ret;
        ret = QMessageBox::warning(this, tr("Application"),
                     tr("The document has been modified.\n"
                        "Do you want to save your changes?"),
                     QMessageBox::Save | QMessageBox::Discard | QMessageBox::Cancel);
        if (ret == QMessageBox::Save)
            return save();
        else if (ret == QMessageBox::Cancel)
            return false;
    }
    return true;
}

void MainWindow::loadFile(const QString &fileName)
{
    QFile file(fileName);
    if (!file.open(QFile::ReadOnly | QFile::Text)) {
        QMessageBox::warning(this, tr("Application"),
                             tr("Cannot read file %1:\n%2.")
                             .arg(fileName)
                             .arg(file.errorString()));
        return;
    }

    QTextStream in(&file);
//#ifndef QT_NO_CURSOR
    QApplication::setOverrideCursor(Qt::WaitCursor);
//#endif
    textEdit->setPlainText(in.readAll());
//#ifndef QT_NO_CURSOR
    QApplication::restoreOverrideCursor();
//#endif

    setCurrentFile(fileName);
    statusBar()->showMessage(tr("File loaded"), 2000);
}

bool MainWindow::saveFile(const QString &fileName)
{
    QFile file(fileName);
    if (!file.open(QFile::WriteOnly | QFile::Text)) {
        QMessageBox::warning(this, tr("Application"),
                             tr("Cannot write file %1:\n%2.")
                             .arg(fileName)
                             .arg(file.errorString()));
        return false;
    }

    QTextStream out(&file);
//#ifndef QT_NO_CURSOR
    QApplication::setOverrideCursor(Qt::WaitCursor);
//#endif
	QString text = "This is a test string";
	text += "\n";
    out << text;
	text = "This is another test string";
	text += "\n";
    out << text;
//#ifndef QT_NO_CURSOR
    QApplication::restoreOverrideCursor();
//#endif

    return true;
}

void MainWindow::setCurrentFile(const QString &fileName)
{
    curFile = fileName;
    textEdit->document()->setModified(false);
    setWindowModified(false);

    QString shownName = curFile;
    if (curFile.isEmpty())
        shownName = "untitled.txt";
    setWindowFilePath(shownName);
}

QString MainWindow::strippedName(const QString &fullFileName)
{
    return QFileInfo(fullFileName).fileName();
}

void MainWindow::on_radioButton_oxygen_clicked()
{

}

void MainWindow::on_radioButton_glucose_clicked(bool checked)
{

}
*/

