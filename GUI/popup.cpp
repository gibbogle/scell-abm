#include "mainwindow.h"
#include "QMessageBox"
#include "QFile"
#include <QDebug>

#include "plotwin.h"

#define NPLOT 100

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::setupPopup()
{
    connect(pushButton_radSF_1,SIGNAL(clicked()),this,SLOT(pushButton_clicked()));
    connect(pushButton_radSF_2,SIGNAL(clicked()),this,SLOT(pushButton_clicked()));
    connect(pushButton_glucose_0,SIGNAL(clicked()),this,SLOT(pushButton_clicked()));
    connect(pushButton_drugKF_0,SIGNAL(clicked()),this,SLOT(pushButton_clicked()));
    connect(pushButton_drugSF_0,SIGNAL(clicked()),this,SLOT(pushButton_clicked()));
    connect(pushButton_complementarySF_0,SIGNAL(clicked()),this,SLOT(pushButton_clicked()));
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::pushButton_clicked()
{

    int HMI_SCALE = 1;
    int cellType;
    QString title, plotType, plotName;

    QObject *senderObj = sender(); // This will give Sender object
    QString senderObjName = senderObj->objectName();
    plotwin = new PlotWin(this);
    QWidget *cw = plotwin->centralWidget();
    QFrame *plotFrame = cw->findChild<QFrame *>("plotFrame");
    QCustomPlot *popup_plot = new QCustomPlot(plotFrame);
    popup_plot->setObjectName("popup_plot");
    popup_plot->setGeometry(QRect(5*HMI_SCALE, 5*HMI_SCALE, 660*HMI_SCALE, 380*HMI_SCALE));

    // generate data:
    // Extract plot type and cell type from senderObjName
    QStringList list = senderObjName.split("_");
    if (list.size() < 3) return;
    if (list.size() == 3) {
        plotType = list[1];
        cellType = list[2].toInt();
        if (plotType == "radSF") {
            plotName = "Survival Fraction cell type " + list[2];
            plotwin->setWindowTitle(plotName);
            double C_O2;
            double maxdose = 50;
            C_O2 = 0;
            QVector<double> x0(NPLOT), y0(NPLOT); // initialize with entries 0..100
            makeSFPlot(list[2], C_O2, maxdose, &x0, &y0);
            // create graph and assign data to it:
            popup_plot->addGraph();
            popup_plot->graph(0)->setData(x0, y0);
            popup_plot->graph(0)->setPen(QPen(Qt::red));
            popup_plot->graph(0)->setName("O2 = 0%");

            C_O2 = 0.18;
            QVector<double> x1(NPLOT), y1(NPLOT); // initialize with entries 0..100
            makeSFPlot(list[2], C_O2, maxdose, &x1, &y1);
            // create graph and assign data to it:
            popup_plot->addGraph();
            popup_plot->graph(1)->setData(x1, y1);
            popup_plot->graph(1)->setPen(QPen(Qt::blue));
            popup_plot->graph(1)->setName("O2 = 20%");

            // give the axes some labels:
            popup_plot->xAxis->setLabel("Dose (Gy)");
            popup_plot->yAxis->setLabel("SF");
            // set axes ranges
            popup_plot->xAxis->setAutoTickStep(false);
            popup_plot->xAxis->setTickStep(5);
            popup_plot->xAxis->setRange(0, maxdose);
            popup_plot->yAxis->setRange(1.e-5, 1);
            popup_plot->yAxis->setScaleType(QCPAxis::stLogarithmic);
            popup_plot->yAxis->setScaleLogBase(10);
            popup_plot->yAxis->setNumberFormat("eb"); // e = exponential, b = beautiful decimal powers
            popup_plot->yAxis->setNumberPrecision(0); // makes sure "1*10^4" is displayed only as "10^4"
            popup_plot->legend->setVisible(true);
        } else if (plotType == "glucose") {
            plotName = "Medium Glucose Depletion";
            plotwin->setWindowTitle(plotName);
            double ndays;
            double C_O2;
            C_O2 = 0.18;
            QVector<double> x0(NPLOT), y0(NPLOT); // initialize with entries 0..100
            makeGlucosePlot(&ndays, &x0, &y0);
            // create graph and assign data to it:
            popup_plot->addGraph();
            popup_plot->graph(0)->setData(x0, y0);
            popup_plot->graph(0)->setPen(QPen(Qt::blue));
            // give the axes some labels:
            popup_plot->xAxis->setLabel("Day");
            popup_plot->yAxis->setLabel("Concentration (mM)");
            // set axes ranges
            popup_plot->xAxis->setAutoTickStep(false);
            popup_plot->xAxis->setTickStep(1);
            popup_plot->xAxis->setRange(0, ndays);
            popup_plot->yAxis->setAutoTickStep(false);
            popup_plot->yAxis->setTickStep(2);
            popup_plot->yAxis->setRange(0, 10);
        } else if (plotType == "drugKF") {
            plotName = "Drug Kill Fraction";
            plotwin->setWindowTitle(plotName);
            QString cellTypeStr, drugTypeStr;
            if (radioButton_drugcelltype_1->isChecked()) {
                cellTypeStr = "CT1";
            } else {
                cellTypeStr = "CT2";
            }
            if (radioButton_drugtype_1->isChecked()) {
                drugTypeStr = "PARENT";
            } else if (radioButton_drugtype_2->isChecked()) {
                drugTypeStr = "METAB1";
            } else {
                drugTypeStr = "METAB2";
            }
            QVector<double> x0(NPLOT), y0(NPLOT);
            double maxdrug;
            x0[0] = 1;
            makeDrugPlot(drugTypeStr, cellTypeStr, &maxdrug, "KF", &x0, &y0);
            if (x0[0] == 1) return; // does not kill
            // create graph and assign data to it:
            popup_plot->addGraph();
            popup_plot->graph(0)->setData(x0, y0);
            popup_plot->graph(0)->setPen(QPen(Qt::blue));
            // give the axes some labels:
            popup_plot->xAxis->setLabel("Drug concentration (mM)");
            popup_plot->yAxis->setLabel("Kill fraction/hour");
            // set axes ranges
            popup_plot->xAxis->setAutoTickStep(false);
            popup_plot->xAxis->setTickStep(maxdrug/5);
            popup_plot->xAxis->setRange(0, maxdrug);
            popup_plot->yAxis->setAutoTickStep(false);
            popup_plot->yAxis->setTickStep(0.2);
            popup_plot->yAxis->setRange(0, 1);
        } else if (plotType == "drugSF") {
            plotName = "Drug Survival Fraction";
            plotwin->setWindowTitle(plotName);
            QString cellTypeStr, drugTypeStr;
            if (radioButton_drugcelltype_1->isChecked()) {
                cellTypeStr = "CT1";
            } else {
                cellTypeStr = "CT2";
            }
            if (radioButton_drugtype_1->isChecked()) {
                drugTypeStr = "PARENT";
            } else if (radioButton_drugtype_2->isChecked()) {
                drugTypeStr = "METAB1";
            } else {
                drugTypeStr = "METAB2";
            }
            QVector<double> x0(NPLOT), y0(NPLOT);
            double maxdrug;
            x0[0] = 1;
            makeDrugPlot(drugTypeStr, cellTypeStr, &maxdrug, "SF", &x0, &y0);
            if (x0[0] == 1) return; // does not kill
            // create graph and assign data to it:
            popup_plot->addGraph();
            popup_plot->graph(0)->setData(x0, y0);
            popup_plot->graph(0)->setPen(QPen(Qt::blue));
            // give the axes some labels:
            popup_plot->xAxis->setLabel("Drug concentration (mM)");
            popup_plot->yAxis->setLabel("Survival fraction/hour");
            // set axes ranges
            popup_plot->xAxis->setAutoTickStep(false);
            popup_plot->xAxis->setTickStep(maxdrug/5);
            popup_plot->xAxis->setRange(0, maxdrug);
            popup_plot->yAxis->setRange(1.e-5, 1);
            popup_plot->yAxis->setScaleType(QCPAxis::stLogarithmic);
            popup_plot->yAxis->setScaleLogBase(10);
            popup_plot->yAxis->setNumberFormat("eb"); // e = exponential, b = beautiful decimal powers
            popup_plot->yAxis->setNumberPrecision(0); // makes sure "1*10^4" is displayed only as "10^4"
        } else if (plotType == "complementarySF") {
            plotName = "Drug + Radiation Survival Fraction";
            plotwin->setWindowTitle(plotName);
            QString cellTypeStr, drugTypeStr;
            if (radioButton_drugcelltype_1->isChecked()) {
                cellTypeStr = "CT1";
            } else {
                cellTypeStr = "CT2";
            }
            if (radioButton_drugtype_1->isChecked()) {
                drugTypeStr = "PARENT";
            } else if (radioButton_drugtype_2->isChecked()) {
                drugTypeStr = "METAB1";
            } else {
                drugTypeStr = "METAB2";
            }
            QVector<double> x0(NPLOT), y0(NPLOT);
            double maxO2;
            x0[0] = 1;
            makeDrugRadiationPlot(drugTypeStr, cellTypeStr, &maxO2, "SF", &x0, &y0);
            if (x0[0] == 1) return; // does not kill
            // create graph and assign data to it:
            popup_plot->addGraph();
            popup_plot->graph(0)->setData(x0, y0);
            popup_plot->graph(0)->setPen(QPen(Qt::blue));
            // give the axes some labels:
            popup_plot->xAxis->setLabel("O2 concentration (mM)");
            popup_plot->yAxis->setLabel("Survival fraction/hour");
            // set axes ranges
            popup_plot->xAxis->setAutoTickStep(false);
            popup_plot->xAxis->setTickStep(maxO2/5);
            popup_plot->xAxis->setRange(0, maxO2);
            popup_plot->yAxis->setRange(1.e-5, 1);
            popup_plot->yAxis->setScaleType(QCPAxis::stLogarithmic);
            popup_plot->yAxis->setScaleLogBase(10);
            popup_plot->yAxis->setNumberFormat("eb"); // e = exponential, b = beautiful decimal powers
            popup_plot->yAxis->setNumberPrecision(0); // makes sure "1*10^4" is displayed only as "10^4"
        }
    }
    plotwin->show();
}

//--------------------------------------------------------------------------------------------------------
// No O2 dependence of glucose consumption
// (NPLOT-1)*nt*dt = ndays*24*60*60 ==> number of time steps nt
//--------------------------------------------------------------------------------------------------------
void MainWindow::makeGlucosePlot(double *ndays, QVector<double> *x, QVector<double> *y)
{
    QLineEdit *line;
    double dt=100;  // sec
    double C, vol_cm3, t, metab, dCdt;
    double MM_C0, max_cell_rate;
    int ncells, Ng, nt;

    // Get parameter values from the GUI fields for glucose
    line = findChild<QLineEdit *>("lineEdit_glucose_ncells");
    ncells = line->text().toInt();
    line = findChild<QLineEdit *>("lineEdit_glucose_ndays");
    *ndays = line->text().toDouble();
    line = findChild<QLineEdit *>("line_MEDIUM_VOLUME");
    vol_cm3 = line->text().toDouble();
    line = findChild<QLineEdit *>("line_GLUCOSE_BDRY_CONC");
    C = line->text().toDouble();
    line = findChild<QLineEdit *>("line_GLUCOSE_MM_KM");
    MM_C0 = line->text().toDouble();
    line = findChild<QLineEdit *>("line_GLUCOSE_HILL_N");
    Ng = line->text().toInt();
    line = findChild<QLineEdit *>("line_GLUCOSE_CONSUMPTION");
    max_cell_rate = line->text().toDouble();

    max_cell_rate *= 1.0e6;     // mol/cell/s -> mumol/cell/s
    nt = ((*ndays)*24*60*60)/((NPLOT-1)*dt);
    (*x)[0] = 0;
    (*y)[0] = C;
    t = 0;
    for (int i=1; i<NPLOT; ++i)
    {
        for (int k=0; k<nt; k++) {
            metab = pow(C,Ng)/(pow(MM_C0,Ng) + pow(C,Ng));
//            qDebug("C: %f Ng: %d MM_C0: %f metab: %f",C,Ng,MM_C0,metab);
            dCdt = (-metab*max_cell_rate)/vol_cm3;	// convert mass rate (mol/s) to concentration rate (mM/s)
            dCdt = ncells*dCdt;
            C = C + dCdt*dt;
            if (C < 0) C = 0;
            t = t + dt;
        }
        (*x)[i] = t/(24*60*60);     // sec -> days
        (*y)[i] = C;
    }
}

//--------------------------------------------------------------------------------------------------------
// Assume that SER = 1 (no drug sensitisation)
//--------------------------------------------------------------------------------------------------------
void MainWindow::makeSFPlot(QString cellTypeStr, double C_O2, double maxdose, QVector<double> *x, QVector<double> *y)
{
    QLineEdit *line;
    double dose;

    // Get parameter values from the GUI fields for celltype
    QString objAlphaName = "line_RADIATION_ALPHA_H_" + cellTypeStr;
    QString objBetaName = "line_RADIATION_BETA_H_" + cellTypeStr;
    QString objOERAlphaName = "line_RADIATION_OER_ALPHA_" + cellTypeStr;
    QString objOERBetaName = "line_RADIATION_OER_BETA_" + cellTypeStr;
    QString objKmName = "line_RADIATION_KM_" + cellTypeStr;
    line = findChild<QLineEdit *>(objAlphaName);
    double LQ_alpha_H = line->text().toDouble();
    line = findChild<QLineEdit *>(objBetaName);
    double LQ_beta_H = line->text().toDouble();
    line = findChild<QLineEdit *>(objOERAlphaName);
    double LQ_OER_am = line->text().toDouble();
    line = findChild<QLineEdit *>(objOERBetaName);
    double LQ_OER_bm = line->text().toDouble();
    line = findChild<QLineEdit *>(objKmName);
    double LQ_K_ms = line->text().toDouble();

    double SER = 1;
    for (int i=0; i<NPLOT; ++i)
    {
        dose = (maxdose*i)/NPLOT;
        double OER_alpha_d = dose*(LQ_OER_am*C_O2 + LQ_K_ms)/(C_O2 + LQ_K_ms);
        double OER_beta_d = dose*(LQ_OER_bm*C_O2 + LQ_K_ms)/(C_O2 + LQ_K_ms);

        OER_alpha_d = OER_alpha_d*SER;
        OER_beta_d = OER_beta_d*SER;

        double expon = LQ_alpha_H*OER_alpha_d + LQ_beta_H*pow(OER_beta_d,2);
        double SF = exp(-expon);

        (*x)[i] = dose;
        (*y)[i] = SF;
    }

//    for (int i=0; i<101; ++i)
//    {
//      (*x)[i] = i/50.0 - 1; // x goes from -1 to 1
//      (*y)[i] = (*x)[i]*(*x)[i];  // let's plot a quadratic function
//    }

}

//--------------------------------------------------------------------------------------------------------
// drugTypeStr = "PARENT", "METAB1", "METAB2"
// cellTypeStr = "CT1", "CT2"
//
// Parameter objectname numbers:
//      Kmet0           0
//      C2              1
//      KO2             2
//      Vmax            3
//      Km              4
//      Klesion         5
//      kill_O2         6
//      kill_drug       7
//      kill_duration   8
//      kill_fraction   9
//      kills           13 --> 15
//      Killmodel       14 --> 16
// +
//      n_O2            13
//      death_prob      14
//
// Note that Kmet0, KO2, kill_duration need to be scaled to time units of sec
//--------------------------------------------------------------------------------------------------------
void MainWindow::makeDrugPlot(QString drugTypeStr, QString cellTypeStr, double *maxdose, QString plotStr, QVector<double> *x, QVector<double> *y)
{
    QLineEdit *line;
    QString objName0, objName;
    int killmodel;
    double C_O2, C2, Kmet0, KO2, n_O2, Ckill_O2, f, T, Ckill, Kd, c, dt;
    double Cdrug, kmet, dMdt, SF;   //, kill_prob;

    objName = "cbox_" + drugTypeStr + "_" + cellTypeStr + "_15";
    QCheckBox *cbox = findChild<QCheckBox *>(objName);
    if (!cbox->isChecked()) {
        LOG_MSG("Does not kill");
        return;     // Does not kill
    }

    objName0 = "line_" + drugTypeStr + "_" + cellTypeStr + "_";
    objName = objName0 + "16";
    line = findChild<QLineEdit *>(objName);
    killmodel = line->text().toInt();
    objName = objName0 + "0";
    line = findChild<QLineEdit *>(objName);
    Kmet0 = line->text().toDouble();
    Kmet0 = Kmet0/60;                       // /min -> /sec
    objName = objName0 + "1";
    line = findChild<QLineEdit *>(objName);
    C2 = line->text().toDouble();
    objName = objName0 + "2";
    line = findChild<QLineEdit *>(objName);
    KO2 = line->text().toDouble();
    KO2 = 1.0e-3*KO2;                       // um -> mM
    objName = objName0 + "13";
    line = findChild<QLineEdit *>(objName);
    n_O2 = line->text().toDouble();
    objName = objName0 + "6";
    line = findChild<QLineEdit *>(objName);
    Ckill_O2 = line->text().toDouble();     // kill_O2
    objName = objName0 + "7";
    line = findChild<QLineEdit *>(objName);
    Ckill = line->text().toDouble();        // kill_drug
    objName = objName0 + "8";
    line = findChild<QLineEdit *>(objName);
    T = line->text().toDouble();            // kill_duration
    T = 60*T;                               // min -> sec
//    printf("Kmet0: %f C2: %f KO2: %f Ckill_O2: %f Ckill: %f T: %f\n",Kmet0,C2,KO2,Ckill_O2,Ckill,T);

    objName = objName0 + "9";
    line = findChild<QLineEdit *>(objName);
    f = line->text().toDouble();            // kill_fraction

//    kmet = (1 - C2 + C2*KO2/(KO2 + Ckill_O2))*Kmet0;
    kmet = (1 - C2 + C2*pow(KO2,n_O2)/(pow(KO2,n_O2) + pow(Ckill_O2,n_O2)))*Kmet0;
    if (killmodel == 1) {
        Kd = -log(1-f)/(T*kmet*Ckill);
    } else if (killmodel == 2) {
        Kd = -log(1-f)/(T*kmet*pow(Ckill,2));
    } else if (killmodel == 3) {
        Kd = -log(1-f)/(T*pow(kmet*Ckill,2));
    } else if (killmodel == 4) {
        Kd = -log(1-f)/(T*Ckill);
    } else if (killmodel == 5) {
        Kd = -log(1-f)/(T*pow(Ckill,2));
    }

    line = findChild<QLineEdit *>("lineEdit_drug_O2");
    C_O2 = line->text().toDouble();
    line = findChild<QLineEdit *>("lineEdit_maxdrugconc");
    *maxdose = line->text().toDouble();

//    dt = 1;    // 1 sec
    for (int i=0; i<NPLOT; i++) {
        Cdrug = (i*(*maxdose)/(NPLOT-1));
//        kmet = (1 - C2 + C2*KO2/(KO2 + C_O2))*Kmet0;
        kmet = (1 - C2 + C2*pow(KO2,n_O2)/(pow(KO2,n_O2) + pow(C_O2,n_O2)))*Kmet0;
        dMdt = kmet*Cdrug;
        if (killmodel == 1) {
            c = Kd*dMdt;
        } else if (killmodel == 2) {
            c = Kd*dMdt*Cdrug;
        } else if (killmodel == 3) {
            c = Kd*pow(dMdt,2);
        } else if (killmodel == 4) {
            c = Kd*Cdrug;
        } else if (killmodel == 5) {
            c = Kd*pow(Cdrug,2);
        }
        SF = exp(-c*3600);      // 3600 sec = 1 hour
//        SF = 1;
//        for (int k=0; k<3600; k++) {      // calculate survival fraction after 3600 kill intervals of 1 sec = 1 hour
//            if (killmodel == 1) {
//                kill_prob = Kd*dMdt*dt;
//            } else if (killmodel == 2) {
//                kill_prob = Kd*dMdt*Cdrug*dt;
//            } else if (killmodel == 3) {
//                kill_prob = Kd*pow(dMdt,2)*dt;
//            } else if (killmodel == 4) {
//                kill_prob = Kd*Cdrug*dt;
//            } else if (killmodel == 5) {
//                kill_prob = Kd*pow(Cdrug,2)*dt;
//            }
//            kill_prob = min(kill_prob,1.0);
//            SF = SF*(1 - kill_prob);
//        }
        (*x)[i] = Cdrug;
        if (plotStr == "KF")
            (*y)[i] = (1 - SF);
        else
            (*y)[i] = SF;
    }
//    qDebug("maxdose: %8.4f  SF: %12.6f",*maxdose,SF);
/*
killmodel = dp%kill_model(ityp,im)		// could use %drugclass to separate kill modes

            drug(idrug)%Kmet0(ictyp,im) = drug(idrug)%Kmet0(ictyp,im)/60					! /min -> /sec
            drug(idrug)%KO2(ictyp,im) = 1.0e-3*drug(idrug)%KO2(ictyp,im)					! um -> mM
            drug(idrug)%kill_duration(ictyp,im) = 60*drug(idrug)%kill_duration(ictyp,im)	! min -> sec

                Ckill_O2 = drug(idrug)%kill_O2(ictyp,im)
                f = drug(idrug)%kill_fraction(ictyp,im)
                T = drug(idrug)%kill_duration(ictyp,im)
                Ckill = drug(idrug)%kill_drug(ictyp,im)
                kmet = (1 - C2 + C2*KO2/(KO2 + Ckill_O2))*Kmet0
                if (kill_model == 1) then
                    Kd = -log(1-f)/(T*kmet*Ckill)
                elseif (kill_model == 2) then
                    Kd = -log(1-f)/(T*kmet*Ckill**2)
                elseif (kill_model == 3) then
                    Kd = -log(1-f)/(T*(kmet*Ckill)**2)
                elseif (kill_model == 4) then
                    Kd = -log(1-f)/(T*Ckill)
                elseif (kill_model == 5) then
                    Kd = -log(1-f)/(T*Ckill**2)
                endif

kmet = (1 - dp%C2(ityp,im) + dp%C2(ityp,im)*dp%KO2(ityp,im)/(dp%KO2(ityp,im) + C_O2))*dp%Kmet0(ityp,im)
dMdt = kmet*cell_list(kcell)%conc(ichemo + im)
if (killmodel == 1) then
    kill_prob = kill_prob + Kd*dMdt*dt
elseif (killmodel == 2) then
    kill_prob = kill_prob + Kd*dMdt*cell_list(kcell)%conc(ichemo + im)*dt
elseif (killmodel == 3) then
    kill_prob = kill_prob + Kd*dMdt**2*dt
elseif (killmodel == 4) then
    kill_prob = kill_prob + Kd*cell_list(kcell)%conc(ichemo + im)*dt
elseif (killmodel == 5) then
    kill_prob = kill_prob + Kd*(cell_list(kcell)%conc(ichemo + im)**2)*dt
endif
*/
}

//--------------------------------------------------------------------------------------------------------
// The net survival fraction is the product of the survival fractions from the radiation dose and from the drug.
//--------------------------------------------------------------------------------------------------------
void MainWindow::makeDrugRadiationPlot(QString drugTypeStr, QString cellTypeStr, double *maxO2, QString plotStr, QVector<double> *x, QVector<double> *y)
{
    QLineEdit *line;
    QString objName0, objName, cellTypeNum;
    int killmodel;
    double C_O2, C2, Kmet0, KO2, n_O2, Ckill_O2, f, T, Ckill, Kd, c;
    double Cdrug, rad_dose, kmet, dMdt, SF_drug, SF_rad, SF;   //, kill_prob;

    objName = "cbox_" + drugTypeStr + "_" + cellTypeStr + "_15";
    QCheckBox *cbox = findChild<QCheckBox *>(objName);
    if (!cbox->isChecked()) {
        LOG_MSG("Does not kill");
        return;     // Does not kill
    }
//    qDebug("makeDrugRadiationPlot");
    cellTypeNum = cellTypeStr.mid(2,1);
    objName0 = "line_" + drugTypeStr + "_" + cellTypeStr + "_";
    objName = objName0 + "16";
    line = findChild<QLineEdit *>(objName);
    killmodel = line->text().toInt();
    objName = objName0 + "0";
    line = findChild<QLineEdit *>(objName);
    Kmet0 = line->text().toDouble();
    Kmet0 = Kmet0/60;                       // /min -> /sec
    objName = objName0 + "1";
    line = findChild<QLineEdit *>(objName);
    C2 = line->text().toDouble();
    objName = objName0 + "2";
    line = findChild<QLineEdit *>(objName);
    KO2 = line->text().toDouble();
    KO2 = 1.0e-3*KO2;                       // um -> mM
    objName = objName0 + "13";
    line = findChild<QLineEdit *>(objName);
    n_O2 = line->text().toDouble();
    objName = objName0 + "6";
    line = findChild<QLineEdit *>(objName);
    Ckill_O2 = line->text().toDouble();     // kill_O2
    objName = objName0 + "7";
    line = findChild<QLineEdit *>(objName);
    Ckill = line->text().toDouble();        // kill_drug
    objName = objName0 + "8";
    line = findChild<QLineEdit *>(objName);
    T = line->text().toDouble();            // kill_duration
    T = 60*T;                               // min -> sec
//    qDebug("Kmet0: %f C2: %f KO2: %f Ckill_O2: %f Ckill: %f T: %f",Kmet0,C2,KO2,Ckill_O2,Ckill,T);

    objName = objName0 + "9";
    line = findChild<QLineEdit *>(objName);
    f = line->text().toDouble();            // kill_fraction

    kmet = (1 - C2 + C2*pow(KO2,n_O2)/(pow(KO2,n_O2) + pow(Ckill_O2,n_O2)))*Kmet0;
    if (killmodel == 1) {
        Kd = -log(1-f)/(T*kmet*Ckill);
    } else if (killmodel == 2) {
        Kd = -log(1-f)/(T*kmet*pow(Ckill,2));
    } else if (killmodel == 3) {
        Kd = -log(1-f)/(T*pow(kmet*Ckill,2));
    } else if (killmodel == 4) {
        Kd = -log(1-f)/(T*Ckill);
    } else if (killmodel == 5) {
        Kd = -log(1-f)/(T*pow(Ckill,2));
    }

    line = findChild<QLineEdit *>("lineEdit_drug_O2");
    *maxO2 = line->text().toDouble();
    line = findChild<QLineEdit *>("lineEdit_maxdrugconc");
    Cdrug = line->text().toDouble();

    QString objAlphaName = "line_RADIATION_ALPHA_H_" + cellTypeNum;
    QString objBetaName = "line_RADIATION_BETA_H_" + cellTypeNum;
    QString objOERAlphaName = "line_RADIATION_OER_ALPHA_" + cellTypeNum;
    QString objOERBetaName = "line_RADIATION_OER_BETA_" + cellTypeNum;
    QString objKmName = "line_RADIATION_KM_" + cellTypeNum;
    line = findChild<QLineEdit *>(objAlphaName);
    double LQ_alpha_H = line->text().toDouble();
    line = findChild<QLineEdit *>(objBetaName);
    double LQ_beta_H = line->text().toDouble();
    line = findChild<QLineEdit *>(objOERAlphaName);
    double LQ_OER_am = line->text().toDouble();
    line = findChild<QLineEdit *>(objOERBetaName);
    double LQ_OER_bm = line->text().toDouble();
    line = findChild<QLineEdit *>(objKmName);
    double LQ_K_ms = line->text().toDouble();
    double SER = 1;

    line = findChild<QLineEdit *>("lineEdit_radiationdose");
    rad_dose = line->text().toDouble();

    for (int i=0; i<NPLOT; i++) {
        C_O2 = (i*(*maxO2)/(NPLOT-1));
        kmet = (1 - C2 + C2*pow(KO2,n_O2)/(pow(KO2,n_O2) + pow(C_O2,n_O2)))*Kmet0;
        dMdt = kmet*Cdrug;
        if (killmodel == 1) {
            c = Kd*dMdt;
        } else if (killmodel == 2) {
            c = Kd*dMdt*Cdrug;
        } else if (killmodel == 3) {
            c = Kd*pow(dMdt,2);
        } else if (killmodel == 4) {
            c = Kd*Cdrug;
        } else if (killmodel == 5) {
            c = Kd*pow(Cdrug,2);
        }
        SF_drug = exp(-c*3600);      // 3600 sec = 1 hour

        double OER_alpha_d = rad_dose*(LQ_OER_am*C_O2 + LQ_K_ms)/(C_O2 + LQ_K_ms);
        double OER_beta_d = rad_dose*(LQ_OER_bm*C_O2 + LQ_K_ms)/(C_O2 + LQ_K_ms);

        OER_alpha_d = OER_alpha_d*SER;
        OER_beta_d = OER_beta_d*SER;

        double expon = LQ_alpha_H*OER_alpha_d + LQ_beta_H*pow(OER_beta_d,2);
        SF_rad = exp(-expon);

        SF = SF_rad*SF_drug;
        (*x)[i] = C_O2;
        if (plotStr == "KF")
            (*y)[i] = (1 - SF);
        else
            (*y)[i] = SF;
    }
}
