// drug.cpp
#include "mainwindow.h"
#include "QMessageBox"
#include "QFile"
#include <QDebug>

#include "drug.h"

DRUG_STR drug[NDRUGS];

char *DRUG_param_name[] = {"diff_coef","medium_diff","cell_diff_in","cell_diff_out","halflife"};
char *KILL_param_name[] = {"Kmet0","C2","KO2","Vmax","Km","Klesion","expt_O2_conc","expt_drug_conc","expt_duration",
                           "expt_kill_fraction","SER_max","SER_Km","SER_KO2","n_O2","death_prob","kills","expt_kill_model","sensitises"};


//--------------------------------------------------------------------------------------------------------
// From the drug index idrug (= TPZ_DRUG, DNB_DRUG, ...) and the drug name drugname, the default
// drug parameter file name is constructed.
// The file is opened and the contents are read into the corresponding DRUG_STR drug[idrug]
//--------------------------------------------------------------------------------------------------------
void MainWindow::readDrugParams(int idrug, QString fileName)
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
    QString line;
    int ival;
    double dval;

    line = in.readLine();
    QStringList data = line.split(" ",QString::SkipEmptyParts);
    drug[idrug].classname = data[0];
    for (int kset=0; kset<3; kset++) {     // 0 = parent, 1 = metab_1, 2 = metab_2
        for (int i=0; i<6; i++) {
            line = in.readLine();
            QStringList data = line.split(" ",QString::SkipEmptyParts);
            if (i > 0) dval = data[0].toDouble();
            if (i==0) {
                drug[idrug].param[kset].name = data[0];
//                qDebug() << drug[idrug].param[kset].name;
            } else {
                drug[idrug].param[kset].dparam[i-1] = dval;
                drug[idrug].param[kset].info[i-1] = (line.mid(data[0].size())).trimmed();
            }

        }
        for (int ictyp=0; ictyp<NCELLTYPES; ictyp++) {
            for (int i=0; i<NDKILLPARAMS+NIKILLPARAMS; i++) {
                line = in.readLine();
                QStringList data = line.split(" ",QString::SkipEmptyParts);
                if (i < NDKILLPARAMS) {
                    dval = data[0].toDouble();
                    drug[idrug].param[kset].kill[ictyp].dparam[i] = dval;
                } else {
                    ival = data[0].toInt();
                    drug[idrug].param[kset].kill[ictyp].iparam[i-NDKILLPARAMS] = ival;
                }
                drug[idrug].param[kset].kill[ictyp].info[i] = (line.mid(data[0].size())).trimmed();
//                if (kset == 0) {
//                    LOG_QMSG(drug[idrug].param[kset].kill[ictyp].info[i]);
//                }
            }
        }
    }
}

//--------------------------------------------------------------------------------------------------------
// Note: Need to ensure that the drug use checkboxes are kept up-to-date with the protocol specification.
// I.e. if a drug is requested in the protocol the checkbox must be set and the drug selected.
//--------------------------------------------------------------------------------------------------------
void MainWindow::writeDrugParams(QTextStream *out, int idrug)
{
    QString line;
    int nch;

    line = drug[idrug].classname;
    nch = line.length();
    for (int k=0; k<max(16-nch,1); k++)
        line += " ";
    line += "CLASS_NAME";
    *out << line + "\n";

    for (int kset=0; kset<3; kset++) {     // 0 = parent, 1 = metab_1, 2 = metab_2
        line = drug[idrug].param[kset].name;
        nch = line.length();
        for (int k=0; k<max(16-nch,1); k++)
            line += " ";
        if (kset == 0)
            line += "DRUG_NAME";
        else
            line += "METABOLITE";
        *out << line + "\n";

        for (int i=0; i<NDPARAMS; i++) {
            line = QString::number(drug[idrug].param[kset].dparam[i]);
            nch = line.length();
            for (int k=0; k<max(16-nch,1); k++)
                line += " ";
//            line += DRUG_param_name[i];
            line += drug[idrug].param[kset].info[i];
            *out << line + "\n";
        }

        for (int ictyp=0; ictyp<NCELLTYPES; ictyp++) {
            for (int i=0; i<NDKILLPARAMS; i++) {
                line = QString::number(drug[idrug].param[kset].kill[ictyp].dparam[i]);
                nch = line.length();
                for (int k=0; k<max(16-nch,1); k++)
                    line += " ";
//                line += KILL_param_name[i];
                line += drug[idrug].param[kset].kill[ictyp].info[i];
                *out << line + "\n";
            }
            for (int i=0; i<NIKILLPARAMS; i++) {
                line = QString::number(drug[idrug].param[kset].kill[ictyp].iparam[i]);
                nch = line.length();
                for (int k=0; k<max(16-nch,1); k++)
                    line += " ";
//                line += KILL_param_name[i+NDKILLPARAMS];
                line += drug[idrug].param[kset].kill[ictyp].info[i+NDKILLPARAMS];
                *out << line + "\n";
            }
        }
    }
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::on_pushButton_savedrugdata_clicked()
{
    int idrug;

    if (radioButton_drugA->isChecked()) {
        idrug = 0;
    } else if (radioButton_drugB->isChecked()) {
        idrug = 1;
    } else {
        return;
    }

    const QString fileName = QFileDialog::getSaveFileName(this, tr("Select Drug File"), ".", tr("Drug Data Files (*.drugdata)"));
    if (fileName.compare("") != 0) {
        LOG_MSG("Selected drug file:");
        LOG_QMSG(fileName);
        QFile file(fileName);
        if (!file.open(QFile::WriteOnly | QFile::Text)) {
            QMessageBox::warning(this, tr("Application"),
                                 tr("Cannot write file %1:\n%2.")
                                 .arg(fileName)
                                 .arg(file.errorString()));
            LOG_MSG("File open failed");
            return;
        }
        QTextStream out(&file);

        writeDrugParams(&out,idrug);
        file.close();
    }
}


//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::populateDrugTable(int idrug)
{
    QString line;
    bool flag;
    QString basestr, ctypstr, numstr, objname;

    LOG_QMSG("populateDrugTable");
    for (int kset=0; kset<3; kset++) {     // 0 = parent, 1 = metab_1, 2 = metab_2
        if (kset == 0) {
            basestr = "PARENT_";
        } else if (kset == 1) {
            basestr = "METAB1_";
        } else {
            basestr = "METAB2_";
        }
        line = drug[idrug].param[kset].name;

        for (int i=0; i<NDPARAMS; i++) {
            line = QString::number(drug[idrug].param[kset].dparam[i]);
            // Need to generate the name of the corresponding lineEdit widget
            numstr = QString::number(i);
            objname = "line_" + basestr + numstr;
            QLineEdit* qline = this->findChild<QLineEdit*>(objname);
            qline->setText(line);
        }

        for (int ictyp=0; ictyp<NCELLTYPES; ictyp++) {
            ctypstr = "CT" + QString::number(ictyp+1) + "_";
            for (int i=0; i<NDKILLPARAMS; i++) {
                line = QString::number(drug[idrug].param[kset].kill[ictyp].dparam[i]);
                numstr = QString::number(i);
                objname = "line_" + basestr + ctypstr + numstr;
                QLineEdit* qline = this->findChild<QLineEdit*>(objname);
                qline->setText(line);
            }
            for (int ii=0; ii<NIKILLPARAMS; ii++) {
                int i = ii + NDKILLPARAMS;
                if (i == KILL_kills || i == KILL_sensitises) {
                   flag = (drug[idrug].param[kset].kill[ictyp].iparam[ii] == 1);
                   numstr = QString::number(i);
                   objname = "cbox_" + basestr + ctypstr + numstr;
                   QCheckBox* qbox = this->findChild<QCheckBox*>(objname);
                   qbox->setChecked(flag);
                } else if (i == KILL_expt_kill_model) {
                    line = QString::number(drug[idrug].param[kset].kill[ictyp].iparam[ii]);
                    numstr = QString::number(i);
                    objname = "line_" + basestr + ctypstr + numstr;
                    QLineEdit* qline = this->findChild<QLineEdit*>(objname);
                    qline->setText(line);
                }
            }
        }
    }
    paramSaved = false;
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::on_buttonGroup_drug_buttonClicked(QAbstractButton* button)
{
    int idrug;
    QString drugname;
    LOG_MSG("on_buttonGroup_drug_buttonClicked");
//    int id = buttonGroup_drug->checkedId();
    QRadioButton *rb = (QRadioButton *)button;
//    drugname = rb->text();
//    LOG_QMSG("drugname: " + drugname);
//    if (id == 1) {
//        drugname =
//    }
//    if (rb->objectName().contains("drugA"))
//        idrug = 0;
//    else if (rb->objectName().contains("drugB"))
//        idrug = 1;
    if (radioButton_drugA->isChecked())
        populateDrugTable(0);
    if (radioButton_drugB->isChecked())
        populateDrugTable(1);
}


//--------------------------------------------------------------------------------------------------------
// One of the drug parameters has changed, need to copy the change to the array drug[]
//--------------------------------------------------------------------------------------------------------
void MainWindow::changeDrugParam(QObject *w)
{
    int idrug, kset, ict, ictyp, i, ii;
    // get idrug
    if (radioButton_drugA->isChecked()) {
        idrug = 0;
    } else if (radioButton_drugB->isChecked()) {
        idrug = 1;
    } else {
        return;
    }
    QString name = w->objectName();
    // get number i
    QStringList data = name.split("_");
    int n = data.size();
    QString numstr = data[n-1];
    i = numstr.toInt();
    ii = i - NDKILLPARAMS;
    if (name.contains("PARENT")) {
        kset = 0;
    } else if (name.contains("METAB1")) {
        kset = 1;
    } else if (name.contains("METAB2")) {
        kset = 2;
    }
    if (name.contains("_CT1")) {            // cell type 1  data
        ict = 1;
    } else if (name.contains("_CT2")) {     // cell type 2  data
        ict = 2;
    } else {                                // basic data
        ict = 0;
    }
    if (ict == 0) {
        QLineEdit *lineEdit = (QLineEdit *)w;
        drug[idrug].param[kset].dparam[i] = lineEdit->text().toDouble();
        return;
    }
    ictyp = ict-1;
    if (name.contains("line_")) {
        QLineEdit *lineEdit = (QLineEdit *)w;
        if (i == KILL_expt_kill_model)
            drug[idrug].param[kset].kill[ictyp].iparam[ii] = lineEdit->text().toInt();
        else
            drug[idrug].param[kset].kill[ictyp].dparam[i] = lineEdit->text().toDouble();
    } else if (name.contains("cbox_")) {
//        qDebug()<<"drug checkbox changed: " + w->objectName() + "\n";
        QCheckBox *cbox = (QCheckBox *)w;
        int val;
        if (cbox->isChecked())
            val = 1;
        else
            val = 0;
        drug[idrug].param[kset].kill[ictyp].iparam[ii] = val;
    }
}


//--------------------------------------------------------------------------------------------------------
// Make QLists of file names postfixed by ".drugdata"
//--------------------------------------------------------------------------------------------------------
void MainWindow::makeDrugFileLists()
{
    QDir myDir(".");
    QStringList filters;

    Drug_FilesList.clear();
    filters << "*.drugdata";
    myDir.setNameFilters(filters);
    Drug_FilesList = myDir.entryList(filters);
    LOG_MSG("Drug files");
    for (int i=0; i<Drug_FilesList.size(); i++) {
        LOG_QMSG(Drug_FilesList[i]);
    }
}

//--------------------------------------------------------------------------------------------------------
// This may be superceded by a widget for selecting a file ("TPZ_..." or "DNB_..." or ...) then deducing
// the drug name from the file name.  This will ensure that only a drug with a parameter file can be chosen.
//--------------------------------------------------------------------------------------------------------
void MainWindow::initDrugComboBoxes()
{
    for (int i=0; i<Drug_FilesList.size(); i++) {
        comb_DRUG_A->addItem(Drug_FilesList[i]);
        comb_DRUG_B->addItem(Drug_FilesList[i]);
    }

}

//--------------------------------------------------------------------------------------------------------
// Read the drug section in the input file
//--------------------------------------------------------------------------------------------------------
void MainWindow::readDrugData(QTextStream *in)
{
    QString line;

    line = in->readLine();
    QStringList data = line.split(" ",QString::SkipEmptyParts);
    int ndrugs = data[0].toInt();
    if (ndrugs == 0) return;
    makeDrugFileLists();
    for (int idrug=0; idrug<ndrugs; idrug++) {
        line = in->readLine();  // CLASS_NAME
        line = in->readLine();  // DRUG_NAME
        QStringList data = line.split(" ",QString::SkipEmptyParts);
        QString drugname = data[0];
        // choose a drugdata file
        bool success;
        QString fileName;
        for (;;) {
            bool ok;
            QString item = QInputDialog::getItem(this, tr("QInputDialog::getItem()"),
                                                     "Select a drug file for: "+drugname, Drug_FilesList, 0, false, &ok);
            if (ok && !item.isEmpty()) {
                fileName = item;
                readDrugParams(idrug, fileName);
                if (fileName.contains(drugname)) {
                    success = true;
                    break;
                } else {
                    QMessageBox::StandardButton reply;
                    reply = QMessageBox::critical(this, tr("QMessageBox::critical()"),
                                                  "This is not a file for drug: "+drugname,
                                                    QMessageBox::Abort | QMessageBox::Retry);
                    if (reply == QMessageBox::Abort) {
                        success = false;
                        break;
                    }
                }
            }
        }
        if (success) {
            if (idrug == 0) {
                // Make fileName selected in comb_DRUG_A
                int k = Drug_FilesList.indexOf(fileName);
                comb_DRUG_A->setCurrentIndex(k);
//                text_DRUG_A_NAME->setText(drugname);
                cbox_USE_DRUG_A->setChecked(true);
            } else if (idrug == 1) {
                // Make fileName selected in comb_DRUG_B
                int k = Drug_FilesList.indexOf(fileName);
                comb_DRUG_B->setCurrentIndex(k);
//                text_DRUG_B_NAME->setText(drugname);
                cbox_USE_DRUG_B->setChecked(true);
            }
        }
        int nskip = 3*(6 + 2*(NDKILLPARAMS + NIKILLPARAMS)) - 1;
        for (int k=0; k<nskip; k++) {
            line = in->readLine();
        }
    }
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::makeHeaderText(QString *header, bool interact)
{
    QString name = qgetenv("USER");
    if (name.isEmpty())
        name = qgetenv("USERNAME");
    QString datestr, text;
    datestr = QDate::currentDate().toString("dd/MM/yy");
    text = "User: " + name;
    if (cbox_USE_DRUG_A->isChecked()) {
        QString drugname = comb_DRUG_A->currentText();
        int i = drugname.indexOf('.');
        drugname = drugname.left(i);
        text += " DrugA: " + drugname;
    }
    if (cbox_USE_DRUG_B->isChecked()) {
        QString drugname = comb_DRUG_B->currentText();
        int i = drugname.indexOf('.');
        drugname = drugname.left(i);
        text += " DrugB: " + drugname;
    }
    if (!interact) {
        text.replace(" ","_");
        *header = datestr + " " + text;
        return;
    }
    bool ok;
    text = QInputDialog::getText(this, tr("QInputDialog::getText()"),
                                         "Header:", QLineEdit::Normal,
                                         text, &ok);
    text.replace(" ","_");
    *header = datestr + " " + text;
}

