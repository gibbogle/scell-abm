#include "mainwindow.h"
#include "QMessageBox"
#include "QFile"
#include <QDebug>

void MainWindow::SetupProtocol()
{
    QStringList tableHeader;

    tableWidget->setRowCount(64);
    tableWidget->setColumnCount(7);
    tableHeader<<"Hour"<<"Drug"<<"Duration"<<"Volume"<<"Conc"<<"Radiation"<<"Medium vol";
    tableWidget->setHorizontalHeaderLabels(tableHeader);

    connect(tableWidget, SIGNAL(cellChanged(int, int)  ),this, SLOT(ProtocolChanged(int, int)) );

}

void MainWindow::LoadProtocol(QString fileName)
{
    int nTimes;
    QTableWidget *table = tableWidget;
    table->clearContents();
    QFile file(fileName);
    if (!file.open(QFile::ReadOnly | QFile::Text)) {
        QMessageBox::warning(this, tr("Application"),
                             tr("Cannot read file %1:\n%2.")
                             .arg(fileName)
                             .arg(file.errorString()));
        return;
    }
    qDebug() << "Opened " << fileName;
    QTextStream in(&file);
    QString line;
    // Skip lines
    for (;;) {
        line = in.readLine();
//        qDebug() << line;
        if (line.compare("PROTOCOL")==0) break;
        if (in.atEnd()) {
            QMessageBox::warning(this, "LoadProtocol", "No protocol lines in the input file");
            return;
        }
    }
//    line = in.readLine();
//    qDebug() << line;
    line = in.readLine();
//    qDebug() << line;
    nTimes = line.toInt();
    for (int row=0; row<nTimes; row++) {
        QString mode = in.readLine();
//        qDebug() << mode;
        if (mode.compare("DRUG") == 0) {
            QString drug = in.readLine();
//            qDebug() << drug;
            setField(table, row, 1, drug);
            QString hour = in.readLine();
//            qDebug() << hour;
            setField(table, row, 0, hour);
            QString duration = in.readLine();
//            qDebug() << duration;
            setField(table, row, 2, duration);
            QString volume = in.readLine();
//            qDebug() << volume;
            setField(table, row, 3, volume);
            QString conc = in.readLine();
//            qDebug() << conc;
            setField(table, row, 4, conc);
        } else if (mode.compare("RADIATION") == 0) {
            QString hour = in.readLine();
//            qDebug() << hour;
            setField(table, row, 0, hour);
            QString dose = in.readLine();
//            qDebug() << dose;
            setField(table, row, 5, dose);
        } else if (mode.compare("MEDIUM") == 0) {
            QString hour = in.readLine();
//            qDebug() << hour;
            setField(table, row, 0, hour);
            QString volume = in.readLine();
//            qDebug() << volume;
            setField(table, row, 6, volume);
        }
    }
    paramSaved = false;
}

void MainWindow::SaveProtocol(QTextStream *out, int ndrugs)
{
    int nTimes, eventType, kevents;
    int err;
    QString entry, hour;
    QString drugEntry, radiationEntry, mediumEntry;
    QTableWidgetItem *item;
    QTableWidget *table = tableWidget;
    QMessageBox msgBox;
    int idrug=1, iradiation=2,imedium=3;

    int row = tableWidget->currentRow();
    int col = tableWidget->currentColumn();
    tableWidget->setCurrentCell(row+1,col);
    /*
    QFile file(fileName);
    if (!file.open(QFile::WriteOnly | QFile::Text | QIODevice::Append)) {
        QMessageBox::warning(this, tr("Application"),
                             tr("Cannot write file %1:\n%2.")
                             .arg(fileName)
                             .arg(file.errorString()));
        qDebug() << "File open failed";
        return;
    }
    QTextStream out(&file);
    */
    *out << "PROTOCOL\n";
    nTimes = 0;
    if (ndrugs == 0) {
        for (int row=0; row<tableWidget->rowCount(); row++) {
            item = tableWidget->item(row,5);
            if (item != 0) {
                entry = item->text();
                if (entry.compare("")) {    // true if <>
                    nTimes++;
                }
            }
            item = tableWidget->item(row,6);
            if (item != 0) {
                entry = item->text();
                if (entry.compare("")) {    // true if <>
                    nTimes++;
                }
            }
        }
        if (nTimes == 0) {      // nothing in the protocol
            *out << 0 << "\n";
            return;
        }
    }
    nTimes = 0;
    for (int row=0; row<tableWidget->rowCount(); row++) {
        item = tableWidget->item(row,0);
        if (item != 0) {
            entry = item->text();
            if (entry.compare("")) {    // true if <>
                nTimes++;
            }
        }
    }
    *out << nTimes << "\n";
    for (int row=0; row<tableWidget->rowCount(); row++) {
        item = tableWidget->item(row,0);
        if (item != 0) {
            entry = item->text();
            if (entry.compare("")) {    // true if <>
                hour = entry;
                kevents = 0;
                err = getField(table,row,1,&drugEntry);
                if (drugEntry.compare("")) {   // Entry in DRUG column
                    eventType = idrug;
                    kevents++;
                }
                err = getField(table,row,5,&radiationEntry);
                if (radiationEntry.compare("")) {   // Entry in RADIATION column
                    eventType = iradiation;
                    kevents++;
                }
                err = getField(table,row,6,&mediumEntry);
                if (mediumEntry.compare("")) {   // Entry in MEDIUM column
                    eventType = imedium;
                    kevents++;
                }
                if (kevents == 0) {
                    QString msg = "No DRUG, RADIATION or MEDIUM data for event at hour: " + hour;
                    msgBox.setText(msg);
                    msgBox.exec();
//                    file.close();
                    return;
                }
                if (kevents > 1) {
                    QString msg = "More than one event at hour: " + hour;
                    msgBox.setText(msg);
                    msgBox.exec();
//                    file.close();
                    return;
                }
                if (eventType == idrug) {
                    *out << "DRUG" << "\n";
                    *out << drugEntry << "\n";
                    *out << hour << "\n";
                    err = getField(table,row,2,&entry);
                    if (entry.compare("")) {   // Entry in Duration column
                        *out << entry << "\n";
                    } else {
                        msgBox.setText("Missing entry in Duration column");
                        msgBox.exec();
                        qDebug() << "Missing entry in Duration column";
//                        file.close();
                        return;
                    }
                    err = getField(table,row,3,&entry);
                    if (entry.compare("")) {   // Entry in Volume column
                        *out << entry << "\n";
                    } else {
                        msgBox.setText("Missing entry in Volume column");
                        msgBox.exec();
                        qDebug() << "Missing entry in Volume column";
//                        file.close();
                        return;
                    }
                    err = getField(table,row,4,&entry);
                    if (entry.compare("")) {   // Entry in Conc column
                        *out << entry << "\n";
                    } else {
                        msgBox.setText("Missing entry in Conc column");
                        msgBox.exec();
                        qDebug() << "Missing entry in Conc column";
//                        file.close();
                        return;
                    }
                    continue;
                }
//                err = getField(table,row,5,&entry);
//                if (entry.compare("")) {   // Entry in RADIATION column
                else if (eventType == iradiation) {
                    *out << "RADIATION" << "\n";
                    *out << hour << "\n";
                    *out << radiationEntry << "\n";
                    continue;
                }
//                err = getField(table,row,6,&entry);
//                if (entry.compare("")) {   // Entry in MEDIUM column
                else if (eventType == imedium) {
                    *out << "MEDIUM" << "\n";
                    *out << hour << "\n";
                    *out << mediumEntry << "\n";
                    continue;
                }
//                qDebug() << "No DRUG, RADIATION or MEDIUM data for event at hour: " << hour;
//                return;
            }
        }
    }
//    file.close();
}

bool MainWindow::ProtocolUsesDrug()
{
    QTableWidgetItem *item;

    for (int row=0; row<tableWidget->rowCount(); row++) {
        item = tableWidget->item(row,1);
        if (item != 0) {
            if (item->text().compare("")) {    // true if <>
                return true;
            }
        }
    }
    return false;
}

int MainWindow::getField(QTableWidget *table, int row, int col, QString *entry)
{
    int err=0;
    QTableWidgetItem *item;
    item = table->item(row,col);
    if (item != 0) {
        *entry = item->text();
    } else {
        *entry = "";
    }
    return err;
}

void MainWindow::setField(QTableWidget *table, int row, int col, QString entry)
{
    QTableWidgetItem *item = new QTableWidgetItem;
    item->setText(entry);
    table->setItem(row,col,item);
}

void MainWindow::ProtocolChanged(int row, int col)
{
    paramSaved = false;
    if (col == 1) {
        QTableWidgetItem *item;
        item = tableWidget->item(row,col);
        QString drugname = item->text();
        LOG_QMSG(drugname);
        if (drugname == "") return;
        if (drugname != text_DRUG_A_NAME->text() && drugname != text_DRUG_B_NAME->text()) {
            LOG_QMSG("drug " + drugname + " is not selected");
            QMessageBox::warning(this, tr("Protocol"),
                                 tr("This drug has not been selected: %1\n %2")
                                 .arg(drugname)
                                 .arg("Select on the Treatment screen."));
            item->setText("");
        }
    }
}
