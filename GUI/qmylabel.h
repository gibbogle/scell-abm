#ifndef QMYLABEL_H
#define QMYLABEL_H

#include <qlabel.h>
#include <qgroupbox.h>
#include <QMouseEvent>

#include "log.h"
LOG_USE();

class QMyLabel: public QLabel 
{
    Q_OBJECT

public:
	QMyLabel(QWidget *parent = 0);
	
signals:
	void labelClicked(QString text);

private:

	void mousePressEvent (QMouseEvent *event);

};

class QMyGroupBox: public QGroupBox
{
    Q_OBJECT

public:
    QMyGroupBox(QWidget *parent = 0);

signals:
    void groupBoxClicked(QString text);

private:

    void mousePressEvent (QMouseEvent *event);

};

#endif

