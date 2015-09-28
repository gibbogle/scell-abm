#include "qmycheckbox.h"
#include "qmylabel.h"
#include "params.h"
#include "qdebug.h"

extern Params *parm;

QMyLabel::QMyLabel(QWidget *parent) : QLabel(parent)
{}
//--------------------------------------------------------------------------------------------------------
// Redefines mousePressEvent for QMyLabel, which extends QLabel.
// This is used to display info about a model parameter,
// or to display the info attached to an infolabel.
//--------------------------------------------------------------------------------------------------------
void QMyLabel::mousePressEvent (QMouseEvent *event) {
	event->accept();
    QString text;
    QString objName = objectName();
    if (objName.contains("infolabel")) {
        QString tag = objName.mid(10);
        parm->infoLabelInfo(tag,&text);
        if (text != "") {
            emit labelClicked(text);
        }
        return;
    }
    QString sname = objName.mid(6);
	// Find which label_ sent the signal, and read its text
	for (int k=0; k<parm->nParams; k++) {
		PARAM_SET param = parm->get_param(k);
		if (sname.compare(param.tag) == 0)
			text = param.text;
	}
    emit labelClicked(text);
};

QMyCheckBox::QMyCheckBox(QWidget *parent) : QCheckBox(parent)
{}
//--------------------------------------------------------------------------------------------------------
// Redefines mousePressEvent for QMyCheckBox, which extends QCheckBox.  This is used to display info about
// a model parameter.
//--------------------------------------------------------------------------------------------------------
void QMyCheckBox::mousePressEvent (QMouseEvent *event) {
    event->accept();
    QString text;
    if (objectName().contains("cbox_")) {
        QString sname = objectName().mid(5);
        // Find which cbox_ sent the signal, and read its text
        for (int k=0; k<parm->nParams; k++) {
            PARAM_SET param = parm->get_param(k);
            if (sname.compare(param.tag) == 0) {
                text = param.text;
            }
        }
    } else if (objectName().contains("checkBox_")) {
        text = this->description;
    }
    if (event->button() == Qt::LeftButton) {
        this->toggle();
    }
    emit checkBoxClicked(text);
};

QMyGroupBox::QMyGroupBox(QWidget *parent) : QGroupBox(parent)
{}
//--------------------------------------------------------------------------------------------------------
// Redefines mousePressEvent for QMyGroupBox, which extends QGroupBox.
// This is used to save a plot on right-button click.
//--------------------------------------------------------------------------------------------------------
void QMyGroupBox::mousePressEvent (QMouseEvent *event) {
    event->accept();
    if (event->button() == Qt::RightButton) {
        QString sname = objectName().mid(9);
        // Find which groupBox_ sent the signal
        LOG_QMSG("groupBox clicked: " + sname);
        emit groupBoxClicked(sname);
    }
};
