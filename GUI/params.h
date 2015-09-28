#ifndef PARAMS_H
#define PARAMS_H

struct infoStruct {
    QString tag;
    QString info;
};
typedef infoStruct INFOSTRUCT;

struct param_set {
	QString tag;
	double value;
	double minvalue;
	double maxvalue;
	QString label;
	QString text;
};
typedef param_set PARAM_SET;

class Params 
{
	PARAM_SET *workingParameterList;
    INFOSTRUCT *workingInfolabelList;
//    INFOSTRUCT *workingInfocheckboxList;

public:

	Params();
	~Params();
	PARAM_SET get_param(int);
	int nParams;
	void set_value(int, double);
	void set_label(int, QString);

    int nInfolabel;
    void infoLabelInfo(QString, QString *);
    void get_labeltag(int, QString *);

//    int nInfocheckbox;
//    void infoCheckboxInfo(QString, QString *);
//    void get_checkboxtag(int, QString *);
};

#endif
