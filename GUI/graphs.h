#ifndef GRAPHS_H
#define GRAPHS_H

#include "global.h"

#define maxGraphs 16
#define TS_TYPE 0
#define PROF_TYPE 1
#define DIST_TYPE 2

struct graph_set {
	QString tag;
	QString title;
	QString yAxisTitle;
    QString description;
	int dataIndex;		// this must be consistent with the ordering of summary_data[]
	bool active;		// false for dummy graphs
	double maxValue;
	double scaling;		// multiplying factor for scaling of summary_data
    double yscale;
    int type;           // 0 = time-series, 1 = profile, 2 = distribution
};

typedef graph_set GRAPH_SET;

class Graphs
{
	GRAPH_SET *graphList;

public:

	Graphs();
	~Graphs();
    GRAPH_SET *tsGraphs;
    GRAPH_SET get_graph(int);
    int n_tsGraphs;
	int nGraphs;
	int get_dataIndex(int);
	QString get_tag(int);
	QString get_title(int);
	QString get_yAxisTitle(int);
    QString get_description(int);
	double get_maxValue(int);
	double get_scaling(int);
    double get_xscale(double);
    double get_yscale(int);
    int get_type(int);
	bool isActive(int);
    bool isTimeseries(int);
    bool isProfile(int);
    bool isDistribution(int);
    void set_maxValue(int, double);
    void makeGraphList(int);

};

#endif
