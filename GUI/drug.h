#ifndef DRUG_H
#define DRUG_H

#define NCELLTYPES 2
#define NDRUGS 2
#define TPZ_DRUG 0
#define DNB_DRUG 1
#define DRUG_A 0
#define DRUG_B 1

#define KILL_Kmet0 0
#define KILL_C2 1
#define KILL_KO2 2
#define KILL_Vmax 3
#define KILL_Km 4
#define KILL_Klesion 5
#define KILL_expt_O2_conc 6
#define KILL_expt_drug_conc 7
#define KILL_expt_duration 8
#define KILL_expt_kill_fraction 9
#define KILL_SER_max0 10
#define KILL_SER_Km 11
#define KILL_SER_KO2 12
#define KILL_kills 13
#define KILL_expt_kill_model 14
#define KILL_sensitises 15

#define NDPARAMS 5
#define NDKILLPARAMS 13
#define NIKILLPARAMS 3

//extern char *DRUG_param_name[];
//extern char *KILL_param_name[];

struct kill_params {
//    double Kmet0;
//    double C2;
//    double KO2;
//    double Vmax;
//    double Km;
//    double Klesion;
//    bool kills;
//    double expt_O2_conc;
//    double expt_drug_conc;
//    double expt_duration;
//    double expt_kill_fraction;
//    int expt_kill_model;
//    bool sensitises;
//    double SER_max0;
//    double SER_Km;
//    double SER_KO2;
    QString info[NDKILLPARAMS+NIKILLPARAMS];
    double dparam[NDKILLPARAMS];
    int iparam[NIKILLPARAMS];
};
typedef kill_params KILL_PARAMS;

struct drug_param_set {
    QString name;
//    double diff_coef;
//    double medium_diff;
//    double cell_diff_in;
//    double cell_diff_out;
//    double halflife;
    QString info[NDPARAMS];
    double dparam[NDPARAMS];
    KILL_PARAMS kill[NCELLTYPES];
};
typedef drug_param_set DRUG_PARAM_SET;

struct drug_str {
    QString classname;
    DRUG_PARAM_SET param[3];
};
typedef drug_str DRUG_STR;

extern DRUG_STR drug[NDRUGS];

#endif // DRUG_H
