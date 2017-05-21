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
#define KILL_NO2 13
#define KILL_PROB 14
#define KILL_Kd 15

#define KILL_kills 16
#define KILL_expt_kill_model 17
#define KILL_sensitises 18

#define NDPARAMS 5
#define NDKILLPARAMS 16
#define NIKILLPARAMS 3

struct kill_params {
    QString info[NDKILLPARAMS+NIKILLPARAMS];
    double dparam[NDKILLPARAMS];
    int iparam[NIKILLPARAMS];
};
typedef kill_params KILL_PARAMS;

struct drug_param_set {
    QString name;
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
