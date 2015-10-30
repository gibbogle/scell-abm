#include "global.h"

namespace Global
{
    QString GUI_build_version;
    QString DLL_build_version;

    int MAX_CHEMO;
    int N_EXTRA;
    int NX, NY, NZ;
    double DELTA_T;
    double DELTA_X;
    double dfraction;
    int nt_vtk;
    int istep;
    bool leftb;

    int nvars_used;
    int nfieldvars_used;
    int GUI_to_DLL_index[32];
    int DLL_to_GUI_index[32];
    QString var_string[32];

    double *FACS_data=NULL;
    int nFACS_cells=0;
    int nFACS_dim=0;

    double *histo_data=NULL;
    double *histo_data_log=NULL;
    int nhisto_bins;
    int nhisto_dim=0;
    double histo_vmin[3*32];
    double histo_vmax[3*32];
    double histo_vmin_log[3*32];
    double histo_vmax_log[3*32];
    int histo_celltype=0;

    int summaryData[100];
    int i_hypoxia_cutoff;
    int i_growth_cutoff;

    double concData[4000];
    int conc_nvars;
    int conc_nc;
    double conc_dx;
    QString casename;

    double volProb[100];
    int vol_nv = 20;
    double vol_v0;
    double vol_dv;
    double oxyProb[100];
    int oxy_nv = 20;
    double oxy_v0;
    double oxy_dv;

//    double distData[4000];
//    bool dist_used[20];
    int dist_nv;
    DIST_SET distParams[20];

//    int cell_list[N_CELLINFO*MAX_CELLS];
    int ncell_list;
    CELL_DATA cell_list[MAX_CELLS];
    double blobcentre[3];

//    double *profile_x[20];
//    double *profile_y[20];
//    int profile_n[20];

    bool showingVTK;
    bool recordingVTK;
    bool showingFACS;
    bool recordingFACS;
    bool showingField;
    bool recordingField;


} // namespace Global
