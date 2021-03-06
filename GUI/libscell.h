#ifndef LIBBCELL32_H
#define LIBBCELL32_H

#ifdef __cplusplus
extern "C" {
#endif
//
//
void get_dll_build_version(char *, int *);
void execute(int *,char *, int *,char *, int *, double *);
void simulate_step(int *);
void terminate_run(int *);
void get_dimensions(int *, int *, int *, int *, double *, int *, int *, bool *, double *, double *);
void get_summary(double *, int *, int *);
void get_concdata(int *, int *, double *, double *);
void get_ic_concdata(int *, int *, double *, double *);
//void get_distdata(int *, DIST_SET *);
void get_volprob(int *, double *, double *, double*);
void get_oxyprob(int *, double *, double *, double *);
void get_nfacs(int *);
void get_facs(double *, double *, double *, int);
void get_histo(int, double *, double *, double *, double *, double *, double *, int);

void get_constituents(int *, int *, int *, char *, int *);

//void get_profile_oxygen(double *, double *, int *);

void test_array(int *, int *, double **);   // To test accessing a Fortran 2D array.
                                            // This returns a pointer to the array, which must be 'allocatable, target'.
//
//
#ifdef __cplusplus
}
#endif

#endif // LIBBCELL32_H
