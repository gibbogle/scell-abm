#include <qstring.h>
#include "params.h"

Params::Params()
{
    PARAM_SET params[] = {

/*
{"BC_AVIDITY_MEDIAN", 1.0, 0.1, 10.0,
"BCR avidity median parameter",
"BCR avidity has a lognormal distribution, described by the median and shape parameters.\n\
(BCR stimulation rate is proportional to the product of BC avidity and antigen load.)"},

{"BC_AVIDITY_SHAPE", 1.1, 1.01, 3.0,
"BCR avidity shape parameter",
"BCR avidity has a lognormal distribution, described by the median and shape parameters.\n\
The shape value must be greater than 1, and values close to 1 give distributions that are close to normal."},
*/

{"GUI_VERSION_NAME", 0, 0, 0,
 "GUI0.00",
 "GUI program version number."},

{"DLL_VERSION_NAME", 0, 0, 0,
 "DLL0.00",
 "DLL version number."},

{"NX", 33, 0, 0,
"Fine grid size",
"Dimension of the fine grid (number of grid pts in X,Y and Z directions).  Must = 1 + multiple of 8."},

{"INITIAL_COUNT", 1000, 0, 0,
"Initial number of tumour cells",
"Initial number of tumour cells"},

{"DIVIDE_TIME_1_MEDIAN", 24, 0, 0,
"Median (h)",
"The time taken for tumour cell division has a lognormal distribution, described by the median and shape parameters. \n\
[hours]"},

{"DIVIDE_TIME_1_SHAPE", 1.2, 0, 0,
"Shape parameter",
"The time taken for tumour cell division has a lognormal distribution, described by the median and shape parameters."},

{"DIVIDE_TIME_2_MEDIAN", 24, 0, 0,
"Division time median parameter",
"The time taken for tumour cell division has a lognormal distribution, described by the median and shape parameters. \n\
[hours]"},

{"DIVIDE_TIME_2_SHAPE", 1.2, 0, 0,
"Division time shape parameter",
"The time taken for tumour cell division has a lognormal distribution, described by the median and shape parameters."},

{"V_DEPENDENT_GROWTH_RATE", 0, 0, 1,
"V-dependent growth rate",
"The growth rate of a cell is proportional to the volume."},

{"RANDOMISE_INITIAL_V", 1, 0, 1,
"Randomise initial cell volumes",
"The volumes of the initial cell population are randomised."},

{"NDAYS", 10.0, 0.0, 30.0,
"Number of days",
"Length of the simulation.\n\
[days]"},

{"DELTA_T", 600, 0, 0,
"Time step",
"Length of main time step, for cell death, division, etc.  Should be a divisor of 3600. \n\
[mins]"},

{"NXB", 35, 0, 0,
"Coarse grid size",
"Dimension of the coarse grid (number of grid pts in X,Y and Z directions).  Grid spacing is 4 times fine grid spacing.  Must be odd"},

{"DELTA_X", 30, 0, 0,
"Fine grid spacing (um)",
"Grid-cell size in um.  Constituent transport and consumption/production is computed on this grid."},

{"A_SEPARATION", 1.0, 0, 0,
"Separation force factor",
"During mitosis the two capped spheres are effectively connected by a nonlinear spring. \n\
The length of the spring s is determined by the mitosis level, and if the centre-centre distance is d \n\
the contribution to the force of repulsion between the spheres is a_separation*(s-d)^3."},

{"A_FORCE", 1., 0, 0,
"Repulsion force factor 'a'",
"The cell-cell force is a function of x = distance/(sum of radii): x = d/(R1+R2). \n\
The force function is: F(x) = a/((x-x0)(x1-x)) + b, where x0 and x1 are the locations of the bounding asymptotes. \n\
The parameter 'b' is calculated by setting the minimum value of F(x) (which occurs at x = (x0+x1)/2 ) equal to -c, this is the maximum attraction force. \n\
The function F(x) is zero at two points, xc1 and xc2.  There is a hysteresis loop for x > xc1: for two cells not in contact, the force is zero for x > xc1. \n\
After contact is made the force is non-zero until x > xc2 - this is the effect of cell-cell adhesion."},

{"C_FORCE", 0.2, 0, 0,
"Attraction force factor 'c'",
"The cell-cell force is a function of x = distance/(sum of radii): x = d/(R1+R2). \n\
The force function is: F(x) = a/((x-x0)(x1-x)) + b, where x0 and x1 are the locations of the bounding asymptotes. \n\
The parameter 'b' is calculated by setting the minimum value of F(x) (which occurs at x = (x0+x1)/2 ) equal to -c, this is the maximum attraction force. \n\
The function F(x) is zero at two points, xc1 and xc2.  There is a hysteresis loop for x > xc1: for two cells not in contact, the force is zero for x > xc1. \n\
After contact is made the force is non-zero until x > xc2 - this is the effect of cell-cell adhesion."},

{"X0_FORCE", 0.7, 0, 0,
"Left asymptote 'x0'",
"The cell-cell force is a function of x = distance/(sum of radii): x = d/(R1+R2). \n\
The force function is: F(x) = a/((x-x0)(x1-x)) + b, where x0 and x1 are the locations of the bounding asymptotes. \n\
The parameter 'b' is calculated by setting the minimum value of F(x) (which occurs at x = (x0+x1)/2 ) equal to -c, this is the maximum attraction force. \n\
The function F(x) is zero at two points, xc1 and xc2.  There is a hysteresis loop for x > xc1: for two cells not in contact, the force is zero for x > xc1. \n\
After contact is made the force is non-zero until x > xc2 - this is the effect of cell-cell adhesion."},

{"X1_FORCE", 1.3, 0, 0,
"Right asymptote 'x1'",
"The cell-cell force is a function of x = distance/(sum of radii): x = d/(R1+R2). \n\
The force function is: F(x) = a/((x-x0)(x1-x)) + b, where x0 and x1 are the locations of the bounding asymptotes. \n\
The parameter 'b' is calculated by setting the minimum value of F(x) (which occurs at x = (x0+x1)/2 ) equal to -c, this is the maximum attraction force. \n\
The function F(x) is zero at two points, xc1 and xc2.  There is a hysteresis loop for x > xc1: for two cells not in contact, the force is zero for x > xc1. \n\
After contact is made the force is non-zero until x > xc2 - this is the effect of cell-cell adhesion."},

{"KDRAG", 5, 0, 0,
 "Drag factor",
 "Displacement = dt*F/drag"},

{"FRANDOM",0.02, 0, 0,
 "Random force factor",
 "Magnitude of random additive force"},

{"NT_CONC", 1, 0, 0,
"Number of ODE solver sub-steps.",
"The number of subdivisions of the major time step, for the ODE diffusion-reaction solver."},

{"NMM3", 500000, 0, 0,
"Cells/cubic mm",
"Number of cells per cubic mm of non-necrotic tumour."},

{"FLUID_FRACTION", 0.5, 0, 0,
"Fluid fraction",
"Fraction of non-necrotic tumour that is extracellular fluid."},

{"MEDIUM_VOLUME", 0.074, 0, 0,
"Medium volume",
"Volume of the medium in which the spheroid is growing."},

{"UNSTIRRED_LAYER", 0.001, 0, 0,
"Unstirred layer width",
"Thickness of the unstirred layer around the spheroid (cm)."},

{"VDIVIDE0", 1.6, 0, 0,
"Nominal divide volume",
"Nominal multiple of normal cell volume at which division occurs."},

{"DVDIVIDE", 0.3, 0, 0,
"Divide volume variation",
"Variation (+/-) about nominal divide volume multiple."},

{"MM_THRESHOLD", 0.1, 0, 0,
"Michaelis-Menten O2 threshold",
"O2 concentration at which the 'soft-landing' adjustment to the Michaelis-Menten function kicks in.\n\
[uM]"},

{"ANOXIA_THRESHOLD", 0.15, 0, 0,
"Anoxia threshold",
"A cell begins to experience anoxia leading to cell death at the O2 concentration given by this threshold value."},

{"ANOXIA_TAG_TIME", 3.0, 0, 0,
"Anoxic time limit",
"Length of time under hypoxia (O2 < anoxic threshold) after which a cell is tagged to die of anoxia.\n\
[h]"},

{"ANOXIA_DEATH_TIME", 3.0, 0, 0,
"Anoxic death delay time",
"Time taken for a cell to die after it is tagged to die of anoxia.\n\
[h]"},

{"TEST_CASE", 0, 0, 0,
"Test case #",
"Number of the test case to run.  The default value of 0 is for a normal run"},

{"SEED1", 1234, 0, 0,
"First RNG seed",
"The random number generator is seeded by a pair of integers.  Changing the seed generates a different Monte Carlo realization."},

{"SEED2", 5678, 0, 0,
"Second RNG seed",
"The random number generator is seeded by a pair of integers.  Changing the seed generates a different Monte Carlo realization."},

{"NCPU", 4, 1, 8,
"Number of CPUs",
"Number of CPUs to use for the simulation."},

{"NCELLTYPES", 2, 0, 0,
"Number of cell types",
"Maximum number of cell types in the spheroid.  The initial percentage of each type must be specified"},

{"CELLPERCENT_1", 100, 0, 100,
"Percentage of cell type 1",
"Percentage of cell type 1"},

//{"CELLDISPLAY_1", 1, 0, 1,
//"Display cell type 1",
//"Display cell type 1"},

{"CELLPERCENT_2", 0, 0, 100,
"Percentage of cell type 2",
"Percentage of cell type 2"},

//{"CELLDISPLAY_2", 1, 0, 1,
//"Display cell type 2",
//"Display cell type 2"},

{"NT_ANIMATION", 1, 0, 0,
 "Animation interval (timesteps)",
 "Interval between animation screen updates (timesteps).  One timestep = 15 sec."},

{"SHOW_PROGENY", 0, 0, 0,
 "Show descendents of cell #",
 "All the descendents of cell with the specified ID are highlighted.  (0 = no selection)"},

{"USE_OXYGEN", 1, 0, 1,
"Use Oxygen?",
"Oxygen is simulated"},

{"OXYGEN_DIFF_COEF", 2.0e-5, 0, 0,
 "Spheroid diffusion coeff",
 "Constituent diffusion coefficient in the spheroid"},

{"OXYGEN_MEDIUM_DIFF", 2.5e-5, 0, 0,
 "Medium diffusion coeff",
 "Constituent diffusion coefficient in the medium"},

{"OXYGEN_CELL_DIFF", 600, 0, 0,
 "Membrane diff constant",
 "Cell membrane diffusion constant Kd"},

{"OXYGEN_BDRY_CONC", 0.18, 0, 0,
 "Boundary concentration",
 "Constituent concentration in the medium"},

{"OXYGEN_CONSTANT", 0, 0, 1,
 "Constant concentration",
 "Extracellular concentration to be held constant everywhere at the specified boundary value"},

{"OXYGEN_CONSUMPTION", 6.25e-17, 0, 0,
 "Max consumption rate",
 "Maximum rate of consumption of the constituent"},

{"OXYGEN_MM_KM", 1.33, 0, 0,
 "Michaelis-Menten Km",
 "Michaelis-Menten Km (uM)"},

{"OXYGEN_HILL_N", 1, 1, 2,
 "Hill function N",
 "Oxygen uptake rate Hill function N"},

{"USE_GLUCOSE", 1, 0, 1,
"Use Glucose?",
"Glucose is simulated"},

{"GLUCOSE_DIFF_COEF", 3.0e-6, 0, 0,
 "Spheroid diffusion coeff",
 "GLUCOSE diffusion coefficient"},

{"GLUCOSE_MEDIUM_DIFF", 6.0e-6, 0, 0,
 "Medium diffusion coeff",
 "Constituent diffusion coefficient in the medium"},

{"GLUCOSE_CELL_DIFF", 100, 0, 0,
 "Membrane diff constant",
 "Cell membrane diffusion coefficient Kd"},

{"GLUCOSE_BDRY_CONC", 9.0, 0, 0,
 "Boundary concentration",
 "GLUCOSE boundary concentration"},

{"GLUCOSE_CONSTANT", 0, 0, 1,
 "Constant concentration",
 "Extracellular concentration to be held constant everywhere at the specified boundary value"},

{"GLUCOSE_CONSUMPTION", 3.8e-17, 0, 0,
 "Max consumption rate",
 "GLUCOSE consumption rate"},

{"GLUCOSE_MM_KM", 1.33, 0, 0,
 "Michaelis-Menten Km",
 "Michaelis-Menten Km (uM)"},

{"GLUCOSE_HILL_N", 1, 1, 2,
 "Hill function N",
 "Glucose uptake rate Hill function N"},

{"USE_TRACER", 0, 0, 1,
"Use Tracer?",
"Tracer is simulated"},

{"TRACER_DIFF_COEF", 6.0e-7, 0, 0,
 "Spheroid diffusion coeff",
 "TRACER diffusion coefficient"},

{"TRACER_MEDIUM_DIFF", 6.0e-6, 0, 0,
 "Medium diffusion coeff",
 "Constituent diffusion coefficient in the medium"},

{"TRACER_CELL_DIFF", 20, 0, 0,
 "Membrane diff constant",
 "Cell membrane diffusion coefficient Kd"},

{"TRACER_BDRY_CONC", 1.0, 0, 0,
 "Boundary concentration",
 "TRACER boundary concentration"},

{"TRACER_CONSTANT", 1, 0, 1,
 "Constant concentration",
 "Extracellular concentration to be held constant everywhere at the specified boundary value"},

{"TRACER_CONSUMPTION", 0, 0, 0,
 "Consumption rate",
 "TRACER consumption rate"},

{"TRACER_MM_KM", 0, 0, 0,
 "Michaelis-Menten Km",
 "Michaelis-Menten Km (uM)"},

{"TRACER_HILL_N", 0, 0, 2,
 "Hill function N",
 "Tracer uptake rate Hill function N"},

//==========================
// Radiotherapy parameters
//==========================

{"RADIATION_ALPHA_H_1", 0.0473, 0, 0,
"Alpha (hypoxia)",
"alpha for irradiation of cells under anoxia (zero oxygen)"},

{"RADIATION_BETA_H_1", 0.0017, 0, 0,
"Beta (hypoxia)",
"beta for irradiation of cells under anoxia (zero oxygen)"},

{"RADIATION_OER_ALPHA_1", 2.5, 0, 0,
"OER alpha",
"Maximum oxygen enhancement ratio for alpha component of radiosensitivity "},

{"RADIATION_OER_BETA_1", 3.0, 0, 0,
"OER beta",
"Maximum oxygen enhancement ratio for beta component of radiosensitivity"},

{"RADIATION_KM_1", 4.3e-3, 0, 0,
"Km for radiosensitivity",
"Oxygen concentration for half maximal radiosensitivity relative to hypoxic cell exposure"},

{"RADIATION_DEATH_PROB_1", 1.0, 0, 0,
"Death prob",
"Probability of death at mitosis for a cell tagged for damage by radiation"},

{"RADIATION_ALPHA_H_2", 0.0473, 0, 0,
"Alpha (hypoxia)",
"alpha for irradiation of cells under anoxia (zero oxygen)"},

{"RADIATION_BETA_H_2", 0.0017, 0, 0,
"Beta (hypoxia)",
"beta for irradiation of cells under anoxia (zero oxygen)"},

{"RADIATION_OER_ALPHA_2", 2.5, 0, 0,
"OER alpha",
"Maximum oxygen enhancement ratio for alpha component of radiosensitivity "},

{"RADIATION_OER_BETA_2", 3.0, 0, 0,
"OER beta",
"Maximum oxygen enhancement ratio for beta component of radiosensitivity"},

{"RADIATION_KM_2", 4.3e-3, 0, 0,
"Km for radiosensitivity",
"Oxygen concentration for half maximal radiosensitivity relative to hypoxic cell exposure"},

{"RADIATION_DEATH_PROB_2", 1.0, 0, 0,
"Death prob",
"Probability of death at mitosis for a cell tagged for damage by radiation"},

{"HYPOXIA_1", 0.1, 0, 0,
"Hypoxia threshold 1",
"Hypoxia threshold 1"},

{"HYPOXIA_2", 1.0, 0, 0,
"Hypoxia threshold 2",
"Hypoxia threshold 2"},

{"HYPOXIA_3", 4.0, 0, 0,
"Hypoxia threshold 3",
"Hypoxia threshold 3"},

{"GROWTH_FRACTION_1", 0.25, 0, 0,
"Growth fraction threshold 1",
"Growth fraction threshold 1"},

{"GROWTH_FRACTION_2", 0.1, 0, 0,
"Growth fraction threshold 2",
"Growth fraction threshold 2"},

{"GROWTH_FRACTION_3", 0.01, 0, 0,
"Growth fraction threshold 3",
"Growth fraction threshold 3"},

{"SPCRAD", 200.0, 0, 0,
"Spectral radius",
"Spectral radius value used by RKC solver"},

{"USE_EXTRA", 0, 0, 1,
"Use extra conc",
"Use extracellular O2 concentration to determine cell death"},

{"USE_RELAX", 1, 0, 1,
"Use O2 relaxation solver",
"Use over- and under-relaxation to solve reaction-diffusion for oxygen"},

{"USE_PAR_RELAX", 1, 0, 1,
"Use parallel O2 relaxation solver",
"Use over- and under-relaxation to solve reaction-diffusion for oxygen, with parallelized over-relaxation"},

{"FD_SOLVER_1", 0, 0, 1,
"Use FD solver?",
"Use the FD solver in the far field"},

//{"USE_RADIATION", 1, 0, 1,
//"Use radiation?",
//"Treatment with radiation"},

//{"USE_TREATMENT_FILE", 0, 0, 1,
//"Use treatment file?",
//"Treatment programme is specified in the treatment file"},

//{"TREATMENT_FILE_NAME", 0, 0, 0,
//"treatment.data",
//"The treatment file contains data describing the drug and radiation dosing schedule"},

{"USE_DROP", 0, 0, 1,
"Account for drop deformation",
"Account for drop deformation when it is released to sit at the bottom of the well"},

{"NDROP", 1000, 0, 0,
"Dropping cell count",
"Number of cells in the spheroid when it is dropped."},

{"DROP_ALPHA", 0.4, 0, 0,
"Contact_diameter/diameter",
"Drop parameter alpha = initial (surface contact diameter)/(blob diameter).  Must be < 1."},

{"DROP_BETA", 0.6, 0, 0,
"Height/diameter",
"Drop parameter beta = initial (blob height)/(blob diameter).  Must be < 1."},

    {"SAVE_PROFILE_DATA",0,0,1,
     "Save profile data",
     "Save data for profile plots as a specified interval"},

    {"SAVE_PROFILE_DATA_FILE_NAME",0,0,0,
     "profile_data",
     "Base file name for saving profile data"},

    {"SAVE_PROFILE_DATA_INTERVAL",0,0,0,
     "Interval",
     "Time interval for saving profile data"},

    {"SAVE_PROFILE_DATA_NUMBER",1,0,0,
     "Number",
     "Number of times to save profile data"},

// This is the end of the parameters that are actually read by the DLL
// Entries after this point are QMyLabel dummies, to enable display of explanatory info  - no input data is transmitted,
// followed by the list of time-series and profile plots selected for this run.

{"DUMMY_HYPOXIA_THRESHOLD", 0, 0, 0,
"Hypoxia threshold",
"Select the intracellular O2 level below which the cell is counted as hypoxic"},

{"DUMMY_GROWTH_FRACTION", 0, 0, 0,
"Growth fraction",
"Select the threshold fraction of average growth rate (i.e. with no nutrient limits) used to count slow-growing cells"},

// Time-series plots
    {"nlive",                     1, 0,1,"","Number of live cells"},
    {"nanoxiadead",               1, 0,1,"","Total number of cells that have been killed by anoxia"},
    {"ndrugAdead",                1, 0,1,"","Total number of cells that have been killed by drugA"},
    {"ndrugBdead",                0, 0,1,"","Total number of cells that have been killed by drugB"},
    {"nradiationdead",            1, 0,1,"","Total number of cells that have been killed by radiation"},
    {"nanoxiatagged",             1, 0,1,"","Current number of cells tagged to die by anoxia"},
    {"ndrugAtagged",              0, 0,1,"","Current number of cells tagged to die by drugA"},
    {"ndrugBtagged",              0, 0,1,"","Current number of cells tagged to die by drugB"},
    {"nradiationtagged",          1, 0,1,"","Current number of cells tagged to die by radiation"},
    {"diameter",                  1, 0,1,"","Spheroid diameter (um)"},
    {"volume",                    0, 0,1,"","Spheroid volume (mm3)"},
    {"hypoxicfraction",           1, 0,1,"","Fraction of cells with oxygen level below the specified threshold for hypoxia"},
    {"growthfraction",            1, 0,1,"","Percentage of cells that are growing at a rate less than the specified fraction of the mean growth rate with no nutrient limits"},
    {"necroticfraction",          1, 0,1,"","Percentage of the spheroid that is necrotic = (number of vacant sites)/(number of sites taken up by the spheroid)"},
    {"platingefficiency",         0, 0,1,"","Percentage of live cells that are viable"},
    {"mediumoxygen",              0, 0,1,"","Average concentration of oxygen in the medium (far-field)"},
    {"mediumglucose",             0, 0,1,"","Average concentration of glucose in the medium (far-field)"},
    {"mediumdrugA",               0, 0,1,"","Average concentration of drug A in the medium (far-field)"},
    {"mediumdrugB",               0, 0,1,"","Average concentration of drug B in the medium (far-field)"},
    {"cellspermm3",               0, 0,1,"","Number of cells per mm3 in the blob"},

// Profile plots
    {"MULTI",                     1, 0,1,"","Selected constituent on a line through the blob centre"},
    {"CFSE",                      0, 0,1,"","Extracellular CFSE concentration on a line through the blob centre"},
    {"Oxygen",                    0, 0,1,"","Extracellular oxygen concentration on a line through the blob centre"},
    {"Glucose",                   0, 0,1,"","Extracellular glucose concentration on a line through the blob centre"},
    {"Tracer",                    0, 0,1,"","Extracellular tracer concentration on a line through the blob centre"},
    {"Drug_A",                    0, 0,1,"","Extracellular drug A concentration on a line through the blob centre"},
    {"Drug_A_metab1",             0, 0,1,"","Extracellular drug A metabolite 1 concentration on a line through the blob centre"},
    {"Drug_A_metab2",             0, 0,1,"","Extracellular drug A metabolite 2 concentration on a line through the blob centre"},
    {"Drug_B",                    0, 0,1,"","Extracellular drug Bconcentration on a line through the blob centre"},
    {"Drug_B_metab1",             0, 0,1,"","Extracellular drug B metabolite 1 concentration on a line through the blob centre"},
    {"Drug_B_metab2",             0, 0,1,"","Extracellular drug B metabolite 2 concentration on a line through the blob centre"},
    {"growthrate",                1, 0,1,"","Cell growth rate on a line through the blob centre"},
    {"cellvolume",                1, 0,1,"","Cell volume fraction on a line through the blob centre"},
    {"O2byvolume",                0, 0,1,"","Cell volume fraction on a line through the blob centre"},
// Distribution plots
//    {"Oxygen",                    0, 0,1,"","Probability distribution of extracellular oxygen concentration"},
//    {"cellvolume",                0, 0,1,"","Probability distribution of cell volume fraction"}

};
	nParams = sizeof(params)/sizeof(PARAM_SET);
	workingParameterList = new PARAM_SET[nParams];
	for (int i=0; i<nParams; i++) {
		workingParameterList[i] = params[i];
	}
}


PARAM_SET Params::get_param(int k)
{
	return workingParameterList[k];
}

void Params::set_value(int k, double v)
{
	workingParameterList[k].value = v;
}

void Params::set_label(int k, QString str)
{
	workingParameterList[k].label = str;
}
