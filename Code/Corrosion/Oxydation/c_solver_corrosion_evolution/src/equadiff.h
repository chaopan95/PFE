#include "c_solver_corrosion_evolution_parameters.h"

double fit_exp(long t);
void time_integration_RK2_Seyeux_2010(double h,int n,double t[],double y[]);
void time_integration_RK2_Leistner_2012(double h,int n,double t[],double y[],double J_v0[]);
void time_integration_RK2_Voyshnis_2014(double h,int n,double t[],double y[],double J_v0[], double J_MCr[], double J_ICr[], double J_H[]);
char* history_synthesis(int n, double t[], double y[], int step);
void affiche(int n,double t[],double y[]);
void sauver(char *nom,int n,double t[],double y[]);
void write_output(char *nom,int n,double t[],double y[],double J_v0[], double J_H[]);

double f_( double x);
double f_Seyeux_2010(double t, double y);
double f_Leistner_2012(double x_0, double y);
double approx_sol(double x_0, double y);

void set_user_parameters(Parameters my_parameters);
void recopie_param_xi(double *x, int *opt);
void recopie_xi_param(double *x);

double J_vO_Seyeux_2010( double x); // oxygen vacancy flux
double J_vO_Leistner_2012( double x); // oxygen vacancy flux
double J_MCr_Leistner_2012( double x); // Chrome ion flux through metallic positions
double J_ICr_Leistner_2012( double x); // Chrome ion flux through interstitial positions
double J_dissol_Leistner_2012( double x); // film dissolution rate
double dC_H_total(double J_ICr, double J_vO, double dt); // total concentration of H

void sensitivity_analysis(double h,int n,double t[],double y[],double J_v0[],char *nom);
void load_params(double *param);
