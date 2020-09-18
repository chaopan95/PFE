//*******************************************
// Boubakar Diawara
// Laboratoire de Physico-chimie des surfaces
//*******************************************

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "equadiff.h"
#include "c_solver_corrosion_evolution_parameters.h"

// I - Constants
//double Omega=1; // ??
//double Nv=1; // Avogadro number (must be equal to 6.022e23 ?)
double Omega=2.92e-3; // oxide molar volume in m³/mol
double Nv=6.022e23; // Avogadro number (must be equal to 6.022e23 ?)
double F=9.65e4; // Faraday constant
double R=8.314472; // Perfect gas constant in J/mol/K
double e_charge=1.602e-19; // electronic charge in C
double kB=1.38e-23; // Boltzmann's constant in J/K
double rho_oxide=5.22; // oxide density in g/cm³
double M_oxide=151.9904; // oxide molar mass in g/mol
double stoich_Cr=2.; // Cr stoichiometry in oxide

// II - Variables
// 1) end-user level
double T=603.; // Temperature in degree K
double pH=7.2; // water chemistry
double x_Cr=0.3; // Cr mole fraction
double x_Fe=0.3; // Fe mole fraction
double x_Ni=0.3; // Ni mole fraction
// 2) expert-user level
double alpha=0.5; // interface polarisability -between 0 and 1)
double DV=0.5; // potential change with respect to a stationnary state 
double F_0f=10.; // potential drop in the film
double F_0mf=0.2; // potential drop at metal-film interface
double F_0fs=0.2; // potential drop at film-solution interface
double CtotM_mf=4.14e28; // total cation concentration in oxide at the metal-film interface 
double CmFe_mf =40.6e14; // Fe3+ concentration in oxide at the metal-film interface
double CMNi_mf =40.6e14; // Ni2+ concentration in oxide at the metal-film interface
double CtotI =2.07e28; // Interstice concentration in Cr2O3 oxide
double CVCr_mf =40.6e14; // Cr3+ vacancy concentration in oxide at the metal-film interface
double CtotM_fs=4.14e28; // total cation concentration at oxide surface (oxide-solution interface)
double CmFe_fs =40.6e7; // Fe3+ concentration at oxide surface (oxide-solution interface)
double CMNi_fs =40.6e7; // Ni2+ concentration at oxide surface (oxide-solution interface)
double CVCr_fs =40.6e20; // Cr3+ vacancy concentration at oxide surface (oxide-solution interface)
// lacunes d'oxygene
double q_vO=2.; // electronic charge of oxygen vacancy
double D_vO=1e-19; // diffusion coefficient of oxygen vacancies
// J_mCr
double DG1=100000.; // Gibbs energy of formation of reaction 1 (Cr_M -> Cr3+_ox + 3e- + 3/2 V_o with V_o = oxygen vacancy)
double DG8=-100000.; // Gibbs energy of formation of reaction 8 (V_o + H2O -> 2 H+ + O_ox)
double DG2=100000.; // Gibbs energy of formation of reaction 2 (Cr_M + V_Cr -> Cr3+_ox + 3e)
double DG4=-100000.; // Gibbs energy of formation of reaction 4 (Fe_M + V_Cr -> Fe3+_ox + 3e)
double DG6=100000.; // Gibbs energy of formation of reaction 6 (Ni_M + 2/3 V_Cr -> Ni2+_ox + 2e-)
double DG9=-100000.; // Gibbs energy of formation of reaction 9 (Cr3+_ox + 3/2 H2O -> 3/2 O_ox + Cr3+_ox + V_Cr)
double DG11=100000.; // Gibbs energy of formation of reaction 11 (Fe3+_ox -> V_Cr + Fe3+_aq)
double DG13=-100000.; // Gibbs energy of formation of reaction 14 (M_Ni -> 2/3 V_Cr + Ni2+_aq)
// J_ICr
double DG3=-100000.; // Gibbs energy of formation of reaction 3 (Cr_M + V_I -> I_Cr + 3e)
double DG5=-100000.; // Gibbs energy of formation of reaction 5 (Fe_M + V_I -> I_Fe + 3e)
double DG7=-100000.; // Gibbs energy of formation of reaction 7 (Ni_M + V_I -> I_Ni + 2e)
double DG10=-100000.; // Gibbs energy of formation of reaction 10 (I_Cr + 3/2 H2O -> 3/2 O_ox + Cr3+_ox + 3H+ + V_I)
double DG12=-100000.; // Gibbs energy of formation of reaction 12 (I_Fe -> V_I + Fe3+_aq)
double DG14=-100000.; // Gibbs energy of formation of reaction 14 (I_Ni -> V_I + Ni2+_aq)
// Cr intersticiels
double q_lCr=3.; // electronic charge of Cr3+ cation
double D_lCr=1e-21; // diffusion coefficient of Cr3+ cation
double D_ICr=1.e-21; // diffusion coefficient of Cr3+ cation
// Cr reseau
double q_mCr=3.; //
double D_mCr; //
double decay_length=2.e-7; //length caracterising the influence zone of the potential
//dissolution
double  charge_number = 3.; //number of electrons transferred during dissolution reaction
double  dissol_order = 1.;
double  dissol_preexp = 0.;
double  dissol_Ea = 100000.;
// parametrs of Voyshnis for H
double gamma_H = 1.;

int i;

void set_user_parameters(Parameters my_parameters){
  double ratio;
  double K1, P1, P2, P3, P4, P5, P6, P, Q, Sdiss, Rdiss;
  double C_VCr_mf, C_MFe_mf, C_MNi_mf, C_MCr_mf;
  double C_VI_mf, C_IFe_mf, C_INi_mf, C_ICr_mf;
  double eterm_1, eterm_2, eterm_3;
  double H_plus;

  T = my_parameters.temperature_in_K;
  pH = my_parameters.pH_temperature;
  x_Cr = my_parameters.x_Cr;
  x_Fe = my_parameters.x_Fe;
  x_Ni = my_parameters.x_Ni;

  alpha = my_parameters.alpha;
  DV = my_parameters.DV;
  F_0f = my_parameters.F_0f;
  F_0mf = my_parameters.F_0mf;
  F_0fs = my_parameters.F_0fs;

  D_vO = my_parameters.D_vO;
  DG1 = my_parameters.DG1;
  DG8 = my_parameters.DG8;
  DG2 = my_parameters.DG2;
  DG4 = my_parameters.DG4;
  DG6 = my_parameters.DG6;
  DG9 = my_parameters.DG9;
  DG11 = my_parameters.DG11;
  DG13 = my_parameters.DG13;
  DG10 = my_parameters.DG10;
  DG12 = my_parameters.DG12;
  DG14 = my_parameters.DG14;
  DG3 = my_parameters.DG3;
  DG5 = my_parameters.DG5;
  DG7 = my_parameters.DG7;

  D_mCr = my_parameters.D_mCr;
  D_ICr = my_parameters.D_ICr;

  CtotM_mf = my_parameters.CtotM_mf;
  CtotI = my_parameters.CtotI;

  decay_length = my_parameters.decay_length;

  charge_number = my_parameters.charge_number;

  dissol_order = my_parameters.dissol_order;
  dissol_preexp = my_parameters.dissol_preexp;
  dissol_Ea = my_parameters.dissol_Ea;
// parameters of Voyshnis for H
  gamma_H = my_parameters.gamma_H;
  ratio = (F_0f + DV*(1.-alpha)) / F_0f ;
  printf("ratio caracterising the square root shape of the curve (when close to 1.) = %f\n", ratio);

  
/******************************** CALCULATE LUMPED PARAMETERS **********************************/
  /* JvO */
  K1 = (q_vO*F*D_vO)/(R*T) * exp(((-1.*2./3.*DG1) + (2.*F*F_0mf))/(R*T) + 2./3.*log(x_Cr)) ;
  //K1 = (q_vO*F*D_vO)/R/T ;
  P1 = K1 * (F_0f+DV);
  P2 = -1. *K1 * DV * alpha ;

  /* JmCr */
  C_VCr_mf = exp((DG2-3.*F*F_0mf)/(R*T) - log(x_Cr));
  C_MFe_mf = exp((-1.*DG4+3.*F*F_0mf)/(R*T)) * x_Fe * C_VCr_mf;
  C_MNi_mf = exp((-1.*DG6+2.*F*F_0mf)/(R*T)) * x_Ni * pow(C_VCr_mf,(2./3.));
  C_MCr_mf = CtotM_mf*Omega/Nv - (1.*C_MNi_mf + 1.*C_MFe_mf + C_VCr_mf);
  P3 = ((q_mCr*F*D_mCr)/(R*T)) * C_MCr_mf * (F_0f+DV);
  P4 = -1.*((q_mCr*F*D_mCr)/(R*T)) * C_MCr_mf * alpha * DV;

  /* JICr */
  eterm_1 = x_Cr*exp((-1.*DG3+3.*F*F_0mf)/R/T);
  eterm_2 = x_Fe*exp((-1.*DG5+3.*F*F_0mf)/R/T);
  eterm_3 = x_Ni*exp((-1.*DG7+2.*F*F_0mf)/R/T);
  C_VI_mf = CtotI*Omega/Nv / (1. + eterm_1 + eterm_2 + eterm_3);
  C_ICr_mf = eterm_1 * C_VI_mf;
  P5 = ((q_mCr*F*D_ICr)/(R*T)) * C_ICr_mf * (F_0f+DV);
  P6 = -1.*((q_mCr*F*D_ICr)/(R*T)) * C_ICr_mf * alpha * DV;

  P = P1 + P3 + P5;
  Q = P2 + P4 + P6;

  /* Jdissol */
  H_plus = exp(-1.*2.303*pH);
  Sdiss = charge_number*e_charge/kB/T * alpha * DV;
  Rdiss = dissol_preexp * exp((-1.*dissol_Ea/kB/Nv/T + charge_number*e_charge*F_0fs/kB/T));


  printf("K1 = %.12e\n", K1);
  printf("P1 = %.12e\n", P1);
  printf("P2 = %.12e\n", P2);
  printf("P3 = %.12e\n", P3);
  printf("P4 = %.12e\n", P4);
  printf("P5 = %.12e\n", P5);
  printf("P6 = %.12e\n", P6);
  printf("P = %.12e\n", P);
  printf("Q = %.12e\n", Q);
  printf("Rdiss = %.12e\n", Rdiss);
  printf("Sdiss = %.12e\n", Sdiss);

  return;
}

void recopie_param_xi(double *x, int *opt){
  // Si opt[i]= 1 x[i] est optimisé
  // Si opt[i]= 0 x[i] n'est pas optimisé

  // lacunes d'oxygene
  x[0]=F_0fs;		opt[0]=1;
  x[1]=F_0f;		opt[1]=1;
  x[2]=F_0mf;		opt[2]=1;
  x[3]=alpha;		opt[3]=1;

  x[4]=DG1;		opt[4]=1;
  x[5]=D_vO;		opt[5]=1;
  x[6]=DG8; 		opt[6]=1;

  x[7]=q_vO;		opt[7]=0;

  x[8]=CmFe_mf;	opt[8]=0;
  x[9]=CMNi_mf;	opt[9]=0;
  x[10]=CVCr_mf;	opt[10]=0;
  x[11]=CtotM_fs;	opt[11]=0;
  x[12]=CmFe_fs;	opt[12]=0;
  x[13]=CMNi_fs;	opt[13]=0;
  x[14]=CVCr_fs;	opt[14]=0;
  // Cr intersticiels
  x[15]= q_lCr;	opt[15]=0;
  x[16]= D_lCr;	opt[16]=0;
  x[21]=CtotM_mf;	opt[17]=0;

  // Cr reseau
  x[22]= q_mCr;	opt[18]=0;
  x[23]= D_mCr;	opt[19]=0;
}

/************************************/
void recopie_xi_param (double *x){
  F_0fs = x[0];
  F_0f  = x[1];
  F_0mf = x[2];

  // lacunes d'oxygene
  alpha= x[3];
  DG1 = x[4];
  D_vO = x[5];
  DG8  = x[6] ;
  q_vO = x[7] ;

  CmFe_mf  = x[8];
  CMNi_mf  = x[9];
  CVCr_mf  = x[10];
  CtotM_fs = x[11];
  CmFe_fs  = x[12];
  CMNi_fs  = x[13];
  CVCr_fs  = x[14];

  // Cr intersticiels
  q_lCr    = x[15];
  D_lCr    = x[16];
  CtotM_mf = x[17];

  // Cr reseau
  q_mCr    = x[18];
  decay_length = x[19];
  DV       = x[20];
/*
  x_Cr     = x[21];
  x_Fe     = x[22];
  x_Ni     = x[23];
  T        = x[24];
  pH       = x[25];
  DG2      = x[26];
  DG4      = x[27];
  DG6      = x[28];
  DG9      = x[29];
  DG11     = x[30];
  DG13     = x[31];
  DG14     = x[32];
  D_mCr    = x[33];
  dissol_order  = x[34];
  dissol_preexp = x[35];
  dissol_Ea     = x[36];
  DG10     = x[37];
  DG12     = x[38];
  DG3      = x[39];
  DG5      = x[40];
  DG7      = x[41];
  D_ICr    = x[42];*/
}

double f_Seyeux_2010(double t, double x){
  double val;

  val = (Omega/Nv) * (J_vO_Seyeux_2010(x));
  return val;
}

double f_Leistner_2012(double x_0, double x){
  double val, f_approx;
  double j_diss;
  int num_intervals;
  // test on x_0, initial value of thickness because x could be to big because
  // of RK2 integration scheme
  if (x_0 > 1e-10) {
    j_diss = J_dissol_Leistner_2012(x);
  } else {
    j_diss = 0.;
  }
  val = (Omega/Nv) * (J_vO_Leistner_2012(x) + J_MCr_Leistner_2012(x) + J_ICr_Leistner_2012(x) - j_diss);

  /* comparison with approximated solution with lumped parameters */
  f_approx = approx_sol(x_0,x);
  /*i = i+1;
  num_intervals = 36000;
  if((num_intervals/10)%i==0) {
  printf("val=%e, approx=%e\n", val, f_approx);
  }*/
 /*****************************************************************/

  if (val < 0.) {
    printf("val < 0. => x_0=%e, x=%e, JvO=%e, JmCr=%e, Jdiss=%e\n", x_0, x, J_vO_Leistner_2012(x), J_MCr_Leistner_2012(x), j_diss);
  }
  return val;
}

double approx_sol(double x_0, double x){
  double f_approx;
  double K1, P1, P2, P3, P4, P5, P6, ftest, P, Q, Rdiss, Sdiss, H_plus, j_diss;
  double C_VCr_mf, C_MFe_mf, C_MNi_mf, C_MCr_mf;
  double C_VI_mf, C_IFe_mf, C_INi_mf, C_ICr_mf;
  double eterm_1, eterm_2, eterm_3;

/* JvO */
  K1 = (q_vO*F*D_vO)/(R*T) * exp(((-1.*2./3.*DG1) + (2.*F*F_0mf))/(R*T) + 2./3.*log(x_Cr)) ;
  P1 = K1 * (F_0f+DV);
  P2 = -1. *K1 * DV * alpha ;

  /* JmCr */
  C_VCr_mf = exp((DG2-3.*F*F_0mf)/(R*T) - log(x_Cr));
  C_MFe_mf = exp((-1.*DG4+3.*F*F_0mf)/(R*T)) * x_Fe * C_VCr_mf;
  C_MNi_mf = exp((-1.*DG6+2.*F*F_0mf)/(R*T)) * x_Ni * pow(C_VCr_mf,(2./3.));
  C_MCr_mf = CtotM_mf*Omega/Nv - (1.*C_MNi_mf + 1.*C_MFe_mf + C_VCr_mf);
  P3 = ((q_mCr*F*D_mCr)/(R*T)) * C_MCr_mf * (F_0f+DV);
  P4 = -1.*((q_mCr*F*D_mCr)/(R*T)) * C_MCr_mf * alpha * DV;

  /* JICr */
  eterm_1 = x_Cr*exp((-1.*DG3+3.*F*F_0mf)/R/T);
  eterm_2 = x_Fe*exp((-1.*DG5+3.*F*F_0mf)/R/T);
  eterm_3 = x_Ni*exp((-1.*DG7+2.*F*F_0mf)/R/T);
  C_VI_mf = CtotI*Omega/Nv / (1. + eterm_1 + eterm_2 + eterm_3);
  C_ICr_mf = eterm_1 * C_VI_mf;
  P5 = ((q_mCr*F*D_ICr)/(R*T)) * C_ICr_mf * (F_0f+DV);
  P6 = -1.*((q_mCr*F*D_ICr)/(R*T)) * C_ICr_mf * alpha * DV;

  P = P1 + P3 + P5;
  Q = P2 + P4 + P6;

  /* Jdissol */
  H_plus = exp(-1.*2.303*pH);
  Sdiss = charge_number*e_charge/kB/T * alpha * DV;
  Rdiss = dissol_preexp * H_plus * exp((-1.*dissol_Ea/kB/Nv/T + charge_number*e_charge*F_0fs/kB/T));
  if (x_0 > 1e-10) {
    j_diss =  (Omega/Nv) * Rdiss*exp(Sdiss*f_(x));
  } else {
    j_diss = 0.;
  }

  f_approx = (P + Q*f_(x))/x - j_diss ;

  return f_approx;
}

double J_dissol_Leistner_2012(double x){
  double H_plus, F_fs;
  double val, val2;

  H_plus = exp(-2.303*pH);
  F_fs = F_0fs + alpha * DV * exp(-1.*(x/decay_length));

  val = dissol_preexp * exp(-1.*(dissol_Ea/kB/Nv/T - charge_number*e_charge*Nv*F_fs/kB/Nv/T)) * pow(H_plus,dissol_order); //1/m2/s
  return val;
}

double J_vO_Seyeux_2010(double x){
  double val;
  double term_1, term_2, term_3, term_4, term_5;

  /* initial expression of val
  val=(Nv/Omega)*q_vO*(F*(F_0f+DV*(1-alpha*f_(x)))/(R*T*x))*D_vO*
    (exp(((-2./3)*DG1+2*F*F_0mf)/(R*T) +1.5*log(x_Cr)) *exp(q_vO*F*(F_0f+DV*(1-alpha*f_(x)))/(R*T))
     -exp(((DG8-2*F*(F_0fs+DV*alpha*f_(x)))/(R*T))-4.606*pH))/
    (exp(q_vO*F*(F_0f+DV*(1-alpha*f_(x)))/(R*T))-1);
  */

  // new expression reviewed with Antoine Seyeux

  term_1 = (q_vO*F) * (F_0f + DV*(1.-alpha*f_(x))) / (R*T*x);

  term_2 = exp(((-2./3.)*DG1+2.*F*F_0mf)/(R*T) +1.5*log(x_Cr)) ;

  term_3 = exp(q_vO*F*(F_0f+DV*(1.-alpha*f_(x))) / (R*T) );

  term_4 = -exp(((DG8-2.*F*(F_0fs+DV*alpha*f_(x))) / (R*T) )-4.606*pH);

  term_5 = exp(q_vO*F*(F_0f+DV*(1.-alpha*f_(x))) / (R*T) ) - 1.;

  val = (Nv/Omega)* term_1 * D_vO * (term_2 * term_3- term_4)/term_5 ;

  return val;
}

double J_vO_Leistner_2012(double x){
  double result;
  double term_1, term_2, term_3, term_4, term_5;
  double ratio;

  term_1 = (q_vO*F) * (F_0f + DV*((1.-alpha*f_(x)))) / (R*T*x);
  term_2 = exp(((-2./3.)*DG1+2.*F*F_0mf)/(R*T) +(2./3.)*log(x_Cr)) ;
  term_3 = exp(q_vO*F*(F_0f+DV*((1.-alpha*f_(x)))) / (R*T) );
  term_4 = exp(((DG8-2.*F*(F_0fs+DV*alpha*f_(x))) / (R*T) )-4.606*pH);
  term_5 = exp(q_vO*F*(F_0f+DV*((1.-alpha*f_(x)))) / (R*T) ) - 1.;

  result = (Nv/Omega)* term_1 * D_vO * (term_2 * term_3- term_4)/term_5 ;

  return result;
}

double J_MCr_Leistner_2012(double x){
  double val;
  double term_1, term_2, term_3, term_4, term_5;
  double C_VCr_mf, C_MFe_mf, C_MNi_mf, C_MCr_mf;
  double C_VCr_fs, C_MFe_fs, C_MNi_fs, C_MCr_fs, F_fs, C_Fe3aq, C_Ni2aq;
  int num_intervals;

  /* calculate concentrations at m/f interface */
  C_VCr_mf = exp((DG2-3.*F*F_0mf)/(R*T) - log(x_Cr));
  C_MFe_mf = exp((-1.*DG4+3.*F*F_0mf)/(R*T)) * x_Fe * C_VCr_mf;
  C_MNi_mf = exp((-1.*DG6+2.*F*F_0mf)/(R*T)) * x_Ni * pow(C_VCr_mf,(2./3.));
  C_MCr_mf = CtotM_mf*Omega/Nv - (1.*C_MNi_mf + 1.*C_MFe_mf + C_VCr_mf);

  /* calculate concentrations at f/s interface */
  F_fs = F_0fs + alpha * DV * f_(x);
  C_VCr_fs = exp((-1.*DG9+3.*F*F_fs)/(R*T));

  C_Ni2aq = rho_oxide*stoich_Cr/M_oxide/1.e-6*x_Ni/x_Cr; // [mol/m^3]
  C_Fe3aq = rho_oxide*stoich_Cr/M_oxide/1.e-6*x_Fe/x_Cr; // [mol/m^3]
  C_Ni2aq = C_Ni2aq*Omega; // [_]
  C_Fe3aq = C_Fe3aq*Omega; // [_]

  C_MFe_fs = exp((DG11-3.*F*F_fs)/(R*T)) * C_VCr_fs * C_Fe3aq;
  C_MNi_fs = exp((DG13-2.*F*F_fs)/(R*T)) * pow(C_VCr_fs,(2./3.)) * C_Ni2aq;
  C_MCr_fs = CtotM_fs*Omega/Nv - (1.*C_MNi_fs + 1.*C_MFe_fs + C_VCr_fs);

 /******  printout: aq concentrations *****/
 /* i = i+1;
  num_intervals = 36000;
  if((num_intervals/10)%i==0) {
    printf("C_Ni2aq=%f, C_Fe3aq=%f\n", C_Ni2aq, C_Fe3aq);
  }*/


  /****** printouts *****/
  /*i = i+1;
  num_intervals = 36000;
  if((num_intervals/10)%i==0) {
    printf("i=%i\n", i);
    printf("C_VCr_fs/C_VCr_mf=%e\n", C_VCr_fs/C_VCr_mf);
    printf("C_MCr_mf/C_MCr_fs=%e\n", C_MCr_mf/C_MCr_fs);
    printf("C_vO_mf/C_vO_fs=%e\n", (exp(((-2./3.)*DG1+2.*F*F_0mf)/(R*T) +(2./3.)*log(x_Cr))) / (exp(((DG8-2.*F*(F_0fs+DV*alpha*f_(x))) / (R*T) )-4.606*pH)));
    printf("exp(q_mCr*F*phi_f(x) / (R*T) )=%e\n", exp(q_mCr*F*(F_0f+DV*((1.-alpha*f_(x)))) / (R*T) ));
    printf("exp(q_vO*F*phi_f(x) / (R*T) )=%e\n", exp(q_vO*F*(F_0f+DV*((1.-alpha*f_(x)))) / (R*T) ));
    printf("C_vO_mf=%e\n", (exp(((-2./3.)*DG1+2.*F*F_0mf)/(R*T) +(2./3.)*log(x_Cr))));
    printf("C_vO_fs=%e\n", (exp(((DG8-2.*F*(F_0fs+DV*alpha*f_(x))) / (R*T) )-4.606*pH)));
    printf("C_VCr_mf=%e\n", C_VCr_mf);
    printf("C_VCr_fs=%e\n", C_VCr_fs);
    printf("C_MCr_mf=%e\n", C_MCr_mf);
    printf("CtotM_mf*Omega/Nv=%e, C_VCr_mf=%e, C_MNi_mf=%e, C_MFe_mf=%e\n", CtotM_mf*Omega/Nv, C_VCr_mf, C_MNi_mf, C_MFe_mf);
    printf("C_MCr_fs=%e\n", C_MCr_fs);
    printf("C_Ni_aq=%e\n", C_Ni2aq);
    printf("C_Fe_aq=%e\n", C_Fe3aq);
    printf("CtotM_fs*Omega/Nv=%e, C_VCr_fs=%e, C_MNi_fs=%e, C_MFe_fs=%e\n", CtotM_fs*Omega/Nv, C_VCr_fs, C_MNi_fs, C_MFe_fs);
    printf("-1*DG9=%e, 3*F*F_fs=%e, R*T=%e\n", -1*DG9, 3*F*F_fs, R*T);
    printf("\n");
    printf("i=%i\n", i);
    printf("Ea/F_fs=%e, F_fs=%e, Ea/kB/Nv/T=%e, 3*e*Nv/kB/Nv/T=%e, exp()=%e\n", dissol_Ea/F_fs, F_fs, dissol_Ea/kB/Nv/T, charge_number*e_charge*Nv/kB/Nv/T, exp(-1*(dissol_Ea/kB/Nv/T - charge_number*e_charge*Nv*F_fs/kB/Nv/T)));
    printf("\n");
  }*/

  /* calculate flux across film */
  term_1 = (q_mCr*F) * (F_0f + DV*((1.-alpha*f_(x)))) / (R*T*x);

  term_3 = exp(q_mCr*F*(F_0f+DV*((1.-alpha*f_(x)))) / (R*T) );

  term_5 = exp(q_mCr*F*(F_0f+DV*((1.-alpha*f_(x)))) / (R*T) ) - 1.;

  val = (Nv/Omega)* term_1 * D_mCr * (C_MCr_mf * term_3 - C_MCr_fs)/term_5 ;

 /* printf("(Nv/Omega)* term_1 * D_vO=%e\n", (Nv/Omega)* term_1 * D_vO);
  printf("term3=%e\n", term_3);
  printf("term5=%e\n", term_5);*/

  return val;
}


double J_ICr_Leistner_2012(double x){
  double val;
  double C_VI_mf, C_IFe_mf, C_INi_mf, C_ICr_mf;
  double C_VI_fs, C_IFe_fs, C_INi_fs, C_ICr_fs, F_fs, C_Fe3aq, C_Ni2aq;
  double term_1, term_2, term_3, term_4, term_5;
  double eterm_1, eterm_2, eterm_3, eterm_4, eterm_5, eterm_6;
  int num_intervals;

  /* calculate concentrations at m/f interface */
  eterm_1 = x_Cr*exp((-1.*DG3+3.*F*F_0mf)/R/T);
  eterm_2 = x_Fe*exp((-1.*DG5+3.*F*F_0mf)/R/T);
  eterm_3 = x_Ni*exp((-1.*DG7+2.*F*F_0mf)/R/T);

  C_VI_mf = CtotI*Omega/Nv / (1. + eterm_1 + eterm_2 + eterm_3);
  C_ICr_mf = eterm_1 * C_VI_mf;
  C_IFe_mf = eterm_2 * C_VI_mf;
  C_INi_mf = eterm_3 * C_VI_mf;

  /* calculate concentrations at f/s interface */
  F_fs = F_0fs + alpha * DV * f_(x);

  C_Ni2aq = rho_oxide*stoich_Cr/M_oxide/1.e-6*x_Ni/x_Cr; // [mol/m^3]
  C_Fe3aq = rho_oxide*stoich_Cr/M_oxide/1.e-6*x_Fe/x_Cr; // [mol/m^3]
  C_Ni2aq = C_Ni2aq*Omega; // [_]
  C_Fe3aq = C_Fe3aq*Omega; // [_]

  eterm_4 = exp((DG10-3.*F*F_fs)/R/T - 6.909*pH);
  eterm_5 = C_Fe3aq*exp((DG12-3.*F*F_fs)/R/T);
  eterm_6 = C_Ni2aq*exp((DG14-2.*F*F_fs)/R/T);

  C_VI_fs = CtotI*Omega/Nv / (1. + eterm_4 + eterm_5 + eterm_6);
  C_ICr_fs = eterm_4 * C_VI_fs;
  C_IFe_fs = eterm_5 * C_VI_fs;
  C_INi_fs = eterm_6 * C_VI_fs;

  /* calculate flux across film */
  term_1 = (q_mCr*F) * (F_0f + DV*((1.-alpha*f_(x)))) / (R*T*x);

  term_3 = exp(q_mCr*F*(F_0f+DV*((1.-alpha*f_(x)))) / (R*T) );

  term_5 = exp(q_mCr*F*(F_0f+DV*((1.-alpha*f_(x)))) / (R*T) ) - 1.;

  val = (Nv/Omega)* term_1 * D_ICr * (C_ICr_mf * term_3 - C_ICr_fs)/term_5 ;

  /****** printouts *****/
  /*i = i+1;
  num_intervals = 36000;
  if((num_intervals/10)%i==0) {
    printf("i=%i\n", i);
    printf("C_ICr_mf/C_ICr_fs=%e\n", C_ICr_mf/C_ICr_fs);
    printf("exp(q_mCr*F*phi_f(x) / (R*T) )=%e\n", exp(q_mCr*F*(F_0f+DV*((1.-alpha*f_(x)))) / (R*T) ));
    printf("exp(q_mCr*F*phi_f(x) / (R*T) )=%e\n", exp(q_mCr*F*(F_0f+DV*((1.-alpha*f_(x)))) / (R*T) ));
    printf("exp(q_vO*F*phi_f(x) / (R*T) )=%e\n", exp(q_vO*F*(F_0f+DV*((1.-alpha*f_(x)))) / (R*T) ));
    printf("C_ICr_mf=%e\n", C_ICr_mf);
    printf("CtotI*Omega/Nv=%e, eterm_1=%e, eterm_2=%e, eterm_3=%e\n", CtotI*Omega/Nv, eterm_1, eterm_2, eterm_3);
    printf("C_ICr_fs=%e\n", C_ICr_fs);
    printf("CtotI*Omega/Nv=%e, eterm_4=%e, eterm_5=%e, eterm_6=%e\n", CtotI*Omega/Nv, eterm_4, eterm_5, eterm_6);
  }*/

  return val;
}

void load_params(double* param){
  param[0]=F_0fs;
  param[1]=F_0f;
  param[2]=F_0mf;
  param[3]=alpha;
  param[4]=DG1;
  param[5]=D_vO;
  param[6]=DG8;
  param[7]=q_vO;
  param[8]=CmFe_mf;
  param[9]=CMNi_mf;
  param[10]=CVCr_mf;
  param[11]=CtotM_fs;
  param[12]=CmFe_fs;
  param[13]=CMNi_fs;
  param[14]=CVCr_fs;
  param[15]= q_lCr;
  param[16]= D_lCr;
  param[17]=CtotM_mf;
  param[18]= q_mCr;
  param[19]= decay_length;
  param[20]= DV;
  /*  param[21]= x_Cr;
  param[22]= x_Fe;
  param[23]= x_Ni;
  param[24]= T;
  param[25]= pH;
  param[26]=DG2;  
  param[27]=DG4;  
  param[28]=DG6;  
  param[29]=DG9;  
  param[30]=DG11;  
  param[31]=DG13;  
  param[32]=DG14;  
  param[33]=D_mCr;
  param[34]=dissol_order;  
  param[35]=dissol_preexp;  
  param[36]=dissol_Ea;  
*/
}

void sensitivity_analysis(double h,int n,double t[],double y[],double J_v0[],char* nom){
  int nn= 20; // Number of parameters - 1
  int i, j, k;
  double param[nn];
  double y_old[n];
  double y_new[n];
  double sens[n];
  char param_names[nn];
  char *sensitivity;
  char temp_name[100] = "";

  load_params(param);

  time_integration_RK2_Leistner_2012(h, n, t, y, J_v0);
  for(i=0;i<=n; i++){
  y_old[i] = y[i];
  }

  sensitivity=strncat(temp_name,nom,strlen(nom)-7);
  sensitivity = strcat(temp_name,".sensitivity_analysis");
  FILE *out;
  out=fopen(sensitivity,"w");
  for(i=0;i<=nn; i++){
     param[i] = param[i]*1.1;
     recopie_xi_param (param); //increase the value of parameter i by 10%

      time_integration_RK2_Leistner_2012(h, n, t, y, J_v0);

      param[i] = param[i]/1.1; //set all parameters back to original values
      fprintf(out, "%i", i);
      fprintf(out, "\t");
      for(j=0;j<=n; j++){
          y_new[j] = y[j];
          //sensitivity of film width with respect to parameter i, at time j:
          sens[j] = (param[i]/y_old[j]) * (y_new[j] - y_old[j])/(param[i]*1.1 - param[i]);
          if (j%360==0){
              fprintf(out, "%.3e", sens[j]);
              fprintf(out, "\t");
          }
      }
      fprintf(out, "\n");
  }
  fclose(out);
}

double f_(double x){
  double ratio;

  ratio = x/decay_length;
  return exp(- ratio);
}

double fit_exp(long t){
  double val,tt;

  tt=t/3600.0; //conversion en heures
  val=13.6*log(0.079*tt+1);

  return(val*1e-9); //conversion en Angstrom
}

void time_integration_RK2_Seyeux_2010(double h, int n, double t[], double y[]){
  double k1;
  int i;
  double yi_plus_one_half;

  //printf("time integration with RK2 method\n");
  for(i=0;i<=n; i++){
    k1 = h*f_Seyeux_2010(t[i], y[i]);
    // first estimation of intermediate time step value
    yi_plus_one_half = y[i] + k1/2;
    y[i+1] = y[i] + h*(f_Seyeux_2010(t[i]+h/2, yi_plus_one_half));
    t[i+1] = t[i] + h;
    if (yi_plus_one_half > y[i+1]){
      y[i+1] = y[i] + h*(f_Seyeux_2010(t[i]+h/2, y[i]));
      // second estimation of intermediate time step value
      yi_plus_one_half = 0.5 * (y[i] + y[i+1]);
      y[i+1] = y[i] + h*(f_Seyeux_2010(t[i]+h/2, yi_plus_one_half));
    }
  }
}

void time_integration_RK2_Leistner_2012(double h, int n, double t[], double y[], double J_v0[]){
  double k1;
  int i;
  double yi_plus_one_half;

  //printf("time integration with RK2 method\n");
  J_v0[0]=J_vO_Leistner_2012(y[0]);
  for(i=0;i<=n; i++){
    k1 = h*f_Leistner_2012(y[i], y[i]);
    // first estimation of intermediate time step value
    yi_plus_one_half = y[i] + k1/2;
    y[i+1] = y[i] + h*(f_Leistner_2012(y[i], yi_plus_one_half));
    t[i+1] = t[i] + h;
    if (yi_plus_one_half > y[i+1]){
      y[i+1] = y[i] + h*(f_Leistner_2012(y[i], y[i]));
      // second estimation of intermediate time step value
      yi_plus_one_half = 0.5 * (y[i] + y[i+1]);
      y[i+1] = y[i] + h*(f_Leistner_2012(y[i], yi_plus_one_half));
    }
    J_v0[i+1] = J_vO_Leistner_2012(y[i+1]);
  }
}

void time_integration_RK2_Voyshnis_2014(double h, int n, double t[], double y[], double J_v0[], double J_MCr[], double J_ICr[], double J_H[]){
  double k1;
  int i;
  double yi_plus_one_half;

  printf("time integration with RK2 method of Voyshnis 2014 model\n");
  J_v0[0] = J_vO_Leistner_2012(y[0]);
  J_ICr[0] = J_ICr_Leistner_2012(y[0]);
  J_MCr[0] = J_MCr_Leistner_2012(y[0]);
  J_H[0] = 0.;

  for(i = 0 ; i <= n ; i++){
    k1 = h * f_Leistner_2012(y[i], y[i]);
    // first estimation of intermediate time step value
    yi_plus_one_half = y[i] + k1/2;
    y[i+1] = y[i] + h * (f_Leistner_2012(y[i], yi_plus_one_half));
    t[i+1] = t[i] + h;

    if (yi_plus_one_half > y[i+1]){
      y[i+1] = y[i] + h * (f_Leistner_2012(y[i], y[i]));
      // second estimation of intermediate time step value
      yi_plus_one_half = 0.5 * (y[i] + y[i+1]);
      y[i+1] = y[i] + h * (f_Leistner_2012(y[i], yi_plus_one_half));
    }

    J_v0[i+1] = 2. * J_vO_Leistner_2012(y[i+1]);
    J_ICr[i+1] = 3. * J_ICr_Leistner_2012(y[i+1]);
    J_MCr[i+1] = 3. * J_MCr_Leistner_2012(y[i+1]);
    J_H[i+1] = gamma_H * (J_v0[i+1] + J_ICr[i+1] + J_MCr[i+1]);
/*
    C_H[i+1] = gamma_H*(C_H[i] + dC_H_total(J_ICr, J_v0[i+1], h));
    J_H[i+1] = (C_H[i+1] - C_H[i]) / h;
*/
    printf("t = %f, J_H = %e\n", t[i+1], J_H[i+1]);
  }
}

char* history_synthesis(int n, double t[], double y[], int step){
  long i;
  char *result, *string;

  string = (char *)malloc(sizeof(char)*(1000));
  result = (char *)malloc(sizeof(char)*(100000));
  strcpy(result, "t[s], X_calc[m], X_exp[m]\n");
  for(i=0;i<=n; i++){
    if (i%step==0){
      sprintf(string, "%.12e, %.12e, %.12e\n", t[i],y[i],fit_exp(t[i]));
      strcat(result, string);
    }
  }
  free(string);

  return result;
}

void affiche(int n, double t[], double y[]){
  long i;

  printf("t[s], X_calc[m], X_exp[m]\n");
  for(i=0;i<=n; i++){
    if (i%3600==0){
      printf("%.12e, %.12e, %.12e\n", t[i],y[i],fit_exp(t[i]));
    }
  }
  printf("\n");
}

//*******************************************************************
// Sauvegarde le résultat dans un fichier texte les valeurs calculées
// avec les valeurs attendues d'aprÃšs le fit experimental
//*******************************************************************
void sauver(char *nom,int n, double t[],double y[]){
  long i;
  FILE *out;

  out=fopen(nom,"w");

  fprintf(out, "t[s], X_calc[m], X_exp[m]\n");
  for(i=0;i<=n; i++){
    if (i%3600==0){
      fprintf(out, "%.12e, %.12e, %.12e\n", t[i],y[i],fit_exp(t[i]));
    }
  }
  fprintf(out, "\n");
  fclose(out);
}

void write_output(char *nom,int n,double t[],double y[],double J_v0[], double J_H[]){
  long i;
  FILE *out;

/*************************************************************/
/* write output for non regression tests */
  out=fopen(nom,"w");

  fprintf(out, "t[s], X_calc[m], X_exp[m], J_v0[1/m2/s], J_H[1/m2/s]\n");
  for(i=0;i<=n; i++){
    if (i%360==0){
      fprintf(out, "%.12e, %.12e, %.12e, %.12e, %.12e\n", t[i],y[i],fit_exp(t[i]), J_v0[i], J_H[i]);
    }
  }
  fprintf(out, "\n");
  fclose(out);

/*************************************************************/
/* write output for graphing */
  char *graph_fig;
  char temp_name[100] = "";

  graph_fig=strncat(temp_name,nom,strlen(nom)+7);
  graph_fig = strcat(temp_name,".graphs");
  
  out=fopen(graph_fig,"w");

  fprintf(out, "t[s], X_calc[m], X_exp[m], J_vO[1/m2/s], J_vO_x[1/m/s], f(X_calc)], J_dissol[1/m2/s], J_vO-J_dissol, J_MCr, J_ICr, J_H[/m2.s)]\n");
  for(i=0;i<=n; i++){
    if (i%360==0){
      fprintf(out, "%.12e, %.12e, %.12e, %.12e, %.12e, %.12e, %.12e, %.12e, %.12e, %.12e, %.12e\n", t[i], y[i], fit_exp(t[i]), J_vO_Leistner_2012(y[i]), y[i]*J_v0[i], f_(y[i]), J_dissol_Leistner_2012(y[i]), J_vO_Leistner_2012(y[i])-J_dissol_Leistner_2012(y[i]), J_MCr_Leistner_2012(y[i]), J_ICr_Leistner_2012(y[i]), J_H[i]);
      printf("J_H[i] = %f\n", J_H[i]);
    }
  }
  fprintf(out, "\n");
  fclose(out);
}

void write_potential(char *nom,int n,double y[],double t[],double d,double F_0mf,double F_0f,double F_0fs,double alpha,double DV){
  long i;
  char *potential;
  char temp_name[100] = "";
  double DF_f, DF_fs, V, F_f, F_fs, V0;

  potential=strncat(temp_name,nom,strlen(nom)-7);
  potential = strcat(temp_name,".potential");
  FILE *out;
  out=fopen(potential,"w");

  fprintf(out, "$t$[s], $\\phi^{\\circ}_{f/s}(d_{b})$, $\\phi_{f}(d_{b})$, $\\phi_{f/s}(x)$, $\\phi_{f}(x)$, $V(d_b)$, $V(x)$\n");
  for(i=0;i<=n; i++){
    if (i%360==0){
      DF_fs = alpha * DV * exp(-1.*(y[i]/d));
      DF_f = DV *(1.-alpha *  exp(-1.*(y[i]/d)));
      F_fs = F_0fs + alpha * DV * exp(-1.*(y[i]/d));
      F_f = F_0f + DV *(1.-alpha *  exp(-1.*(y[i]/d)));
      V = F_0mf + (F_0f + DF_f) + (F_0fs + DF_fs);
      V0 = F_0mf + F_0f + F_0fs;
      fprintf(out, "%.12e, %.12e, %.12e, %.12e, %.12e, %.12e, %.12e\n", t[i], F_0fs, F_0f, F_fs, F_f, V0, V);
    }
  }
  fprintf(out, "\n");
  fclose(out);
}
