#ifndef CONFIG_H
#define CONFIG_H

struct Parameters{
  double temperature_in_K;
  double pH_temperature;
  double x_Cr;
  double x_Fe;
  double x_Ni;
  double time_in_seconds;
  char *save_history;
  char *output_file_name;
  char *model;

  double alpha; // interface polarisability -between 0 and 1)
  double DV; // potential change with respect to a stationnary state
  double F_0f; // potential drop in the film
  double F_0mf; // potential drop at metal-film interface
  double F_0fs; // potential drop at film-solution interface
  /* double CtotM_mf=40.6e28; // total cation concentration in oxide at the metal-film interface  */
  /* double CmFe_mf =40.6e14; // Fe3+ concentration in oxide at the metal-film interface */
  /* double CMNi_mf =40.6e14; // Ni2+ concentration in oxide at the metal-film interface */
  /* double CVCr_mf =40.6e14; // Cr3+ vacancy concentration in oxide at the metal-film interface */
  /* double CtotM_fs=40.6e28; // total cation concentration at oxide surface (oxide-solution interface) */
  /* double CmFe_fs =40.6e7; // Fe3+ concentration at oxide surface (oxide-solution interface) */
  /* double CMNi_fs =40.6e7; // Ni2+ concentration at oxide surface (oxide-solution interface) */
  /* double CVCr_fs =40.6e20; // Cr3+ vacancy concentration at oxide surface (oxide-solution interface) */
  /* // lacunes d'oxygene */
  /* double q_vO=+2; // electronic charge */
  double D_mCr;  // diffusion coefficient of chromium ions
  double D_ICr;  // diffusion coefficient of chromium ions in interstices
  double D_vO; // diffusion coefficient of oxygen vacancies
  double DG1; // Gibbs enrgy of formation of reaction 1 (Cr_M -> Cr3+_ox + 3e- + 3/2 V_o with V_o = oxygen vacancy)
  double DG8; // Gibbs enrgy of formation of reaction 8 (V_o + H2O -> 2 H+ + O_ox) */
  double DG2; // Gibbs energy of formation of reaction 2 (Cr_M + V_Cr -> Cr3+_ox + 3e)
  double DG4; // Gibbs energy of formation of reaction 4 (Fe_M + V_Cr -> Fe3+_ox + 3e)
  double DG6; // Gibbs energy of formation of reaction 6 (Ni_M + 2/3 V_Cr -> Ni2+_ox + 2e-)
  double DG9; // Gibbs energy of formation of reaction 9 (Cr3+_ox + 3/2 H2O -> 3/2 O_ox + Cr3+_ox + V_Cr)
  double DG11; // Gibbs energy of formation of reaction 11 (Fe3+_ox -> V_Cr + Fe3+_aq)
  double DG13; // Gibbs energy of formation of reaction 13 (M_Ni -> 2/3 V_Cr + Ni2+_aq)
  double DG3; // Gibbs energy of formation of reaction 3 (Cr_M + V_I -> I_Cr + 3e)
  double DG5; // Gibbs energy of formation of reaction 5 (Fe_M + V_I -> I_Fe + 3e)
  double DG7; // Gibbs energy of formation of reaction 7 (Ni_M + V_I -> I_Ni + 2e)
  double DG10; // Gibbs energy of formation of reaction 10 (I_Cr + 3/2 H2O -> 3/2 O_ox + Cr3+_ox + 3H+ + V_I)
  double DG12; // Gibbs energy of formation of reaction 12 (I_Fe -> V_I + Fe3+_aq)
  double DG14; // Gibbs energy of formation of reaction 14 (I_Ni -> V_I + Ni2+_aq)

  /* // Cr intersticiels */
  /* double q_lCr=+3; // electronic charge of Cr3+ cation */
  /* double D_lCr=1e-21; // diffusion coefficient of Cr3+ cation */

  /* // Cr reseau */
  /* double q_mCr=+3; // */
  /* double D_mCr=1e-21; // */

  double CtotM_mf;
  double CtotI;

  double decay_length;
  double charge_number;

  double dissol_order;
  double dissol_preexp;
  double dissol_Ea;
// parameters of Voyshnis 2014 for H
  double gamma_H;
};

typedef struct Parameters Parameters;

void get_parameters();
void print_parameters();

#endif
