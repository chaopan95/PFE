#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "c_solver_corrosion_evolution_functions.h"
#include "equadiff.h"
#include "c_solver_corrosion_evolution_optimisation_functions.h"

void optimisation(){
  int nn= 19;   // Nombre total de paramètres dans le modèle de croissance
  int NP = 6;   // Les NP+1 premiers parametres sont optimisees
  int ALGO = 2; // ALGO=1 : Plus grande pente
                // ALGO=2 : Fletcher-Reeves
                // ALGO=3 : Polak-Ribiere
  char *algo[] = {"","Plus grande pente","Fletcher-Reeves","Polak-Ribiere"};
  int MAXITER = 10;    // Nombre maximum d'itérations
  double H = 1e-4;     // Pas des différences finies
  double GTOL = 1e-20; // Tolérance sur le gradient
  double ETOL = 1e-19; // Tolérance sur l'erreur
  int nb_heures = 10;  // Intervalle de temps de la simulation
  long NINTG = 3600*(nb_heures); // Nombre de points d'intégration en s

  // variables
  int *opt;
  double h;
  double *x, *z;
  double *dx, *ddx; // Vecteur gradient et auxiliaire
  double t[3600*(nb_heures+1)], y[3600*(nb_heures+1)];

  // print parameters
  printf("Nombre de paramètres           = %d\n", NP);
  printf("Algorithme de descente         = %s\n", algo[ALGO]);
  printf("Nombre maximum d'itérations    = %d\n", MAXITER);
  printf("Pas des différences finies     = %.1e\n", H);
  printf("Tolérance sur le gradient      = %.1e\n", GTOL);
  printf("Tolérance sur l'erreur         = %.1e\n", ETOL);
  printf("Nombre de points d'intégration = %d\n", (int)(NINTG));
  printf("\n");

  // allocation des tableaux utilisés pour la résolution
  opt = ivector(nn+1);
  x = dvector(nn+1);
  z = dvector(nn+1);
  dx = dvector(nn+1);
  ddx = dvector(nn+1);

  //*********************************************************************
  // On commence par résoudre l'equation diff. avec les valeurs de départ
  // des paramètres et on sauvegarde le resultat dans un fichier
  // texte avec les valeurs attendues d'après le fit experimental
  //*********************************************************************
  t[0] = 0; y[0] = 1e-19; h = 1;
  time_integration_RK2_Seyeux_2010(h, NINTG, t, y);
  affiche( NINTG, t, y);
  sauver("c_solver_corrosion_evolution_test_1.output", NINTG, t, y);

  //**************************************************************
  // On calcule et on affiche la valeur de l'erreur et du gradient
  // pour les valeurs initiales   des paramètres
  //**************************************************************
  recopie_param_xi(x, opt);
  affiche_param(x, NP);
  /* test de erreur() */
  printf("Valeurs de l'erreur avec les paramètres de départ e=%e\n",erreur(NINTG, h, x, t, y));
  /* test de gradient() */
  printf("Valeurs du gradient avec les paramètres de départ g=%e\n",gradient(H, NP, NINTG, h, x, dx, t, y));
  affiche_vect(dx,NP);

  //****************************************************************
  // On optimise les paramètres avec la méthode du gradient conjugué
  // et on sauvegarde le resultat dans un fichier  texte avec les
  // valeurs attendues d'après le fit experimental
  //****************************************************************
  gradient_conjugue(NP, ALGO, MAXITER, H, GTOL, ETOL, NINTG, x, dx, ddx, t, y);
  affiche_param(x, NP);
  sauver("c_solver_corrosion_evolution_optimisation_1.output", NINTG, t, y);
}

//*************************************************
// Optimisation par la méthode du gradient conjugué
//*************************************************
double gradient_conjugue(int NP, int ALGO, int MAXITER,
			 double H, double GTOL, double ETOL,
			 long NINTG,
			 double *x, double *dx, double *ddx,
			 double t[], double y[]){
  int i, n;
  double e, g, h=1., b;
  double e_old;

  printf("debut gradient_conjugue\n");
  g = gradient(NP, H, NINTG, h, x, dx, t, y); // La première direction, c'est le gradient [dEx, dEy, dEz], g est le carré de sa norme
  for (n=1;n<=MAXITER;n++){
    e = descendre(NP, NINTG, h, x, dx, t, y); // Récupère le minimum d'énergie dans la direction [dEx, dEy, dEz]
    for (i=0;i<=NP;i++) {// Copie l'ancienne direction [dEx, dEy, dEz] dans [dGx, dGy, dGz]
      ddx[i] = dx[i];
    }

    h = gradient(NP, H, NINTG, h, x, dx, t, y); // Calcule le nouveau gradient [dEx, dEy, dEz] et le carré de sa norme */
    printf("\nn = %5d, e = %.10e, g = %e\n",n,e,sqrt(h/(NP))); // Affiche [itération, énergie, gradient]
    affiche_param(x, NP);
    if (sqrt(h/(NP))<GTOL) break; // Compare le gradient actuel au gradient terminal souhaité (GTOL)
    if(n>1) if(fabs(e_old-e)<ETOL){printf("Energie ne varie plus\n");	break;} // Pour Fletcher-Reeves on a tout : h = dE*dE

    if (ALGO==3) {// Pour Polak-Ribiere il faut calculer h = (dE-dG)*dG
      h=0;
      for (i=0;i<=NP;i++) h+= (dx[i]-ddx[i])*dx[i];
    }
    b=h/g; // Facteur de pondération entre ancienne direction et nouveau gradient
    g=0; // Il reste à calculer la nouvelle direction et le carré de sa norme

    for (i=0;i<=NP;i++) {
      dx[i]+=b*ddx[i];
      g+=dx[i]*dx[i];
    }
    e_old=e;
  }
  printf("fin gradient_conjugue\n");

  return e;
}

//*****************************
// Calcul numérique du gradient
//*****************************
double gradient(int NP, double H, long NINTG, double h,
		double *x, double *dx,
		double t[], double y[]){
  int i;
  double I1,I2;
  double g;
  //  double t;
  double t_temp;

  for (i=0;i<=NP;i++)
    {
      t_temp=x[i];
      x[i]=t_temp*(1+H);
      I1=erreur(NINTG, h, x, t, y);

      x[i]=t_temp*(1-H);
      I2=erreur(NINTG, h, x, t, y);

      x[i]=t_temp;
      dx[i]=(I1-I2)/(2*H);
    }
  g=0;
  for (i=0;i<=NP;i++) g+=dx[i]*dx[i];

  return(g);
}

double erreur(long NINTG, double h, double *x,
	      double t[], double y[]){
  int i;
  double I;

  recopie_xi_param(x);
  t[0]=0; y[0]=1e-19; h=1; //k=1;
  time_integration_RK2_Seyeux_2010(h, NINTG, t, y);

  // Calcul de l'integrale par la méthode des trapèzes
  I=f1(1, y)+f1(NINTG, y);
  for(i=2;i<NINTG; i++){
    for(i=3600;i<NINTG; i+=3600){
      {
        I=I+(2+2*(i%2))*f1(i, y);
      }
    }
  }
  I=(h/3)*I;
  return(I);
}

double f1(int i, double y[]){
  return((y[i]-fit_exp(i))*(y[i]-fit_exp(i))*1e9);
}

double descendre(int NP, long NINTG, double h, double *x, double *dx,
		 double t[], double y[]){
  double R0 = 1e-6; // pas de départ de la descente
  double e, f, r;

  e = erreur(NINTG, h, x, t, y); // La direction de descente est déjà dans [dEx, dEy, dEz]
  r=R0; // Recherche r, le pas de départ, en commençant à R0
  f=e;
  // Tant que l'énergie augmente, il faut réduire le pas
  do{
    progresser(NP, r, x, dx); // Avance
    e = erreur(NINTG, h, x, t, y);
    if (e>f){
        progresser(NP, -r, x, dx); // Recule
        r/=2; // Divise le pas par 2
    }
  } while (e>f);
  if (r<R0) printf("r<R0 (%.1e) : %e\n",R0,r); // Informe que le pas initial est inférieur à R0

  // Le pas de départ est trouvé
  // Tant que l'énergie baisse on avance
  do{
    f = e;
    progresser(NP, r, x, dx); // Avance
    e = erreur(NINTG, h, x, t, y);
    //printf("descente normal  e=%.15e  r=%e\n",e,r);
    //if (e<f) r*=R; // Si l'énergie baisse on augmente le pas
    if (e<f) r*=2; // Si l'énergie baisse on augmente le pas
  } while (e<f);
  progresser(NP, -r, x, dx); // On est allé trop loin, il faut reculer

  return f; // Retourne la meilleure énergie (c'est f, pas e)
}

void progresser(int NP , double r, double *x, double *dx){
  int i;

  for (i=0;i<=NP;i++){
    x[i]-=r*dx[i]*x[i];
  }
}
