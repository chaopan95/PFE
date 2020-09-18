#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "c_solver_corrosion_evolution_functions.h"

/**
 * \file c_solver_corrosion_evolution_functions.c
 * \date March 2012
 * \version 0.5
 * \author B. Diawara, C. Toulemonde
 * \brief c_solver_corrosion_evolution C component that computes corrosion thickness layer growth.\n
 * \n
 */

//**********************
// Fonctions d'affichage
//**********************

//**********************************************
// Attends l'appui sur une touche pour continuer
//**********************************************

void pause_clavier()
{
  printf("Appuyez sur une touche pour continuer \n");
  getchar();
}

//*********************************************************
// Affiche la valeur des param�tres en cours d'optimisation
//*********************************************************
void affiche_param(double *x, int NP)
{
  printf("Parametres :");
  //affiche_vect_scientifique(x,NP);
  affiche_vect(x,NP);
  printf("]\n");
}

void affiche_vect(double *x, int n)
{
  int i;

  for (i=0;i<=n;i++){
    //printf("%.12f",x[i]);
    printf("%e ",x[i]);
    if (i<n) printf(", ");
  }
  printf("\n");
}

void affiche_vect_scientifique(double *x, int n)
{
  int i;
  for (i=1;i<=n/2;i++){
    printf("alpha%d=%e ",i,x[2*(i-1)+2]);
    printf("d%d=%e \n",i,x[2*(i-1)+1]);
  }
  printf("\n");
}

//***********************************************
// Fonctions d'allocation dynamique de la m�moire
//***********************************************

//******************************************
// Alloue un vecteur de r�els de type double
//*******************************************

double *dvector(int n)
{
  return (double *) malloc(n*sizeof(double));
}

//****************************
// Alloue un vecteur d'entiers
//****************************

int *ivector(int n)
{
  return (int *) malloc(n*sizeof(int));
}

//*******************************************
// Alloue une matrice de r�els de type double
//*******************************************
double **dmatrix(int row, int col)
{
  int i;
  double **m;

  m=(double **) malloc(row*sizeof(double*));
  if (m) {
    m[0]=(double *) malloc(col*row*sizeof(double));
    if (!m[0]) return 0L;
    for (i=1;i<row;i++) m[i]=m[i-1]+col;
  }
  return m;
}
