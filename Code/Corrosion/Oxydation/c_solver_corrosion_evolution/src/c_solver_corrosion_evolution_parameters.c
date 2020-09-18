#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dictionary.h>
#include <iniparser.h>
#include <MAP_utils_parameters.h>
#include "c_solver_corrosion_evolution_parameters.h"

void get_parameters(Parameters *my_parameters, char *input_file_name){

  char *temp_string;
  dictionary *ini;
  temp_string = (char*)malloc(1000*sizeof(char));

  ini = iniparser_load(input_file_name);

  if(ini==NULL){
    printf("cannot parse file: %s\n\n", input_file_name);
    exit(0);
  }else{
    my_parameters->model = get_string_parameter(ini, "c_solver_corrosion_evolution:model");
    my_parameters->temperature_in_K = get_float_parameter(ini, "c_solver_corrosion_evolution:temperature_in_K");
    my_parameters->pH_temperature = get_float_parameter(ini, "c_solver_corrosion_evolution:pH_temperature");
    my_parameters->x_Cr = get_float_parameter(ini, "c_solver_corrosion_evolution:x_Cr");
    my_parameters->x_Fe = get_float_parameter(ini, "c_solver_corrosion_evolution:x_Fe");
    my_parameters->x_Ni = get_float_parameter(ini, "c_solver_corrosion_evolution:x_Ni");
    my_parameters->time_in_seconds = get_float_parameter(ini, "c_solver_corrosion_evolution:time_in_seconds");
    my_parameters->save_history = get_string_parameter(ini, "c_solver_corrosion_evolution:save_history");
    my_parameters->output_file_name = get_string_parameter(ini, "c_solver_corrosion_evolution:output_file_name");

    temp_string = iniparser_getstring(ini, "c_solver_corrosion_evolution:alpha", NULL);
    if(temp_string != NULL){
      my_parameters->alpha = atof(temp_string);
    } else { my_parameters->alpha = 0.5; }

    temp_string = iniparser_getstring(ini, "c_solver_corrosion_evolution:DV", NULL);
    if(temp_string != NULL){
      my_parameters->DV = atof(temp_string);
    } else { my_parameters->DV = 0.5; }

    temp_string = iniparser_getstring(ini, "c_solver_corrosion_evolution:F_0f", NULL);
    if(temp_string != NULL){
      my_parameters->F_0f = atof(temp_string);
    } else { my_parameters->F_0f = 10.; }

    temp_string = iniparser_getstring(ini, "c_solver_corrosion_evolution:F_0mf", NULL);
    if(temp_string != NULL){
      my_parameters->F_0mf = atof(temp_string);
    } else { my_parameters->F_0mf = 0.2; }

    temp_string = iniparser_getstring(ini, "c_solver_corrosion_evolution:F_0fs", NULL);
    if(temp_string != NULL){
      my_parameters->F_0fs = atof(temp_string);
    } else { my_parameters->F_0fs = 0.2; }

    temp_string = iniparser_getstring(ini, "c_solver_corrosion_evolution:DG8", NULL);
    if(temp_string != NULL){
      my_parameters->DG8 = atof(temp_string);
    } else { my_parameters->DG8 = -100000.; } 

    temp_string = iniparser_getstring(ini, "c_solver_corrosion_evolution:D_vO", NULL);
    if(temp_string != NULL){
      my_parameters->D_vO = atof(temp_string);
    } else { my_parameters->D_vO = 1.e-19; }

    temp_string = iniparser_getstring(ini, "c_solver_corrosion_evolution:D_mCr", NULL);
    if(temp_string != NULL){
      my_parameters->D_mCr = atof(temp_string);
    } else { my_parameters->D_mCr = 0.; }

    temp_string = iniparser_getstring(ini, "c_solver_corrosion_evolution:D_ICr", NULL);
    if(temp_string != NULL){
      my_parameters->D_ICr = atof(temp_string);
    } else { my_parameters->D_ICr = 0.; }

    temp_string = iniparser_getstring(ini, "c_solver_corrosion_evolution:DG1", NULL);
    if(temp_string != NULL){
      my_parameters->DG1 = atof(temp_string);
    } else { my_parameters->DG1 = -100000.; }

    temp_string = iniparser_getstring(ini, "c_solver_corrosion_evolution:DV", NULL);
    if(temp_string != NULL){
      my_parameters->DV = atof(temp_string);
    } else { my_parameters->DV = 0.5; }

    temp_string = iniparser_getstring(ini, "c_solver_corrosion_evolution:decay_length", NULL);
    if(temp_string != NULL){
      my_parameters->decay_length = atof(temp_string);
    } else { my_parameters->decay_length = 1.; }

    temp_string = iniparser_getstring(ini, "c_solver_corrosion_evolution:charge_number", NULL);
    if(temp_string != NULL){
      my_parameters->charge_number = atof(temp_string);
    } else { my_parameters->charge_number = 3.; }

    temp_string = iniparser_getstring(ini, "c_solver_corrosion_evolution:dissol_order", NULL);
    if(temp_string != NULL){
      my_parameters->dissol_order = atof(temp_string);
    } else { my_parameters->dissol_order = 1.; }

    temp_string = iniparser_getstring(ini, "c_solver_corrosion_evolution:dissol_preexp", NULL);
    if(temp_string != NULL){
      my_parameters->dissol_preexp = atof(temp_string);
    } else { my_parameters->dissol_preexp = 0.; }

    temp_string = iniparser_getstring(ini, "c_solver_corrosion_evolution:dissol_Ea", NULL);
    if(temp_string != NULL){
      my_parameters->dissol_Ea = atof(temp_string);
    } else { my_parameters->dissol_Ea = 100000.; }

    temp_string = iniparser_getstring(ini, "c_solver_corrosion_evolution:DG2", NULL);
    if(temp_string != NULL){
      my_parameters->DG2 = atof(temp_string);
    } else { my_parameters->DG2 = -100000.; }

    temp_string = iniparser_getstring(ini, "c_solver_corrosion_evolution:DG4", NULL);
    if(temp_string != NULL){
      my_parameters->DG4 = atof(temp_string);
    } else { my_parameters->DG4 = -100000.; }

    temp_string = iniparser_getstring(ini, "c_solver_corrosion_evolution:DG6", NULL);
    if(temp_string != NULL){
      my_parameters->DG6 = atof(temp_string);
    } else { my_parameters->DG6 = -100000.; }

    temp_string = iniparser_getstring(ini, "c_solver_corrosion_evolution:DG9", NULL);
    if(temp_string != NULL){
      my_parameters->DG9 = atof(temp_string);
    } else { my_parameters->DG9 = -120000.; }

    temp_string = iniparser_getstring(ini, "c_solver_corrosion_evolution:DG11", NULL);
    if(temp_string != NULL){
      my_parameters->DG11 = atof(temp_string);
    } else { my_parameters->DG11 = -100000.; }

    temp_string = iniparser_getstring(ini, "c_solver_corrosion_evolution:DG13", NULL);
    if(temp_string != NULL){
      my_parameters->DG13 = atof(temp_string);
    } else { my_parameters->DG13 = -100000.; }

    temp_string = iniparser_getstring(ini, "c_solver_corrosion_evolution:DG10", NULL);
    if(temp_string != NULL){
      my_parameters->DG10 = atof(temp_string);
    } else { my_parameters->DG10 = -100000.; }

    temp_string = iniparser_getstring(ini, "c_solver_corrosion_evolution:DG12", NULL);
    if(temp_string != NULL){
      my_parameters->DG12 = atof(temp_string);
    } else { my_parameters->DG12 = -100000.; }

    temp_string = iniparser_getstring(ini, "c_solver_corrosion_evolution:DG14", NULL);
    if(temp_string != NULL){
      my_parameters->DG14 = atof(temp_string);
    } else { my_parameters->DG14 = -100000.; }

    temp_string = iniparser_getstring(ini, "c_solver_corrosion_evolution:DG3", NULL);
    if(temp_string != NULL){
      my_parameters->DG3 = atof(temp_string);
    } else { my_parameters->DG3 = -100000.; }

    temp_string = iniparser_getstring(ini, "c_solver_corrosion_evolution:DG5", NULL);
    if(temp_string != NULL){
      my_parameters->DG5 = atof(temp_string);
    } else { my_parameters->DG5 = -100000.; }

    temp_string = iniparser_getstring(ini, "c_solver_corrosion_evolution:DG7", NULL);
    if(temp_string != NULL){
      my_parameters->DG7 = atof(temp_string);
    } else { my_parameters->DG7 = -100000.; }

    temp_string = iniparser_getstring(ini, "c_solver_corrosion_evolution:CtotM_mf", NULL);
    if(temp_string != NULL){
      my_parameters->CtotM_mf = atof(temp_string);
    } else { my_parameters->CtotM_mf = 4.14e28; }

    temp_string = iniparser_getstring(ini, "c_solver_corrosion_evolution:CtotI", NULL);
    if(temp_string != NULL){
      my_parameters->CtotI = atof(temp_string);
    } else { my_parameters->CtotI = 2.07e28; }

//  parameters of Voyshnis 2014 for H
    temp_string = iniparser_getstring(ini,
"c_solver_corrosion_evolution:gamma_H", NULL);
    if(temp_string != NULL){
      my_parameters->gamma_H = atof(temp_string);
    } else { my_parameters->gamma_H = 1.; }
  }
  free(ini);
  return;
}

void check_parameters(const Parameters *my_parameters){

 int error=0;

 if(my_parameters->temperature_in_K <= 0.){
   printf("Parameter temperature_in_K must be > 0 \n");
   error++;
 }
 if(my_parameters->pH_temperature <= 0.){
   printf("Parameter pH_temperature must be > 0 \n");
   error++;
 }
 if(my_parameters->x_Cr <= 0.){
   printf("Parameter x_Cr must be > 0 \n");
   error++;
 }
 if(my_parameters->time_in_seconds <= 0.){
   printf("Parameter time_in_seconds must be > 0 \n");
   error++;
// parameters of Voyshnis for H
 }
 if(my_parameters->gamma_H <= 0.){
   printf("Parameter gamma_H must be [0.,1.] \n");
   error++;
 }
 if((strstr(my_parameters->model, (char*)"Seyeux_2010")==NULL)&&(strstr(my_parameters->model, (char*)"Leistner_2012")==NULL)&&(strstr(my_parameters->model, (char*)"Voyshnis_2014")==NULL)){
   printf("model does not belong to [Seyeux_2010, Leistner_2012, Voyshnis_2014]\n");
   error++;
 }
 if(error != 0) exit(0);
}

void print_parameters(const Parameters *my_parameters){
  printf("model = %s\n", my_parameters->model);
  printf("temperature_in_K = %f\n", my_parameters->temperature_in_K);
  printf("pH_temperature = %f\n", my_parameters->pH_temperature);
  printf("x_Cr = %f\n", my_parameters->x_Cr);
  printf("x_Fe = %f\n", my_parameters->x_Fe);
  printf("x_Ni = %f\n", my_parameters->x_Ni);
  printf("time_in_seconds = %f\n", my_parameters->time_in_seconds);
  printf("save_history = %s\n", my_parameters->save_history);
  printf("output_file_name = %s\n", my_parameters->output_file_name);
  printf("alpha = %f\n", my_parameters->alpha);
  printf("DV = %f\n", my_parameters->DV);
  printf("F_0f = %f\n", my_parameters->F_0f);
  printf("F_0mf = %f\n", my_parameters->F_0mf);
  printf("F_0fs = %f\n", my_parameters->F_0fs);
  printf("D_vO = %.3e\n", my_parameters->D_vO);
  printf("D_mCr = %.3e\n", my_parameters->D_mCr);
  printf("D_ICr = %.3e\n", my_parameters->D_ICr);
  printf("DG1 = %f\n", my_parameters->DG1);
  printf("DG8 = %f\n", my_parameters->DG8);
  printf("DG2 = %f\n", my_parameters->DG2);
  printf("DG4 = %f\n", my_parameters->DG4);
  printf("DG6 = %f\n", my_parameters->DG6);
  printf("DG9 = %f\n", my_parameters->DG9);
  printf("DG11 = %f\n", my_parameters->DG11);
  printf("DG13 = %f\n", my_parameters->DG13);
  printf("DG10 = %f\n", my_parameters->DG10);
  printf("DG12 = %f\n", my_parameters->DG12);
  printf("DG14 = %f\n", my_parameters->DG14);
  printf("DG3 = %f\n", my_parameters->DG3);
  printf("DG5 = %f\n", my_parameters->DG5);
  printf("DG7 = %f\n", my_parameters->DG7);
  printf("CtotM_mf = %.3e\n", my_parameters->CtotM_mf);
  printf("CtotI = %.3e\n", my_parameters->CtotI);
  printf("decay_length = %.3e\n", my_parameters->decay_length);
  printf("charge number = %.3e\n", my_parameters->charge_number);
  printf("order, dissolution = %f\n", my_parameters->dissol_order);
  printf("preexponential factor, dissolution = %.3e\n", my_parameters->dissol_preexp);
  printf("activation energy, dissolution = %.3e\n", my_parameters->dissol_Ea);
  printf("gamma_H = %.3e\n", my_parameters->gamma_H);
  printf("\n");
}
