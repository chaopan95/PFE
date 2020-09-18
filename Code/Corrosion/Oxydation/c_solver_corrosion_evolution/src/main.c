#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "c_solver_corrosion_evolution_parameters.h"
#include <dictionary.h>
#include <iniparser.h>
#include <argtable2.h>
#include <print_doc_string.h>

#include "c_solver_corrosion_evolution_functions.h"
#include "equadiff.h"
#include "c_solver_corrosion_evolution_optimisation_functions.h"

/**
 * \file main.c
 * \brief c_solver_corrosion_evolution is a c MAP component that computes the evolution of oxyde thickness.
 * see MAP online documentation for a more complete information.
 * authors : B. Diawara, A. Seyeux, C. Toulemonde, K. Leistner
 * History :
 *       2010 : creation of the model by Antoine Seyeux (Laboratoire de Physico-chimie des surfaces)
 *       2011 : implementation of the model into a c code by Boubakar Diawara (Laboratoire de Physico-chimie des surfaces)
 *       2012 : integration into MAP by Charles Toulemonde, add new growth and dissolution model during Kirsten Leistner post-doc period
 *       2014 : prepare code for integration of Svetlana Voyshnis PhD thesis
 */


int main(int argc, char **argv){

  struct arg_file *inifile = arg_file0("f",NULL,"file.input",        "configuration file (default is \"-\")");
  struct arg_lit  *verbose = arg_lit0("v","verbose,debug",            "verbose messages");
  struct arg_lit  *help    = arg_lit0("h","help",                    "print this help and exit");
  struct arg_end  *end     = arg_end(20);
  void* argtable[] = {verbose,help,inifile,end};
  int nerrors;
  int exitcode=0;
  FILE *stream=NULL;
  char *code_name=(char*)"c_solver_corrosion_evolution";

  if(arg_nullcheck(argtable) != 0){
    printf("%s: insufficient memory\n",code_name);
    exitcode=1;
    goto exit;
  }

  inifile->filename[0]="-";
  nerrors = arg_parse(argc,argv,argtable);

  // print usage
  if(help->count > 0){
    print_doc_string((char*)"c_solver_corrosion_evolution");
    printf("Usage: %s", code_name);
    arg_print_syntax(stdout,argtable,"\n");
    arg_print_glossary(stdout,argtable,"  %-25s %s\n");
    exitcode=0;
    goto exit;
  }

  if(nerrors > 0){
    arg_print_errors(stdout,end,code_name);
    printf("Try '%s --help' for more information.\n",code_name);
    exitcode=1;
    goto exit;
  }

  // launch Bob Diawara first optimisation functions when no arguments are given
  if(argc==1){
    optimisation();
    exitcode=0;
    goto exit;
  }

  // log into file
  if((verbose->count==0)&&(argc>1)){
    if((stream = freopen("c_solver_corrosion_evolution.log", "a", stdout)) == NULL){
      fprintf(stderr, "unable to redirect stdout to file c_solver_corrosion_evolution.log\n");
    }
  }

  // variables
  double h;
  char *line, *option, *input_file_name, *string;
  int index;
  Parameters my_parameters;

  // read parameters
  line = (char*)malloc(1000*sizeof(char));
  option = (char*)malloc(1000*sizeof(char));
  input_file_name = (char*)malloc(1000*sizeof(char));
  string = (char*)malloc(1000*sizeof(char));
  line = "********************************************************************************";

  printf("%s\ncode : %s\n\nBEGIN\n\n", line, code_name);
  printf("%d arguments found for code '%s' :\n", argc-1, argv[0]);
  for (index=1; index<argc; index++) {
    printf("- argument(%d) = '%s'\n", index, argv[index]);
  }
  printf("\n");

  if(inifile->count>0){
    printf("- option '-f' has been found\n\n");
    printf("- reading input file %s :\n",*inifile->filename);
    get_parameters(&my_parameters,*inifile->filename);
    check_parameters(&my_parameters);

    print_parameters(&my_parameters);
    set_user_parameters(my_parameters);

    long NINTG = 36000;
    double t[NINTG+1], y[NINTG+1], J_v0[NINTG+1], J_MCr[NINTG+1], J_ICr[NINTG+1], J_H[NINTG+1];
    t[0] = 0; y[0] = 1e-12;
    h = my_parameters.time_in_seconds/(1.*NINTG);

    if(strstr(my_parameters.model, (char*)"Seyeux_2010")!=NULL){
      time_integration_RK2_Seyeux_2010(h, NINTG, t, y);
      if (strstr(my_parameters.save_history, "True")!=NULL){
        printf("%s", history_synthesis(NINTG, t, y, 3600));
        sauver( my_parameters.output_file_name, NINTG, t, y);
      } else {
      printf("y(t = %.6e) = \n%.6e\n", t[NINTG], y[NINTG]);
      }
    }
    if(strstr(my_parameters.model, (char*)"Leistner_2012")!=NULL){
      time_integration_RK2_Leistner_2012(h, NINTG, t, y, J_v0);
      write_potential( my_parameters.output_file_name, NINTG, y, t, my_parameters.decay_length, my_parameters.F_0mf, my_parameters.F_0f, my_parameters.F_0fs, my_parameters.alpha, my_parameters.DV);
      if (strstr(my_parameters.save_history, "True")!=NULL){
        printf("%s", history_synthesis(NINTG, t, y, 3600));
        write_output( my_parameters.output_file_name, NINTG, t, y, J_v0, J_H);
      } else {
        printf("y(t = %.6e) = \n%.6e\n", t[NINTG], y[NINTG], J_v0[NINTG]);
      }
      sensitivity_analysis(h, NINTG, t, y, J_v0, my_parameters.output_file_name);
    }
    if(strstr(my_parameters.model, (char*)"Voyshnis_2014")!=NULL){
      time_integration_RK2_Voyshnis_2014(h, NINTG, t, y, J_v0, J_MCr, J_ICr, J_H);
      write_potential( my_parameters.output_file_name, NINTG, y, t, my_parameters.decay_length, my_parameters.F_0mf, my_parameters.F_0f, my_parameters.F_0fs, my_parameters.alpha, my_parameters.DV);
      if (strstr(my_parameters.save_history, "True")!=NULL){
        printf("%s", history_synthesis(NINTG, t, y, 3600));
        write_output( my_parameters.output_file_name, NINTG, t, y, J_v0, J_H);
      } else {
        printf("y(t = %.6e) = \n%.6e\n", t[NINTG], y[NINTG], J_v0[NINTG]);
      }
      sensitivity_analysis(h, NINTG, t, y, J_v0, my_parameters.output_file_name);
    }
    exitcode=0;
    goto exit;
  }else{
    printf("Option '-f' has not been found !!!\nSTOP\n");
    exitcode=1;
    goto exit;
  }

  exit:
    if(stream!=NULL){
      fclose(stream);
      freopen ("/dev/tty", "w", stdout);
    }
    arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));
    return exitcode;
}
