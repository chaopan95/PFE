void optimisation();
double f1(int i, double y[]);
double erreur(long NINTG, double h, double *x, double t[], double y[]);
double gradient(int NP, double H, long NINTG, double h,
		double *x, double *dx,
		double t[], double y[]);
double gradient_conjugue(int NP, int ALGO, int MAXITER,
			 double H, double GTOL, double ETOL,
			 long NINTG,
			 double *x, double *dx, double *ddx,
			 double t[], double y[]);
double gradient_conjugue_cense(void);
double descendre(int NP, long NINTG, double h, double *x, double *dx,
		 double t[], double y[]);
void progresser(int NP , double r, double *x, double *dx);
