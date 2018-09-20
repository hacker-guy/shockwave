/***************************************************************************
 *
 *   File        : tasks.h
 *   Student Id  : <758397>
 *   Name        : <JUSTIN BUGEJA>
 *
 ***************************************************************************/

#ifndef TASKS_H

typedef struct q2_eqpar_t q2_eqpar_t;
typedef struct q3_matrixpar_t q3_matrixpar_t;

/* Question 2 */
void shockwave(const char* q2_file);

double newton_raph(double b, double mval, double theta);
double shock_der(double B, double m, double theta);
double shock_funct(double B, double m, double theta);

double rad2deg(double rad);
double deg2rad(double deg);
double cot(double theta);

/* Question 3 */
void linalgbsys(const char* q3_file);

void gauss_sol(q3_matrixpar_t* matrix, int n, double *heateq, int Nt, int t);

/* Question 5 */
void interp(const char* q5_file, const double xo);
int cmpfnc (const void * a, const void * b);

/* Question 6 */
void heateqn(const char* q6_file);

double explicit_int(double* fn, int Nt, int Nx, double deltat, double deltax,
		int t, int x, double mu);
double explicit_int_ve(double* fn, int Nt, int Nx, double deltat, double deltax,
		int t, int x, double mu);
void implicit_int(double* fn, int Nt, int Nx, double deltat, double deltax,
		int t, int x, double mu);
double RHS(double* fn, int x, int t, int Nx, int Nt, double deltax, double mu);

#endif
