/***************************************************************************
 *
 *   File        : tasks.c
 *   Student Id  : <758397>
 *   Name        : <JUSTIN BUGEJA>
 *
 ***************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <string.h>
#include <assert.h>
#include "tasks.h"

#define MAXITR 100
#define LARGE_NUM 2147483646
#define LAGNUM 3
#define BUFF 5
#define NVALS 1
#define EPSILON 1.0e-12
#define MAXDEG 90
#define Y 1.4
#define XMAX 1.0
#define TMAX 2.0
#define NTPRINT 99
#define M_PI 3.14159265359
#define OUTPUTFILE1 "out_shock.csv"
#define OUTPUTFILE2 "out_linalsys.csv"
#define OUTPUTFILE3 "out_interp.csv"
#define OUTPUTFILE4 "out_heateqn_explicit_fe.csv"
#define OUTPUTFILE5 "out_heateqn_explicit_ve.csv"
#define OUTPUTFILE6 "out_heateqn_implicit_fe.csv"

/* structure to find which three points to use for lagrange */
typedef struct {

	double x;
	double fx;
	double dist;

} dist_t;

/* parameters for Q2 */
struct q2_eqpar_t {

	double m;
	double theta;
	double betaL;
	double betaU;
	double gamma;

};

/* parameters for a diagonal matrix */
struct q3_matrixpar_t {

	double *a;
	double *b;
	double *c;
	double *q;
	double *x;

};

/* Question 2 */
void shockwave(const char* q2_file) {
	FILE *file;
	/* open file */
	if ((file = fopen(q2_file, "r")) == NULL) {
		printf("ERROR - file1 could not be opened\n\n");
		exit(EXIT_FAILURE);
	}

	q2_eqpar_t eq;

	/* read in file to equation */
	fscanf(file, "%*s");
	fscanf(file, "%lf, %lf, %lf, %lf, %lf", &eq.m, &eq.theta, &eq.betaL,
			&eq.betaU, &eq.gamma);

	/* initialising */
	int msize = NVALS;
	double *m = malloc(msize * sizeof(double));

	/* to skip the next line */
	fscanf(file, "%*s");

	/* read each line */
	int n = 0;
	while (fscanf(file, "%lf%*s\n", &m[n]) == 1) {
		n++;
		/* reallocate memory */
		if (n >= msize) {
			msize *= 2;
			double* temp = realloc(m, msize * sizeof(double));
			if (!temp) {
				printf("realloc failed\n");
				exit(EXIT_FAILURE);
			}
			m = temp;
		}
	}
	fclose(file);

	/* PRINT OUTPUT TO FILE*/
	FILE *fp = fopen(OUTPUTFILE1, "w");
	if (fp == NULL) {
		printf("Error opening file!\n");
		exit(EXIT_FAILURE);
	}

	/* print header */
	fprintf(fp, "M,theta,beta_lower,beta_upper\n");

	/* differing values of m */
	for (int i = 0; i < n; i++) {
		float mval = m[i];
		/* differing values of theta */
		for (int deg = 0; deg <= MAXDEG; deg++) {
			double theta = deg2rad(deg);

			/* inital guess */
			double bLi, bL = deg2rad(eq.betaL);
			/* calculating B_lower through Newton-Raphson */
			bLi = rad2deg(newton_raph(bL, mval, theta));

			/* inital guess */
			double bUi, bU = deg2rad(eq.betaU);
			/* calculating B_upper through Newton-Raphson */
			bUi = rad2deg(newton_raph(bU, mval, theta));

			/* stop when b is not physically realisable */
			if (bLi < deg || bUi < deg || bLi > MAXDEG || bUi > MAXDEG) {
				break;
			}
			/* print output */
			fprintf(fp, "%.4lf,%d,%.4lf,%.4lf\n", mval, deg, bLi, bUi);
		}
	}
	fclose(fp);
	free(m);

}

/* Newton Raphson method for matrix solution */
double newton_raph(double b, double mval, double theta) {
	double fn, bi;
	int i;
	/* Iteratively solve */
	for (i = 0; i <= MAXITR; i++) {
		fn = shock_funct(b, mval, theta) / shock_der(b, mval, theta);
		bi = b - fn;
		/* stop when solution is found */
		if (fabs(fn) < EPSILON) {
			break;
		}
		b = bi;
	}
	/* if solution could not be found */
	if (i > MAXITR) {
		return -1;
	}
	return bi;
}

/* derivative of the shock function */
double shock_der(double B, double m, double theta) {

	/* naive implementation of derivative */
	double b1 = B + EPSILON;
	double b2 = B - EPSILON;
	double f1 = shock_funct(b1, m, theta);
	double f2 = shock_funct(b2, m, theta);

	return (f2 - f1) / (b2 - b1);
}

/* rad to degrees */
double rad2deg(double rad) {
	return rad * 180 / M_PI;
}

/* rad to degrees */
double deg2rad(double deg) {
	return deg * M_PI / 180;
}

/* shock function */
double shock_funct(double B, double m, double theta) {
	return (2 * cot(B)
			* ((pow(m, 2) * pow(sin(B), 2) - 1)
					/ (pow(m, 2) * (Y + cos(2 * B)) + 2)) - tan(theta));
}

/* cotangent function */
double cot(double theta) {
	return cos(theta) / sin(theta);
}

/* Question 3*/
void linalgbsys(const char* q3_file) {
	FILE *file;
	if ((file = fopen(q3_file, "r")) == NULL) {
		printf("ERROR - file2 could not be opened\n\n");
		exit(EXIT_FAILURE);
	}

	int size = NVALS;

	/* initialising matrix */
	q3_matrixpar_t matrix;
	matrix.a = malloc(size * sizeof(double));
	matrix.b = malloc(size * sizeof(double));
	matrix.c = malloc(size * sizeof(double));
	matrix.q = malloc(size * sizeof(double));
	if (!matrix.a || !matrix.b || !matrix.c || !matrix.q) {
		printf("malloc failed\n");
		exit(EXIT_FAILURE);
	}

	/* skip header line */
	fscanf(file, "%*s");

	int n = 0;
	/* read through matrix */
	while (fscanf(file, "%lf, %lf, %lf, %lf", &matrix.a[n], &matrix.b[n],
			&matrix.c[n], &matrix.q[n]) == 4) {
		/*
		 printf("%lf, %lf, %lf, %lf\n n = %d\n", matrix.a[n], matrix.b[n],
		 matrix.c[n], matrix.q[n], n);
		 */
		n++;
		/* reallocating memory */
		if (n >= size) {
			size *= 2;
			double* temp1 = realloc(matrix.a, size * sizeof(double));
			double* temp2 = realloc(matrix.b, size * sizeof(double));
			double* temp3 = realloc(matrix.c, size * sizeof(double));
			double* temp4 = realloc(matrix.q, size * sizeof(double));

			if (!temp1 || !temp2 || !temp3 || !temp4) {
				printf("realloc failed\n");
				exit(EXIT_FAILURE);
			}
			matrix.a = temp1;
			matrix.b = temp2;
			matrix.c = temp3;
			matrix.q = temp4;

		}
	}
	fclose(file);
	/* initialise x vector*/
	matrix.x = malloc(size * sizeof(double));
	if (!matrix.x) {
		printf("malloc failed\n");
		exit(EXIT_FAILURE);
	}
	/* solve through gauss */
	gauss_sol(&matrix, n, NULL, 0, 0);

	/* PRINT OUTPUT TO FILE*/
	FILE *fp = fopen(OUTPUTFILE2, "w");
	if (fp == NULL) {
		printf("Error opening file!\n");
		exit(EXIT_FAILURE);
	}

	/* print header */
	fprintf(fp, "x\n");
	/* print x vals */
	for (int i = 0; i < n; i++) {
		fprintf(fp, "%.4lf\n", matrix.x[i]);
	}

	/* cleanup */
	fclose(fp);
	free(matrix.a);
	free(matrix.b);
	free(matrix.c);
	free(matrix.q);
	free(matrix.x);

}

/* Gauss solution for a tri-diagonal matrix*/
void gauss_sol(q3_matrixpar_t* matrix, int n, double *heateq, int Nt, int t) {
	int i, j;
	double mgauss_a[n];
	double mgauss_q[n];

	/* iterating forwards */
	for (i = 0; i < n; i++) {

		if (i == 0) {
			mgauss_a[i] = matrix->a[i];
			mgauss_q[i] = matrix->q[i];
		} else {
			mgauss_a[i] = matrix->a[i]
					- (matrix->c[i] * matrix->b[i - 1]) / mgauss_a[i - 1];
			mgauss_q[i] = matrix->q[i]
					- (matrix->c[i] * mgauss_q[i - 1]) / mgauss_a[i - 1];
		}
	}
	int num = n - 1;

	/* iterating backwards */
	for (j = num; j >= 0; j--) {
		/* solving into heat eqn (used for Q6) */
		if (heateq != NULL) {
			if (j == num) {
				heateq[(Nt) * j + t] = mgauss_q[j] / mgauss_a[j];
			} else {
				heateq[(Nt) * j + t] = (mgauss_q[j]
						- matrix->b[j] * matrix->x[j + 1]) / mgauss_a[j];
			}
		}
		if (j == num) {
			matrix->x[j] = mgauss_q[j] / mgauss_a[j];
		} else {
			matrix->x[j] = (mgauss_q[j] - matrix->b[j] * matrix->x[j + 1])
					/ mgauss_a[j];
		}
	}
}

/* Question 5*/
void interp(const char* q5_file, const double xo) {
	FILE* file;

	/* read in file */
	if ((file = fopen(q5_file, "r")) == NULL) {
		printf("ERROR - file3 could not be opened\n\n");
		exit(EXIT_FAILURE);
	}

	/* skip header */
	fscanf(file, "%*s");

	/* initialise values */
	int size = NVALS;
	double *x = malloc(size * sizeof(double));
	double *fx = malloc(size * sizeof(double));
	if (!x || !fx) {
		printf("malloc failed\n");
		exit(EXIT_FAILURE);
	}

	int n = 0;
	/* read in values */
	while ((fscanf(file, "%lf, %lf", &x[n], &fx[n])) == 2) {
		n++;
		/* reallocating memory */
		if (n >= size) {
			size *= 2;
			double *temp1 = realloc(x, size * sizeof(double));
			double *temp2 = realloc(fx, size * sizeof(double));

			if (!temp1 || !temp2) {
				printf("realloc failed\n");
				exit(EXIT_FAILURE);
			}
			x = temp1;
			fx = temp2;
		}
	}
	fclose(file);

	/* initialising distance struct */
	dist_t xdist[n];

	int p;

	/* setting the distance struct */
	for (p = 0; p < n; p++) {
		xdist[p].x = x[p];
		xdist[p].fx = fx[p];
		xdist[p].dist = fabs(x[p] - xo);
	}

	/* sort to find which three points should be used for lagrange */
	qsort(xdist, n, 3 * sizeof(double), cmpfnc);

	double lagsum = 0;
	double l[LAGNUM];

	/* lagrange formula */
	int i, j;
	for (i = 0; i < LAGNUM; i++) {
		l[i] = 1.0;
		for (j = 0; j < LAGNUM; j++) {
			if (i == j) {
				continue;
			}
			l[i] *= (xo - xdist[j].x) / (xdist[i].x - xdist[j].x);
		}
		lagsum += l[i] * xdist[i].fx;
	}

//	double b0 = fx[1];
//	double b1 = (fx[2] - b0) / (x[2] - x[1]);
//	double b2 = (((fx[3] - fx[2]) / (x[3] - x[2])) - b1) / (x[3] - x[1]);

	/* initialising values for tri-diagonal matrix */
	q3_matrixpar_t matrix;
	matrix.a = malloc(n * sizeof(double));
	matrix.b = malloc(n * sizeof(double));
	matrix.c = malloc(n * sizeof(double));
	matrix.q = malloc(n * sizeof(double));
	matrix.x = malloc(n * sizeof(double));
	if (!matrix.a || !matrix.b || !matrix.c || !matrix.q || !matrix.x) {
		printf("malloc failed\n");
		exit(EXIT_FAILURE);
	}

	/* function variables */
	double a[n];
	double b[n];
	double c[n];
	double d[n];
	double h[n];

	/* natural spline, set double der to 0 at end points */
	c[0] = 0;
	c[n - 1] = 0;
	/* initialising */
	for (i = 0; i < n; i++) {
		a[i] = fx[i];
		b[i] = 0;
		c[i] = 0;
		d[i] = 0;
		h[i] = x[i + 1] - x[i];
	}
	for (i = 0; i < n; i++) {
		/* fill triangular matrix */
		if (i == 0) {
			matrix.a[i] = 1;
			matrix.b[i] = 0;
			matrix.c[i] = 0;
			matrix.q[i] = 0;
		} else if (i == n - 1) {
			matrix.a[i] = 1;
			matrix.b[i] = 0;
			matrix.c[i] = 0;
			matrix.q[i] = 0;
		} else {
			/* spline formula */
			matrix.a[i] = 2 * (h[i - 1] + h[i]);
			matrix.b[i] = h[i];
			matrix.c[i] = h[i - 1];
			matrix.q[i] = (3 / h[i]) * (a[i + 1] - a[i])
					+ (3 / h[i - 1]) * (a[i - 1] - a[i]);
		}
	}

	/* solving triangular matrix */
	gauss_sol(&matrix, n, NULL, 0, 0);

	int m = 0;
	for (i = 1; i < n; i++) {

		/* setting c values */
		c[i] = matrix.x[i];

		/* to find interval to evaluate xo */
		if (x[i] < xo) {
			m = i;
		}
	}

	/* solving for the remaining variables */
	for (i = 0; i < n; i++) {
		b[i] = (1 / h[i]) * (a[i + 1] - a[i])
				- (h[i] / 3) * (2 * c[i] + c[i + 1]);
		d[i] = (c[i + 1] - c[i]) / (3 * h[i]);
	}

	/* cubic spline formula */
	double cubspl = a[m] + b[m] * (xo - x[m]) + c[m] * pow((xo - x[m]), 2)
			+ d[m] * pow((xo - x[m]), 3);

	/* PRINT OUTPUT TO FILE*/
	FILE *fp = fopen(OUTPUTFILE3, "w");
	if (fp == NULL) {
		printf("Error opening file!\n");
		exit(EXIT_FAILURE);
	}

	/* print output */
	fprintf(fp, "lagrange\n");
	fprintf(fp, "%.4lf\n", lagsum);
	fprintf(fp, "cubic\n");
	fprintf(fp, "%.4lf\n", cubspl);

	/* cleanup */
	fclose(fp);
	free(x);
	free(fx);
	free(matrix.a);
	free(matrix.b);
	free(matrix.c);
	free(matrix.q);
	free(matrix.x);

}

/* function to sort x based on distance from xo */
int cmpfnc(const void * a, const void * b) {
	dist_t ab = *(dist_t*) a;
	dist_t bb = *(dist_t*) b;
	if (ab.dist > bb.dist)
		return 1;
	else if (ab.dist < bb.dist)
		return -1;
	else
		return 0;
}

/* Question 6 */
void heateqn(const char* q6_file) {
	/* read in file */
	FILE *file;
	if ((file = fopen(q6_file, "r")) == NULL) {
		printf("ERROR - file4 could not be opened\n\n");
		exit(EXIT_FAILURE);
	}

	/* initialising */
	double mu = 0.0;
	int Nx = 0, Nt = 0, i = 0, j = 0, k = 0;

	/* skip header */
	fscanf(file, "%*s");
	fscanf(file, "%lf, %d, %d", &mu, &Nx, &Nt);

	/* close file */
	fclose(file);

	/* initialise arrays to store explicit fe, explicit ve, and implicit*/
	double *fn_e = calloc((Nx + 1) * Nt, sizeof(double));
	if (!fn_e) {
		printf("malloc failed\n");
		exit(EXIT_FAILURE);
	}
	double *fn_e_ve = calloc((Nx + 1) * Nt, sizeof(double));
	if (!fn_e_ve) {
		printf("malloc failed\n");
		exit(EXIT_FAILURE);
	}
	double *fn_i = calloc((Nx + 1) * Nt, sizeof(double));
	if (!fn_i) {
		printf("malloc failed\n");
		exit(EXIT_FAILURE);
	}
	/* delta x and delta t*/
	double deltax = XMAX / Nx;
	double deltat = TMAX / Nt;

	/* setting initial conditions */
	for (i = 0; i <= Nx; i++) {
		double x = 1.0 * (XMAX * i) / Nx;
		if (x >= 0.125 && x <= 0.35) {
			fn_e[(Nt) * i + 0] = 0.5 * (1 - cos(8 * M_PI * (x - 0.125)));
			fn_e_ve[(Nt) * i + 0] = 0.5 * (1 - cos(8 * M_PI * (x - 0.125)));
			fn_i[(Nt) * i + 0] = 0.5 * (1 - cos(8 * M_PI * (x - 0.125)));
		} else {
			fn_e[(Nt) * i + 0] = 0;
			fn_e_ve[(Nt) * i + 0] = 0;
			fn_i[(Nt) * i + 0] = 0;
		}
	}

	/* use explicit Euler to find remaining f */
	for (k = 1; k < Nt; k++) {
		for (j = 0; j <= Nx; j++) {
			fn_e[(Nt) * j + k] = explicit_int(fn_e, Nt, Nx, deltat, deltax, j,
					k, mu);
			//		printf("%.4lf, %.4lf\n", 1.0 * (XMAX * j) / Nx, fn_e[(Nt) * j + k]);
			fn_e_ve[(Nt) * j + k] = explicit_int_ve(fn_e, Nt, Nx, deltat,
					deltax, j, k, mu);
		}
		/* use implicit Euler to find remaining f */
		implicit_int(fn_i, Nt, Nx, deltat, deltax, j, k, mu);
	}

	/* PRINT OUTPUT TO FILE - explicit fixed ends*/
	FILE *fp1 = fopen(OUTPUTFILE4, "w");
	/* PRINT OUTPUT TO FILE - explicit variable ends*/
	FILE *fp2 = fopen(OUTPUTFILE5, "w");
	/* PRINT OUTPUT TO FILE - implicit fixed ends*/
	FILE *fp3 = fopen(OUTPUTFILE6, "w");
	if (!fp1 || !fp2 || !fp3) {
		printf("Error opening file!\n");
		exit(EXIT_FAILURE);
	}

	/* print header */
	fprintf(fp1, "x,f(x)\n");
	fprintf(fp2, "x,f(x)\n");
	fprintf(fp3, "x,f(x)\n");

	/* iteratively print values at 100th timestep i.e. t = 99 */
	for (i = 0; i <= Nx; i++) {
		fprintf(fp1, "%.4lf, %.4lf\n", 1.0 * (XMAX * i) / Nx,
				fn_e[(Nt) * i + NTPRINT]);
		fprintf(fp2, "%.4lf, %.4lf\n", 1.0 * (XMAX * i) / Nx,
				fn_e_ve[(Nt) * i + NTPRINT]);
		fprintf(fp3, "%.4lf, %.4lf\n", 1.0 * (XMAX * i) / Nx,
				fn_i[(Nt) * i + NTPRINT]);
	}

	/* cleanup */
	fclose(fp1);
	fclose(fp2);
	fclose(fp3);
	free(fn_e);
	free(fn_e_ve);
	free(fn_i);

}

/* explicit solution with fixed ends */
double explicit_int(double* fn, int Nt, int Nx, double deltat, double deltax,
		int x, int t, double mu) {
	t = t - 1;
	double fval = fn[(Nt) * x + t];
	/* fixed end implementation */
	if (x == 0 || x == Nx) {
		return fval;
	}
	return fval + deltat * RHS(fn, x, t, Nx, Nt, deltax, mu);
}

/* explicit solution with variable ends*/
double explicit_int_ve(double* fn, int Nt, int Nx, double deltat, double deltax,
		int x, int t, double mu) {
	t = t - 1;
	double fval = fn[(Nt) * x + t];
	return fval + deltat * RHS(fn, x, t, Nx, Nt, deltax, mu);
}

/* implicit solution with fixed ends */
void implicit_int(double* fn, int Nt, int Nx, double deltat, double deltax,
		int x, int t, double mu) {

	double a, b;
	/* tri-diagonal values */
	a = 1 + (2 * deltat * mu) / pow(deltax, 2);
	b = -(deltat * mu) / pow(deltax, 2);

	q3_matrixpar_t matrix;
	matrix.a = malloc(Nx * sizeof(double));
	matrix.b = malloc(Nx * sizeof(double));
	matrix.c = malloc(Nx * sizeof(double));
	matrix.q = malloc(Nx * sizeof(double));
	matrix.x = malloc(Nx * sizeof(double));
	if (!matrix.a || !matrix.b || !matrix.c || !matrix.q || !matrix.x) {
		printf("malloc failed\n");
		exit(EXIT_FAILURE);
	}

	int i;
	for (i = 0; i < Nx; i++) {
		/* fill triangular matrix */
		if (i == 0) {
			matrix.a[i] = 1;
			matrix.b[i] = 0;
			matrix.c[i] = 0;
			matrix.q[i] = fn[(Nt) * i + (t - 1)];
		} else if (i == Nx - 1) {
			matrix.a[i] = 1;
			matrix.b[i] = 0;
			matrix.c[i] = 0;
			matrix.q[i] = fn[(Nt) * i + (t - 1)];
		} else {
			matrix.a[i] = a;
			matrix.b[i] = b;
			matrix.c[i] = b;
			matrix.q[i] = fn[(Nt) * i + (t - 1)];
		}
	}
	/* solve tri-diagonal matrix to fill fn with values */
	gauss_sol(&matrix, Nx, fn, Nt, t);

	/* cleanup */
	free(matrix.a);
	free(matrix.b);
	free(matrix.c);
	free(matrix.q);
	free(matrix.x);

}

/* RHS solution */
double RHS(double* fn, int x, int t, int Nx, int Nt, double deltax, double mu) {
	/* variable ends formula*/
	if (x == 0) {
		return mu
				* ((fn[(Nt) * (x) + t]) - 2 * (fn[(Nt) * (x + 1) + t])
						+ (fn[(Nt) * (x + 2) + t])) / pow(deltax, 2);
	}
	/* variable ends formula */
	if (x == Nx) {
		return mu
				* ((fn[(Nt) * (x) + t]) - 2 * (fn[(Nt) * (x - 1) + t])
						+ (fn[(Nt) * (x - 2) + t])) / pow(deltax, 2);
	}
	/* general sol formula */
	return mu
			* ((fn[(Nt) * (x + 1) + t]) - 2 * (fn[(Nt) * (x) + t])
					+ (fn[(Nt) * (x - 1) + t])) / pow(deltax, 2);
}

