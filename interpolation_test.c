//This is to test interpolation efficiency

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

//Random test array
#define TEST 1000000
double test[TEST];
double fres[TEST];
double lres[TEST];

//Linear Interpolation array
double * linx;
double * liny;
int linmax;

//Function direct calculation
double f(double);

//Second Derivative
double fpp(double);

//Function linear interpolate
double flin(double);

//Create linear interpolation array
void create_int(double, double, double);

//Destroy linear interpolation array
void free_int();

//Populate test array with random samples
void populate(double, double);

// MAIN FUNCTION
int main() {
	srand(time(0));
	
	double a, b, lambda, error, maxerr, errtot;
	clock_t fstart, fend, lstart, lend;
	double ftime, ltime;
	
	a = 0.5;
	b = 10.0;
	lambda = 0.1;
	
	populate(a, b);
	create_int(a, b, lambda);
	
	fstart = clock();
	for (int i=0; i < TEST; i++) {
		fres[i] = f(test[i]);
	}
	fend = clock();
	
	lstart = clock();
	for (int i=0; i < TEST; i++) {
		lres[i] = flin(test[i]);
	}
	lend = clock();
	
	error = 0;
	maxerr = 0;
	errtot = 0;
	for (int i=0; i < TEST; i++) {
		error = fabs(fres[i] - lres[i]);
		errtot += error;
		if (error > maxerr) {
			maxerr = error;
		}
	}
	
	ftime = 1000.0 * (fend - fstart) / CLOCKS_PER_SEC;
	ltime = 1000.0 * (lend - lstart) / CLOCKS_PER_SEC;
	
	printf("\n%d Sample Points Run\n", TEST);
	printf("%d Interpolation Points Created\n\n", linmax);
	printf("  f Total Time: %lf ms\n", ftime);
	printf("  l Total Time: %lf ms\n", ltime);
	printf("Maximum Error:  %e\n", maxerr);
	printf("Average Error:  %e\n\n", errtot / (double) TEST);
	
	free_int();

  return 0;
}

//Populate test array with random samples
void populate(double a, double b) {
	for (int i=0; i < TEST; i++) {
		test[i] = a + (b-a) * rand()/((double)RAND_MAX + 1);
	}
	
	return;
}

//Function direct calculation
double f(double r) {
  /*double r12, r6;
  
  r12 = r*r*r*r * r*r*r*r * r*r*r*r;
  r6 = r*r*r * r*r*r;
  
  return 1.0/r12 - 1.0/r6;*/
  
  return sqrt(r);
}

//Second Derivative
double fpp(double r) {
	double sr;
	
	sr = sqrt(r);
	
	/*double r14, r8;
	
	r14 = r*r*r*r*r * r*r*r*r*r * r*r*r*r;
	r8  = r*r*r*r*r * r*r*r;*/
	
	return -1 / (4 * sr * sr * sr);
}

//Function linear interpolate
double flin(double x) {
	double y;
	
	if (x == linx[0]) {
		return liny[0];
	} else {
		for (int i = 1; i < linmax; i++) {
			if (x > linx[i]) {
				y = liny[i-1];
				y += (liny[i] - liny[i-1]) * (x - linx[i-1]) / (linx[i] - linx[i-1]);
				return y;
			}
		}
	}
	
	return liny[linmax - 1];
}

//Create linear interpolation array
void create_int(double a, double b, double lambda) {
	int N;
	double x;
	
	linmax = 1000;
	
	//guess an initial allocation for interpolation array
	linx = (double *) malloc(linmax * sizeof(double));
	liny = (double *) malloc(linmax * sizeof(double));
	
	
	N = 0;
	x = b;
	do {
		linx[N] = x;
		liny[N] = f(x);
		
		/*
		printf("x[%d] = %lf\n", N, x);
		printf("fpp(x) = %lf\n", fpp(x));
		printf("|fpp(x)| = %lf\n", fabs(fpp(x)));
		printf("sqrt(|fpp(s)|) = %lf\n", sqrt(fabs(fpp(x))));
		printf("lambda / sqrt(...) = %lf\n", lambda / sqrt(fabs(fpp(x))));
		*/
		
		N++;
		x -= lambda / sqrt(fabs(fpp(x)));
		
		if (N >= linmax) {
			linmax += 100;
			linx = (double *) realloc(linx, linmax * sizeof(double));
			liny = (double *) realloc(liny, linmax * sizeof(double));
		}
	} while (x > a);
	
	linx[N] = x;
	liny[N] = f(x);
	N++;
	
	if (N < linmax) {
		linmax = N;
		linx = (double *) realloc(linx, linmax * sizeof(double));
		liny = (double *) realloc(liny, linmax * sizeof(double));
	}
	
	return;
}


//Destroy linear interpolation array
void free_int() {
	free(linx);
	free(liny);
	
	return;
}

