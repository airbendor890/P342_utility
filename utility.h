#include<stdio.h>


/*A4*/
void mat_mul(int n,double arr1[n][n],double arr2[n][n]);
double lu_dcmpsn(int n,double a[n][n],int * P);
void lu_partial_pivot(int n,double arr[n][n],int * P,int r);
void forward_back_sub(int n,double a[n][n],int *P,double *b);

/*A5*/
double root_newton_raphson(double(*func)(double),double x0,FILE*);
double root_false_position(double (*func)(double),double x1,double x2,FILE*);
double root_bisection(double(*func)(double),double x1,double x2,FILE *);
int brac(double (*func)(double),double* x1,double* x2);

double func(double);
double D1(double(*func)(double),double);
double D2(double(*func)(double),double);

void drive_laguerre();
double laguerre(double * arr,int n,double x0,int * iter);
double poly(double* arr,int n,double x0);
double D1_poly(double* arr,int n,double x0);
double D2_poly(double * arr,int n,double x0);

/*A6*/
