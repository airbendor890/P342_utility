#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"utility.h"

/*--------------------------------A4-----------------------------------------*/

//-----------------------------------------------------------------------
// function for partial pivot used in lu dcmpsn
void lu_partial_pivot(int n,double arr[n][n],int* P,int r){
		int r1=r+1;
		while(arr[r][r]==0 && r1<n){
			if(abs(arr[r1][r])>abs(arr[r][r])){
				double temp;
				for(int c=0;c<n;c++){
					temp=arr[r1][c];
					arr[r1][c]=arr[r][c];
					arr[r][c]=temp;
				}

				int dum;
				dum=P[r1];
				P[r1]=P[r];
				P[r]=dum;
			}
			r1++;
		}



}
//-----------------------------------------------------------------------
//lu dcmpsn function
double lu_dcmpsn(int n,double a[n][n],int *P){
	/* P will hold the permutation infromation that will occur during 
	partial pivoting*/

	double sum=0.0;
	double det=1.0;
	for(int j=0;j<n;j++){	/*loop over columns of crouts method*/
		for(int i=0;i<=j;i++){
			sum=a[i][j];
			for(int k=0;k<i;k++) sum=sum-a[i][k]*a[k][j];
			a[i][j]=sum;
			
		}
		if(a[j][j]==0) lu_partial_pivot(n,a,P,j);

		for(int i=j+1;i<n;i++){
			sum=a[i][j];
			for(int k=0;k<j;k++) sum=sum-a[i][k]*a[k][j];
			//still if diag element is 0 its a singular matrix
			//return 0  
			if(a[j][j]==0.0){	
				printf("singular matrix !!");
				return 0;
			}
			else
			a[i][j]=sum/a[j][j];
		}

	}

	//its not singular ..return its determinant
	for(int i=0;i<n;i++)det=det*a[i][i];

	return det;

}	

//-------------------------------------------------------------------

void forward_back_sub(int n,double a[n][n],int *P,double * b){
	/* first make the permutation matrix from P*/
	int perm[n][n];
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			if(j==P[i])
				perm[i][j]=1;
			else
				perm[i][j]=0;

		}
	}

	/* now permute the vector b*/
	double c[n];
	for(int i=0;i<n;i++){
		double sum=0.0;
		for(int j=0;j<n;j++)
			sum+=perm[i][j]*b[j];
		c[i]=sum;
	}
	for(int i=0;i<n;i++)b[i]=c[i];

	
	double sum=0.0;
	/*forward substitution Ly=b */
	for(int i=0;i<n;i++){
		sum=b[i];
		for(int j=0;j<i;j++){
			sum=sum-a[i][j]*b[j];
		}
		b[i]=sum;
	}

	/* back subtitution Ux=y */

	for(int i=n-1;i>=0;i--){
		sum=b[i];
		for(int j=i+1;j<n;j++){
			sum=sum-a[i][j]*b[j];
		}

		b[i]=sum/a[i][i];
	}

}

//------------------------------------------------------------
void mat_mul(int n,double arr1[n][n],double arr2[n][n]){

	double mul[n][n];
	for(int i=0;i<n;i++)
	for(int j=0;j<n;j++)mul[i][j]=0;

	for(int i=0;i<n;i++)
	for(int j=0;j<n;j++)
	for(int k=0;k<n;k++)
	mul[i][j]+=arr1[i][k]*arr2[k][j];

	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++)
		printf("%lf  ",mul[i][j]);
		printf("\n\n");
	}
}


/*----------------------A5-------------------------------------------------------*/


//--------------------------------------------------------------------------------
/*root bracketing function*/
int brac(double (*func)(double),double* x1,double* x2){
	double f1,f2;
	double FACTOR=1.5;
	
	//if(*x1=*x2)break;/*bad range*/
	f1=(*func)(*x1);
	f2=(*func)(*x2);
	for(int j=1;j<=50;j++){
		if(f1*f2<0.0) return 1;	/*change of sign...bracketed*/		
		if(fabs(f1)<fabs(f2))
			f1=(*func)(*x1+=FACTOR*(*x1-*x2));
		else
			f2=(*func)(*x2+=FACTOR*(*x2-*x1));
	}

	return 0;	/*not bracketed even after 50 trys*/
}


//---------------------------------------------------------------------------------
/*root bisection method function*/
double root_bisection(double (*func)(double),double x1,double x2,FILE* fptr)
{	/*assuming input is such that x2>x1*/
	/* fptr will write the iterations to file*/
	/*it must be passed with file opened*/
	double xmid;	
	int JMAX=200;
	double epsilon=0.000001;

	for(int j=1;j<=JMAX;j++){
		xmid=0.5*(x1+x2);

		if((*func)(x1)*(*func)(xmid)<0)
			x2=xmid;
		else			
			x1=xmid;
		
		printf("\n%d	%lf	%lf	%lf\n",j,x1,x2,fabs(x1-x2));
		fprintf(fptr,"%d	%lf\n",j,fabs(x1-x2));
		/*if converged to desirable accuracy*/
		if(fabs(x1-x2)<epsilon){fclose(fptr);return xmid;}

	}

}

//------------------------------------------------------------------------------------
/*false position root finding method function*/
double root_false_position(double (*func)(double),double x1,double x2,FILE * fptr)
{
	/*assumption is that x1 and x2 are passsed after bracketing*/
	/*fptr will write the data to opened file*/
	/*assumption is that file is already open before passing*/
	double epsilon=0.000001;
	int JMAX=200;
	double c_k,c_k1;

	c_k1=x2-((x2-x1)*(*func)(x2))/((*func)(x2)-(*func)(x1));
	if((*func)(x1)*(*func)(c_k1)<0)
		x2=c_k1;
	else
		x1=c_k1;


	for(int j=1;j<=JMAX;j++){
		c_k=c_k1;
		c_k1=x2-((x2-x1)*(*func)(x2))/((*func)(x2)-(*func)(x1));
		
		if((*func)(x1)*(*func)(c_k1	)<0)
			x2=c_k1;
		else
			x1=c_k1;

		printf("\n%d	%lf	%lf	%lf\n",j,c_k,c_k1,fabs(c_k1-c_k));
		fprintf(fptr,"%d	%lf\n",j,fabs(c_k1-c_k));

		if(fabs(c_k1-c_k)<epsilon) return c_k1;
	}

}

//--------------------------------------------------------------------------------------
/*newton raphson root finding method function*/
double root_newton_raphson(double(*func)(double),double x0,FILE * fptr){
	/*x0 is the gussed root*/
	/*fptr will write convergence data to file*/
	/*file pointer must be passed already linked to a file*/
	double epsilon=0.000001;
	int JMAX=40;
	double x_n,x_n1;

	x_n1=x0;
	for(int j=1;j<=JMAX;j++){
		x_n=x_n1;
		x_n1=x_n-((*func)(x_n)/D1(func,x_n));	/*first derivative*/
		fprintf(fptr,"%d	%lf\n",j,fabs(x_n1-x_n));
		printf("\n%d	%lf	%lf	%lf\n",j,x_n,x_n1,fabs(x_n1-x_n));
		if(fabs(x_n1-x_n)<epsilon) return x_n;
	}

}

//-------------------------------------------------------------------------------------
/*1st derivative for one variable function*/
double D1(double(*func)(double),double x0){
	double h=0.001;
	
	double d=((*func)(x0+h)-(*func)(x0-h))/(2*h);
	return d;
}

//------------------------------------------------------------------------------------

/*2nd derivative for  one variavble function*/

double D2(double(*func)(double),double x0){
	double h=0.0001;
	double d=((*func)(x0+h)+(*func)(x0-h)-2*(*func)(x0))/(h*h);

	return d;

}

//------------------------------------------------------------------------------------
void drive_laguerre(){
	int deg;
	FILE* file;
	printf("\nLaguerre and synthetic division method for real roots of polynomial of degree n\n");
	printf("\nEnter degree here and store coefficients of polynomial in a text file(poly.txt) in descending power:\n");
	scanf("%d",&deg);
	double arr[deg+1];
	double roots[deg];/*store roots*/
	int iter;
	file=fopen("poly.txt","r");
	for(int i=0;i<=deg;i++)
		fscanf(file,"%lf",&(arr[i]));
	fclose(file);	

	printf("\nroot		iterations taken\n");
	for(int j=0;j<deg;j++){
		double x0=5.00;
		roots[j]=laguerre(arr,deg-j,x0,&iter);
		printf("\n%lf		%d\n",roots[j],iter);
	
		/*perform deflation*/
		double r=0.0;
		for(int jj=0;jj<=deg-j;jj++){
			arr[jj]+=r*(roots[j]);
			r=arr[jj];

		}
		/*end deflation*/		
	}


}

//------------------------------------------------------------------------------------
/*function: laguerre method to find a real root of a polynomial*/

double laguerre(double * arr,int n,double x0,int * iter){
	/*array input n degree polynomial,guess root x0*/
	/*no of iteration taken in *iter */
	int MAXIT=50;
	double epsilon=1.0e-6;
	if(fabs(poly(arr,n,x0))<epsilon){*iter=0;return x0;}
	double G,H,a_k,a_k1;
	
	a_k1=x0;
	for(int j=1;j<MAXIT;j++){
		double a,d1,d2;
		a_k=a_k1;
		if(fabs(poly(arr,n,a_k))<epsilon){*iter=j-1;return a_k;}

		G=D1_poly(arr,n,a_k)/poly(arr,n,a_k);
		H=pow(G,2)-(D2_poly(arr,n,a_k)/poly(arr,n,a_k));
		d1=G+sqrt((n-1)*(n*H-pow(G,2)));
		d2=G-sqrt((n-1)*(n*H-pow(G,2)));
	
		if(fabs(d1)>fabs(d2))
			a=n/d1;
		else
			a=n/d2;
				
		a_k1=a_k-a;
		//printf("\n%lf	%lf	%lf\n",a,d1,d2);
		if(fabs(a_k1-a_k)<epsilon){*iter=j; return a_k;}
		if(fabs(poly(arr,n,a_k1))<epsilon){*iter=j;return a_k1;}

	}
	
	*iter=MAXIT;/*failed to find root within MAXIT*/

}

//------------------------------------------------------------------------------------
/*function:polynomial value at x0(array input)*/
double poly(double* arr,int n,double x0){
	double value=0.0;
	for(int i=0;i<=n;i++)
		value+=arr[i]*pow(x0,n-i);
	return value;
}

//------------------------------------------------------------------------------------
/*function: first derivarive of a polynomial function At x0(array input)*/
double D1_poly(double* arr,int n,double x0){
	double value=0.0;
	for(int i=0;i<n;i++)
		value+=arr[i]*(n-i)*pow(x0,n-i-1);
//	printf("\n%lf\n",value);
	return value;
}

//------------------------------------------------------------------------------------
/*function: second derivative of a polynomial function at x0(array input)*/
double D2_poly(double * arr,int n,double x0){
	double value=0.0;
	for(int i=0;i<n-1;i++)
		value+=arr[i]*(n-i)*(n-i-1)*pow(x0,n-i-2);

	return value;
}

//-----------------------------------------------------------------------------------
