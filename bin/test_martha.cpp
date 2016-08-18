/*
Kyriakos Hadjiyiannakou
Extract Axial form factor
Mimic Thomas code
*/

/*
  Contents
-> global scope
1) defintions of operators overload
2) macros for defining constants
3) forward declaration of functions

-> main body
1) Greeting
2) check passing inputs right
3) fill filenames and print paths
4) calculate and print number of configurations
5) read configuration's number
6) create filename names of threep and twop
7) allocate memory & read data
8) calculate number of bins
9) flip the sign of momenta. Because q=p_final-p_initial, while in formulas: p_initial
10) aver twop with momenta giving same q^2
11) create averages over bins for threep and twop
12) calculate naive effective mass
13) calculate naive ratio
14) binning for jackknife analysis
15) jackknife effective and ratio
16) mean bin values
17) jacknife errors
18) fit to extract plateaus
19) jacknife error to fitted values
20) extract quantities naive
21) extract quantities for each bin
22) find jackknife errors for form factors
23) calculate data for ratio to print out
24) print results

-> function
1) fit_ratio_plato
2) fit_EffMass_plato
3) extractGAGP
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex>
#include <string>
#include <fstream>
typedef std::complex<double> Complex;
#include <lapacke.h>
#include <time.h>


typedef struct{
  int mom;
  int mu;
} Combinations;


Complex operator*(double a,Complex b){
  Complex res;
  res.real() = a * b.real();
  res.imag() = a * b.imag();
  return res;
}

Complex operator*(Complex a,double b){
  Complex res;
  res.real() = a.real() * b;
  res.imag() = a.imag() * b;
  return res;
}

Complex operator/(Complex a,double b){
  Complex res;
  res.real() = a.real() / b;
  res.imag() = a.imag() / b;
  return res;
}
int main(){
  double A[2*72];
  //      make weighted least square solution:
  FILE* ptr_file;
  ptr_file = fopen("table.dat","r");

  for(int i = 0 ; i < 72 ; i++)
    fscanf(ptr_file,"%lf %lf",&(A[0*72+i]) , &(A[1*72+i]) );      
      // ready for singular value decomposition

      
      int info;
      int lwork;
      double wkopt;
      double* work=NULL;
      int m=72;
      int n=2;
      int lda=m;
      int ldu=m;
      int ldvt=n;

      double *u=(double*)malloc(ldu*m*sizeof(double));
      double *w=(double*)malloc(n*sizeof(double));
      double *vt=(double*)malloc(ldvt*n*sizeof(double));
      double *oneOverW =(double*)calloc(n*n,sizeof(double));

      double x[2][72];

      lwork = -1;
      dgesvd_((char*) "All",(char*) "All", &m, &n, A, &lda, w, u, &ldu, vt, &ldvt, &wkopt, &lwork, &info );
      lwork = (int)wkopt;
      work = (double*)malloc( lwork*sizeof(double) );
      // Compute SVD 
      dgesvd_((char*) "All",(char*) "All", &m, &n, A, &lda, w, u, &ldu, vt, &ldvt, work, &lwork, &info );
      // Check for convergence 
      if( info > 0 ) {
	printf( "The algorithm computing SVD failed to converge.\n" );
	exit( 1 );
      }     
      
      oneOverW[0*2+0] = 1./w[0];
      oneOverW[1*2+1] = 1./w[1];

      for(int alpha = 0 ; alpha < 2 ; alpha++)
	for(int beta = 0 ; beta < 2 ; beta++)
	  for(int gamma = 0 ; gamma < 2 ; gamma++)
	    for(int delta = 0 ; delta < 72 ; delta++){
	      x[alpha][delta] += vt[alpha*2+beta] * oneOverW[beta*2+gamma] * u[gamma*m+delta]; 
	    }

      for(int i = 0 ; i < 72 ; i++)
	printf("%d %+f %+f\n",i,x[0][i],x[1][i]);

}

