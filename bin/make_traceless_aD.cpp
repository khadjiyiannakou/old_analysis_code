#include <stdio.h>
#include <stdlib.h>
#include <complex>

typedef std::complex<double> Complex;

#define TSINK 12
#define NMOM 257

int main(int argc , char *argv[]){
  if(argc != 3){
    fprintf(stderr,"error wrong number of input\n");
    exit(EXIT_FAILURE);
  }

  FILE *ptr_in, *ptr_out;

  ptr_in = fopen(argv[1],"r");
  ptr_out = fopen(argv[2],"w");

  if(ptr_in == NULL){
    fprintf(stderr,"error open file for reading\n");
    exit(EXIT_FAILURE);
  }

  if(ptr_out == NULL){
    fprintf(stderr,"error open file for writting\n");
    exit(EXIT_FAILURE);
  }

  int dummy;
  Complex *threep = (Complex*)malloc((TSINK+1)*4*4*NMOM*sizeof(Complex));
  Complex *threep_traceless = (Complex*)malloc((TSINK+1)*4*4*NMOM*sizeof(Complex));

  int momList[3][NMOM];

  for(int it=0; it < TSINK+1 ; it++)                                                                                                                                        
    for(int mu = 0 ; mu < 4 ; mu++)
      for(int nu = 0 ; nu <= mu ; nu++)
	for(int imom = 0 ; imom < NMOM ; imom++)
	  int returnValue = fscanf(ptr_in,"%d %d %d %d %lf %lf %d %d",&dummy,&(momList[0][imom]), &(momList[1][imom]), &(momList[2][imom]), &(threep[it*4*4*NMOM + mu*4*NMOM + nu*NMOM + imom].real()), &(threep[it*4*4*NMOM + mu*4*NMOM + nu*NMOM + imom].imag()), &dummy, &dummy );

  for(int it=0; it < TSINK+1 ; it++)
    for(int mu = 0 ; mu < 4 ; mu++)
      for(int nu = 0 ; nu <= mu ; nu++)
        for(int imom = 0 ; imom < NMOM ; imom++){

	  if( mu == 0 && nu == 0){
	    threep_traceless[it*4*4*NMOM + 0*4*NMOM + 0*NMOM + imom].real() = 0.75 * threep[it*4*4*NMOM + 0*4*NMOM + 0*NMOM + imom].real() - 0.25* ( threep[it*4*4*NMOM + 1*4*NMOM + 1*NMOM + imom].real() + threep[it*4*4*NMOM + 2*4*NMOM + 2*NMOM + imom].real() + threep[it*4*4*NMOM + 3*4*NMOM + 3*NMOM + imom].real()); 
	    threep_traceless[it*4*4*NMOM + 0*4*NMOM + 0*NMOM + imom].imag() = 0.75 * threep[it*4*4*NMOM + 0*4*NMOM + 0*NMOM + imom].imag() - 0.25* ( threep[it*4*4*NMOM + 1*4*NMOM + 1*NMOM + imom].imag() + threep[it*4*4*NMOM + 2*4*NMOM + 2*NMOM + imom].imag() + threep[it*4*4*NMOM + 3*4*NMOM + 3*NMOM + imom].imag()); 
	  }
	  else if( mu == 1 && nu == 1){
	    threep_traceless[it*4*4*NMOM + 1*4*NMOM + 1*NMOM + imom].real() = 0.75 * threep[it*4*4*NMOM + 1*4*NMOM + 1*NMOM + imom].real() - 0.25* ( threep[it*4*4*NMOM + 0*4*NMOM + 0*NMOM + imom].real() + threep[it*4*4*NMOM + 2*4*NMOM + 2*NMOM + imom].real() + threep[it*4*4*NMOM + 3*4*NMOM + 3*NMOM + imom].real()); 
	    threep_traceless[it*4*4*NMOM + 1*4*NMOM + 1*NMOM + imom].imag() = 0.75 * threep[it*4*4*NMOM + 1*4*NMOM + 1*NMOM + imom].imag() - 0.25* ( threep[it*4*4*NMOM + 0*4*NMOM + 0*NMOM + imom].imag() + threep[it*4*4*NMOM + 2*4*NMOM + 2*NMOM + imom].imag() + threep[it*4*4*NMOM + 3*4*NMOM + 3*NMOM + imom].imag()); 
	  }
	  else if( mu == 2 && nu == 2){
	    threep_traceless[it*4*4*NMOM + 2*4*NMOM + 2*NMOM + imom].real() = 0.75 * threep[it*4*4*NMOM + 2*4*NMOM + 2*NMOM + imom].real() - 0.25* ( threep[it*4*4*NMOM + 0*4*NMOM + 0*NMOM + imom].real() + threep[it*4*4*NMOM + 1*4*NMOM + 1*NMOM + imom].real() + threep[it*4*4*NMOM + 3*4*NMOM + 3*NMOM + imom].real()); 
	    threep_traceless[it*4*4*NMOM + 2*4*NMOM + 2*NMOM + imom].imag() = 0.75 * threep[it*4*4*NMOM + 2*4*NMOM + 2*NMOM + imom].imag() - 0.25* ( threep[it*4*4*NMOM + 0*4*NMOM + 0*NMOM + imom].imag() + threep[it*4*4*NMOM + 1*4*NMOM + 1*NMOM + imom].imag() + threep[it*4*4*NMOM + 3*4*NMOM + 3*NMOM + imom].imag()); 
	  }
	  else if( mu == 3 && nu == 3){
	    threep_traceless[it*4*4*NMOM + 3*4*NMOM + 3*NMOM + imom].real() = 0.75 * threep[it*4*4*NMOM + 3*4*NMOM + 3*NMOM + imom].real() - 0.25* ( threep[it*4*4*NMOM + 0*4*NMOM + 0*NMOM + imom].real() + threep[it*4*4*NMOM + 1*4*NMOM + 1*NMOM + imom].real() + threep[it*4*4*NMOM + 2*4*NMOM + 2*NMOM + imom].real()); 
	    threep_traceless[it*4*4*NMOM + 3*4*NMOM + 3*NMOM + imom].imag() = 0.75 * threep[it*4*4*NMOM + 3*4*NMOM + 3*NMOM + imom].imag() - 0.25* ( threep[it*4*4*NMOM + 0*4*NMOM + 0*NMOM + imom].imag() + threep[it*4*4*NMOM + 1*4*NMOM + 1*NMOM + imom].imag() + threep[it*4*4*NMOM + 2*4*NMOM + 2*NMOM + imom].imag()); 
	  }
	  else{
	    threep_traceless[it*4*4*NMOM + mu*4*NMOM + nu*NMOM + imom].real() = threep[it*4*4*NMOM + mu*4*NMOM + nu*NMOM + imom].real();
	    threep_traceless[it*4*4*NMOM + mu*4*NMOM + nu*NMOM + imom].imag() = threep[it*4*4*NMOM + mu*4*NMOM + nu*NMOM + imom].imag();
	  }

	}
  
  
  for(int it=0; it < TSINK+1 ; it++)
    for(int mu = 0 ; mu < 4 ; mu++)
      for(int nu = 0 ; nu <= mu ; nu++)
        for(int imom = 0 ; imom < NMOM ; imom++)
	  fprintf(ptr_out,"%d %+d %+d %+d %+e %+e %d %d\n",it,momList[0][imom], momList[1][imom], momList[2][imom], threep_traceless[it*4*4*NMOM + mu*4*NMOM + nu*NMOM + imom].imag(),  -threep_traceless[it*4*4*NMOM + mu*4*NMOM + nu*NMOM + imom].real(), mu, nu );

  

  return 0;
}
