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

#define L 32
#define T 64
#define Z_A 0.757
#define BINSIZE 1
#define TSINK 12
#define LOW2PT 6
#define HIGH2PT 17
#define LOW3PT 2
#define HIGH3PT 10
#define AINV (0.197/0.0863)
#define AMNPHYS (0.9382720/AINV)
#define FITMAX 4
#define NMOM 257
#define MAXMOMSQ 17 // ATTENTION 7,15 missing
//#define NMOM 93
//#define MAXMOMSQ 9

#define PI 3.14159265359
#define MAX_STRING 257

// forward declaration
Complex fit_ratio_plato(Complex* ratio, Complex error_ratio[TSINK+1][4][NMOM], int fit_plato_low, int fit_plato_high, int mu, int imom);
double fit_EffMass_plato(double MEff[], double dMEff[],int fit_plato_low, int fit_plato_high);
void extractGAGP(double *GA,double *GP,Complex *RR,Complex *dRR, double *Ep,int momList[][NMOM], int p2[], int numberMomPerQ2[] );


int main(int argc, char *argv[]){

  clock_t time;

  // Greeting //
  printf("This program calculate Axial Form Factor\n");
  printf("Right inputs order must be\n");
  printf("(1) Executable , (2) Trajectory list, (3) Threep base name, (4) Twop base name, (5) Output name\n\n ");
  //////////////////


  // check passing inputs right //
  if(argc != 5){
    fprintf(stderr,"Error: Wrong number of input files \n");
    exit(EXIT_FAILURE);
  }
  //////////////////////


  // fill filenames and print paths //
  char filename_Traj[MAX_STRING];
  char filename_Threep[MAX_STRING];
  char filename_Twop[MAX_STRING];
  char filename_Output[MAX_STRING];

  //output files
  char filename_R_mom0[MAX_STRING];
  char filename_R_mom1[MAX_STRING];


  strcpy(filename_Traj,argv[1]);
  strcpy(filename_Threep,argv[2]);
  strcpy(filename_Twop,argv[3]);
  strcpy(filename_Output,argv[4]);

  printf("Got filename for trajectories : %s\n",filename_Traj);
  printf("Got filename for threep : %s\n",filename_Threep);
  printf("Got filename for twop : %s\n",filename_Twop);
  printf("Got filename for output : %s\n\n",filename_Output);

  sprintf(filename_R_mom0,"%s_%s.dat",filename_Output,"R_mom0");
  sprintf(filename_R_mom1,"%s_%s.dat",filename_Output,"R_mom1");

  FILE *out_R_mom0 = NULL, *out_R_mom1 = NULL;
  out_R_mom0 = fopen(filename_R_mom0,"w");
  out_R_mom1 = fopen(filename_R_mom1,"w");

  if(out_R_mom0 == NULL || out_R_mom1 == NULL ){ fprintf(stderr,"Error open files for writting\n"); exit(EXIT_FAILURE);}

  ///////////////////////////////////////

  // calculate and print number of configurations //
  int numLines = 0;
  std::ifstream in(filename_Traj);
  while ( in.good() )
    {
      std::string line;
      std::getline(in, line);
      ++numLines;
    }

  int Nconfs = numLines-1;
  printf("Attention trajectories file must have new line at the end to calculate right number of lines\n");
  printf("Filename with trajectories has %d configurations\n\n",Nconfs);
  //////////////////////////////////////////////////


  // read configuration's number// 
  FILE *ptr_traj;
  ptr_traj=fopen(filename_Traj,"r");

  if(ptr_traj==NULL){ fprintf(stderr,"Error open trajectories file for reading\n");exit(EXIT_FAILURE); }

  char **confString;
  confString = (char**)malloc(Nconfs*sizeof(char*));
  for(int iconf=0;iconf<Nconfs;iconf++) confString[iconf] = (char*)malloc(MAX_STRING*sizeof(char));
  int returnValue;
  for(int iconf=0;iconf<Nconfs;iconf++){
    returnValue=fscanf(ptr_traj,"%s",confString[iconf]);
  }
  fclose(ptr_traj);
  //////////////////////////////////////////////////

  // create filename names of threep and twop //

  char **threep_names;
  char **twop_names;

  threep_names = (char**)malloc(Nconfs*sizeof(char*));
  twop_names = (char**)malloc(Nconfs*sizeof(char*));

  for(int iconf=0;iconf<Nconfs;iconf++){
    threep_names[iconf] = (char*)malloc(MAX_STRING*sizeof(char));
    twop_names[iconf] = (char*)malloc(MAX_STRING*sizeof(char));
  }

  for(int iconf=0;iconf<Nconfs;iconf++){
    sprintf(threep_names[iconf],"%s.%s_isov",filename_Threep,confString[iconf]);
    sprintf(twop_names[iconf],"%s.%s",filename_Twop,confString[iconf]);
  }

  ///////////////////////////////////////////////


  // allocate memory & read data //
  Complex *threep = (Complex*)malloc(Nconfs*(TSINK+1)*4*NMOM*sizeof(Complex));
  Complex *twop = (Complex*)malloc(Nconfs*T*NMOM*sizeof(Complex));
  int momList[3][NMOM];
  int dummy;

  if(threep == NULL || twop == NULL){ fprintf(stderr,"Error: Out of memory\n"); exit(EXIT_FAILURE);}
  FILE *ptr_threep,*ptr_twop;

  time = clock();
  for(int iconf = 0 ; iconf < Nconfs ; iconf++){
    ptr_threep = fopen(threep_names[iconf],"r");
    ptr_twop = fopen(twop_names[iconf],"r");
    if(ptr_threep == NULL || ptr_twop == NULL){ fprintf(stderr,"Error: open files for reading\n"); exit(EXIT_FAILURE);}
    // read threep


    for(int it=0; it < TSINK+1 ; it++)
      for(int mu = 0 ; mu < 4 ; mu++)
	for(int imom = 0 ; imom < NMOM ; imom++)
	  returnValue = fscanf(ptr_threep,"%d %d %d %d %lf %lf %d",&dummy,&(momList[0][imom]), &(momList[1][imom]), &(momList[2][imom]), &(threep[iconf*(TSINK+1)*4*NMOM + it*4*NMOM + mu*NMOM + imom].real()), &(threep[iconf*(TSINK+1)*4*NMOM + it*4*NMOM + mu*NMOM + imom].imag()), &dummy );


    for(int it=0; it < T ; it++)
      for(int imom = 0 ; imom < NMOM ; imom++)
	returnValue = fscanf(ptr_twop,"%d %d %d %d %lf %lf",&dummy,&(momList[0][imom]), &(momList[1][imom]), &(momList[2][imom]),&(twop[iconf*T*NMOM + it*NMOM + imom].real()), &(twop[iconf*T*NMOM + it*NMOM + imom].imag()) );


    fclose(ptr_threep);
    fclose(ptr_twop);

  }

  time = clock() - time;
  fprintf(stdout,"Finished reading data in %f seconds\n",(float)time/CLOCKS_PER_SEC);
  //////////////////
  time = clock();

  // calculate number of bins //
  int Nbins = Nconfs/BINSIZE;

  printf("\nNumber of bins is: %d\n",Nbins);
  //////////////////////

  // flip the sign of momenta. Because q=p_final-p_initial, while in formulas: p_initial


  for(int imom = 0 ; imom < NMOM; imom++){
    momList[0][imom] *= -1;
    momList[1][imom] *= -1;
    momList[2][imom] *= -1;
  }

  ////////////////////////

  // aver twop with momenta giving same q^2 //
  int p2[NMOM];
  int numberMomPerQ2[MAXMOMSQ] = {};
  Complex *twop_q2 = (Complex*)calloc(Nconfs*T*MAXMOMSQ,sizeof(Complex));

  for(int imom=0;imom<NMOM;imom++)
    p2[imom]=momList[0][imom]*momList[0][imom] + momList[1][imom]*momList[1][imom] + momList[2][imom]*momList[2][imom];

  for(int imom=0;imom<NMOM;imom++)
    numberMomPerQ2[p2[imom]]++;

  for(int iconf = 0 ; iconf < Nconfs ; iconf++)
    for(int it =0 ; it < T ; it++)
      for(int imom = 0 ; imom < NMOM ; imom++){
	twop_q2[iconf*T*MAXMOMSQ + it*MAXMOMSQ + p2[imom] ] = twop_q2[iconf*T*MAXMOMSQ + it*MAXMOMSQ + p2[imom] ] + twop[iconf*T*NMOM+it*NMOM+imom];
	  }


  for(int iconf = 0 ; iconf < Nconfs ; iconf++)
    for(int it =0 ; it < T ; it++)
      for(int imom2 = 0 ; imom2 < MAXMOMSQ ; imom2++){
	if(numberMomPerQ2[imom2] != 0) twop_q2[iconf*T*MAXMOMSQ + it*MAXMOMSQ + imom2 ].real() /= numberMomPerQ2[imom2];
	if(numberMomPerQ2[imom2] != 0) twop_q2[iconf*T*MAXMOMSQ + it*MAXMOMSQ + imom2 ].imag() /= numberMomPerQ2[imom2];
      }

  /////////////////////////////////////////


  // create averages over bins for threep and twop //
  Complex *threep_mean = (Complex*)calloc((TSINK+1)*4*NMOM,sizeof(Complex));
  Complex *twop_q2_mean = (Complex*)calloc(T*MAXMOMSQ,sizeof(Complex));


  for(int it = 0 ; it < TSINK+1 ; it++)
    for(int mu = 0 ; mu < 4 ; mu++)
      for(int imom = 0 ; imom < NMOM ; imom++){
	for(int iconf = 0 ; iconf < Nconfs ; iconf++){
	  threep_mean[it*4*NMOM + mu*NMOM + imom] = threep_mean[it*4*NMOM + mu*NMOM + imom] + threep[iconf*(TSINK+1)*4*NMOM + it*4*NMOM + mu*NMOM + imom];
	    }
	threep_mean[it*4*NMOM + mu*NMOM + imom].real() /= Nconfs;
	threep_mean[it*4*NMOM + mu*NMOM + imom].imag() /= Nconfs;
      }


  for(int it = 0 ; it < T ; it++)
    for(int imom2 = 0 ; imom2 < MAXMOMSQ ; imom2++){
      for(int iconf = 0 ; iconf < Nconfs ; iconf++){
	twop_q2_mean[it*MAXMOMSQ+imom2] = twop_q2_mean[it*MAXMOMSQ+imom2] + twop_q2[iconf*T*MAXMOMSQ + it*MAXMOMSQ + imom2 ];
      }
      twop_q2_mean[it*MAXMOMSQ+imom2].real() /= Nconfs;
      twop_q2_mean[it*MAXMOMSQ+imom2].imag() /= Nconfs;
    }

  ////////////////////////////

  // calculate naive effective mass //
  double MEff_mean[T] = {};
  for(int it = 1 ; it < T ; it++){
    MEff_mean[it] = log(twop_q2_mean[it*MAXMOMSQ+0].real() / twop_q2_mean[(it+1)*MAXMOMSQ+0].real());
  }

  //////////////////////////////////////



  // calculate naive ratio //
  Complex *RR_mean = (Complex*)malloc((TSINK+1)*4*NMOM*sizeof(Complex));

  Complex squareRoot;
  double subRoot;
  for(int it=0;it < TSINK+1;it++)
    for(int mu = 0 ; mu < 4 ; mu++)
      for(int imom = 0 ; imom < NMOM ; imom++){


	subRoot = (twop_q2_mean[(TSINK-it)*MAXMOMSQ+p2[imom]].real()  * twop_q2_mean[(it)*MAXMOMSQ+0].real() *  twop_q2_mean[TSINK*MAXMOMSQ+0].real() ) /(twop_q2_mean[(TSINK-it)*MAXMOMSQ+0].real()  * twop_q2_mean[(it)*MAXMOMSQ+p2[imom]].real() *  twop_q2_mean[TSINK*MAXMOMSQ+p2[imom]].real() );
	squareRoot = std::sqrt( (Complex)  subRoot  );

	RR_mean[it*4*NMOM+mu*NMOM+imom] = ( threep_mean[it*4*NMOM + mu*NMOM + imom] / twop_q2_mean[TSINK*MAXMOMSQ+0].real() )* squareRoot;
      }

  //////////////////////////////


  // binning for jackknife analysis //
  Complex *threep_b = (Complex*)calloc(Nbins*(TSINK+1)*4*NMOM,sizeof(Complex));
  Complex *twop_q2_b = (Complex*)calloc(Nbins*T*MAXMOMSQ,sizeof(Complex));
  int listDiscardConfs[BINSIZE];
  bool checkFlag;

  for(int ibin = 0 ; ibin < Nbins ; ibin++){

    int istart = ibin*BINSIZE;
    for(int i = 0; i < BINSIZE ; i++)listDiscardConfs[i]=istart+i;


    for(int iconf = 0 ; iconf < Nconfs ; iconf++){
      checkFlag=true;

      for(int icheck = 0 ; icheck < BINSIZE ; icheck++)
	if( listDiscardConfs[icheck] == iconf ) checkFlag = false;
	 
      if(checkFlag == true){
	for(int it = 0 ; it < TSINK+1 ; it++)
	  for(int mu = 0 ; mu < 4 ; mu++)
	    for(int imom = 0 ; imom < NMOM ; imom++){

	      threep_b[ibin*(TSINK+1)*4*NMOM + it*4*NMOM + mu*NMOM + imom] = threep_b[ibin*(TSINK+1)*4*NMOM + it*4*NMOM + mu*NMOM + imom] + threep[iconf*(TSINK+1)*4*NMOM + it*4*NMOM + mu*NMOM + imom];

	    }

	for(int it = 0 ; it < T ; it++)
	  for(int imom2 = 0 ; imom2 < MAXMOMSQ; imom2++)
	    twop_q2_b[ibin*T*MAXMOMSQ+it*MAXMOMSQ+imom2] =  twop_q2_b[ibin*T*MAXMOMSQ+it*MAXMOMSQ+imom2] + twop_q2[iconf*T*MAXMOMSQ + it*MAXMOMSQ + imom2 ];
      }
   
      
    }

    for(int it = 0 ; it < TSINK+1 ; it++)
      for(int mu = 0 ; mu < 4 ; mu++)
	for(int imom = 0 ; imom < NMOM ; imom++){
	  threep_b[ibin*(TSINK+1)*4*NMOM + it*4*NMOM + mu*NMOM + imom].real() /= (Nconfs-BINSIZE);
	  threep_b[ibin*(TSINK+1)*4*NMOM + it*4*NMOM + mu*NMOM + imom].imag() /= (Nconfs-BINSIZE);
	}

	for(int it = 0 ; it < T ; it++)
	  for(int imom2 = 0 ; imom2 < MAXMOMSQ; imom2++){
	    twop_q2_b[ibin*T*MAXMOMSQ+it*MAXMOMSQ+imom2].real() /= (Nconfs-BINSIZE);
	    twop_q2_b[ibin*T*MAXMOMSQ+it*MAXMOMSQ+imom2].imag() /= (Nconfs-BINSIZE);
	  }

  }

  ///////////////////////////////////////


  // jackknife effective and ratio //
  double *MEff_b = (double*)calloc(Nbins*T,sizeof(double));
  Complex *RR_b = (Complex*)calloc(Nbins*(TSINK+1)*4*NMOM,sizeof(Complex));

  for(int ibin = 0 ; ibin < Nbins ; ibin++)
    for(int it = 1 ; it < T ; it++){
      MEff_b[ibin*T+it] = log(twop_q2_b[ibin*T*MAXMOMSQ+it*MAXMOMSQ+0].real() / twop_q2_b[ibin*T*MAXMOMSQ+(it+1)*MAXMOMSQ+0].real());
    }


  for(int ibin = 0 ; ibin < Nbins ; ibin++)
    for(int it=0;it < TSINK+1;it++)
      for(int mu = 0 ; mu < 4 ; mu++)
	for(int imom = 0 ; imom < NMOM ; imom++){


	  subRoot = (twop_q2_b[ibin*T*MAXMOMSQ+(TSINK-it)*MAXMOMSQ+p2[imom]].real()  * twop_q2_b[ibin*T*MAXMOMSQ+(it)*MAXMOMSQ+0].real() *  twop_q2_b[ibin*T*MAXMOMSQ+TSINK*MAXMOMSQ+0].real() ) /(twop_q2_b[ibin*T*MAXMOMSQ+(TSINK-it)*MAXMOMSQ+0].real()  * twop_q2_b[ibin*T*MAXMOMSQ+(it)*MAXMOMSQ+p2[imom]].real() *  twop_q2_b[ibin*T*MAXMOMSQ+TSINK*MAXMOMSQ+p2[imom]].real() );
	  squareRoot = std::sqrt( (Complex)  subRoot  );

	  RR_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM+mu*NMOM+imom] = ( threep_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM + mu*NMOM + imom] / twop_q2_b[ibin*T*MAXMOMSQ+TSINK*MAXMOMSQ+0].real() )* squareRoot;
      
	}

  ///////////////////////////////////// 

  // mean bin values //
  double MEff_bmean[T]={};
  Complex RR_bmean[TSINK+1][4][NMOM]={};

  for(int it = 0 ; it < T ; it++){
    for(int ibin = 0 ; ibin < Nbins ; ibin++){
      MEff_bmean[it] += MEff_b[ibin*T+it];
    }
    MEff_bmean[it] /= Nbins;
  }

  for(int it = 0 ; it < TSINK+1 ; it++)
    for(int mu = 0 ; mu < 4 ; mu++)
      for(int imom = 0 ; imom < NMOM ; imom++){
	for(int ibin = 0 ; ibin < Nbins ; ibin++){
	  RR_bmean[it][mu][imom] = RR_bmean[it][mu][imom] + RR_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM+mu*NMOM+imom];
	}
	RR_bmean[it][mu][imom].real() /= Nbins;
	RR_bmean[it][mu][imom].imag() /= Nbins;

      }

  ////////////////////////////

  // jacknife errors //
  Complex temp_dthreep = (Complex) {0,0};
  Complex temp_dtwop =(Complex) {0,0};
  double temp_dMEff = 0;
  Complex temp_dRR = (Complex) {0,0};
  
  Complex dthreep[TSINK+1][4][NMOM];
  Complex dtwop[T][MAXMOMSQ];
  double dMEff[T];
  Complex dRR[TSINK+1][4][NMOM];


  for(int it = 0 ; it < TSINK+1 ; it++)
    for(int mu = 0 ; mu < 4 ; mu++)
      for(int imom = 0 ; imom < NMOM ; imom++){

	for(int ibin = 0 ; ibin < Nbins ; ibin++){
	  temp_dthreep.real() = temp_dthreep.real() + ( threep_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM + mu*NMOM + imom].real() - threep_mean[it*4*NMOM + mu*NMOM + imom].real() )*(threep_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM + mu*NMOM + imom].real() - threep_mean[it*4*NMOM + mu*NMOM + imom].real() ); 
	  temp_dthreep.imag() = temp_dthreep.imag() + ( threep_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM + mu*NMOM + imom].imag() - threep_mean[it*4*NMOM + mu*NMOM + imom].imag() )*(threep_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM + mu*NMOM + imom].imag() - threep_mean[it*4*NMOM + mu*NMOM + imom].imag() ); 



	  temp_dRR.real() = temp_dRR.real() + ( RR_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM + mu*NMOM + imom].real() - RR_bmean[it][mu][imom].real() )*(RR_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM + mu*NMOM + imom].real() - RR_bmean[it][mu][imom].real() ); 
	  temp_dRR.imag() = temp_dRR.imag() + ( RR_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM + mu*NMOM + imom].imag() - RR_bmean[it][mu][imom].imag() )*(RR_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM + mu*NMOM + imom].imag() - RR_bmean[it][mu][imom].imag() ); 


	}

	dthreep[it][mu][imom].real() = sqrt( (Nbins-1)/( (double) Nbins ) * temp_dthreep.real()  ); 
	dthreep[it][mu][imom].imag() = sqrt( (Nbins-1)/( (double) Nbins ) * temp_dthreep.imag()  ); 

	dRR[it][mu][imom].real() = sqrt( (Nbins-1)/( (double) Nbins ) * temp_dRR.real()  ); 
	dRR[it][mu][imom].imag() = sqrt( (Nbins-1)/( (double) Nbins ) * temp_dRR.imag()  ); 

	temp_dthreep = (Complex) {0,0};
	temp_dRR = (Complex) {0,0};

      }

  for(int it = 0 ; it < T ; it++)
    for(int imom2=0;imom2 < MAXMOMSQ;imom2++){
      for(int ibin =0 ; ibin < Nbins ; ibin++){
	temp_dtwop.real() = temp_dtwop.real() + (twop_q2_b[ibin*T*MAXMOMSQ+it*MAXMOMSQ+imom2].real() - twop_q2_mean[it*MAXMOMSQ+imom2].real() )*( twop_q2_b[ibin*T*MAXMOMSQ+it*MAXMOMSQ+imom2].real() - twop_q2_mean[it*MAXMOMSQ+imom2].real());
	temp_dtwop.imag() = temp_dtwop.imag() + (twop_q2_b[ibin*T*MAXMOMSQ+it*MAXMOMSQ+imom2].imag() - twop_q2_mean[it*MAXMOMSQ+imom2].imag() )*( twop_q2_b[ibin*T*MAXMOMSQ+it*MAXMOMSQ+imom2].imag() - twop_q2_mean[it*MAXMOMSQ+imom2].imag());

      }

      dtwop[it][imom2].real() = sqrt((Nbins-1)/( (double) Nbins ) *  temp_dtwop.real());
      dtwop[it][imom2].imag() = sqrt((Nbins-1)/( (double) Nbins ) *  temp_dtwop.imag());
      temp_dtwop = (Complex) {0,0};
    }

  for(int it = 0 ; it < T ; it++){
    for(int ibin =0 ; ibin < Nbins ; ibin++){
      temp_dMEff = temp_dMEff + (MEff_b[ibin*T+it] - MEff_bmean[it] )*(MEff_b[ibin*T+it] - MEff_bmean[it] );
    }
    dMEff[it] = sqrt((Nbins-1)/( (double) Nbins ) * temp_dMEff);
    temp_dMEff=0;
  }

  /////////////////////////////////


  // fit to extract plateaus //


  // 1) fit to naive data
  Complex *RRplateau_mean = (Complex*)calloc(4*NMOM,sizeof(Complex));
  Complex *RRplateau_b = (Complex*)calloc(Nbins*4*NMOM,sizeof(Complex));
  Complex *RRplateau_bmean = (Complex*)calloc(4*NMOM,sizeof(Complex));



  for(int mu = 0 ; mu < 4 ; mu++)
    for(int imom = 0 ; imom < NMOM ; imom++){
      RRplateau_mean[mu*NMOM+imom] =  fit_ratio_plato(RR_mean, dRR, LOW3PT, HIGH3PT, mu, imom);
    }
  // 2) fit to jacknife data

  for(int ibin = 0 ; ibin < Nbins ; ibin++)
    for(int mu = 0 ; mu < 4 ; mu++)
      for(int imom = 0 ; imom < NMOM ; imom++){
	RRplateau_b[ibin*4*NMOM+mu*NMOM+imom] = fit_ratio_plato(RR_b+ibin*(TSINK+1)*4*NMOM, dRR,LOW3PT, HIGH3PT, mu, imom);
      }    

  for(int mu = 0 ; mu < 4 ; mu++)
    for(int imom = 0 ; imom < NMOM ; imom++){
      
      for(int ibin = 0 ; ibin < Nbins ; ibin++){
	RRplateau_bmean[mu*NMOM+imom] = RRplateau_bmean[mu*NMOM+imom] + RRplateau_b[ibin*4*NMOM+mu*NMOM+imom];
      }
      RRplateau_bmean[mu*NMOM+imom].real() /= Nbins;
      RRplateau_bmean[mu*NMOM+imom].imag() /= Nbins;
    }

  
    
  
  ////////////////////////////////////////

  // jacknife error to fitted values //
  Complex temp_dRRplateau = (Complex) {0.,0.};
  Complex *dRRplateau = (Complex*)calloc(4*NMOM,sizeof(Complex));


    for(int mu = 0 ; mu < 4 ; mu++)
      for(int imom = 0 ; imom < NMOM ; imom++){

	for(int ibin = 0 ; ibin < Nbins ; ibin++){

	  temp_dRRplateau.real() = temp_dRRplateau.real() + (RRplateau_b[ibin*4*NMOM+mu*NMOM+imom].real() - RRplateau_bmean[mu*NMOM+imom].real())*(RRplateau_b[ibin*4*NMOM+mu*NMOM+imom].real() - RRplateau_bmean[mu*NMOM+imom].real());

	  temp_dRRplateau.imag() = temp_dRRplateau.imag() + (RRplateau_b[ibin*4*NMOM+mu*NMOM+imom].imag() - RRplateau_bmean[mu*NMOM+imom].imag())*(RRplateau_b[ibin*4*NMOM+mu*NMOM+imom].imag() - RRplateau_bmean[mu*NMOM+imom].imag());
	  
	}
	
	dRRplateau[mu*NMOM+imom].real() = sqrt( (Nbins-1)/( (double) Nbins ) * temp_dRRplateau.real()  );
	dRRplateau[mu*NMOM+imom].imag() = sqrt( (Nbins-1)/( (double) Nbins ) * temp_dRRplateau.imag()  );

	temp_dRRplateau = (Complex) {0,0};

      }

  ///////////////////////////////

    // extract quantities naive //
    double mass_mean;
    double Ep_mean[MAXMOMSQ];
    mass_mean=fit_EffMass_plato(MEff_mean, dMEff,LOW2PT, HIGH2PT);

    for(int imom = 0 ; imom < NMOM ; imom++)
      Ep_mean[p2[imom]] = sqrt(mass_mean*mass_mean + p2[imom]*(2.*PI/L)*(2.*PI/L) );

    
    double GA_mean[MAXMOMSQ] = {};
    double GP_mean[MAXMOMSQ] = {};

    extractGAGP(GA_mean,GP_mean,RRplateau_mean,dRRplateau,Ep_mean,momList,p2,numberMomPerQ2);


    // extract quantities for each bin //
    double *mass_b = (double*)calloc(Nbins,sizeof(double));
    double *Ep_b = (double*)calloc(Nbins*MAXMOMSQ,sizeof(double));
    double *GA_b = (double*)calloc(Nbins*MAXMOMSQ,sizeof(double));
    double *GP_b = (double*)calloc(Nbins*MAXMOMSQ,sizeof(double));

    double *Ep_bmean = (double*)calloc(MAXMOMSQ,sizeof(double));
    double *GA_bmean = (double*)calloc(MAXMOMSQ,sizeof(double));
    double *GP_bmean = (double*)calloc(MAXMOMSQ,sizeof(double));

    for(int ibin = 0 ; ibin < Nbins ; ibin++)
      mass_b[ibin] = fit_EffMass_plato(MEff_b+ibin*T, dMEff,LOW2PT, HIGH2PT);

    for(int ibin = 0 ; ibin < Nbins ; ibin++)
      for(int imom = 0 ; imom < NMOM ; imom++)
	Ep_b[ibin*MAXMOMSQ+p2[imom]] = sqrt(mass_b[ibin]*mass_b[ibin] + p2[imom]*(2.*PI/L)*(2.*PI/L) );

    for(int ibin = 0 ; ibin < Nbins ; ibin++){
      extractGAGP(GA_b + ibin*MAXMOMSQ,GP_b + ibin*MAXMOMSQ,RRplateau_b + ibin*4*NMOM ,dRRplateau,Ep_b + ibin*MAXMOMSQ,momList,p2,numberMomPerQ2);
    }    
    ////////////////////////////////////

    // find jackknife errors for form factors //
    double temp_Ep_bmean=0.;
    double temp_GA_bmean=0.;
    double temp_GP_bmean=0.;

    double *dGA = (double*)calloc(MAXMOMSQ,sizeof(double));
    double *dGP = (double*)calloc(MAXMOMSQ,sizeof(double));
    double *dEp = (double*)calloc(MAXMOMSQ,sizeof(double));

    for(int imom2 = 0 ; imom2 < MAXMOMSQ ; imom2++){
      for(int ibin = 0 ; ibin < Nbins ; ibin++){
	Ep_bmean[imom2] += Ep_b[ibin*MAXMOMSQ+imom2];
	GA_bmean[imom2] += GA_b[ibin*MAXMOMSQ+imom2];
	GP_bmean[imom2] += GP_b[ibin*MAXMOMSQ+imom2];
      }
      Ep_bmean[imom2] /= Nbins;
      GA_bmean[imom2] /= Nbins;
      GP_bmean[imom2] /= Nbins;
    }

    for(int imom2 = 0 ; imom2 < MAXMOMSQ ; imom2++){
      for(int ibin = 0 ; ibin < Nbins ; ibin++){
	temp_Ep_bmean += (Ep_bmean[imom2] - Ep_b[ibin*MAXMOMSQ+imom2])*(Ep_bmean[imom2] - Ep_b[ibin*MAXMOMSQ+imom2]); 
	temp_GA_bmean += (GA_bmean[imom2] - GA_b[ibin*MAXMOMSQ+imom2])*(GA_bmean[imom2] - GA_b[ibin*MAXMOMSQ+imom2]); 
	temp_GP_bmean += (GP_bmean[imom2] - GP_b[ibin*MAXMOMSQ+imom2])*(GP_bmean[imom2] - GP_b[ibin*MAXMOMSQ+imom2]); 
      }
      dGA[imom2] = sqrt( (Nbins-1)/( (double) Nbins ) * temp_GA_bmean  );
      dGP[imom2] = sqrt( (Nbins-1)/( (double) Nbins ) * temp_GP_bmean  );
      dEp[imom2] = sqrt( (Nbins-1)/( (double) Nbins ) * temp_Ep_bmean  );

      temp_Ep_bmean=0.;
      temp_GA_bmean=0.;
      temp_GP_bmean=0.;
    }
    ////////////////////////////////////////////

    // calculate data for ratio to print out //

    //printOut R zero momentum

    for(int it = 0 ; it < TSINK+1 ; it++)
      fprintf(out_R_mom0,"%d %+.8f %+.8f\n",it, RR_bmean[it][0][0].real() , dRR[it][0][0].real());


    /*
    // printOut RE momentum 1
    if(MAXMOMSQ > 1){
      double C_kinematics;
      double coeff[12];
      double *value = (double*)malloc(Nbins*(TSINK+1)*12*sizeof(double));
      double *dirValue = (double*)calloc(Nbins*(TSINK+1),sizeof(double));
      double dirMeanValue[TSINK+1]={0};

      C_kinematics = sqrt( ( 2.*Ep_mean[0]*Ep_mean[0] ) / ( Ep_mean[1] * ( Ep_mean[1]+Ep_mean[0] ) ) );
      for(int i = 0 ; i < 6 ; i++)
	coeff[i] = ( C_kinematics * (Ep_mean[1]+Ep_mean[0]) )/(2.*Ep_mean[0]);
      for(int i = 6 ; i < 8 ; i++)
	coeff[i] = ( C_kinematics/(2.*Ep_mean[0])) * (-momList[0][i-5]*(2.*PI/L));
      for(int i = 8 ; i < 10 ; i++)
	coeff[i] = ( C_kinematics/(2.*Ep_mean[0])) * (-momList[1][i-5]*(2.*PI/L));
      for(int i = 10 ; i < 12 ; i++)
	coeff[i] = ( C_kinematics/(2.*Ep_mean[0])) * (-momList[2][i-5]*(2.*PI/L));
      
      for(int i = 0 ; i < 6 ; i++) 
	for(int it = 0 ; it < TSINK+1 ; it++)
	  for(int ibin = 0; ibin < Nbins; ibin++)
	    value[ibin*(TSINK+1)*12+it*12+i] = Z_A*RRE_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM+0*NMOM+(i+1)].real();

	for(int it = 0 ; it < TSINK+1 ; it++)
	  for(int ibin = 0; ibin < Nbins; ibin++){
	    value[ibin*(TSINK+1)*12+it*12+6] = Z_A*RRE_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM+1*NMOM+1].imag();
	    value[ibin*(TSINK+1)*12+it*12+7] = Z_A*RRE_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM+1*NMOM+2].imag();
	    value[ibin*(TSINK+1)*12+it*12+8] = Z_A*RRE_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM+2*NMOM+3].imag();
	    value[ibin*(TSINK+1)*12+it*12+9] = Z_A*RRE_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM+2*NMOM+4].imag();
	    value[ibin*(TSINK+1)*12+it*12+10] = Z_A*RRE_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM+3*NMOM+5].imag();
	    value[ibin*(TSINK+1)*12+it*12+11] = Z_A*RRE_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM+3*NMOM+6].imag();
	  }

	for(int it = 0 ; it < TSINK+1 ; it++)
          for(int ibin = 0; ibin < Nbins; ibin++){
	    for(int i = 0 ; i < 12 ; i++)
	      dirValue[ibin*(TSINK+1)+it] += value[ibin*(TSINK+1)*12+it*12+i]/coeff[i];
	    dirValue[ibin*(TSINK+1)+it] /= 12;
	  }

	for(int it = 0 ; it < TSINK+1 ; it++){
	  for(int ibin = 0; ibin < Nbins; ibin++)
	    dirMeanValue[it] += dirValue[ibin*(TSINK+1)+it];
	  dirMeanValue[it] /= Nbins;
	}

	double tempValue=0;
	double dirMeanValue_err[TSINK+1];

	for(int it = 0 ; it < TSINK+1 ; it++){
	  for(int ibin = 0 ; ibin < Nbins ; ibin++)
	    tempValue += (dirMeanValue[it] - dirValue[ibin*(TSINK+1)+it])*(dirMeanValue[it] - dirValue[ibin*(TSINK+1)+it]);
	  dirMeanValue_err[it] = sqrt( (Nbins-1)/( (double) Nbins ) * tempValue  );
	  tempValue=0.;
	}

	for(int it = 0 ; it < TSINK+1 ; it++)
	  fprintf(out_RE_mom1,"%d %+.8f %+.8f\n",it, dirMeanValue[it] , dirMeanValue_err[it]);


	// printOut RM momentum 1
	// warning : this works only for the existing convention for the order of momenta
	coeff[0] = ( C_kinematics/(2.*Ep_mean[0])) * (1*(2.*PI/L));	
	coeff[1] = ( C_kinematics/(2.*Ep_mean[0])) * (-1*(2.*PI/L));	
	coeff[2] = ( C_kinematics/(2.*Ep_mean[0])) * (-1*(2.*PI/L));	
	coeff[3] = ( C_kinematics/(2.*Ep_mean[0])) * (1*(2.*PI/L));	
	coeff[4] = ( C_kinematics/(2.*Ep_mean[0])) * (-1*(2.*PI/L));	
	coeff[5] = ( C_kinematics/(2.*Ep_mean[0])) * (+1*(2.*PI/L));	
	coeff[6] = ( C_kinematics/(2.*Ep_mean[0])) * (+1*(2.*PI/L));	
	coeff[7] = ( C_kinematics/(2.*Ep_mean[0])) * (-1*(2.*PI/L));	
	coeff[8] = ( C_kinematics/(2.*Ep_mean[0])) * (1*(2.*PI/L));	
	coeff[9] = ( C_kinematics/(2.*Ep_mean[0])) * (-1*(2.*PI/L));	
	coeff[10] = ( C_kinematics/(2.*Ep_mean[0])) * (-1*(2.*PI/L));	
	coeff[11] = ( C_kinematics/(2.*Ep_mean[0])) * (1*(2.*PI/L));	

	for(int it = 0 ; it < TSINK+1 ; it++)
	  for(int ibin = 0; ibin < Nbins; ibin++){	
	    value[ibin*(TSINK+1)*12+it*12+0] = Z_A*RRM_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM+1*NMOM+3].real();
	    value[ibin*(TSINK+1)*12+it*12+1] = Z_A*RRM_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM+1*NMOM+4].real();
	    value[ibin*(TSINK+1)*12+it*12+2] = Z_A*RRM_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM+1*NMOM+5].real();
	    value[ibin*(TSINK+1)*12+it*12+3] = Z_A*RRM_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM+1*NMOM+6].real();
	    value[ibin*(TSINK+1)*12+it*12+4] = Z_A*RRM_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM+2*NMOM+1].real();
	    value[ibin*(TSINK+1)*12+it*12+5] = Z_A*RRM_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM+2*NMOM+2].real();
	    value[ibin*(TSINK+1)*12+it*12+6] = Z_A*RRM_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM+2*NMOM+5].real();
	    value[ibin*(TSINK+1)*12+it*12+7] = Z_A*RRM_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM+2*NMOM+6].real();
	    value[ibin*(TSINK+1)*12+it*12+8] = Z_A*RRM_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM+3*NMOM+1].real();
	    value[ibin*(TSINK+1)*12+it*12+9] = Z_A*RRM_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM+3*NMOM+2].real();
	    value[ibin*(TSINK+1)*12+it*12+10] = Z_A*RRM_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM+3*NMOM+3].real();
	    value[ibin*(TSINK+1)*12+it*12+11] = Z_A*RRM_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM+3*NMOM+4].real();
	  }
	memset(dirValue,0,Nbins*(TSINK+1)*sizeof(double));
	memset(dirMeanValue,0,(TSINK+1)*sizeof(double));

	for(int it = 0 ; it < TSINK+1 ; it++)
          for(int ibin = 0; ibin < Nbins; ibin++){
	    for(int i = 0 ; i < 12 ; i++)
	      dirValue[ibin*(TSINK+1)+it] += value[ibin*(TSINK+1)*12+it*12+i]/coeff[i];
	    dirValue[ibin*(TSINK+1)+it] /= 12;
	  }

	for(int it = 0 ; it < TSINK+1 ; it++){
	  for(int ibin = 0; ibin < Nbins; ibin++)
	    dirMeanValue[it] += dirValue[ibin*(TSINK+1)+it];
	  dirMeanValue[it] /= Nbins;
	}
	tempValue=0.;
	for(int it = 0 ; it < TSINK+1 ; it++){
	  for(int ibin = 0 ; ibin < Nbins ; ibin++)
	    tempValue += (dirMeanValue[it] - dirValue[ibin*(TSINK+1)+it])*(dirMeanValue[it] - dirValue[ibin*(TSINK+1)+it]);
	  dirMeanValue_err[it] = sqrt( (Nbins-1)/( (double) Nbins ) * tempValue  );
	  tempValue=0.;
	}

	for(int it = 0 ; it < TSINK+1 ; it++)
	  fprintf(out_RM_mom1,"%d %+.8f %+.8f\n",it, dirMeanValue[it] , dirMeanValue_err[it]);


	free(value);
	free(dirValue);

    }
    */



    // print results //
    double q2GeV[MAXMOMSQ];

    for(int imom2 = 0 ; imom2 < MAXMOMSQ; imom2++){
      if(numberMomPerQ2[imom2] != 0)
	q2GeV[imom2] = (2.*Ep_mean[imom2]*mass_mean - 2.*mass_mean*mass_mean)*AINV*AINV;
      else
	q2GeV[imom2] = log(-1);
    }
    time=clock()-time;
    printf("Analysis took %f seconds\n",(float)time/CLOCKS_PER_SEC);
    printf("\n\n RESULTS \n\n");
    
    printf("\nGA FORM FACTOR\n");
    for(int imom2 = 0 ; imom2 < MAXMOMSQ; imom2++)
      printf("%f %f %f\n",q2GeV[imom2],GA_bmean[imom2],dGA[imom2]);

    printf("\nGP FORM FACTOR\n");
    for(int imom2 = 0 ; imom2 < MAXMOMSQ; imom2++)
      printf("%f %f %f\n",q2GeV[imom2],GP_bmean[imom2],dGP[imom2]);

    //////////////////////////


    // free memory //
    for(int iconf=0;iconf<Nconfs;iconf++){
      free(confString[iconf]);
      free(threep_names[iconf]);
      free(twop_names[iconf]);
    }
    free(confString);
    free(threep_names);
    free(twop_names);
    free(threep);
    free(twop);
    free(RR_mean);
    free(twop_q2);
    free(threep_mean);
    free(twop_q2_mean);
    free(threep_b);
    free(twop_q2_b);
    free(MEff_b);
    free(RR_b);
    free(RRplateau_mean);
    free(RRplateau_b);
    free(RRplateau_bmean);
    free(dRRplateau);
    free(mass_b);
    free(Ep_b);
    free(GA_b);
    free(GP_b);
    free(Ep_bmean);
    free(GA_bmean);
    free(GP_bmean);
    free(dGA);
    free(dGP);
    free(dEp);

    ///////////////////////

  return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////   FUNCTIONS ///////////////////////////////////////////////////////////////////////////////////////////////

Complex fit_ratio_plato(Complex* ratio, Complex error_ratio[TSINK+1][4][NMOM], int fit_plato_low, int fit_plato_high, int mu, int imom){
  int fitrange;
  Complex yi[fit_plato_high - fit_plato_low + 1];
  Complex sigmai[fit_plato_high - fit_plato_low + 1];
  Complex S , Sy;
  Complex par_constant;

  fitrange = fit_plato_high - fit_plato_low + 1;

  for(int it = 0; it < fitrange ; it++){
    yi[it] = ratio[(fit_plato_low+it)*4*NMOM+mu*NMOM+imom] ;
    sigmai[it] = error_ratio[fit_plato_low+it][mu][imom] ;

  }

  S =(Complex) {0.,0.};
  Sy =(Complex) {0.,0.};

  for(int it = 0 ; it < fitrange ; it++){
    Sy.real() += yi[it].real() / ( sigmai[it].real() * sigmai[it].real() );
    S.real() += 1. / ( sigmai[it].real() * sigmai[it].real() );

    Sy.imag() += yi[it].imag() / ( sigmai[it].imag() * sigmai[it].imag() );
    S.imag() += 1. / ( sigmai[it].imag() * sigmai[it].imag() );
  }

  par_constant.real() = Sy.real() / S.real() ;
  par_constant.imag() = Sy.imag() / S.imag() ;

  return(par_constant);
}

double fit_EffMass_plato(double MEff[], double dMEff[],int fit_plato_low, int fit_plato_high){

  int fitrange;
  double yi[fit_plato_high - fit_plato_low + 1];
  double sigmai[fit_plato_high - fit_plato_low + 1];
  double S , Sy;
  double par_constant;

  fitrange = fit_plato_high - fit_plato_low + 1;

  for(int it = 0; it < fitrange ; it++){
    yi[it] = MEff[fit_plato_low+it] ;
    sigmai[it] = dMEff[fit_plato_low+it] ;

  }

  S =0.;
  Sy =0.;

  for(int it = 0 ; it < fitrange ; it++){
    Sy += yi[it] / ( sigmai[it] * sigmai[it] );
    S += 1. / ( sigmai[it] * sigmai[it] );
  }

  par_constant = Sy / S ;

  return(par_constant);

}

void extractGAGP(double *GA,double *GP,Complex *RR,Complex *dRR, double *Ep,int momList[][NMOM], int p2[], int numberMomPerQ2[] ){

  double C; // kinematics term
  Complex *CoeffGA[MAXMOMSQ] = {NULL};
  Complex *CoeffGP[MAXMOMSQ] = {NULL};
  int startMom = 0;

  double *b[MAXMOMSQ] = {NULL};
  double *db[MAXMOMSQ] = {NULL};
  double *A[MAXMOMSQ] = {NULL}; // important because of linkage with fortran matrices rows and columns are stored different in memory
  double *u[MAXMOMSQ] = {NULL};
  double *w[MAXMOMSQ] = {NULL};
  double *vt[MAXMOMSQ] = {NULL};
  double *oneOverW[MAXMOMSQ] = {NULL};
  double *x[MAXMOMSQ] = {NULL};

  for(int imom2 = 0 ; imom2 < MAXMOMSQ ; imom2++){
    GA[imom2]=0.;
    GP[imom2]=0.;

    if(numberMomPerQ2[imom2] != 0){
      C = sqrt( ( 2.*Ep[0]*Ep[0] ) / ( Ep[imom2] * ( Ep[imom2]+Ep[0] ) ) );    
      
      CoeffGA[imom2] = (Complex*)calloc(numberMomPerQ2[imom2]*3,sizeof(Complex));
      CoeffGP[imom2] = (Complex*)calloc(numberMomPerQ2[imom2]*3,sizeof(Complex));

      b[imom2] = (double*)calloc(3*numberMomPerQ2[imom2],sizeof(double));
      db[imom2] = (double*)calloc(3*numberMomPerQ2[imom2],sizeof(double));
      A[imom2] = (double*)calloc(2*3*numberMomPerQ2[imom2],sizeof(double)); // dimension (# measurements) x 2 // two form factors

      for(int imom = 0 ; imom < NMOM ; imom++){
	int q2 = momList[0][imom]*momList[0][imom] + momList[1][imom]*momList[1][imom] + momList[2][imom]*momList[2][imom];
	if(q2 == imom2){
	  startMom = imom; // get the position which start each momentum square
	  break;
	}

      }
      
      for(int mu = 1 ; mu < 4 ; mu++)
	for(int i = 0 ; i < numberMomPerQ2[imom2] ; i++){
	  CoeffGA[imom2][(mu-1)*numberMomPerQ2[imom2] + i ].real() = -(C*(Ep[imom2]+Ep[0]))/(2.*Ep[0]); 
	  CoeffGP[imom2][(mu-1)*numberMomPerQ2[imom2] + i].real() = ( C/(4.*Ep[0]*Ep[0]) ) * ( momList[mu-1][startMom+i] *(momList[0][startMom+i] + momList[1][startMom+i] + momList[2][startMom+i])) * (2.*PI/L) * (2.*PI/L);  
	 	  
	}      

      for(int mu = 1 ; mu < 4 ; mu++)
	for(int i = 0 ; i < numberMomPerQ2[imom2] ; i++){

	  b[imom2][(mu-1)*numberMomPerQ2[imom2] + i ] = imag(RR[mu*NMOM + startMom + i]);        // sigma is in the imaginary part
	  db[imom2][(mu-1)*numberMomPerQ2[imom2] + i ] = imag(dRR[mu*NMOM + startMom + i]);
	  A[imom2][0*3*numberMomPerQ2[imom2]+(mu-1)*numberMomPerQ2[imom2] + i] = real(CoeffGA[imom2][(mu-1)*numberMomPerQ2[imom2] + i ]);
	  A[imom2][1*3*numberMomPerQ2[imom2]+(mu-1)*numberMomPerQ2[imom2] + i] = real(CoeffGP[imom2][(mu-1)*numberMomPerQ2[imom2] + i ]);
	}	

      /*
      if(imom2 == 12){
	
	for(int mu = 1 ; mu < 4 ; mu++)
	  for(int i = 0 ; i < numberMomPerQ2[imom2] ; i++){
	    printf("%f %f\n",b[imom2][(mu-1)*numberMomPerQ2[imom2] + i],db[imom2][(mu-1)*numberMomPerQ2[imom2] + i]);
	  }

      }
      */
      //      make weighted least square solution: 
      for(int mu = 1 ; mu < 4 ; mu++)
	for(int i = 0 ; i < numberMomPerQ2[imom2] ; i++){
	  A[imom2][0*3*numberMomPerQ2[imom2]+(mu-1)*numberMomPerQ2[imom2] + i] = A[imom2][0*3*numberMomPerQ2[imom2]+(mu-1)*numberMomPerQ2[imom2] + i]/ fabs(db[imom2][(mu-1)*numberMomPerQ2[imom2] + i ]);
	  A[imom2][1*3*numberMomPerQ2[imom2]+(mu-1)*numberMomPerQ2[imom2] + i] = A[imom2][1*3*numberMomPerQ2[imom2]+(mu-1)*numberMomPerQ2[imom2] + i]/ fabs(db[imom2][(mu-1)*numberMomPerQ2[imom2] + i ]);
	  b[imom2][(mu-1)*numberMomPerQ2[imom2] + i] = b[imom2][(mu-1)*numberMomPerQ2[imom2] + i] /fabs(db[imom2][(mu-1)*numberMomPerQ2[imom2] + i]);
	}
      
      // ready for singular value decomposition

      
      int info;
      int lwork;
      double wkopt;
      double* work=NULL;
      int m=numberMomPerQ2[imom2]*3;
      int n=2;
      int lda=m;
      int ldu=m;
      int ldvt=n;

      u[imom2]=(double*)malloc(ldu*m*sizeof(double));
      w[imom2]=(double*)malloc(n*sizeof(double));
      vt[imom2]=(double*)malloc(ldvt*n*sizeof(double));
      oneOverW[imom2] =(double*)calloc(n*n,sizeof(double));
      x[imom2]=(double*)calloc(n,sizeof(double));

      lwork = -1;
      dgesvd_((char*) "All",(char*) "All", &m, &n, A[imom2], &lda, w[imom2], u[imom2], &ldu, vt[imom2], &ldvt, &wkopt, &lwork, &info );
      lwork = (int)wkopt;
      work = (double*)malloc( lwork*sizeof(double) );
      // Compute SVD 
      dgesvd_((char*) "All",(char*) "All", &m, &n, A[imom2], &lda, w[imom2], u[imom2], &ldu, vt[imom2], &ldvt, work, &lwork, &info );
      // Check for convergence 
      if( info > 0 ) {
	printf( "The algorithm computing SVD failed to converge.\n" );
	exit( 1 );
      }     
      
      oneOverW[imom2][0*2+0] = 1./w[imom2][0];
      oneOverW[imom2][1*2+1] = 1./w[imom2][1];
      if(imom2 == 0) oneOverW[imom2][1*2+1] = 0.;

      for(int alpha = 0 ; alpha < 2 ; alpha++)
	for(int beta = 0 ; beta < 2 ; beta++)
	  for(int gamma = 0 ; gamma < 2 ; gamma++)
	    for(int delta = 0 ; delta < numberMomPerQ2[imom2]*3 ; delta++){
	      x[imom2][alpha] += vt[imom2][alpha*2+beta] * oneOverW[imom2][beta*2+gamma] * u[imom2][gamma*m+delta] * b[imom2][delta]; 
	    }


      GA[imom2] = x[imom2][0] * Z_A;
      GP[imom2] = (imom2 == 0) ? log(-1) : x[imom2][1] * Z_A;
    }
    else{
            GA[imom2] = log(-1); // assing nan
            GP[imom2] = log(-1); // assing nan
    }



    free(CoeffGA[imom2]);
    free(CoeffGP[imom2]);
    free(b[imom2]);
    free(db[imom2]);
    free(A[imom2]);
    free(u[imom2]);
    free(w[imom2]);
    free(vt[imom2]);
    free(oneOverW[imom2]);
    free(x[imom2]);

    }

  }

