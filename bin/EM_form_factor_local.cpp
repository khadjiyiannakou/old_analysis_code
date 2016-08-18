/*
Kyriakos Hadjiyiannakou
Extract Electromagnetic form factor
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
3) extractGE
4) extractGM

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
#define Z_V 0.611
#define BINSIZE 2
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
#define PI 3.14159265359
#define MAX_STRING 257

// forward declaration
Complex fit_ratio_plato(Complex* ratio, Complex error_ratio[TSINK+1][4][NMOM], int fit_plato_low, int fit_plato_high, int mu, int imom);
double fit_EffMass_plato(double MEff[], double dMEff[],int fit_plato_low, int fit_plato_high);
void extractGE(double *GE,Complex *RRE,Complex *dRRE, double *Ep,int momList[][NMOM], int p2[], int numberMomPerQ2[] );
void extractGM(double *GM,Complex *RRM,Complex *dRRM, double *Ep,int momList[][NMOM], int p2[], int numberMomPerQ2[] );

int main(int argc, char *argv[]){

  clock_t time;

  // Greeting //
  printf("This program calculate Electromagnetic Form Factor\n");
  printf("Right inputs order must be\n");
  printf("(1) Executable , (2) Trajectory list, (3) Threep Electric base name, (4) Threep Magnetic base name, (5) Twop base name, (6) Output name\n\n ");
  //////////////////


  // check passing inputs right //
  if(argc != 6){
    fprintf(stderr,"Error: Wrong number of input files \n");
    exit(EXIT_FAILURE);
  }
  //////////////////////


  // fill filenames and print paths //
  char filename_Traj[MAX_STRING];
  char filename_Threep_Electric[MAX_STRING];
  char filename_Threep_Magnetic[MAX_STRING];
  char filename_Twop[MAX_STRING];
  char filename_Output[MAX_STRING];

  //output files
  char filename_RE_mom0[MAX_STRING];
  char filename_RE_mom1[MAX_STRING];
  char filename_RM_mom1[MAX_STRING];

  strcpy(filename_Traj,argv[1]);
  strcpy(filename_Threep_Electric,argv[2]);
  strcpy(filename_Threep_Magnetic,argv[3]);
  strcpy(filename_Twop,argv[4]);
  strcpy(filename_Output,argv[5]);

  printf("Got filename for trajectories : %s\n",filename_Traj);
  printf("Got filename for threep Electric : %s\n",filename_Threep_Electric);
  printf("Got filename for threep Magnetic : %s\n",filename_Threep_Magnetic);
  printf("Got filename for twop : %s\n",filename_Twop);
  printf("Got filename for output : %s\n\n",filename_Output);

  sprintf(filename_RE_mom0,"%s_%s.dat",filename_Output,"RE_mom0");
  sprintf(filename_RE_mom1,"%s_%s.dat",filename_Output,"RE_mom1");
  sprintf(filename_RM_mom1,"%s_%s.dat",filename_Output,"RM_mom1");

  FILE *out_RE_mom0 = NULL, *out_RE_mom1 = NULL, *out_RM_mom1 = NULL;
  out_RE_mom0 = fopen(filename_RE_mom0,"w");
  out_RE_mom1 = fopen(filename_RE_mom1,"w");
  out_RM_mom1 = fopen(filename_RM_mom1,"w");
  if(out_RE_mom0 == NULL || out_RE_mom1 == NULL || out_RM_mom1 == NULL){ fprintf(stderr,"Error open files for writting\n"); exit(EXIT_FAILURE);}

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

  char **threep_Electric_names;
  char **threep_Magnetic_names;
  char **twop_names;

  threep_Electric_names = (char**)malloc(Nconfs*sizeof(char*));
  threep_Magnetic_names = (char**)malloc(Nconfs*sizeof(char*));
  twop_names = (char**)malloc(Nconfs*sizeof(char*));

  for(int iconf=0;iconf<Nconfs;iconf++){
    threep_Electric_names[iconf] = (char*)malloc(MAX_STRING*sizeof(char));
    threep_Magnetic_names[iconf] = (char*)malloc(MAX_STRING*sizeof(char));
    twop_names[iconf] = (char*)malloc(MAX_STRING*sizeof(char));
  }

  for(int iconf=0;iconf<Nconfs;iconf++){
    sprintf(threep_Electric_names[iconf],"%s.%s_isov",filename_Threep_Electric,confString[iconf]);
    sprintf(threep_Magnetic_names[iconf],"%s.%s_isov",filename_Threep_Magnetic,confString[iconf]);
    sprintf(twop_names[iconf],"%s.%s",filename_Twop,confString[iconf]);
  }

  ///////////////////////////////////////////////


  // allocate memory & read data //
  Complex *threep_Electric = (Complex*)malloc(Nconfs*(TSINK+1)*4*NMOM*sizeof(Complex));
  Complex *threep_Magnetic = (Complex*)malloc(Nconfs*(TSINK+1)*4*NMOM*sizeof(Complex));
  Complex *twop = (Complex*)malloc(Nconfs*T*NMOM*sizeof(Complex));
  int momList[3][NMOM];
  int dummy;

  if(threep_Electric == NULL || threep_Magnetic == NULL || twop == NULL){ fprintf(stderr,"Error: Out of memory\n"); exit(EXIT_FAILURE);}
  FILE *ptr_threep_electric,*ptr_threep_magnetic,*ptr_twop;

  time = clock();
  for(int iconf = 0 ; iconf < Nconfs ; iconf++){
    ptr_threep_electric = fopen(threep_Electric_names[iconf],"r");
    ptr_threep_magnetic = fopen(threep_Magnetic_names[iconf],"r");
    ptr_twop = fopen(twop_names[iconf],"r");
    if(ptr_threep_electric == NULL || ptr_threep_magnetic == NULL || ptr_twop == NULL){ fprintf(stderr,"Error: open files for reading\n"); exit(EXIT_FAILURE);}
    // read threep


    for(int it=0; it < TSINK+1 ; it++)
      for(int mu = 0 ; mu < 4 ; mu++)
	for(int imom = 0 ; imom < NMOM ; imom++)
	  returnValue = fscanf(ptr_threep_electric,"%d %d %d %d %lf %lf %d",&dummy,&(momList[0][imom]), &(momList[1][imom]), &(momList[2][imom]), &(threep_Electric[iconf*(TSINK+1)*4*NMOM + it*4*NMOM + mu*NMOM + imom].real()), &(threep_Electric[iconf*(TSINK+1)*4*NMOM + it*4*NMOM + mu*NMOM + imom].imag()), &dummy );

    for(int it=0; it < TSINK+1 ; it++)
      for(int mu = 0 ; mu < 4 ; mu++)
	for(int imom = 0 ; imom < NMOM ; imom++)
	  returnValue = fscanf(ptr_threep_magnetic,"%d %d %d %d %lf %lf %d",&dummy,&(momList[0][imom]), &(momList[1][imom]), &(momList[2][imom]), &(threep_Magnetic[iconf*(TSINK+1)*4*NMOM + it*4*NMOM + mu*NMOM + imom].real()), &(threep_Magnetic[iconf*(TSINK+1)*4*NMOM + it*4*NMOM + mu*NMOM + imom].imag()), &dummy );

    for(int it=0; it < T ; it++)
      for(int imom = 0 ; imom < NMOM ; imom++)
	returnValue = fscanf(ptr_twop,"%d %d %d %d %lf %lf",&dummy,&(momList[0][imom]), &(momList[1][imom]), &(momList[2][imom]),&(twop[iconf*T*NMOM + it*NMOM + imom].real()), &(twop[iconf*T*NMOM + it*NMOM + imom].imag()) );


    fclose(ptr_threep_electric);
    fclose(ptr_threep_magnetic);
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
  Complex *threep_Electric_mean = (Complex*)calloc((TSINK+1)*4*NMOM,sizeof(Complex));
  Complex *threep_Magnetic_mean = (Complex*)calloc((TSINK+1)*4*NMOM,sizeof(Complex));
  Complex *twop_q2_mean = (Complex*)calloc(T*MAXMOMSQ,sizeof(Complex));

  for(int it = 0 ; it < TSINK+1 ; it++)
    for(int mu = 0 ; mu < 4 ; mu++)
      for(int imom = 0 ; imom < NMOM ; imom++){
	for(int iconf = 0 ; iconf < Nconfs ; iconf++){
	  threep_Electric_mean[it*4*NMOM + mu*NMOM + imom] = threep_Electric_mean[it*4*NMOM + mu*NMOM + imom] + threep_Electric[iconf*(TSINK+1)*4*NMOM + it*4*NMOM + mu*NMOM + imom];
	  threep_Magnetic_mean[it*4*NMOM + mu*NMOM + imom] = threep_Magnetic_mean[it*4*NMOM + mu*NMOM + imom] + threep_Magnetic[iconf*(TSINK+1)*4*NMOM + it*4*NMOM + mu*NMOM + imom];
	    }
	threep_Electric_mean[it*4*NMOM + mu*NMOM + imom].real() /= Nconfs;
	threep_Electric_mean[it*4*NMOM + mu*NMOM + imom].imag() /= Nconfs;

	threep_Magnetic_mean[it*4*NMOM + mu*NMOM + imom].real() /= Nconfs;
	threep_Magnetic_mean[it*4*NMOM + mu*NMOM + imom].imag() /= Nconfs;
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
  Complex *RRE_mean = (Complex*)malloc((TSINK+1)*4*NMOM*sizeof(Complex));
  Complex *RRM_mean = (Complex*)malloc((TSINK+1)*4*NMOM*sizeof(Complex));
  Complex squareRoot;
  double subRoot;
  for(int it=0;it < TSINK+1;it++)
    for(int mu = 0 ; mu < 4 ; mu++)
      for(int imom = 0 ; imom < NMOM ; imom++){


	subRoot = (twop_q2_mean[(TSINK-it)*MAXMOMSQ+p2[imom]].real()  * twop_q2_mean[(it)*MAXMOMSQ+0].real() *  twop_q2_mean[TSINK*MAXMOMSQ+0].real() ) /(twop_q2_mean[(TSINK-it)*MAXMOMSQ+0].real()  * twop_q2_mean[(it)*MAXMOMSQ+p2[imom]].real() *  twop_q2_mean[TSINK*MAXMOMSQ+p2[imom]].real() );
	squareRoot = std::sqrt( (Complex)  subRoot  );

	RRE_mean[it*4*NMOM+mu*NMOM+imom] = ( threep_Electric_mean[it*4*NMOM + mu*NMOM + imom] / twop_q2_mean[TSINK*MAXMOMSQ+0].real() )* squareRoot;
	RRM_mean[it*4*NMOM+mu*NMOM+imom] = ( threep_Magnetic_mean[it*4*NMOM + mu*NMOM + imom] / twop_q2_mean[TSINK*MAXMOMSQ+0].real() )* squareRoot;
	

      }

  //////////////////////////////


  // binning for jackknife analysis //
  Complex *threepE_b = (Complex*)calloc(Nbins*(TSINK+1)*4*NMOM,sizeof(Complex));
  Complex *threepM_b = (Complex*)calloc(Nbins*(TSINK+1)*4*NMOM,sizeof(Complex));
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

	      threepE_b[ibin*(TSINK+1)*4*NMOM + it*4*NMOM + mu*NMOM + imom] = threepE_b[ibin*(TSINK+1)*4*NMOM + it*4*NMOM + mu*NMOM + imom] + threep_Electric[iconf*(TSINK+1)*4*NMOM + it*4*NMOM + mu*NMOM + imom];
	      threepM_b[ibin*(TSINK+1)*4*NMOM + it*4*NMOM + mu*NMOM + imom] = threepM_b[ibin*(TSINK+1)*4*NMOM + it*4*NMOM + mu*NMOM + imom] + threep_Magnetic[iconf*(TSINK+1)*4*NMOM + it*4*NMOM + mu*NMOM + imom];
	    }

	for(int it = 0 ; it < T ; it++)
	  for(int imom2 = 0 ; imom2 < MAXMOMSQ; imom2++)
	    twop_q2_b[ibin*T*MAXMOMSQ+it*MAXMOMSQ+imom2] =  twop_q2_b[ibin*T*MAXMOMSQ+it*MAXMOMSQ+imom2] + twop_q2[iconf*T*MAXMOMSQ + it*MAXMOMSQ + imom2 ];
      }
   
      
    }

    for(int it = 0 ; it < TSINK+1 ; it++)
      for(int mu = 0 ; mu < 4 ; mu++)
	for(int imom = 0 ; imom < NMOM ; imom++){
	  threepE_b[ibin*(TSINK+1)*4*NMOM + it*4*NMOM + mu*NMOM + imom].real() /= (Nconfs-BINSIZE);
	  threepE_b[ibin*(TSINK+1)*4*NMOM + it*4*NMOM + mu*NMOM + imom].imag() /= (Nconfs-BINSIZE);
	  threepM_b[ibin*(TSINK+1)*4*NMOM + it*4*NMOM + mu*NMOM + imom].real() /= (Nconfs-BINSIZE);
	  threepM_b[ibin*(TSINK+1)*4*NMOM + it*4*NMOM + mu*NMOM + imom].imag() /= (Nconfs-BINSIZE);
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
  Complex *RRE_b = (Complex*)calloc(Nbins*(TSINK+1)*4*NMOM,sizeof(Complex));
  Complex *RRM_b = (Complex*)calloc(Nbins*(TSINK+1)*4*NMOM,sizeof(Complex));

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

	  RRE_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM+mu*NMOM+imom] = ( threepE_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM + mu*NMOM + imom] / twop_q2_b[ibin*T*MAXMOMSQ+TSINK*MAXMOMSQ+0].real() )* squareRoot;
	  RRM_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM+mu*NMOM+imom] = ( threepM_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM + mu*NMOM + imom] / twop_q2_b[ibin*T*MAXMOMSQ+TSINK*MAXMOMSQ+0].real() )* squareRoot;
      
	}

  ///////////////////////////////////// 

  // mean bin values //
  double MEff_bmean[T]={};
  Complex RRE_bmean[TSINK+1][4][NMOM]={};
  Complex RRM_bmean[TSINK+1][4][NMOM]={};

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
	  RRE_bmean[it][mu][imom] = RRE_bmean[it][mu][imom] + RRE_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM+mu*NMOM+imom];
	  RRM_bmean[it][mu][imom] = RRM_bmean[it][mu][imom] + RRM_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM+mu*NMOM+imom];
	}
	RRE_bmean[it][mu][imom].real() /= Nbins;
	RRE_bmean[it][mu][imom].imag() /= Nbins;

	RRM_bmean[it][mu][imom].real() /= Nbins;
	RRM_bmean[it][mu][imom].imag() /= Nbins;
      }

  ////////////////////////////

  // jacknife errors //
  Complex temp_dthreepE = (Complex) {0,0};
  Complex temp_dthreepM = (Complex) {0,0};
  Complex temp_dtwop =(Complex) {0,0};
  double temp_dMEff = 0;
  Complex temp_dRRE = (Complex) {0,0};
  Complex temp_dRRM = (Complex) {0,0};
  
  Complex dthreepE[TSINK+1][4][NMOM];
  Complex dthreepM[TSINK+1][4][NMOM];
  Complex dtwop[T][MAXMOMSQ];
  double dMEff[T];
  Complex dRRE[TSINK+1][4][NMOM];
  Complex dRRM[TSINK+1][4][NMOM];

  for(int it = 0 ; it < TSINK+1 ; it++)
    for(int mu = 0 ; mu < 4 ; mu++)
      for(int imom = 0 ; imom < NMOM ; imom++){

	for(int ibin = 0 ; ibin < Nbins ; ibin++){
	  temp_dthreepE.real() = temp_dthreepE.real() + ( threepE_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM + mu*NMOM + imom].real() - threep_Electric_mean[it*4*NMOM + mu*NMOM + imom].real() )*(threepE_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM + mu*NMOM + imom].real() - threep_Electric_mean[it*4*NMOM + mu*NMOM + imom].real() ); 
	  temp_dthreepE.imag() = temp_dthreepE.imag() + ( threepE_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM + mu*NMOM + imom].imag() - threep_Electric_mean[it*4*NMOM + mu*NMOM + imom].imag() )*(threepE_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM + mu*NMOM + imom].imag() - threep_Electric_mean[it*4*NMOM + mu*NMOM + imom].imag() ); 

	  temp_dthreepM.real() = temp_dthreepM.real() + ( threepM_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM + mu*NMOM + imom].real() - threep_Magnetic_mean[it*4*NMOM + mu*NMOM + imom].real() )*(threepM_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM + mu*NMOM + imom].real() - threep_Magnetic_mean[it*4*NMOM + mu*NMOM + imom].real() ); 
	  temp_dthreepM.imag() = temp_dthreepM.imag() + ( threepM_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM + mu*NMOM + imom].imag() - threep_Magnetic_mean[it*4*NMOM + mu*NMOM + imom].imag() )*(threepM_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM + mu*NMOM + imom].imag() - threep_Magnetic_mean[it*4*NMOM + mu*NMOM + imom].imag() ); 


	  temp_dRRE.real() = temp_dRRE.real() + ( RRE_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM + mu*NMOM + imom].real() - RRE_bmean[it][mu][imom].real() )*(RRE_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM + mu*NMOM + imom].real() - RRE_bmean[it][mu][imom].real() ); 
	  temp_dRRE.imag() = temp_dRRE.imag() + ( RRE_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM + mu*NMOM + imom].imag() - RRE_bmean[it][mu][imom].imag() )*(RRE_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM + mu*NMOM + imom].imag() - RRE_bmean[it][mu][imom].imag() ); 

	  temp_dRRM.real() = temp_dRRM.real() + ( RRM_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM + mu*NMOM + imom].real() - RRM_bmean[it][mu][imom].real() )*(RRM_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM + mu*NMOM + imom].real() - RRM_bmean[it][mu][imom].real() ); 
	  temp_dRRM.imag() = temp_dRRM.imag() + ( RRM_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM + mu*NMOM + imom].imag() - RRM_bmean[it][mu][imom].imag() )*(RRM_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM + mu*NMOM + imom].imag() - RRM_bmean[it][mu][imom].imag() ); 

	}

	dthreepE[it][mu][imom].real() = sqrt( (Nbins-1)/( (double) Nbins ) * temp_dthreepE.real()  ); 
	dthreepE[it][mu][imom].imag() = sqrt( (Nbins-1)/( (double) Nbins ) * temp_dthreepE.imag()  ); 
	dthreepM[it][mu][imom].real() = sqrt( (Nbins-1)/( (double) Nbins ) * temp_dthreepM.real()  ); 
	dthreepM[it][mu][imom].imag() = sqrt( (Nbins-1)/( (double) Nbins ) * temp_dthreepM.imag()  ); 

	dRRE[it][mu][imom].real() = sqrt( (Nbins-1)/( (double) Nbins ) * temp_dRRE.real()  ); 
	dRRE[it][mu][imom].imag() = sqrt( (Nbins-1)/( (double) Nbins ) * temp_dRRE.imag()  ); 
	dRRM[it][mu][imom].real() = sqrt( (Nbins-1)/( (double) Nbins ) * temp_dRRM.real()  ); 
	dRRM[it][mu][imom].imag() = sqrt( (Nbins-1)/( (double) Nbins ) * temp_dRRM.imag()  ); 

	temp_dthreepE = (Complex) {0,0};
	temp_dthreepM = (Complex) {0,0};
	temp_dRRE = (Complex) {0,0};
	temp_dRRM = (Complex) {0,0};

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
  Complex *RREplateau_mean = (Complex*)calloc(4*NMOM,sizeof(Complex));
  Complex *RRMplateau_mean = (Complex*)calloc(4*NMOM,sizeof(Complex));
  Complex *RREplateau_b = (Complex*)calloc(Nbins*4*NMOM,sizeof(Complex));
  Complex *RRMplateau_b = (Complex*)calloc(Nbins*4*NMOM,sizeof(Complex));
  Complex *RREplateau_bmean = (Complex*)calloc(4*NMOM,sizeof(Complex));
  Complex *RRMplateau_bmean = (Complex*)calloc(4*NMOM,sizeof(Complex));


  for(int mu = 0 ; mu < 4 ; mu++)
    for(int imom = 0 ; imom < NMOM ; imom++){
      RREplateau_mean[mu*NMOM+imom] =  fit_ratio_plato(RRE_mean, dRRE, LOW3PT, HIGH3PT, mu, imom);
      RRMplateau_mean[mu*NMOM+imom] =  fit_ratio_plato(RRM_mean, dRRM, LOW3PT, HIGH3PT, mu, imom);
    }
  // 2) fit to jacknife data

  for(int ibin = 0 ; ibin < Nbins ; ibin++)
    for(int mu = 0 ; mu < 4 ; mu++)
      for(int imom = 0 ; imom < NMOM ; imom++){
	RREplateau_b[ibin*4*NMOM+mu*NMOM+imom] = fit_ratio_plato(RRE_b+ibin*(TSINK+1)*4*NMOM, dRRE,LOW3PT, HIGH3PT, mu, imom);
	RRMplateau_b[ibin*4*NMOM+mu*NMOM+imom] = fit_ratio_plato(RRM_b+ibin*(TSINK+1)*4*NMOM, dRRM,LOW3PT, HIGH3PT, mu, imom);
      }    

  for(int mu = 0 ; mu < 4 ; mu++)
    for(int imom = 0 ; imom < NMOM ; imom++){
      
      for(int ibin = 0 ; ibin < Nbins ; ibin++){
	RREplateau_bmean[mu*NMOM+imom] = RREplateau_bmean[mu*NMOM+imom] + RREplateau_b[ibin*4*NMOM+mu*NMOM+imom];
	RRMplateau_bmean[mu*NMOM+imom] = RRMplateau_bmean[mu*NMOM+imom] + RRMplateau_b[ibin*4*NMOM+mu*NMOM+imom];
      }
      RREplateau_bmean[mu*NMOM+imom].real() /= Nbins;
      RREplateau_bmean[mu*NMOM+imom].imag() /= Nbins;
      RRMplateau_bmean[mu*NMOM+imom].real() /= Nbins;
      RRMplateau_bmean[mu*NMOM+imom].imag() /= Nbins;
    }

  
    
  
  ////////////////////////////////////////

  // jacknife error to fitted values //
  Complex temp_dRREplateau = (Complex) {0.,0.};
  Complex temp_dRRMplateau = (Complex) {0.,0.};
  Complex *dRREplateau = (Complex*)calloc(4*NMOM,sizeof(Complex));
  Complex *dRRMplateau = (Complex*)calloc(4*NMOM,sizeof(Complex));

    for(int mu = 0 ; mu < 4 ; mu++)
      for(int imom = 0 ; imom < NMOM ; imom++){

	for(int ibin = 0 ; ibin < Nbins ; ibin++){

	  temp_dRREplateau.real() = temp_dRREplateau.real() + (RREplateau_b[ibin*4*NMOM+mu*NMOM+imom].real() - RREplateau_bmean[mu*NMOM+imom].real())*(RREplateau_b[ibin*4*NMOM+mu*NMOM+imom].real() - RREplateau_bmean[mu*NMOM+imom].real());

	  temp_dRREplateau.imag() = temp_dRREplateau.imag() + (RREplateau_b[ibin*4*NMOM+mu*NMOM+imom].imag() - RREplateau_bmean[mu*NMOM+imom].imag())*(RREplateau_b[ibin*4*NMOM+mu*NMOM+imom].imag() - RREplateau_bmean[mu*NMOM+imom].imag());

	  temp_dRRMplateau.real() = temp_dRRMplateau.real() + (RRMplateau_b[ibin*4*NMOM+mu*NMOM+imom].real() - RRMplateau_bmean[mu*NMOM+imom].real())*(RRMplateau_b[ibin*4*NMOM+mu*NMOM+imom].real() - RRMplateau_bmean[mu*NMOM+imom].real());

	  temp_dRRMplateau.imag() = temp_dRRMplateau.imag() + (RRMplateau_b[ibin*4*NMOM+mu*NMOM+imom].imag() - RRMplateau_bmean[mu*NMOM+imom].imag())*(RRMplateau_b[ibin*4*NMOM+mu*NMOM+imom].imag() - RRMplateau_bmean[mu*NMOM+imom].imag());

	  
	}
	
	dRREplateau[mu*NMOM+imom].real() = sqrt( (Nbins-1)/( (double) Nbins ) * temp_dRREplateau.real()  );
	dRREplateau[mu*NMOM+imom].imag() = sqrt( (Nbins-1)/( (double) Nbins ) * temp_dRREplateau.imag()  );

	dRRMplateau[mu*NMOM+imom].real() = sqrt( (Nbins-1)/( (double) Nbins ) * temp_dRRMplateau.real()  );
	dRRMplateau[mu*NMOM+imom].imag() = sqrt( (Nbins-1)/( (double) Nbins ) * temp_dRRMplateau.imag()  );

	temp_dRREplateau = (Complex) {0,0};
	temp_dRRMplateau = (Complex) {0,0};

      }

  ///////////////////////////////

    // extract quantities naive //
    double mass_mean;
    double Ep_mean[MAXMOMSQ];
    mass_mean=fit_EffMass_plato(MEff_mean, dMEff,LOW2PT, HIGH2PT);

    for(int imom = 0 ; imom < NMOM ; imom++)
      Ep_mean[p2[imom]] = sqrt(mass_mean*mass_mean + p2[imom]*(2.*PI/L)*(2.*PI/L) );

    
    double GE_mean[MAXMOMSQ] = {};
    double GM_mean[MAXMOMSQ] = {};

    extractGE(GE_mean,RREplateau_mean,dRREplateau,Ep_mean,momList,p2,numberMomPerQ2);
    extractGM(GM_mean,RRMplateau_mean,dRRMplateau,Ep_mean,momList,p2,numberMomPerQ2);
    /////////////////////////////

    // extract quantities for each bin //
    double *mass_b = (double*)calloc(Nbins,sizeof(double));
    double *Ep_b = (double*)calloc(Nbins*MAXMOMSQ,sizeof(double));
    double *GE_b = (double*)calloc(Nbins*MAXMOMSQ,sizeof(double));
    double *GM_b = (double*)calloc(Nbins*MAXMOMSQ,sizeof(double));

    double *Ep_bmean = (double*)calloc(MAXMOMSQ,sizeof(double));
    double *GE_bmean = (double*)calloc(MAXMOMSQ,sizeof(double));
    double *GM_bmean = (double*)calloc(MAXMOMSQ,sizeof(double));

    for(int ibin = 0 ; ibin < Nbins ; ibin++)
      mass_b[ibin] = fit_EffMass_plato(MEff_b+ibin*T, dMEff,LOW2PT, HIGH2PT);

    for(int ibin = 0 ; ibin < Nbins ; ibin++)
      for(int imom = 0 ; imom < NMOM ; imom++)
	Ep_b[ibin*MAXMOMSQ+p2[imom]] = sqrt(mass_b[ibin]*mass_b[ibin] + p2[imom]*(2.*PI/L)*(2.*PI/L) );

    for(int ibin = 0 ; ibin < Nbins ; ibin++){
      extractGE(GE_b + ibin*MAXMOMSQ,RREplateau_b + ibin*4*NMOM ,dRREplateau,Ep_b + ibin*MAXMOMSQ,momList,p2,numberMomPerQ2);
      extractGM(GM_b + ibin*MAXMOMSQ,RRMplateau_b + ibin*4*NMOM ,dRRMplateau,Ep_b + ibin*MAXMOMSQ,momList,p2,numberMomPerQ2);
    }    
    ////////////////////////////////////

    // find jackknife errors for form factors //
    double temp_Ep_bmean=0.;
    double temp_GE_bmean=0.;
    double temp_GM_bmean=0.;

    double *dGE = (double*)calloc(MAXMOMSQ,sizeof(double));
    double *dGM = (double*)calloc(MAXMOMSQ,sizeof(double));
    double *dEp = (double*)calloc(MAXMOMSQ,sizeof(double));

    for(int imom2 = 0 ; imom2 < MAXMOMSQ ; imom2++){
      for(int ibin = 0 ; ibin < Nbins ; ibin++){
	Ep_bmean[imom2] += Ep_b[ibin*MAXMOMSQ+imom2];
	GE_bmean[imom2] += GE_b[ibin*MAXMOMSQ+imom2];
	GM_bmean[imom2] += GM_b[ibin*MAXMOMSQ+imom2];
      }
      Ep_bmean[imom2] /= Nbins;
      GE_bmean[imom2] /= Nbins;
      GM_bmean[imom2] /= Nbins;
    }

    for(int imom2 = 0 ; imom2 < MAXMOMSQ ; imom2++){
      for(int ibin = 0 ; ibin < Nbins ; ibin++){
	temp_Ep_bmean += (Ep_bmean[imom2] - Ep_b[ibin*MAXMOMSQ+imom2])*(Ep_bmean[imom2] - Ep_b[ibin*MAXMOMSQ+imom2]); 
	temp_GE_bmean += (GE_bmean[imom2] - GE_b[ibin*MAXMOMSQ+imom2])*(GE_bmean[imom2] - GE_b[ibin*MAXMOMSQ+imom2]); 
	temp_GM_bmean += (GM_bmean[imom2] - GM_b[ibin*MAXMOMSQ+imom2])*(GM_bmean[imom2] - GM_b[ibin*MAXMOMSQ+imom2]); 
      }
      dGE[imom2] = sqrt( (Nbins-1)/( (double) Nbins ) * temp_GE_bmean  );
      dGM[imom2] = sqrt( (Nbins-1)/( (double) Nbins ) * temp_GM_bmean  );
      dEp[imom2] = sqrt( (Nbins-1)/( (double) Nbins ) * temp_Ep_bmean  );

      temp_Ep_bmean=0.;
      temp_GE_bmean=0.;
      temp_GM_bmean=0.;
    }
    ////////////////////////////////////////////

    // calculate data for ratio to print out //

    //printOut RE zero momentum

    for(int it = 0 ; it < TSINK+1 ; it++)
      fprintf(out_RE_mom0,"%d %+.8f %+.8f\n",it, RRE_bmean[it][0][0].real() , dRRE[it][0][0].real());

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
	    value[ibin*(TSINK+1)*12+it*12+i] = Z_V*RRE_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM+0*NMOM+(i+1)].real();

	for(int it = 0 ; it < TSINK+1 ; it++)
	  for(int ibin = 0; ibin < Nbins; ibin++){
	    value[ibin*(TSINK+1)*12+it*12+6] = Z_V*RRE_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM+1*NMOM+1].imag();
	    value[ibin*(TSINK+1)*12+it*12+7] = Z_V*RRE_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM+1*NMOM+2].imag();
	    value[ibin*(TSINK+1)*12+it*12+8] = Z_V*RRE_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM+2*NMOM+3].imag();
	    value[ibin*(TSINK+1)*12+it*12+9] = Z_V*RRE_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM+2*NMOM+4].imag();
	    value[ibin*(TSINK+1)*12+it*12+10] = Z_V*RRE_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM+3*NMOM+5].imag();
	    value[ibin*(TSINK+1)*12+it*12+11] = Z_V*RRE_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM+3*NMOM+6].imag();
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
	    value[ibin*(TSINK+1)*12+it*12+0] = Z_V*RRM_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM+1*NMOM+3].real();
	    value[ibin*(TSINK+1)*12+it*12+1] = Z_V*RRM_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM+1*NMOM+4].real();
	    value[ibin*(TSINK+1)*12+it*12+2] = Z_V*RRM_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM+1*NMOM+5].real();
	    value[ibin*(TSINK+1)*12+it*12+3] = Z_V*RRM_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM+1*NMOM+6].real();
	    value[ibin*(TSINK+1)*12+it*12+4] = Z_V*RRM_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM+2*NMOM+1].real();
	    value[ibin*(TSINK+1)*12+it*12+5] = Z_V*RRM_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM+2*NMOM+2].real();
	    value[ibin*(TSINK+1)*12+it*12+6] = Z_V*RRM_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM+2*NMOM+5].real();
	    value[ibin*(TSINK+1)*12+it*12+7] = Z_V*RRM_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM+2*NMOM+6].real();
	    value[ibin*(TSINK+1)*12+it*12+8] = Z_V*RRM_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM+3*NMOM+1].real();
	    value[ibin*(TSINK+1)*12+it*12+9] = Z_V*RRM_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM+3*NMOM+2].real();
	    value[ibin*(TSINK+1)*12+it*12+10] = Z_V*RRM_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM+3*NMOM+3].real();
	    value[ibin*(TSINK+1)*12+it*12+11] = Z_V*RRM_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM+3*NMOM+4].real();
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
    
    printf("\nELECTRIC FORM FACTOR\n");
    for(int imom2 = 0 ; imom2 < MAXMOMSQ; imom2++)
      printf("%f %f %f\n",q2GeV[imom2],GE_bmean[imom2],dGE[imom2]);

    printf("\nMAGNETIC FORM FACTOR\n");
    for(int imom2 = 0 ; imom2 < MAXMOMSQ; imom2++)
      printf("%f %f %f\n",q2GeV[imom2],GM_bmean[imom2],dGM[imom2]);

    //////////////////////////


    // free memory //
    for(int iconf=0;iconf<Nconfs;iconf++){
      free(confString[iconf]);
      free(threep_Electric_names[iconf]);
      free(threep_Magnetic_names[iconf]);
      free(twop_names[iconf]);
    }
    free(confString);
    free(threep_Electric_names);
    free(threep_Magnetic_names);
    free(twop_names);
    free(threep_Electric);
    free(threep_Magnetic);
    free(twop);
    free(RRE_mean);
    free(RRM_mean);
    free(twop_q2);
    free(threep_Electric_mean);
    free(threep_Magnetic_mean);
    free(twop_q2_mean);
    free(threepE_b);
    free(threepM_b);
    free(twop_q2_b);
    free(MEff_b);
    free(RRE_b);
    free(RRM_b);
    free(RREplateau_mean);
    free(RRMplateau_mean);
    free(RREplateau_b);
    free(RRMplateau_b);
    free(RREplateau_bmean);
    free(RRMplateau_bmean);
    free(dRREplateau);
    free(dRRMplateau);
    free(mass_b);
    free(Ep_b);
    free(GE_b);
    free(GM_b);
    free(Ep_bmean);
    free(GE_bmean);
    free(GM_bmean);
    free(dGE);
    free(dGM);
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

void extractGE(double *GE,Complex *RRE,Complex *dRRE, double *Ep,int momList[][NMOM], int p2[], int numberMomPerQ2[] ){

  double C; // kinematics term
  Complex *CoeffGe[MAXMOMSQ] = {NULL};
  Combinations *comb[MAXMOMSQ] = {NULL};
  double *b[MAXMOMSQ] = {NULL};
  double *db[MAXMOMSQ] = {NULL};
  double *A[MAXMOMSQ] = {NULL};
  double *u[MAXMOMSQ] = {NULL};
  double *w[MAXMOMSQ] = {NULL};
  double *vt[MAXMOMSQ] = {NULL};

  int counter = 0;
  
  for(int imom2 = 0 ; imom2 < MAXMOMSQ ; imom2++){
    GE[imom2]=0.;

    if(numberMomPerQ2[imom2] != 0){
      C = sqrt( ( 2.*Ep[0]*Ep[0] ) / ( Ep[imom2] * ( Ep[imom2]+Ep[0] ) ) );    

      
      for(int mu = 0 ; mu < 4 ; mu++)
	for(int imom = 0 ; imom < NMOM ; imom++){
	  int q2 = momList[0][imom]*momList[0][imom] + momList[1][imom]*momList[1][imom] + momList[2][imom]*momList[2][imom];
	  if(q2 == imom2){	    
	    if(mu == 0)counter++;        //printf("mom=%d,mu=%d\n",imom,mu);
	    if(mu == 1 && momList[0][imom] != 0)counter++; //printf("mom=%d,mu=%d\n",imom,mu);
	    if(mu == 2 && momList[1][imom] != 0)counter++; //printf("mom=%d,mu=%d\n",imom,mu);
	    if(mu == 3 && momList[2][imom] != 0)counter++; //printf("mom=%d,mu=%d\n",imom,mu);
	  }
	}

      CoeffGe[imom2] = (Complex*)calloc(counter,sizeof(Complex));
      comb[imom2] = (Combinations*)calloc(counter,sizeof(Combinations));
      b[imom2] = (double*)calloc(counter,sizeof(double));
      db[imom2] = (double*)calloc(counter,sizeof(double));
      A[imom2] = (double*)calloc(counter,sizeof(double));

      counter=0;

      for(int mu = 0 ; mu < 4 ; mu++)
	for(int imom = 0 ; imom < NMOM ; imom++){
	  int q2 = momList[0][imom]*momList[0][imom] + momList[1][imom]*momList[1][imom] + momList[2][imom]*momList[2][imom];
	  if(q2 == imom2){	    
	    if(mu == 0){ comb[imom2][counter].mom = imom; comb[imom2][counter].mu = mu;counter++; }
	    if(mu == 1 && momList[0][imom] != 0){ comb[imom2][counter].mom = imom; comb[imom2][counter].mu = mu;counter++; }
	    if(mu == 2 && momList[1][imom] != 0){ comb[imom2][counter].mom = imom; comb[imom2][counter].mu = mu;counter++; }
	    if(mu == 3 && momList[2][imom] != 0){ comb[imom2][counter].mom = imom; comb[imom2][counter].mu = mu;counter++; }
	  }
	}
      
      for(int i = 0 ; i < counter ; i++){
	if(comb[imom2][i].mu == 0){
	  CoeffGe[imom2][i].real() =  (C*(Ep[imom2]+Ep[0]))/(2.*Ep[0]); 
	}
	else{
	  CoeffGe[imom2][i].imag() =  -( C/(2.*Ep[0])) * (momList[comb[imom2][i].mu-1][comb[imom2][i].mom]*(2.*PI/L)) ;
	}
      }

      /*
      if(imom2 == 1){
	for(int i = 0 ; i < counter ; i++)
	  printf("%d %f %f\n",comb[imom2][i].mu,CoeffGe[imom2][i].real(),CoeffGe[imom2][i].imag());
      }
      printf("\n\n");
      */

      // transform to real problem
      for(int i = 0 ; i < counter ; i++){
	int mu = comb[imom2][i].mu;
	int imom = comb[imom2][i].mom;
	if(comb[imom2][i].mu == 0){
	  A[imom2][i] = real(CoeffGe[imom2][i]);
	  b[imom2][i] = real( RRE[mu*NMOM+imom]);
	  db[imom2][i] = real( dRRE[mu*NMOM+imom]);	
	}
	else{
	  A[imom2][i] = real(CoeffGe[imom2][i] * ( (Complex) {0,1.} ) );
	  b[imom2][i] = real( RRE[mu*NMOM+imom] * ( (Complex) {0,1.} ) );
	  db[imom2][i] = real( dRRE[mu*NMOM+imom] * ( (Complex) {0,1.} ) );
	}

      }


      //      make weighted least square solution: 
      
      for(int i = 0 ; i < counter ; i++){
	A[imom2][i] = A[imom2][i]/ fabs(db[imom2][i]);
	b[imom2][i] = b[imom2][i] /fabs(db[imom2][i]);
      }
      
      // ready for singular value decomposition


      int info;
      int lwork;
      double wkopt;
      double* work=NULL;
      int m=counter;
      int n=1;
      int lda=m;
      int ldu=m;
      int ldvt=n;
      u[imom2]=(double*)malloc(ldu*m*sizeof(double));
      w[imom2]=(double*)malloc(n*sizeof(double));
      vt[imom2]=(double*)malloc(ldvt*n*sizeof(double));


      lwork = -1;
      dgesvd_((char*) "All",(char*) "All", &m, &n, &(A[imom2][0]), &lda, w[imom2], u[imom2], &ldu, vt[imom2], &ldvt, &wkopt, &lwork, &info );
      lwork = (int)wkopt;
      work = (double*)malloc( lwork*sizeof(double) );
      // Compute SVD 
      dgesvd_((char*) "All",(char*) "All", &m, &n, &(A[imom2][0]), &lda, w[imom2], u[imom2], &ldu, vt[imom2], &ldvt, work, &lwork, &info );
      // Check for convergence 
      if( info > 0 ) {
	printf( "The algorithm computing SVD failed to converge.\n" );
	exit( 1 );
      }


      for(int i = 0 ; i < counter ; i++)
      	GE[imom2] += vt[imom2][0]*(1./w[imom2][0])*u[imom2][i*counter+0]*b[imom2][i]*Z_V;
      
    }
    else{
            GE[imom2] = log(-1); // assing nan
    }


    counter=0;
    free(CoeffGe[imom2]);
    free(comb[imom2]);  
    free(b[imom2]);
    free(db[imom2]);
    free(A[imom2]);
    free(u[imom2]);
    free(w[imom2]);
    free(vt[imom2]);
  }

}

void extractGM(double *GM,Complex *RRM,Complex *dRRM, double *Ep,int momList[][NMOM], int p2[], int numberMomPerQ2[] ){


  double C; // kinematics term
  Complex *CoeffGe[MAXMOMSQ] = {NULL};
  Combinations *comb[MAXMOMSQ] = {NULL};
  double *b[MAXMOMSQ] = {NULL};
  double *db[MAXMOMSQ] = {NULL};
  double *A[MAXMOMSQ] = {NULL};
  double *u[MAXMOMSQ] = {NULL};
  double *w[MAXMOMSQ] = {NULL};
  double *vt[MAXMOMSQ] ={NULL};

  int counter = 0;
  GM[0] = log(-1); // cant calculate GM(0) so set it nan

  for(int imom2 = 1 ; imom2 < MAXMOMSQ ; imom2++){
    GM[imom2]=0.;

    if(numberMomPerQ2[imom2] != 0){
      C = sqrt( ( 2.*Ep[0]*Ep[0] ) / ( Ep[imom2] * ( Ep[imom2]+Ep[0] ) ) );    

      
      for(int mu = 1 ; mu < 4 ; mu++)
	for(int imom = 0 ; imom < NMOM ; imom++){
	  int q2 = momList[0][imom]*momList[0][imom] + momList[1][imom]*momList[1][imom] + momList[2][imom]*momList[2][imom];
	  if(q2 == imom2){	    
	    if(mu == 1 && ( momList[1][imom] != 0 || momList[2][imom] != 0) ) counter++;// printf("mom=%d,mu=%d\n",imom,mu);
	    if(mu == 2 && ( momList[0][imom] != 0 || momList[2][imom] != 0 ) )counter++;// printf("mom=%d,mu=%d\n",imom,mu);
	    if(mu == 3 && ( momList[0][imom] != 0 || momList[1][imom] != 0 ) )counter++; //printf("mom=%d,mu=%d\n",imom,mu);
	  }
	}


      CoeffGe[imom2] = (Complex*)calloc(counter,sizeof(Complex));
      comb[imom2] = (Combinations*)calloc(counter,sizeof(Combinations));
      b[imom2] = (double*)calloc(counter,sizeof(double));
      db[imom2] = (double*)calloc(counter,sizeof(double));
      A[imom2] = (double*)calloc(counter,sizeof(double));

      counter=0;

      for(int mu = 1 ; mu < 4 ; mu++)
	for(int imom = 0 ; imom < NMOM ; imom++){
	  int q2 = momList[0][imom]*momList[0][imom] + momList[1][imom]*momList[1][imom] + momList[2][imom]*momList[2][imom];
	  if(q2 == imom2){	    

	    if(mu == 1 && ( momList[1][imom] != 0 || momList[2][imom] != 0) ){ comb[imom2][counter].mom = imom; comb[imom2][counter].mu = mu;counter++; }
	    if(mu == 2 && ( momList[0][imom] != 0 || momList[2][imom] != 0) ){ comb[imom2][counter].mom = imom; comb[imom2][counter].mu = mu;counter++; }
	    if(mu == 3 && ( momList[0][imom] != 0 || momList[1][imom] != 0) ){ comb[imom2][counter].mom = imom; comb[imom2][counter].mu = mu;counter++; }
	  }
	}
      
      for(int i = 0 ; i < counter ; i++){
	if(comb[imom2][i].mu == 1) CoeffGe[imom2][i].real() =  -( C/(2.*Ep[0])) * ( (momList[1][comb[imom2][i].mom] - momList[2][comb[imom2][i].mom])*(2.*PI/L)) ;
	if(comb[imom2][i].mu == 2) CoeffGe[imom2][i].real() =  -( C/(2.*Ep[0])) * ( (momList[2][comb[imom2][i].mom] - momList[0][comb[imom2][i].mom])*(2.*PI/L)) ;
	if(comb[imom2][i].mu == 3) CoeffGe[imom2][i].real() =  -( C/(2.*Ep[0])) * ( (momList[0][comb[imom2][i].mom] - momList[1][comb[imom2][i].mom])*(2.*PI/L)) ;
      }


      /*
      if(imom2 == 1){
	for(int i = 0 ; i < counter ; i++)
	  printf("%d %d %f %f\n",comb[imom2][i].mu,comb[imom2][i].mom,CoeffGe[imom2][i].real(),CoeffGe[imom2][i].imag());
      }
      printf("\n\n");
      */

      // transform to real problem // this problem is already real
      for(int i = 0 ; i < counter ; i++){
	int mu = comb[imom2][i].mu;
	int imom = comb[imom2][i].mom;
	A[imom2][i] = real(CoeffGe[imom2][i]);
	b[imom2][i] = real(RRM[mu*NMOM+imom]);
	db[imom2][i] = real(dRRM[mu*NMOM+imom]);
      }

      //      make weighted least square solution: 
      
      for(int i = 0 ; i < counter ; i++){
	A[imom2][i] = A[imom2][i]/ fabs(db[imom2][i]);
	b[imom2][i] = b[imom2][i] /fabs(db[imom2][i]);
      }
      
      // ready for singular value decomposition


      int info;
      int lwork;
      double wkopt;
      double* work=NULL;
      int m=counter;
      int n=1;
      int lda=m;
      int ldu=m;
      int ldvt=n;
      u[imom2]=(double*)malloc(ldu*m*sizeof(double));
      w[imom2]=(double*)malloc(n*sizeof(double));
      vt[imom2]=(double*)malloc(ldvt*n*sizeof(double));


      lwork = -1;
      dgesvd_((char*) "All",(char*) "All", &m, &n, &(A[imom2][0]), &lda, w[imom2], u[imom2], &ldu, vt[imom2], &ldvt, &wkopt, &lwork, &info );
      lwork = (int)wkopt;
      work = (double*)malloc( lwork*sizeof(double) );
      // Compute SVD 
      dgesvd_((char*) "All",(char*) "All", &m, &n, &(A[imom2][0]), &lda, w[imom2], u[imom2], &ldu, vt[imom2], &ldvt, work, &lwork, &info );
      // Check for convergence 
      if( info > 0 ) {
	printf( "The algorithm computing SVD failed to converge.\n" );
	exit( 1 );
      }


      for(int i = 0 ; i < counter ; i++)
      	GM[imom2] += vt[imom2][0]*(1./w[imom2][0])*u[imom2][i*counter+0]*b[imom2][i] * Z_V;


    }
    else{
            GM[imom2] = log(-1); // assing nan
    }


    counter=0;
    free(CoeffGe[imom2]);
    free(comb[imom2]);  
    free(b[imom2]);
    free(db[imom2]);
    free(A[imom2]);
    free(u[imom2]);
    free(w[imom2]);
    free(vt[imom2]);

  }


}
