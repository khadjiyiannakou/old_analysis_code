/*
Kyriakos Hadjiyiannakou
Extract Tensor form factor
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
3) extractGFF
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
#define ZT 0.77

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
Complex fit_ratio_plato2(Complex* ratio, Complex error_ratio[TSINK+1][4*4][NMOM], int fit_plato_low, int fit_plato_high, int mu, int nu, int imom);
double fit_EffMass_plato(double MEff[], double dMEff[],int fit_plato_low, int fit_plato_high);
void extract_tensor(double *A10,double *B10,double *C10,Complex *RR[],Complex *dRR[], double *Ep,int momList[][NMOM], int p2[], int numberMomPerQ2[] );




int main(int argc, char *argv[]){

  clock_t time;

  // Greeting //
  printf("This program calculate Tensor form factors\n");
  printf("Right inputs order must be\n");
  printf("(1) Executable , (2) Trajectory list, (3) Threep tensor type1, (4) Threep tensor type2, (5) Twop base name, (6) Output name\n\n ");
  //////////////////


  // check passing inputs right //
  if(argc != 6){
    fprintf(stderr,"Error: Wrong number of input files \n");
    exit(EXIT_FAILURE);
  }
  //////////////////////


  // fill filenames and print paths //
  char filename_Traj[MAX_STRING];
  char filename_ThreepType1[MAX_STRING];
  char filename_ThreepType2[MAX_STRING];

  char filename_Twop[MAX_STRING];
  char filename_Output[MAX_STRING];

  //output files
  char filename_R_mom0[MAX_STRING];
  char filename_R_mom1[MAX_STRING];

  strcpy(filename_Traj,argv[1]);
  strcpy(filename_ThreepType1,argv[2]);
  strcpy(filename_ThreepType2,argv[3]);
  strcpy(filename_Twop,argv[4]);
  strcpy(filename_Output,argv[5]);

  printf("Got filename for trajectories : %s\n",filename_Traj);
  printf("Got filename for threep type1 : %s\n",filename_ThreepType1);
  printf("Got filename for threep type2 : %s\n",filename_ThreepType2);

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

  char **threep_namesType1, **threep_namesType2;
  char **twop_names;

  threep_namesType1 = (char**)malloc(Nconfs*sizeof(char*));
  threep_namesType2 = (char**)malloc(Nconfs*sizeof(char*));

  twop_names = (char**)malloc(Nconfs*sizeof(char*));

  for(int iconf=0;iconf<Nconfs;iconf++){
    threep_namesType1[iconf] = (char*)malloc(MAX_STRING*sizeof(char));
    threep_namesType2[iconf] = (char*)malloc(MAX_STRING*sizeof(char));
    twop_names[iconf] = (char*)malloc(MAX_STRING*sizeof(char));
  }

  for(int iconf=0;iconf<Nconfs;iconf++){
    sprintf(threep_namesType1[iconf],"%s.%s_isov",filename_ThreepType1,confString[iconf]);
    sprintf(threep_namesType2[iconf],"%s.%s_isov",filename_ThreepType2,confString[iconf]);
    sprintf(twop_names[iconf],"%s.%s",filename_Twop,confString[iconf]);
  }

  ///////////////////////////////////////////////


  // allocate memory & read data //
  Complex *threep[2];
  for(int i = 0 ; i < 2 ;i++)threep[i] = (Complex*)malloc(Nconfs*(TSINK+1)*4*4*NMOM*sizeof(Complex)); 

  Complex *twop = (Complex*)malloc(Nconfs*T*NMOM*sizeof(Complex));
  int momList[3][NMOM];
  int dummy;

  for(int i = 0 ; i < 2 ; i++)
    if(threep[i] == NULL){ fprintf(stderr,"Error: Out of memory\n"); exit(EXIT_FAILURE);}

  if(twop == NULL){ fprintf(stderr,"Error: Out of memory\n"); exit(EXIT_FAILURE);}
  FILE *ptr_threep[2] = {NULL} ,*ptr_twop = NULL;

  time = clock();

  
  for(int iconf = 0 ; iconf < Nconfs ; iconf++){
    ptr_threep[0] = fopen(threep_namesType1[iconf],"r");
    ptr_threep[1] = fopen(threep_namesType2[iconf],"r");

    ptr_twop = fopen(twop_names[iconf],"r");
    if(ptr_threep[0] == NULL || ptr_threep[1] == NULL || ptr_twop == NULL){ fprintf(stderr,"Error: open files for reading\n"); exit(EXIT_FAILURE);}
    // read threep

    for(int i = 0 ; i < 2 ; i++)
      for(int it=0; it < TSINK+1 ; it++)
	for(int mu=0; mu<4; mu++)
	  for(int nu=mu+1; nu<4; nu++)
	    for(int imom = 0 ; imom < NMOM ; imom++)
	      returnValue = fscanf(ptr_threep[i],"%d %d %d %d %lf %lf %d %d",&dummy,&(momList[0][imom]), &(momList[1][imom]), &(momList[2][imom]), &(threep[i][iconf*(TSINK+1)*4*4*NMOM + it*4*4*NMOM + mu*4*NMOM + nu*NMOM + imom].real()), &(threep[i][iconf*(TSINK+1)*4*4*NMOM + it*4*4*NMOM + mu*4*NMOM + nu*NMOM + imom].imag()), &dummy, &dummy );


    for(int it=0; it < T ; it++)
      for(int imom = 0 ; imom < NMOM ; imom++)
	returnValue = fscanf(ptr_twop,"%d %d %d %d %lf %lf",&dummy,&(momList[0][imom]), &(momList[1][imom]), &(momList[2][imom]),&(twop[iconf*T*NMOM + it*NMOM + imom].real()), &(twop[iconf*T*NMOM + it*NMOM + imom].imag()) );


    fclose(ptr_threep[0]);
    fclose(ptr_threep[1]);
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

  // if you use quda for stochastic I use fourier with different sign  so be carefull
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
  Complex *threep_mean[2];
  for(int i = 0 ; i < 2 ; i++) threep_mean[i] = (Complex*)calloc((TSINK+1)*4*4*NMOM,sizeof(Complex));

  Complex *twop_q2_mean = (Complex*)calloc(T*MAXMOMSQ,sizeof(Complex));

  for(int i = 0 ; i < 2 ; i++)
    for(int it = 0 ; it < TSINK+1 ; it++)
      for(int mu=0; mu<4; mu++)
	for(int nu=mu+1; nu<4; nu++)
	  for(int imom = 0 ; imom < NMOM ; imom++){
	    for(int iconf = 0 ; iconf < Nconfs ; iconf++){
	      threep_mean[i][it*4*4*NMOM + mu*4*NMOM  + nu*NMOM + imom] = threep_mean[i][it*4*4*NMOM + mu*4*NMOM  + nu*NMOM + imom] + threep[i][iconf*(TSINK+1)*4*4*NMOM + it*4*4*NMOM + mu*4*NMOM + nu*NMOM + imom];
	    }
	    threep_mean[i][it*4*4*NMOM + mu*4*NMOM  + nu*NMOM + imom].real() /= Nconfs;
	    threep_mean[i][it*4*4*NMOM + mu*4*NMOM  + nu*NMOM + imom].imag() /= Nconfs;
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
  Complex *RR_mean[2];
  for(int i = 0 ; i < 2 ;i++) RR_mean[i] = (Complex*)malloc((TSINK+1)*4*4*NMOM*sizeof(Complex));

  Complex squareRoot;
  double subRoot;

  for(int i = 0 ; i < 2 ;i++)
    for(int it=0;it < TSINK+1;it++)
      for(int mu=0; mu<4; mu++)
	for(int nu=mu+1; nu<4; nu++)
	  for(int imom = 0 ; imom < NMOM ; imom++){

	    subRoot = (twop_q2_mean[(TSINK-it)*MAXMOMSQ+p2[imom]].real()  * twop_q2_mean[(it)*MAXMOMSQ+0].real() *  twop_q2_mean[TSINK*MAXMOMSQ+0].real() ) /(twop_q2_mean[(TSINK-it)*MAXMOMSQ+0].real()  * twop_q2_mean[(it)*MAXMOMSQ+p2[imom]].real() *  twop_q2_mean[TSINK*MAXMOMSQ+p2[imom]].real() );
	    squareRoot = std::sqrt( (Complex)  subRoot  );

	    RR_mean[i][it*4*4*NMOM+mu*4*NMOM+nu*NMOM+imom] = ( threep_mean[i][it*4*4*NMOM+mu*4*NMOM+nu*NMOM+imom] / twop_q2_mean[TSINK*MAXMOMSQ+0].real() )* squareRoot;
	  }

  //////////////////////////////



  // binning for jackknife analysis //
  Complex *threep_b[2];
  for(int i = 0 ; i < 2 ; i++) threep_b[i] = (Complex*)calloc(Nbins*(TSINK+1)*4*4*NMOM,sizeof(Complex));

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

	for(int i = 0 ; i < 2 ; i++)
	  for(int it = 0 ; it < TSINK+1 ; it++)
	    for(int mu=0; mu<4; mu++)
	      for(int nu=mu+1; nu<4; nu++)
		for(int imom = 0 ; imom < NMOM ; imom++){

		  threep_b[i][ibin*(TSINK+1)*4*4*NMOM + it*4*4*NMOM + mu*4*NMOM + nu*NMOM + imom] = threep_b[i][ibin*(TSINK+1)*4*4*NMOM + it*4*4*NMOM + mu*4*NMOM + nu*NMOM + imom] + threep[i][iconf*(TSINK+1)*4*4*NMOM + it*4*4*NMOM + mu*4*NMOM + nu*NMOM + imom];
		  
		}

	for(int it = 0 ; it < T ; it++)
	  for(int imom2 = 0 ; imom2 < MAXMOMSQ; imom2++)
	    twop_q2_b[ibin*T*MAXMOMSQ+it*MAXMOMSQ+imom2] =  twop_q2_b[ibin*T*MAXMOMSQ+it*MAXMOMSQ+imom2] + twop_q2[iconf*T*MAXMOMSQ + it*MAXMOMSQ + imom2 ];
      }
   
      
    }

    for(int i = 0 ; i < 2 ; i++)
      for(int it = 0 ; it < TSINK+1 ; it++)
	for(int mu=0; mu<4; mu++)
	  for(int nu=mu+1; nu<4; nu++)
	    for(int imom = 0 ; imom < NMOM ; imom++){
	      threep_b[i][ibin*(TSINK+1)*4*4*NMOM + it*4*4*NMOM + mu*4*NMOM + nu*NMOM +imom].real() /= (Nconfs-BINSIZE);
	      threep_b[i][ibin*(TSINK+1)*4*4*NMOM + it*4*4*NMOM + mu*4*NMOM + nu*NMOM +imom].imag() /= (Nconfs-BINSIZE);
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
  Complex *RR_b[2];
  for(int i = 0 ; i < 2 ; i++) RR_b[i] = (Complex*)calloc(Nbins*(TSINK+1)*4*4*NMOM,sizeof(Complex));

  for(int ibin = 0 ; ibin < Nbins ; ibin++)
    for(int it = 1 ; it < T ; it++){
      MEff_b[ibin*T+it] = log(twop_q2_b[ibin*T*MAXMOMSQ+it*MAXMOMSQ+0].real() / twop_q2_b[ibin*T*MAXMOMSQ+(it+1)*MAXMOMSQ+0].real());
    }

  for(int i = 0 ; i < 2 ; i++)
    for(int ibin = 0 ; ibin < Nbins ; ibin++)
      for(int it=0;it < TSINK+1;it++)
	for(int mu=0; mu<4; mu++)
	  for(int nu=mu+1; nu<4; nu++)
	    for(int imom = 0 ; imom < NMOM ; imom++){

	      subRoot = (twop_q2_b[ibin*T*MAXMOMSQ+(TSINK-it)*MAXMOMSQ+p2[imom]].real()  * twop_q2_b[ibin*T*MAXMOMSQ+(it)*MAXMOMSQ+0].real() *  twop_q2_b[ibin*T*MAXMOMSQ+TSINK*MAXMOMSQ+0].real() ) /(twop_q2_b[ibin*T*MAXMOMSQ+(TSINK-it)*MAXMOMSQ+0].real()  * twop_q2_b[ibin*T*MAXMOMSQ+(it)*MAXMOMSQ+p2[imom]].real() *  twop_q2_b[ibin*T*MAXMOMSQ+TSINK*MAXMOMSQ+p2[imom]].real() );
	      squareRoot = std::sqrt( (Complex)  subRoot  );

	      RR_b[i][ibin*(TSINK+1)*4*4*NMOM+it*4*4*NMOM+mu*4*NMOM+nu*NMOM+imom] = ( threep_b[i][ibin*(TSINK+1)*4*4*NMOM + it*4*4*NMOM + mu*4*NMOM + nu*NMOM + imom] / twop_q2_b[ibin*T*MAXMOMSQ+TSINK*MAXMOMSQ+0].real() )* squareRoot;

	    }

  ///////////////////////////////////// 



  // mean bin values //
  double MEff_bmean[T]={};
  Complex RR_bmean[2][TSINK+1][4*4][NMOM]={};

  for(int it = 0 ; it < T ; it++){
    for(int ibin = 0 ; ibin < Nbins ; ibin++){
      MEff_bmean[it] += MEff_b[ibin*T+it];
    }
    MEff_bmean[it] /= Nbins;
  }

  for(int i = 0 ; i < 2 ; i++)
    for(int it = 0 ; it < TSINK+1 ; it++)
      for(int mu=0; mu<4; mu++)
	for(int nu=mu+1; nu<4; nu++)
	  for(int imom = 0 ; imom < NMOM ; imom++){
	    for(int ibin = 0 ; ibin < Nbins ; ibin++){
	      RR_bmean[i][it][mu*4+nu][imom] = RR_bmean[i][it][mu*4+nu][imom] + RR_b[i][ibin*(TSINK+1)*4*4*NMOM+it*4*4*NMOM+mu*4*NMOM+nu*NMOM+imom];
	    }
	    RR_bmean[i][it][mu*4+nu][imom].real() /= Nbins;
	    RR_bmean[i][it][mu*4+nu][imom].imag() /= Nbins;
	  }

  ////////////////////////////




  // jacknife errors //
  Complex temp_dthreep = (Complex) {0,0};
  Complex temp_dtwop =(Complex) {0,0};
  double temp_dMEff = 0;
  Complex temp_dRR = (Complex) {0,0};
  
  Complex dthreep[2][TSINK+1][4*4][NMOM];
  Complex dtwop[T][MAXMOMSQ];
  double dMEff[T];
  Complex dRR[2][TSINK+1][4*4][NMOM];


  for(int i = 0 ; i < 2 ; i++)
    for(int it = 0 ; it < TSINK+1 ; it++)
      for(int mu=0; mu<4; mu++)
	for(int nu=mu+1; nu<4; nu++)
	  for(int imom = 0 ; imom < NMOM ; imom++){
	    for(int ibin = 0 ; ibin < Nbins ; ibin++){
	      temp_dthreep.real() = temp_dthreep.real() + ( threep_b[i][ibin*(TSINK+1)*4*4*NMOM+it*4*4*NMOM + mu*4*NMOM + nu*NMOM + imom].real() - threep_mean[i][it*4*4*NMOM + mu*4*NMOM + nu*NMOM + imom].real() )*(threep_b[i][ibin*(TSINK+1)*4*4*NMOM+it*4*4*NMOM + mu*4*NMOM + nu*NMOM + imom].real() - threep_mean[i][it*4*4*NMOM + mu*4*NMOM + nu*NMOM + imom].real() ); 
	      temp_dthreep.imag() = temp_dthreep.imag() + ( threep_b[i][ibin*(TSINK+1)*4*4*NMOM+it*4*4*NMOM + mu*4*NMOM + nu*NMOM + imom].imag() - threep_mean[i][it*4*4*NMOM + mu*4*NMOM + nu*NMOM + imom].imag() )*(threep_b[i][ibin*(TSINK+1)*4*4*NMOM+it*4*4*NMOM + mu*4*NMOM + nu*NMOM + imom].imag() - threep_mean[i][it*4*4*NMOM + mu*4*NMOM + nu*NMOM + imom].imag() ); 
	      
	      	      
	      temp_dRR.real() = temp_dRR.real() + ( RR_b[i][ibin*(TSINK+1)*4*4*NMOM+it*4*4*NMOM + mu*4*NMOM + nu*NMOM + imom].real() - RR_bmean[i][it][mu*4+nu][imom].real() )*(RR_b[i][ibin*(TSINK+1)*4*4*NMOM+it*4*4*NMOM + mu*4*NMOM + nu*NMOM + imom].real() - RR_bmean[i][it][mu*4+nu][imom].real() ); 
	      temp_dRR.imag() = temp_dRR.imag() + ( RR_b[i][ibin*(TSINK+1)*4*4*NMOM+it*4*4*NMOM + mu*4*NMOM + nu*NMOM + imom].imag() - RR_bmean[i][it][mu*4+nu][imom].imag() )*(RR_b[i][ibin*(TSINK+1)*4*4*NMOM+it*4*4*NMOM + mu*4*NMOM + nu*NMOM + imom].imag() - RR_bmean[i][it][mu*4+nu][imom].imag() ); 

	      
	    }
	    
	    dthreep[i][it][mu*4+nu][imom].real() = sqrt( (Nbins-1)/( (double) Nbins ) * temp_dthreep.real()  ); 
	    dthreep[i][it][mu*4+nu][imom].imag() = sqrt( (Nbins-1)/( (double) Nbins ) * temp_dthreep.imag()  ); 
	    
	    dRR[i][it][mu*4+nu][imom].real() = sqrt( (Nbins-1)/( (double) Nbins ) * temp_dRR.real()  ); 
	    dRR[i][it][mu*4+nu][imom].imag() = sqrt( (Nbins-1)/( (double) Nbins ) * temp_dRR.imag()  ); 

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

  // test write
  //  for(int i = 1 ; i < 4 ; i++){

  //    for(int it = 0 ; it <= TSINK ; it++)
  //    printf("%d %+e %+e\n",it,RR_bmean[i][it][i][0].imag(),dRR[i][it][i][0].imag());
  //  printf("\n\n");
  // }


  // fit to extract plateaus //


  // 1) fit to naive data
  Complex *RRplateau_mean[2];
  Complex *RRplateau_b[2]; 
  Complex *RRplateau_bmean[2]; 
  
  for(int i = 0 ; i < 2 ; i++){
    RRplateau_mean[i] = (Complex*)calloc(4*4*NMOM,sizeof(Complex));
    RRplateau_b[i] = (Complex*)calloc(Nbins*4*4*NMOM,sizeof(Complex));
    RRplateau_bmean[i] = (Complex*)calloc(4*4*NMOM,sizeof(Complex));
  }

  for(int i = 0 ; i < 2 ; i++)
    for(int mu=0; mu<4; mu++)
      for(int nu=mu+1; nu<4; nu++)
	for(int imom = 0 ; imom < NMOM ; imom++){
	  RRplateau_mean[i][mu*4*NMOM+nu*NMOM+imom] =  fit_ratio_plato2(RR_mean[i], dRR[i], LOW3PT, HIGH3PT, mu,nu,imom);
	}
  // 2) fit to jacknife data

  for(int i = 0 ; i < 2 ; i++)
    for(int ibin = 0 ; ibin < Nbins ; ibin++)
      for(int mu=0; mu<4; mu++)
	for(int nu=mu+1; nu<4; nu++)
	  for(int imom = 0 ; imom < NMOM ; imom++){
	    RRplateau_b[i][ibin*4*4*NMOM+mu*4*NMOM+nu*NMOM+imom] = fit_ratio_plato2(RR_b[i]+ibin*(TSINK+1)*4*4*NMOM, dRR[i],LOW3PT, HIGH3PT, mu,nu, imom);
	  }    

  for(int i = 0 ; i < 2 ; i++)
    for(int mu=0; mu<4; mu++)
      for(int nu=mu+1; nu<4; nu++)
	for(int imom = 0 ; imom < NMOM ; imom++){
      
	  for(int ibin = 0 ; ibin < Nbins ; ibin++){
	    RRplateau_bmean[i][mu*4*NMOM+nu*NMOM+imom] = RRplateau_bmean[i][mu*4*NMOM+nu*NMOM+imom] + RRplateau_b[i][ibin*4*4*NMOM+mu*4*NMOM+nu*NMOM+imom];
	  }
	  RRplateau_bmean[i][mu*4*NMOM+nu*NMOM+imom].real() /= Nbins;
	  RRplateau_bmean[i][mu*4*NMOM+nu*NMOM+imom].imag() /= Nbins;
	}


  
  
  ////////////////////////////////////////

  // jacknife error to fitted values //
  Complex temp_dRRplateau = (Complex) {0.,0.};
  Complex *dRRplateau[2];

  for(int i = 0 ; i < 2 ; i++)
    dRRplateau[i] = (Complex*)calloc(4*4*NMOM,sizeof(Complex));


  for(int i = 0 ; i < 2 ; i++)
    for(int mu=0; mu<4; mu++)
      for(int nu=mu+1; nu<4; nu++)
	for(int imom = 0 ; imom < NMOM ; imom++){
	  for(int ibin = 0 ; ibin < Nbins ; ibin++){

	    temp_dRRplateau.real() = temp_dRRplateau.real() + (RRplateau_b[i][ibin*4*4*NMOM+mu*4*NMOM+nu*NMOM+imom].real() - RRplateau_bmean[i][mu*4*NMOM+nu*NMOM+imom].real())*(RRplateau_b[i][ibin*4*4*NMOM+mu*4*NMOM+nu*NMOM+imom].real() - RRplateau_bmean[i][mu*4*NMOM+nu*NMOM+imom].real());

	    temp_dRRplateau.imag() = temp_dRRplateau.imag() + (RRplateau_b[i][ibin*4*4*NMOM+mu*4*NMOM+nu*NMOM+imom].imag() - RRplateau_bmean[i][mu*4*NMOM+nu*NMOM+imom].imag())*(RRplateau_b[i][ibin*4*4*NMOM+mu*4*NMOM+nu*NMOM+imom].imag() - RRplateau_bmean[i][mu*4*NMOM+nu*NMOM+imom].imag());
	  
	  }
	
	  dRRplateau[i][mu*4*NMOM+nu*NMOM+imom].real() = sqrt( (Nbins-1)/( (double) Nbins ) * temp_dRRplateau.real()  );
	  dRRplateau[i][mu*4*NMOM+nu*NMOM+imom].imag() = sqrt( (Nbins-1)/( (double) Nbins ) * temp_dRRplateau.imag()  );

	  temp_dRRplateau = (Complex) {0,0};
	  
	}

  ///////////////////////////////


  // extract quantities naive //
  double mass_mean;
  double Ep_mean[MAXMOMSQ];
  mass_mean=fit_EffMass_plato(MEff_mean, dMEff,LOW2PT, HIGH2PT);
  
  for(int imom = 0 ; imom < NMOM ; imom++)
    Ep_mean[p2[imom]] = sqrt(mass_mean*mass_mean + p2[imom]*(2.*PI/L)*(2.*PI/L) );
  
    
  double A10_mean[MAXMOMSQ] = {};
  double B10_mean[MAXMOMSQ] = {};
  double C10_mean[MAXMOMSQ] = {};


  extract_tensor(A10_mean,B10_mean,C10_mean,RRplateau_mean,dRRplateau,Ep_mean,momList,p2,numberMomPerQ2);

    //    for(int imom2 = 0 ; imom2 < MAXMOMSQ ; imom2++)
    //  printf("%d %+e\n",imom2,GA_mean[imom2]);

    // extract quantities for each bin //
  double *mass_b = (double*)calloc(Nbins,sizeof(double));
  double *Ep_b = (double*)calloc(Nbins*MAXMOMSQ,sizeof(double));
  double *A10_b = (double*)calloc(Nbins*MAXMOMSQ,sizeof(double));
  double *B10_b = (double*)calloc(Nbins*MAXMOMSQ,sizeof(double));
  double *C10_b = (double*)calloc(Nbins*MAXMOMSQ,sizeof(double));
  
  double *Ep_bmean = (double*)calloc(MAXMOMSQ,sizeof(double));
  double *A10_bmean = (double*)calloc(MAXMOMSQ,sizeof(double));
  double *B10_bmean = (double*)calloc(MAXMOMSQ,sizeof(double));
  double *C10_bmean = (double*)calloc(MAXMOMSQ,sizeof(double));
  

  for(int ibin = 0 ; ibin < Nbins ; ibin++)
    mass_b[ibin] = fit_EffMass_plato(MEff_b+ibin*T, dMEff,LOW2PT, HIGH2PT);
  
  for(int ibin = 0 ; ibin < Nbins ; ibin++)
      for(int imom = 0 ; imom < NMOM ; imom++)
	Ep_b[ibin*MAXMOMSQ+p2[imom]] = sqrt(mass_b[ibin]*mass_b[ibin] + p2[imom]*(2.*PI/L)*(2.*PI/L) );


    Complex *RR_plateau_test[2];
    for(int i = 0 ; i < 2 ; i++)RR_plateau_test[i] = (Complex*)calloc(4*4*NMOM,sizeof(Complex));

    for(int ibin = 0 ; ibin < Nbins ; ibin++){

      for(int i = 0 ; i < 2 ; i++)
	for(int mu=0; mu<4; mu++)
	  for(int nu=mu+1; nu<4; nu++)
	    for(int imom = 0 ; imom < NMOM ; imom++)
	      RR_plateau_test[i][mu*4*NMOM+nu*NMOM+imom] = RRplateau_b[i][ibin*4*4*NMOM+mu*4*NMOM+nu*NMOM+imom];

      extract_tensor(A10_b + ibin*MAXMOMSQ,B10_b + ibin*MAXMOMSQ,C10_b + ibin*MAXMOMSQ,RR_plateau_test,dRRplateau,Ep_b + ibin*MAXMOMSQ,momList,p2,numberMomPerQ2);
    }    

    
    ////////////////////////////////////



    // find jackknife errors for form factors //
    double temp_Ep_bmean=0.;
    double temp_A10_bmean=0.;
    double temp_B10_bmean=0.;
    double temp_C10_bmean=0.;

    double *dA10 = (double*)calloc(MAXMOMSQ,sizeof(double));
    double *dB10 = (double*)calloc(MAXMOMSQ,sizeof(double));
    double *dC10 = (double*)calloc(MAXMOMSQ,sizeof(double));

    double *dEp = (double*)calloc(MAXMOMSQ,sizeof(double));

    for(int imom2 = 0 ; imom2 < MAXMOMSQ ; imom2++){
      for(int ibin = 0 ; ibin < Nbins ; ibin++){
	Ep_bmean[imom2] += Ep_b[ibin*MAXMOMSQ+imom2];
	A10_bmean[imom2] += A10_b[ibin*MAXMOMSQ+imom2];
	B10_bmean[imom2] += B10_b[ibin*MAXMOMSQ+imom2];
	C10_bmean[imom2] += C10_b[ibin*MAXMOMSQ+imom2];
      }
      Ep_bmean[imom2] /= Nbins;
      A10_bmean[imom2] /= Nbins;
      B10_bmean[imom2] /= Nbins;
      C10_bmean[imom2] /= Nbins;
    }

    for(int imom2 = 0 ; imom2 < MAXMOMSQ ; imom2++){
      for(int ibin = 0 ; ibin < Nbins ; ibin++){
	temp_Ep_bmean += (Ep_bmean[imom2] - Ep_b[ibin*MAXMOMSQ+imom2])*(Ep_bmean[imom2] - Ep_b[ibin*MAXMOMSQ+imom2]); 
	temp_A10_bmean += (A10_bmean[imom2] - A10_b[ibin*MAXMOMSQ+imom2])*(A10_bmean[imom2] - A10_b[ibin*MAXMOMSQ+imom2]); 
	temp_B10_bmean += (B10_bmean[imom2] - B10_b[ibin*MAXMOMSQ+imom2])*(B10_bmean[imom2] - B10_b[ibin*MAXMOMSQ+imom2]); 
	temp_C10_bmean += (C10_bmean[imom2] - C10_b[ibin*MAXMOMSQ+imom2])*(C10_bmean[imom2] - C10_b[ibin*MAXMOMSQ+imom2]); 
      }
      dA10[imom2] = sqrt( (Nbins-1)/( (double) Nbins ) * temp_A10_bmean  );
      dB10[imom2] = sqrt( (Nbins-1)/( (double) Nbins ) * temp_B10_bmean  );
      dC10[imom2] = sqrt( (Nbins-1)/( (double) Nbins ) * temp_C10_bmean  );

      temp_Ep_bmean=0.;
      temp_A10_bmean=0.;
      temp_B10_bmean=0.;
      temp_C10_bmean=0.;
    }
    ////////////////////////////////////////////

    // calculate data for ratio to print out //

    //printOut R zero momentum

    //    for(int it = 0 ; it < TSINK+1 ; it++)
    //  fprintf(out_R_mom0,"%d %+.8f %+.8f\n",it, RR_bmean[it][0][0].real() , dRR[it][0][0].real());





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
    
    printf("\nA10 FORM FACTOR\n");
    for(int imom2 = 0 ; imom2 < MAXMOMSQ; imom2++)
      printf("%f %f %f\n",q2GeV[imom2],A10_bmean[imom2],dA10[imom2]);

    printf("\nB10 FORM FACTOR\n");
    for(int imom2 = 0 ; imom2 < MAXMOMSQ; imom2++)
      printf("%f %f %f\n",q2GeV[imom2],B10_bmean[imom2],dB10[imom2]);

    printf("\nC10 FORM FACTOR\n");
    for(int imom2 = 0 ; imom2 < MAXMOMSQ; imom2++)
      printf("%f %f %f\n",q2GeV[imom2],C10_bmean[imom2],dC10[imom2]);

    //////////////////////////

    for(int iconf=0;iconf<Nconfs;iconf++){
      free(confString[iconf]);
      free(threep_namesType1[iconf]);
      free(threep_namesType2[iconf]);
      free(twop_names[iconf]);
    }

    free(confString);
    free(threep_namesType1);
    free(threep_namesType2);
    free(twop_names);

    for(int i = 0 ; i < 2 ; i++) free(threep[i]);

    free(twop);
    for(int i = 0 ; i < 2 ; i++) free(RR_mean[i]);
    free(twop_q2);
    for(int i = 0 ; i < 2 ; i++) free(threep_mean[i]);
    free(twop_q2_mean);
    for(int i = 0 ; i < 2 ; i++) free(threep_b[i]);
    free(twop_q2_b);
    free(MEff_b);

    for(int i = 0 ; i < 2 ; i++){
      free(RR_b[i]);
      free(RRplateau_mean[i]);
      free(RRplateau_b[i]);
      free(RRplateau_bmean[i]);
      free(dRRplateau[i]);
      free(RR_plateau_test[i]);
    }

    free(mass_b);
    free(Ep_b);
    free(A10_b);
    free(B10_b);
    free(C10_b);
    free(Ep_bmean);
    free(A10_bmean);
    free(B10_bmean);
    free(C10_bmean);
    free(dA10);
    free(dB10);
    free(dC10);
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

Complex fit_ratio_plato2(Complex* ratio, Complex error_ratio[TSINK+1][4*4][NMOM], int fit_plato_low, int fit_plato_high, int mu, int nu, int imom){
  int fitrange;
  Complex yi[fit_plato_high - fit_plato_low + 1];
  Complex sigmai[fit_plato_high - fit_plato_low + 1];
  Complex S , Sy;
  Complex par_constant;

  fitrange = fit_plato_high - fit_plato_low + 1;

  for(int it = 0; it < fitrange ; it++){
    yi[it] = ratio[(fit_plato_low+it)*4*4*NMOM+mu*4*NMOM+nu*NMOM+imom] ;
    sigmai[it] = error_ratio[fit_plato_low+it][mu*4+nu][imom] ;
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

void extract_tensor(double *A10,double *B10,double *C10,Complex *RR[],Complex *dRR[], double *Ep,int momList[][NMOM], int p2[], int numberMomPerQ2[] ){

  double epsilon[4][4][4][4] = {{{{0.}}}};
  epsilon[0][1][2][3] = 1.;
  epsilon[0][1][3][2] = -1.;
  epsilon[0][2][1][3] = -1.;
  epsilon[0][2][3][1] = 1.;
  epsilon[0][3][1][2] = 1.;
  epsilon[0][3][2][1] = -1.;

  epsilon[1][0][2][3] = -1.;
  epsilon[1][0][3][2] = 1.;
  epsilon[1][2][0][3] = 1.;
  epsilon[1][2][3][0] = -1.;
  epsilon[1][3][0][2] = -1.;
  epsilon[1][3][2][0] = 1.;

  epsilon[2][0][1][3] = 1.;
  epsilon[2][0][3][1] = -1.;
  epsilon[2][1][0][3] = -1.;
  epsilon[2][1][3][0] = 1.;
  epsilon[2][3][0][1] = 1.;
  epsilon[2][3][1][0] = -1.;

  epsilon[3][0][1][2] = -1.;
  epsilon[3][0][2][1] = 1.;
  epsilon[3][1][0][2] = 1.;
  epsilon[3][1][2][0] = -1.;
  epsilon[3][2][0][1] = -1.;
  epsilon[3][2][1][0] = 1.;

  double C;
  Complex *CoeffA10[MAXMOMSQ] = {NULL};
  Complex *CoeffB10[MAXMOMSQ] = {NULL};
  Complex *CoeffC10[MAXMOMSQ] = {NULL};

  double *b[MAXMOMSQ] = {NULL};
  double *db[MAXMOMSQ] = {NULL};
  double *A[MAXMOMSQ] = {NULL}; // important because of linkage with fortran matrices rows and columns are stored different in memory
  double *u[MAXMOMSQ] = {NULL};
  double *w[MAXMOMSQ] = {NULL};
  double *vt[MAXMOMSQ] = {NULL};
  double *oneOverW[MAXMOMSQ] = {NULL};
  double *x[MAXMOMSQ] = {NULL};


  for(int imom2 = 0 ; imom2 < MAXMOMSQ ; imom2++){
    A10[imom2] = 0.;
    B10[imom2] = 0.;
    C10[imom2] = 0.;
    int counter = 0;
    if(numberMomPerQ2[imom2] != 0){ 
      C = sqrt( ( 2.*Ep[0]*Ep[0] ) / ( Ep[imom2] * ( Ep[imom2]+Ep[0] ) ) );    

      CoeffA10[imom2] = (Complex*)calloc(numberMomPerQ2[imom2]*9,sizeof(Complex));
      CoeffB10[imom2] = (Complex*)calloc(numberMomPerQ2[imom2]*9,sizeof(Complex));
      CoeffC10[imom2] = (Complex*)calloc(numberMomPerQ2[imom2]*9,sizeof(Complex));

      b[imom2] = (double*)calloc(numberMomPerQ2[imom2]*9,sizeof(double));
      db[imom2] = (double*)calloc(numberMomPerQ2[imom2]*9,sizeof(double));
      A[imom2] = (double*)calloc(3*numberMomPerQ2[imom2]*9,sizeof(double)); // dimension (# measurements) x 2 // two form factors


      // import coefficients for equations ++++++++++++++++++++++++++++++++++++++++++++++++++++++
      for(int imom = 0 ; imom < NMOM ; imom++){
	int q2 = momList[0][imom]*momList[0][imom] + momList[1][imom]*momList[1][imom] + momList[2][imom]*momList[2][imom];                                                  
	if(q2 == imom2){

	  // P^0j_0
	  for(int j = 1 ; j < 4 ; j++){
	    CoeffA10[imom2][counter].real() = -C*( (momList[j-1][imom]*(2.*PI/L))/(2.*Ep[0]) );
	    CoeffB10[imom2][counter].real() = -C*( (momList[j-1][imom]*(2.*PI/L))*(2.*Ep[0])/(4.*Ep[0]*Ep[0]) );
	    CoeffC10[imom2][counter].real() = -C*( 2.*(Ep[imom2]+Ep[0])*(momList[j-1][imom]*(2.*PI/L)) /(4.*Ep[0]*Ep[0])  );
	    counter++;
	  }

	  // P^0j_(1+2+3)
	  for(int j = 1 ; j < 4 ; j++){
	    for(int l = 1 ; l < 4 ; l++) CoeffA10[imom2][counter].real() += C*( ((epsilon[0][j][1][l]+epsilon[0][j][2][l]+epsilon[0][j][3][l])*(momList[l-1][imom]*(2.*PI/L)))/(2.*Ep[0]) );
	    for(int l = 1 ; l < 4 ; l++) CoeffB10[imom2][counter].real() += C*( (((epsilon[0][j][1][l]+epsilon[0][j][2][l]+epsilon[0][j][3][l])*(momList[l-1][imom]*(2.*PI/L)))*(Ep[0] - Ep[imom2] )) / (4.*Ep[0]*Ep[0])  );
	    counter++;
	  }
	  
	  // P^ij_(1+2+3)
	  for(int i = 1 ; i < 3 ; i++)
	    for(int j = i+1 ; j < 4 ; j++){
	      CoeffA10[imom2][counter].real() = C*( (epsilon[i][j][1][0] + epsilon[i][j][2][0] + epsilon[i][j][3][0]) * (Ep[imom2]+Ep[0]) / (2.*Ep[0]) );
	      for(int l = 1 ; l < 4 ; l++) CoeffB10[imom2][counter].real() += C*( ( (epsilon[j][1][0][l]+epsilon[j][2][0][l]+epsilon[j][3][0][l])*(momList[i-1][imom]*(2.*PI/L)) - (epsilon[i][1][0][l]+epsilon[i][2][0][l]+epsilon[i][3][0][l])*(momList[j-1][imom]*(2.*PI/L)) ) * momList[l-1][imom]*(2.*PI/L)  / (4.*Ep[0]*Ep[0])  );
	      counter++;
	    }
	  
	}
	
      }
      // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      counter = 0;
      // import data from 3pf ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      for(int imom = 0 ; imom < NMOM ; imom++){
	int q2 = momList[0][imom]*momList[0][imom] + momList[1][imom]*momList[1][imom] + momList[2][imom]*momList[2][imom];                                                  
	if(q2 == imom2){

	  // P^0j_0
	  for(int j = 1 ; j < 4 ; j++){
	    b[imom2][counter] = imag(RR[0][0*4*NMOM+j*NMOM+imom]);
	    db[imom2][counter] = imag(dRR[0][0*4*NMOM+j*NMOM+imom]);
	    A[imom2][0*numberMomPerQ2[imom2]*9+counter] = real(CoeffA10[imom2][counter]);
	    A[imom2][1*numberMomPerQ2[imom2]*9+counter] = real(CoeffB10[imom2][counter]);
	    A[imom2][2*numberMomPerQ2[imom2]*9+counter] = real(CoeffC10[imom2][counter]);
	    counter++;
	  }

	  
	  // P^0j_(1+2+3)
	  for(int j = 1 ; j < 4 ; j++){
	    b[imom2][counter] = real(RR[1][0*4*NMOM+j*NMOM+imom]);
	    db[imom2][counter] = real(dRR[1][0*4*NMOM+j*NMOM+imom]);
	    A[imom2][0*numberMomPerQ2[imom2]*9+counter] = real(CoeffA10[imom2][counter]);
	    A[imom2][1*numberMomPerQ2[imom2]*9+counter] = real(CoeffB10[imom2][counter]);
	    A[imom2][2*numberMomPerQ2[imom2]*9+counter] = real(CoeffC10[imom2][counter]);
	    counter++;	    
	  }

	  // P^ij_(1+2+3)
	  for(int i = 1 ; i < 3 ; i++)
	    for(int j = i+1 ; j < 4 ; j++){
	      b[imom2][counter] = imag(RR[1][i*4*NMOM+j*NMOM+imom]);
	      db[imom2][counter] = imag(dRR[1][i*4*NMOM+j*NMOM+imom]);
	      A[imom2][0*numberMomPerQ2[imom2]*9+counter] = real(CoeffA10[imom2][counter]);
	      A[imom2][1*numberMomPerQ2[imom2]*9+counter] = real(CoeffB10[imom2][counter]);
	      A[imom2][2*numberMomPerQ2[imom2]*9+counter] = real(CoeffC10[imom2][counter]);
	      counter++;	    	      
	    }

	}      
      }
      
      // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      //      if(imom2 == 1){
      //	for(int i = 0 ; i < counter ; i++)
      //	  printf("%f %f %f\n",A[imom2][0*numberMomPerQ2[imom2]*9+i],A[imom2][1*numberMomPerQ2[imom2]*9+i],A[imom2][2*numberMomPerQ2[imom2]*9+i]);
      // }
      // make weighted least square +++++++++++++++++++++++++++++++++++++++++++++++++++++++
      for(int i = 0 ; i < counter ; i++){
	b[imom2][i] /= fabs(db[imom2][i]);
	A[imom2][0*numberMomPerQ2[imom2]*9+i] /= fabs(db[imom2][i]);
	A[imom2][1*numberMomPerQ2[imom2]*9+i] /= fabs(db[imom2][i]);
	A[imom2][2*numberMomPerQ2[imom2]*9+i] /= fabs(db[imom2][i]);
      }

      // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


      //      if(imom2 == 1){
      //	for(int i = 0 ; i < counter ; i++)
      //	  printf("%d %+f %+f %+f\n",i,CoeffA10[imom2][i].real(),CoeffB10[imom2][i].real(),CoeffC10[imom2][i].real());
      // }

      fflush(stdout);



      counter = 0;

      // ready for singular value decomposition

      int info;
      int lwork;
      double wkopt;
      double* work=NULL;
      int m=numberMomPerQ2[imom2]*9;
      int n=3;
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
      
      oneOverW[imom2][0*3+0] = 1./w[imom2][0];
      oneOverW[imom2][1*3+1] = 1./w[imom2][1];
      oneOverW[imom2][2*3+2] = 1./w[imom2][2];
      if(imom2 == 0) oneOverW[imom2][1*3+1] = 0.;
      if(imom2 == 0) oneOverW[imom2][2*3+2] = 0.;

      for(int alpha = 0 ; alpha < 3 ; alpha++)
	for(int beta = 0 ; beta < 3 ; beta++)
	  for(int gamma = 0 ; gamma < 3 ; gamma++)
	    for(int delta = 0 ; delta < numberMomPerQ2[imom2]*9 ; delta++){
	      x[imom2][alpha] += vt[imom2][alpha*3+beta] * oneOverW[imom2][beta*3+gamma] * u[imom2][gamma*m+delta] * b[imom2][delta]; 
	    }


      A10[imom2] = x[imom2][0]*ZT;
      B10[imom2] = (imom2 == 0) ? log(-1) : x[imom2][1]*ZT;
      C10[imom2] = (imom2 == 0) ? log(-1) : x[imom2][2]*ZT;
	
    }
    else{
      A10[imom2] = log(-1); // assing nan
      B10[imom2] = log(-1); // assing nan
      C10[imom2] = log(-1);
    }


  } // close for over momenta

}


