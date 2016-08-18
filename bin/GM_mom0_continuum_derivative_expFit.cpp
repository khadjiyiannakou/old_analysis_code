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
#include <mpfit.h>

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

#define L 48
#define T 96
#define Z_V 1
#define BINSIZE 1
#define TSINK 14
#define LOW2PT 10
#define HIGH2PT 17
#define LOW3PT 4
#define HIGH3PT 10
#define AINV (0.197/0.093)
#define AMNPHYS (0.9382720/AINV)
#define FITMAX 4
#define NMOM 1
#define MAXMOMSQ 1 // ATTENTION 7,15 missing

//#define NMOM 93
//#define MAXMOMSQ 9
#define POINTSPERTIMESLICE 100

#define PI 3.14159265359
#define MAX_STRING 257

// forward declaration
double fit_ratio_plato(double* ratio, double error_ratio[TSINK+1], int fit_plato_low, int fit_plato_high);

//Complex fit_ratio_plato(Complex* ratio, Complex error_ratio[TSINK+1][4][NMOM], int fit_plato_low, int fit_plato_high, int mu, int imom);
double fit_EffMass_plato(double MEff[], double dMEff[],int fit_plato_low, int fit_plato_high);

void fit_ratio_exp(double* ratio, double error_ratio[TSINK+1], int fit_low, int fit_high, double fit_values[2]);

void ratioGM_mom0(double *GM,Complex (*RRM_b)[3][3],double *Ep_b,int ibin);

int main(int argc, char *argv[]){

  clock_t time;

  // Greeting //
  printf("This program calculate GM at zero momentum transfer using continuum derivative\n");
  printf("Right inputs order must be\n");
  printf("(1) Executable , (2) Trajectory list, (3) Threep G5G1, (4) Threep G5G2, (5) Threep G5G3, (6) Twop base name, (7) Output name\n\n ");
  //////////////////


  // check passing inputs right //
  if(argc != 7){
    fprintf(stderr,"Error: Wrong number of input files \n");
    exit(EXIT_FAILURE);
  }
  //////////////////////


  // fill filenames and print paths //
  char filename_Traj[MAX_STRING];
  char filename_Threep[3][MAX_STRING];
  char filename_Twop[MAX_STRING];
  char filename_Output[MAX_STRING];

  //output files
  char filename_RGM_mom0[MAX_STRING];
  char filename_RGM_mom0_band[MAX_STRING];

  strcpy(filename_Traj,argv[1]);
  strcpy(filename_Threep[0],argv[2]);
  strcpy(filename_Threep[1],argv[3]);
  strcpy(filename_Threep[2],argv[4]);

  strcpy(filename_Twop,argv[5]);
  strcpy(filename_Output,argv[6]);

  printf("Got filename for trajectories : %s\n",filename_Traj);
  printf("Got filename for threep G5G1 : %s\n",filename_Threep[0]);
  printf("Got filename for threep G5G2 : %s\n",filename_Threep[1]);
  printf("Got filename for threep G5G3 : %s\n",filename_Threep[2]);

  printf("Got filename for twop : %s\n",filename_Twop);
  printf("Got filename for output : %s\n\n",filename_Output);

  sprintf(filename_RGM_mom0,"%s_%s.dat",filename_Output,"RGM_mom0");
  sprintf(filename_RGM_mom0_band,"%s_%s.dat",filename_Output,"RGM_mom0_band");

  FILE  *out_RGM_mom0 = NULL;
  FILE  *out_RGM_mom0_band = NULL;

  out_RGM_mom0 = fopen(filename_RGM_mom0,"w");
  if(out_RGM_mom0 == NULL){ fprintf(stderr,"Error open files for writting\n"); exit(EXIT_FAILURE);}

  out_RGM_mom0_band = fopen(filename_RGM_mom0_band,"w");
  if(out_RGM_mom0_band == NULL){ fprintf(stderr,"Error open files for writting\n"); exit(EXIT_FAILURE);}  

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

  char **threep_names[3];
  char **twop_names;

  for(int i = 0 ; i < 3 ; i++)
    threep_names[i] = (char**)malloc(Nconfs*sizeof(char*));

  twop_names = (char**)malloc(Nconfs*sizeof(char*));

  for(int iconf=0;iconf<Nconfs;iconf++){
    for(int i = 0 ; i < 3 ; i++)
      threep_names[i][iconf] = (char*)malloc(MAX_STRING*sizeof(char));
    twop_names[iconf] = (char*)malloc(MAX_STRING*sizeof(char));
  }

  for(int iconf=0;iconf<Nconfs;iconf++){
    for(int i = 0 ; i < 3 ; i++)
      sprintf(threep_names[i][iconf],"%s.%s_isov",filename_Threep[i],confString[iconf]);
    sprintf(twop_names[iconf],"%s.%s",filename_Twop,confString[iconf]);
  }

  ///////////////////////////////////////////////


  // allocate memory & read data //
  Complex (*threep)[3][3] =  (Complex(*)[3][3])malloc(3*3*Nconfs*(TSINK+1)*4*NMOM*sizeof(Complex)); // index for the projector and index for the derivative direction
  Complex *twop = (Complex*)malloc(Nconfs*T*NMOM*sizeof(Complex));
  int momList[3][NMOM];
  int dummy;

  if(threep == NULL || twop == NULL){ fprintf(stderr,"Error: Out of memory\n"); exit(EXIT_FAILURE);}
  FILE *ptr_threep[3],*ptr_twop;

  time = clock();
  for(int iconf = 0 ; iconf < Nconfs ; iconf++){
    ptr_threep[0] = fopen(threep_names[0][iconf],"r");
    ptr_threep[1] = fopen(threep_names[1][iconf],"r");
    ptr_threep[2] = fopen(threep_names[2][iconf],"r");

    ptr_twop = fopen(twop_names[iconf],"r");
    if(ptr_threep[0] == NULL || ptr_threep[1] == NULL || ptr_threep[2] == NULL  || ptr_twop == NULL){ fprintf(stderr,"Error: open files for reading\n"); exit(EXIT_FAILURE);}
    // read threep

    for(int i = 0 ; i < 3 ; i++)
      for(int dqi = 0 ; dqi < 3 ; dqi++)
	for(int it=0; it < TSINK+1 ; it++)
	  for(int mu = 0 ; mu < 4 ; mu++)
	    for(int imom = 0 ; imom < NMOM ; imom++)
	      returnValue = fscanf(ptr_threep[i],"%d %d %d %d %lf %lf %d",&dummy,&(momList[0][imom]), &(momList[1][imom]), &(momList[2][imom]), &(threep[iconf*(TSINK+1)*4*NMOM + it*4*NMOM + mu*NMOM + imom][i][dqi].real()), &(threep[iconf*(TSINK+1)*4*NMOM + it*4*NMOM + mu*NMOM + imom][i][dqi].imag()), &dummy );


    for(int it=0; it < T ; it++)
      for(int imom = 0 ; imom < NMOM ; imom++)
	returnValue = fscanf(ptr_twop,"%d %d %d %d %lf %lf",&dummy,&(momList[0][imom]), &(momList[1][imom]), &(momList[2][imom]),&(twop[iconf*T*NMOM + it*NMOM + imom].real()), &(twop[iconf*T*NMOM + it*NMOM + imom].imag()) );


    for(int i = 0 ; i < 3 ; i++)
      fclose(ptr_threep[i]);
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
  Complex (*threep_mean)[3][3] = (Complex(*)[3][3])calloc(3*3*(TSINK+1)*4*NMOM,sizeof(Complex));
  Complex *twop_q2_mean = (Complex*)calloc(T*MAXMOMSQ,sizeof(Complex));

  for(int i = 0 ; i < 3 ; i++)
    for(int dqi = 0 ; dqi < 3 ; dqi++)
      for(int it = 0 ; it < TSINK+1 ; it++)
	for(int mu = 0 ; mu < 4 ; mu++)
	  for(int imom = 0 ; imom < NMOM ; imom++){
	    for(int iconf = 0 ; iconf < Nconfs ; iconf++){
	      threep_mean[it*4*NMOM + mu*NMOM + imom][i][dqi] = threep_mean[it*4*NMOM + mu*NMOM + imom][i][dqi] + threep[iconf*(TSINK+1)*4*NMOM + it*4*NMOM + mu*NMOM + imom][i][dqi];
	    }
	    
	    threep_mean[it*4*NMOM + mu*NMOM + imom][i][dqi].real() /= Nconfs;
	    threep_mean[it*4*NMOM + mu*NMOM + imom][i][dqi].imag() /= Nconfs;

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

  Complex (*RR_mean)[3][3] = (Complex(*)[3][3])malloc(3*3*(TSINK+1)*4*NMOM*sizeof(Complex));
  Complex squareRoot;
  double subRoot;

  for(int i = 0 ; i < 3 ; i++)
    for(int dqi = 0 ; dqi < 3 ; dqi++)
      for(int it=0;it < TSINK+1;it++)
	for(int mu = 0 ; mu < 4 ; mu++)
	  for(int imom = 0 ; imom < NMOM ; imom++){

	    
	    subRoot = (twop_q2_mean[(TSINK-it)*MAXMOMSQ+p2[imom]].real()  * twop_q2_mean[(it)*MAXMOMSQ+0].real() *  twop_q2_mean[TSINK*MAXMOMSQ+0].real() ) /(twop_q2_mean[(TSINK-it)*MAXMOMSQ+0].real()  * twop_q2_mean[(it)*MAXMOMSQ+p2[imom]].real() *  twop_q2_mean[TSINK*MAXMOMSQ+p2[imom]].real() );
	    squareRoot = std::sqrt( (Complex)  subRoot  );

	    RR_mean[it*4*NMOM+mu*NMOM+imom][i][dqi] = ( threep_mean[it*4*NMOM + mu*NMOM + imom][i][dqi] / twop_q2_mean[TSINK*MAXMOMSQ+0].real() )* squareRoot;
	  }

  //////////////////////////////


  // binning for jackknife analysis //
  Complex (*threep_b)[3][3] = (Complex(*)[3][3])calloc(3*3*Nbins*(TSINK+1)*4*NMOM,sizeof(Complex));
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

	for(int i = 0 ; i < 3 ; i++)
	  for(int dqi = 0 ; dqi < 3 ; dqi++)
	    for(int it = 0 ; it < TSINK+1 ; it++)
	      for(int mu = 0 ; mu < 4 ; mu++)
		for(int imom = 0 ; imom < NMOM ; imom++){
		  threep_b[ibin*(TSINK+1)*4*NMOM + it*4*NMOM + mu*NMOM + imom][i][dqi] = threep_b[ibin*(TSINK+1)*4*NMOM + it*4*NMOM + mu*NMOM + imom][i][dqi] + threep[iconf*(TSINK+1)*4*NMOM + it*4*NMOM + mu*NMOM + imom][i][dqi];
		}

	for(int it = 0 ; it < T ; it++)
	  for(int imom2 = 0 ; imom2 < MAXMOMSQ; imom2++)
	    twop_q2_b[ibin*T*MAXMOMSQ+it*MAXMOMSQ+imom2] =  twop_q2_b[ibin*T*MAXMOMSQ+it*MAXMOMSQ+imom2] + twop_q2[iconf*T*MAXMOMSQ + it*MAXMOMSQ + imom2 ];
      }
   
      
    }

  for(int i = 0 ; i < 3 ; i++)
    for(int dqi = 0 ; dqi < 3 ; dqi++)
      for(int it = 0 ; it < TSINK+1 ; it++)
	for(int mu = 0 ; mu < 4 ; mu++)
	  for(int imom = 0 ; imom < NMOM ; imom++){
	    threep_b[ibin*(TSINK+1)*4*NMOM + it*4*NMOM + mu*NMOM + imom][i][dqi].real() /= (Nconfs-BINSIZE);
	    threep_b[ibin*(TSINK+1)*4*NMOM + it*4*NMOM + mu*NMOM + imom][i][dqi].imag() /= (Nconfs-BINSIZE);
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
  Complex (*RR_b)[3][3] = (Complex(*)[3][3])calloc(3*3*Nbins*(TSINK+1)*4*NMOM,sizeof(Complex));

  for(int ibin = 0 ; ibin < Nbins ; ibin++)
    for(int it = 1 ; it < T ; it++){
      MEff_b[ibin*T+it] = log(twop_q2_b[ibin*T*MAXMOMSQ+it*MAXMOMSQ+0].real() / twop_q2_b[ibin*T*MAXMOMSQ+(it+1)*MAXMOMSQ+0].real());
    }


  for(int i = 0 ; i < 3 ; i++)
    for(int dqi = 0 ; dqi < 3 ; dqi++)
      for(int ibin = 0 ; ibin < Nbins ; ibin++)
	for(int it=0;it < TSINK+1;it++)
	  for(int mu = 0 ; mu < 4 ; mu++)
	    for(int imom = 0 ; imom < NMOM ; imom++){
	      

	      subRoot = (twop_q2_b[ibin*T*MAXMOMSQ+(TSINK-it)*MAXMOMSQ+p2[imom]].real()  * twop_q2_b[ibin*T*MAXMOMSQ+(it)*MAXMOMSQ+0].real() *  twop_q2_b[ibin*T*MAXMOMSQ+TSINK*MAXMOMSQ+0].real() ) /(twop_q2_b[ibin*T*MAXMOMSQ+(TSINK-it)*MAXMOMSQ+0].real()  * twop_q2_b[ibin*T*MAXMOMSQ+(it)*MAXMOMSQ+p2[imom]].real() *  twop_q2_b[ibin*T*MAXMOMSQ+TSINK*MAXMOMSQ+p2[imom]].real() );
	      squareRoot = std::sqrt( (Complex)  subRoot  );

	  RR_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM+mu*NMOM+imom][i][dqi] = ( threep_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM + mu*NMOM + imom][i][dqi] / twop_q2_b[ibin*T*MAXMOMSQ+TSINK*MAXMOMSQ+0].real() )* squareRoot;
      
	}

  ///////////////////////////////////// 

  // mean bin values //
  double MEff_bmean[T]={};
  Complex RR_bmean[TSINK+1][4][NMOM][3][3]={};


  for(int it = 0 ; it < T ; it++){
    for(int ibin = 0 ; ibin < Nbins ; ibin++){
      MEff_bmean[it] += MEff_b[ibin*T+it];
    }
    MEff_bmean[it] /= Nbins;
  }

  double dMEff[T];
  double temp_dMEff = 0.;
  for(int it = 0 ; it < T ; it++){
    for(int ibin =0 ; ibin < Nbins ; ibin++){
      temp_dMEff = temp_dMEff + (MEff_b[ibin*T+it] - MEff_bmean[it] )*(MEff_b[ibin*T+it] - MEff_bmean[it] );
    }
    dMEff[it] = sqrt((Nbins-1)/( (double) Nbins ) * temp_dMEff);
    temp_dMEff=0;
  }

  //  for(int it = 0 ; it < T ; it++)
  //  printf("%d %+e %+e\n",it, MEff_bmean[it], dMEff[it]);


  for(int i = 0 ; i < 3 ; i++)
    for(int dqi = 0 ; dqi < 3 ; dqi++)
      for(int it = 0 ; it < TSINK+1 ; it++)
	for(int mu = 0 ; mu < 4 ; mu++)
	  for(int imom = 0 ; imom < NMOM ; imom++){
	    for(int ibin = 0 ; ibin < Nbins ; ibin++){
	      RR_bmean[it][mu][imom][i][dqi] = RR_bmean[it][mu][imom][i][dqi] + RR_b[ibin*(TSINK+1)*4*NMOM+it*4*NMOM+mu*NMOM+imom][i][dqi];
	    }
	    RR_bmean[it][mu][imom][i][dqi].real() /= Nbins;
	    RR_bmean[it][mu][imom][i][dqi].imag() /= Nbins;

	  }


    // extract quantities for each bin //
    double *mass_b = (double*)calloc(Nbins,sizeof(double));
    double *Ep_b = (double*)calloc(Nbins*MAXMOMSQ,sizeof(double));
    double *RGM_b = (double*)calloc(Nbins*(TSINK+1)*MAXMOMSQ,sizeof(double));

    double *Ep_bmean = (double*)calloc(MAXMOMSQ,sizeof(double));
    double *RGM_bmean = (double*)calloc((TSINK+1)*MAXMOMSQ,sizeof(double));
    double *dRGM = (double*)calloc((TSINK+1)*MAXMOMSQ,sizeof(double));

    for(int ibin = 0 ; ibin < Nbins ; ibin++)
      mass_b[ibin] = fit_EffMass_plato(MEff_b+ibin*T, dMEff,LOW2PT, HIGH2PT);

    for(int ibin = 0 ; ibin < Nbins ; ibin++)
      for(int imom = 0 ; imom < NMOM ; imom++)
	Ep_b[ibin*MAXMOMSQ+p2[imom]] = sqrt(mass_b[ibin]*mass_b[ibin] + p2[imom]*(2.*PI/L)*(2.*PI/L) );

    for(int ibin = 0 ; ibin < Nbins ; ibin++)
      Ep_bmean[0] += Ep_b[0];

    Ep_bmean[0] /= Nbins;

    printf(" mass = %+e\n",Ep_bmean[0]);

    for(int ibin = 0 ; ibin < Nbins ; ibin++){
      ratioGM_mom0(RGM_b ,RR_b,Ep_b ,ibin);
    }    

    ////////////////////////////////////
    
    for(int it = 0 ; it < TSINK+1 ; it++){
      for(int ibin = 0 ; ibin < Nbins ; ibin++)
	RGM_bmean[it] += RGM_b[ibin*(TSINK+1) + it];
      RGM_bmean[it] /= Nbins;
    }
      
    double temp_error=0.;

    for(int it = 0 ; it < TSINK+1 ; it++){
      for(int ibin = 0 ; ibin < Nbins ; ibin++){
	temp_error += (RGM_bmean[it] - RGM_b[ibin*(TSINK+1) + it])*(RGM_bmean[it] - RGM_b[ibin*(TSINK+1) + it]); 
      }
      dRGM[it] = sqrt( (Nbins-1)/( (double) Nbins ) * temp_error  );

      temp_error=0.;
    }

   
    for(int it = 0 ; it < TSINK+1 ; it++)
      fprintf(out_RGM_mom0,"%d %+.8f %+.8f\n",it, RGM_bmean[it] , dRGM[it]);


    double *GM0_b = (double*)calloc(Nbins,sizeof(double));
    double *DE_b = (double*)calloc(Nbins,sizeof(double));

    double GM0_bmean=0.;
    double GM0_error=0.;
    double temp_GM0_error=0.;

    double DE_bmean=0.;
    double DE_error=0.;
    double temp_DE_error=0.;

    double fit_values[2] = {0};
    for(int ibin = 0 ; ibin < Nbins ; ibin++){
      fit_ratio_exp(RGM_b + ibin*(TSINK+1)*MAXMOMSQ, dRGM, LOW3PT, HIGH3PT,fit_values);
      GM0_b[ibin] = fit_values[0];
      DE_b[ibin] = fit_values[1];
    }

    for(int ibin = 0 ; ibin < Nbins ; ibin++)
      GM0_bmean += GM0_b[ibin];
    GM0_bmean /= Nbins;

    for(int ibin = 0 ; ibin < Nbins ; ibin++)
      DE_bmean += DE_b[ibin];
    DE_bmean /= Nbins;


    for(int ibin = 0 ; ibin < Nbins ; ibin++){
      temp_GM0_error += (GM0_bmean - GM0_b[ibin])*(GM0_bmean - GM0_b[ibin]); 
    }
    GM0_error = sqrt( (Nbins-1)/( (double) Nbins ) * temp_GM0_error  );

    for(int ibin = 0 ; ibin < Nbins ; ibin++){
      temp_DE_error += (DE_bmean - DE_b[ibin])*(DE_bmean - DE_b[ibin]); 
    }
    DE_error = sqrt( (Nbins-1)/( (double) Nbins ) * temp_DE_error  );

    fprintf(out_RGM_mom0,"\n[%d %d] %f %f \t %f %f\n",LOW3PT, HIGH3PT,GM0_bmean,GM0_error, DE_bmean, DE_error);


    // calcualte the band for the exponential fit

    double *yy_b = (double*)calloc(Nbins*(HIGH3PT-LOW3PT)*POINTSPERTIMESLICE,sizeof(double));
    double *yy_bmean = (double*)calloc((HIGH3PT-LOW3PT)*POINTSPERTIMESLICE,sizeof(double));
    double *yy_error = (double*)calloc((HIGH3PT-LOW3PT)*POINTSPERTIMESLICE,sizeof(double));
    double *xx = (double*)calloc((HIGH3PT-LOW3PT)*POINTSPERTIMESLICE,sizeof(double));

    for(int i = LOW3PT ; i < HIGH3PT ; i++)
      for(int j = 0 ; j < POINTSPERTIMESLICE ; j++)
	xx[i*POINTSPERTIMESLICE+j] = i + j/((double)POINTSPERTIMESLICE);

    for(int ibin = 0 ; ibin < Nbins ; ibin++)
      for(int i = LOW3PT ; i < HIGH3PT ; i++)
	for(int j = 0 ; j < POINTSPERTIMESLICE ; j++)
	  yy_b[ibin*(HIGH3PT-LOW3PT)*POINTSPERTIMESLICE+i*POINTSPERTIMESLICE+j] = GM0_b[ibin] * exp(-DE_b[ibin]*xx[i*POINTSPERTIMESLICE+j]);
    

    for(int i = LOW3PT ; i < HIGH3PT ; i++)
      for(int j = 0 ; j < POINTSPERTIMESLICE ; j++){
	for(int ibin = 0 ; ibin < Nbins ; ibin++)
	  yy_bmean[i*POINTSPERTIMESLICE+j] += yy_b[ibin*(HIGH3PT-LOW3PT)*POINTSPERTIMESLICE+i*POINTSPERTIMESLICE+j];
	yy_bmean[i*POINTSPERTIMESLICE+j] /= Nbins;
      }


    double temp_yy_err = 0.;

    for(int i = LOW3PT ; i < HIGH3PT ; i++)
      for(int j = 0 ; j < POINTSPERTIMESLICE ; j++){
	temp_yy_err = 0.;
	for(int ibin = 0 ; ibin < Nbins ; ibin++)
	  temp_yy_err += (yy_bmean[i*POINTSPERTIMESLICE+j] - yy_b[ibin*(HIGH3PT-LOW3PT)*POINTSPERTIMESLICE+i*POINTSPERTIMESLICE+j]) * (yy_bmean[i*POINTSPERTIMESLICE+j] - yy_b[ibin*(HIGH3PT-LOW3PT)*POINTSPERTIMESLICE+i*POINTSPERTIMESLICE+j]);
	yy_error[i*POINTSPERTIMESLICE+j] = sqrt( (Nbins-1)/( (double) Nbins ) * temp_yy_err  );
      }

    for(int i = LOW3PT ; i < HIGH3PT ; i++)
      for(int j = 0 ; j < POINTSPERTIMESLICE ; j++){
	fprintf(out_RGM_mom0_band,"%f %f %f\n",xx[i*POINTSPERTIMESLICE+j],yy_bmean[i*POINTSPERTIMESLICE+j]+yy_error[i*POINTSPERTIMESLICE+j], yy_bmean[i*POINTSPERTIMESLICE+j]-yy_error[i*POINTSPERTIMESLICE+j]);
      }



    /*
    for(int ibin = 0 ; ibin < Nbins ; ibin++){
      fitVal_b[ibin] = fit_ratio_plato(RGM_b + ibin*(TSINK+1)*MAXMOMSQ, dRGM, LOW3PT, HIGH3PT);
    }

    for(int ibin = 0 ; ibin < Nbins ; ibin++)
      fitVal_bmean += fitVal_b[ibin];
    fitVal_bmean /= Nbins;

    for(int ibin = 0 ; ibin < Nbins ; ibin++){
      temp_err_fitVal += (fitVal_bmean - fitVal_b[ibin])*(fitVal_bmean - fitVal_b[ibin]);
    }
    fitVal_error = sqrt( (Nbins-1)/( (double) Nbins ) * temp_err_fitVal  );

    fprintf(out_RGM_mom0, "\n Plateau Fit [%d %d] %f %f\n",LOW3PT, HIGH3PT,fitVal_bmean,fitVal_error);
    */

    

    // TODO
    // find dRGM_b and RGM_bmean
    // perform fit to take GM_b from RGM_b
    // find GM_bmean and dGM

    /*
    // find jackknife errors for form factors //
    double temp_Ep_bmean=0.;
    double temp_GM_bmean=0.;

    double *dGM = (double*)calloc(MAXMOMSQ,sizeof(double));
    double *dEp = (double*)calloc(MAXMOMSQ,sizeof(double));

    for(int imom2 = 0 ; imom2 < MAXMOMSQ ; imom2++){
      for(int ibin = 0 ; ibin < Nbins ; ibin++){
	Ep_bmean[imom2] += Ep_b[ibin*MAXMOMSQ+imom2];
	GM_bmean[imom2] += GM_b[ibin*MAXMOMSQ+imom2];
      }
      Ep_bmean[imom2] /= Nbins;
      GM_bmean[imom2] /= Nbins;
    }

    for(int imom2 = 0 ; imom2 < MAXMOMSQ ; imom2++){
      for(int ibin = 0 ; ibin < Nbins ; ibin++){
	temp_Ep_bmean += (Ep_bmean[imom2] - Ep_b[ibin*MAXMOMSQ+imom2])*(Ep_bmean[imom2] - Ep_b[ibin*MAXMOMSQ+imom2]); 
	temp_GM_bmean += (GM_bmean[imom2] - GM_b[ibin*MAXMOMSQ+imom2])*(GM_bmean[imom2] - GM_b[ibin*MAXMOMSQ+imom2]); 
      }
      dGM[imom2] = sqrt( (Nbins-1)/( (double) Nbins ) * temp_GM_bmean  );
      dEp[imom2] = sqrt( (Nbins-1)/( (double) Nbins ) * temp_Ep_bmean  );

      temp_Ep_bmean=0.;
      temp_GM_bmean=0.;
    }
    */
    ////////////////////////////////////////////

    // calculate data for ratio to print out //

    //printOut RE zero momentum

    //    for(int it = 0 ; it < TSINK+1 ; it++)
    //  fprintf(out_RE_mom0,"%d %+.8f %+.8f\n",it, RRE_bmean[it][0][0].real() , dRRE[it][0][0].real());


    // print results //
    /*
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
    

    printf("\nMAGNETIC FORM FACTOR\n");
    for(int imom2 = 0 ; imom2 < MAXMOMSQ; imom2++)
      printf("%f %f %f\n",q2GeV[imom2],GM_bmean[imom2],dGM[imom2]);
    */
    //////////////////////////


    // free memory //
    for(int iconf=0;iconf<Nconfs;iconf++){
      free(confString[iconf]);
      free(threep_names[0][iconf]);
      free(threep_names[1][iconf]);
      free(threep_names[2][iconf]);
      free(twop_names[iconf]);
    }
    free(confString);
    free(threep_names[0]);
    free(threep_names[1]);
    free(threep_names[2]);
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
    //    free(RR_bmean);
    free(mass_b);
    free(RGM_b);
    free(RGM_bmean);
    free(dRGM);
    
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

struct vars_struct {
  double *x;
  double *y;
  double *ey;
};

int expfunc(int m, int n, double *p, double *dy, double **dvec, void *vars)
{
  int i;
  struct vars_struct *v = (struct vars_struct *) vars;
  double *x, *y, *ey, f;

  x = v->x;
  y = v->y;
  ey = v->ey;

  for (i=0; i<m; i++) {
    f = p[0]*exp(-p[1]*x[i]);     /* Linear fit function; note f = a*exp(-b*x) */
    dy[i] = (y[i] - f)/ey[i];
  }

  return 0;
}

void fit_ratio_exp(double* ratio, double error_ratio[TSINK+1], int fit_low, int fit_high, double fit_values[2]){

  struct vars_struct v;
  int status;
  mp_result result;
  double perror[2];
  fit_values[0] = 0.;
  fit_values[1] = 0.;
  memset(&result,0,sizeof(result));       /* Zero results structure */
  result.xerror = perror;

  int numData = fit_high - fit_low + 1;
  double *x = (double*) malloc((TSINK+1)*sizeof(double));
  for(int i = 0 ; i < TSINK+1 ; i++) x[i] = i;

  v.x = x + fit_low;
  v.y = ratio + fit_low;
  v.ey = error_ratio + fit_low;

  status = mpfit(expfunc, numData, 2, fit_values, 0, 0, (void *) &v, &result);
  free(x);
}


double fit_ratio_plato(double* ratio, double error_ratio[TSINK+1], int fit_plato_low, int fit_plato_high){
  int fitrange;
  double yi[fit_plato_high - fit_plato_low + 1];
  double sigmai[fit_plato_high - fit_plato_low + 1];
  double S , Sy;
  double par_constant;

  fitrange = fit_plato_high - fit_plato_low + 1;

  for(int it = 0; it < fitrange ; it++){
    yi[it] = ratio[it] ;
    sigmai[it] = error_ratio[fit_plato_low+it] ;
  }

  S = 0.;
  Sy = 0.;

  for(int it = 0 ; it < fitrange ; it++){
    Sy += yi[it] / ( sigmai[it] * sigmai[it] );
    S += 1. / ( sigmai[it] * sigmai[it] );
  }

  par_constant = Sy / S ;

  return(par_constant);
}
/*
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
*/
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

void  ratioGM_mom0(double *RGM_b ,Complex (*RR_b)[3][3], double *Ep_b, int ibin){

  for(int imom = 0 ; imom < NMOM ; imom++)
    for(int it = 0 ; it < TSINK+1 ; it++){
      RGM_b[ibin*(TSINK+1)*MAXMOMSQ + it*MAXMOMSQ + imom] = (2*Ep_b[ibin*MAXMOMSQ+imom]) * (
									     +RR_b[ibin*(TSINK+1)*4*MAXMOMSQ + it*4*MAXMOMSQ + 1*MAXMOMSQ + imom][1][2].real()
									     -RR_b[ibin*(TSINK+1)*4*MAXMOMSQ + it*4*MAXMOMSQ + 2*MAXMOMSQ + imom][0][2].real()
									     +RR_b[ibin*(TSINK+1)*4*MAXMOMSQ + it*4*MAXMOMSQ + 3*MAXMOMSQ + imom][0][1].real()
									     -RR_b[ibin*(TSINK+1)*4*MAXMOMSQ + it*4*MAXMOMSQ + 1*MAXMOMSQ + imom][2][1].real()
									     +RR_b[ibin*(TSINK+1)*4*MAXMOMSQ + it*4*MAXMOMSQ + 2*MAXMOMSQ + imom][2][0].real()
									     -RR_b[ibin*(TSINK+1)*4*MAXMOMSQ + it*4*MAXMOMSQ + 3*MAXMOMSQ + imom][1][0].real()
									   )/6.;
    }

}
