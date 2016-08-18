QCD=../lib/libqcd.a
PATH_LIB_FIT=../lib/external_libs/cmpfit-1.2
PATH_LIBS=../lib/external_libs/lapack-3.4.0

INC=-I../include -I../lib/external_libs/lapack-3.4.0/lapacke/include -I$(PATH_LIB_FIT)
CC = g++
CFLAGS = -Wall -O3 $(INC)
LIBS = $(PATH_LIBS)/liblapack.a $(PATH_LIBS)/librefblas.a $(PATH_LIBS)/libtmglib.a $(PATH_LIB_FIT)/libmpfit.a -L/usr/local/software/jureca/Stage2/software/Core/GCC/4.9.3/lib64 -lgfortran
EXECUTABLES = EM_form_factor_noether EM_form_factor_local Axial_form_factor Axial_form_factor_Seperate_Projectors Vector_Derivative_form_factor Axial_Derivative_form_factor Tensor_form_factor Scalar_form_factor test_martha GM_mom0_continuum_derivative F3_mom0_continuum_derivative F3_mom0_continuum_derivative_expFit GM_mom0_continuum_derivative_expFit

all: $(EXECUTABLES)

F3_mom0_continuum_derivative_expFit: F3_mom0_continuum_derivative_expFit.o $(QCD)
	$(CC) $(CFLAGS) $^ -o $@ $(LIBS) -lm

F3_mom0_continuum_derivative: F3_mom0_continuum_derivative.o $(QCD)
	$(CC) $(CFLAGS) $^ -o $@ $(LIBS) -lm

GM_mom0_continuum_derivative: GM_mom0_continuum_derivative.o $(QCD)
	$(CC) $(CFLAGS) $^ -o $@ $(LIBS) -lm

GM_mom0_continuum_derivative_expFit: GM_mom0_continuum_derivative_expFit.o $(QCD)
	$(CC) $(CFLAGS) $^ -o $@ $(LIBS) -lm

test_martha: test_martha.o $(QCD)
	$(CC) $(CFLAGS) $^ -o $@ $(LIBS) -lm

EM_form_factor_noether: EM_form_factor_noether.o $(QCD)
	$(CC) $(CFLAGS) $^ -o $@ $(LIBS) -lm

EM_form_factor_local: EM_form_factor_local.o $(QCD)
	$(CC) $(CFLAGS) $^ -o $@ $(LIBS) -lm

Axial_form_factor: Axial_form_factor.o $(QCD)
	$(CC) $(CFLAGS) $^ -o $@ $(LIBS) -lm

Axial_form_factor_Seperate_Projectors: Axial_form_factor_Seperate_Projectors.o $(QCD)
	$(CC) $(CFLAGS) $^ -o $@ $(LIBS) -lm

Vector_Derivative_form_factor: Vector_Derivative_form_factor.o $(QCD)
	$(CC) $(CFLAGS) $^ -o $@ $(LIBS) -lm	

Axial_Derivative_form_factor: Axial_Derivative_form_factor.o $(QCD)
	$(CC) $(CFLAGS) $^ -o $@ $(LIBS) -lm

Tensor_form_factor: Tensor_form_factor.o $(QCD)
	$(CC) $(CFLAGS) $^ -o $@ $(LIBS) -lm

Scalar_form_factor: Scalar_form_factor.o $(QCD)
	$(CC) $(CFLAGS) $^ -o $@ $(LIBS) -lm

%.o: %.cpp 
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	-rm $(EXECUTABLES)
