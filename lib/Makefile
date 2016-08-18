INC= -I../include 
CC=g++
CFLAGS = -Wall -O3 $(INC)
QCD = libqcd.a

QCD_OBJS = nrutil.o svd.o 

QCD_HEADERS =  nrutil.h svd.h

all: $(QCD)

$(QCD): $(QCD_OBJS)
	ar cru $(QCD) $(QCD_OBJS)

clean:
	-rm -f *.o $(QCD)
