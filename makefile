# Makefile created by Pedro for the DissObjt with _PORT_FORTRAN_

#define the C compiler to use
FC = ifort
CPPFLAGS = -P -traditional -D_GLIBCXX_USE_CXX11_ABI=0
CC = gcc
CXX = gcc
CFLAGS = -D_GLIBCXX_USE_CXX11_ABI=0
CXXFLAGS = -D_GLIBCXX_USE_CXX11_ABI=0
LDFLAGS =
AR = ar
ARFLAGS = r
#LIBS = -L/cluster/home/pduarte/models2/ecodynamo/ecolib -lBGCFunctions
INCLUDE = -L/cluster/home/pduarte/pduarte8/EcodynamoLibFortran/src/BGCLibraries
#FFLAGS -mcmodel=large -xHOST #-Nmpi

SRC = /cluster/home/pduarte/pduarte8/EcodynamoLibFortran/src/Main.f90 /cluster/home/pduarte/pduarte8/EcodynamoLibFortran/src/BGCLibraries/BiogeochemicalProcesses.f90

Main.o:
	$(FC) -o Main.exe $(SRC)  
#Dissobjt.o:
#	$(CC) $(CFLAGS) $(INCLUDE) -o ../../ecolib/libdissobjt.so $(SRC)	

#Dissobjt.so:
#	$(CC) -shared -o Dissobjt.so Dissobjt.o
#     clean:
#	rm -f core Dissobjt.o
