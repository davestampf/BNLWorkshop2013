.SUFFIXES : .f90


.f90.o :
	pgf90 -m64 -Mcuda=cc35 -c $^

NFORFLAGS = -m64 -Mcuda=cc35 -Bstatic -Bstatic_pgi 

all: cuda-devices map-noncuda map-cuda-small map-cuda-large

cuda-devices : cuda-devices.o
	pgf90 $(NFORFLAGS) -o cuda-devices $^

map-noncuda : map-noncuda.o
	pgf90 $(NFORFLAGS) -o map-noncuda $^

map-cuda-small : map-cuda-small.o
	pgf90 $(NFORFLAGS) -o map-cuda-small $^

map-cuda-large : map-cuda-large.o
	pgf90 $(NFORFLAGS) -o map-cuda-large $^

clean:
	- rm core
	- rm *.o
	- rm cuda-devices
	- rm map-noncuda
	- rm map-cuda-small
	- rm map-cuda-large
