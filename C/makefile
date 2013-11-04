.SUFFIXES : .cu .f90

.cu.o :
	nvcc $(NVCCFLAGS) -c $^

NVCCFLAGS = -m64 -arch=sm_35


all: cuda-devices map-cuda-small map-noncuda map-cuda-large \
	scan-cuda-small scan-noncuda null reduce-noncuda reduce-cuda \
	transpose-noncuda transpose-cuda-serial transpose-cuda-row \
	transpose-cuda-full transpose-cuda-coalesce 

cuda-devices : cuda-devices.o
	nvcc $(NVCCFLAGS) -o cuda-devices $^

map-noncuda : map-noncuda.o clock.o
	nvcc $(NVCCFLAGS) -o map-noncuda $^

map-cuda-small : map-cuda-kernel-small.o map-cuda-small.o clock.o
	nvcc $(NVCCFLAGS) -o map-cuda-small $^

map-cuda-large : map-cuda-kernel-large.o map-cuda-large.o clock.o
	nvcc $(NVCCFLAGS) -o map-cuda-large $^

null : null-kernel.o map-cuda-large.o clock.o
	nvcc $(NVCCFLAGS) -o null $^

scan-cuda-small: scan-cuda-kernel-small.o scan-cuda-test-small.o clock.o
	nvcc $(NVCCFLAGS) -o scan-cuda-small $^

scan-noncuda: scan-noncuda.o clock.o
	nvcc $(NVCCFLAGS) -o scan-noncuda $^

reduce-noncuda: reduce-noncuda.o clock.o
	nvcc $(NVCCFLAGS) -o reduce-noncuda $^

reduce-cuda: reduce-cuda-large.o reduce-cuda-kernel.o clock.o
	nvcc $(NVCCFLAGS) -o reduce-cuda $^

transpose-noncuda: transpose-noncuda.o clock.o
	nvcc $(NVCCFLAGS) -o transpose-noncuda $^

transpose-cuda-serial: transpose-cuda-serial.o clock.o
	nvcc $(NVCCFLAGS) -o transpose-cuda-serial $^

transpose-cuda-row: transpose-cuda-row.o clock.o
	nvcc $(NVCCFLAGS) -o transpose-cuda-row $^

transpose-cuda-full: transpose-cuda-full.o clock.o
	nvcc $(NVCCFLAGS) -o transpose-cuda-full $^

transpose-cuda-coalesce: transpose-cuda-coalesce.o clock.o
	nvcc $(NVCCFLAGS) -o transpose-cuda-coalesce $^

clean:
	- rm core
	- rm *.o
	- rm cuda-devices
	- rm map-noncuda
	- rm map-cuda-small
	- rm map-cuda-large
	- rm scan-cuda-small
	- rm scan-noncuda
	- rm null
	- rm reduce-noncuda
	- rm reduce-cuda
	- rm transpose-noncuda
	- rm transpose-cuda-serial
	- rm transpose-cuda-row
	- rm transpose-cuda-full
	- rm transpose-cuda-coalesce
