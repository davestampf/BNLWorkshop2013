/*
 * Just how many cuda enabled devices on this machine?
 * Also, what are their properties?
 * 
 * Note - EVERY cuda call returns an error value. While
 * this is vital in real code, it gets in the way of 
 * tutorial code.  I'm showing it here for cudaGetDeviceCount
 * but will omit it for the rest of the tutorial.
 */

#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv) {

	int numberOfDevices;
	cudaError_t err;

	err = cudaGetDeviceCount(&numberOfDevices);
	if (err != cudaSuccess) {
		fprintf(stderr,"fail - cudaGetDeviceCount %d\n",err);
		exit(1);
	}
	printf("Number of cuda devices = %d\n",numberOfDevices);

	/* the cudaDeviceProp struct is fairly large - read about it in the
	   docs. */
	for (int dev = 0; dev < numberOfDevices; dev++) {
		cudaDeviceProp props;
		cudaGetDeviceProperties(&props,dev);
		printf("Device # %d\n",dev);
		printf(" name = %s\n",props.name);
		printf(" version = %d.%d\n",props.major,props.minor);
		printf(" total global memory = %ld\n",props.totalGlobalMem);
		printf(" shared Memory/Block = %ld\n",props.sharedMemPerBlock);
		printf(" registers/block = %d\n",props.regsPerBlock);
		printf(" warp size = %d\n",props.warpSize);
		printf(" Max threads/block = %d\n",props.maxThreadsPerBlock);
		printf(" Max Threads Dim = %d x %d x %d\n",props.maxThreadsDim[0],
			props.maxThreadsDim[1],props.maxThreadsDim[2]);
		printf(" Max Grid Size = %d x %d x %d\n",props.maxGridSize[0],
			props.maxGridSize[1],props.maxGridSize[2]);
		printf(" Multi-processor count = %d\n",props.multiProcessorCount);
		printf(" Max Threads/multiprocessor = %d\n",props.maxThreadsPerMultiProcessor);

	}
	return 0;
}
