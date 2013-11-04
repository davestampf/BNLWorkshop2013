/*
 * A very simple cuda implementation of map
 */

#include <stdio.h>
#include <stdlib.h>

__global__ void map(float* out, float* in, int size); 

void startClock(char*);
void stopClock(char*);
void printClock(char*);

int main(int argc, char** argv) {

        if (argc < 2) {
                printf("Usage: %s #-of-floats\n",argv[0]);
                exit(1);
        }
        long size = atoi(argv[1]);
        printf("size = %d\n",size);
	
	cudaDeviceProp props;
	cudaGetDeviceProperties(&props,0);
	if (size > props.maxThreadsPerBlock*(long)props.maxGridSize[0]) {
		fprintf(stderr,"Max size for the large model is %d\n",
			props.maxThreadsPerBlock*props.maxGridSize[0]);
		exit(1);
	}

	void *d_in;	// device data
	void *d_out;
	float *h_in;	// host data
	float *h_out;

	cudaMalloc(&d_in,size*sizeof(float));
	cudaMalloc(&d_out,size*sizeof(float));
	h_in = (float*) malloc(size*sizeof(float));
	h_out =(float*) malloc(size*sizeof(float));

	for (int i = 0; i < size; i++) {
		h_in[i] = i;
	}

	startClock("copy data to device");	
	cudaMemcpy(d_in,h_in,size*sizeof(float),cudaMemcpyHostToDevice);
	stopClock("copy data to device");	

	startClock("compute");
	
	// use max threads/block and the required # of blocks

	int numBlocks = ceil(1.0*size/props.maxThreadsPerBlock);
	map<<<numBlocks,props.maxThreadsPerBlock>>>((float*) d_out,(float*) d_in,size);
	cudaThreadSynchronize();
	stopClock("compute");
	
	startClock("copy data to host");
	cudaMemcpy(h_out,d_out,size*sizeof(float),cudaMemcpyDeviceToHost);
	stopClock("copy data to host");

	for (int i = 0; i < size; i++) {
		printf("%f -> %f\n",h_in[i],h_out[i]);
	}

	free(h_in);
	free(h_out);
	cudaFree(d_in);
	cudaFree(d_out);

	printClock("copy data to device");
	printClock("compute");
	printClock("copy data to host");
}

