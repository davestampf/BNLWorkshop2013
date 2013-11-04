/*
 * transpose an array - using a cuda device, but no parallelism
 */

#include <stdio.h>
#include <stdlib.h>

__global__ void cudaTransposeSerial(float* out, float *in, int size);
void startClock(char*);
void stopClock(char*);
void printClock(char*);

#define DIM 1024 

int main(int argc, char** argv) {

	float *h_in;
	float *h_out;

	h_in = (float*) malloc(DIM*DIM*sizeof(float));
	h_out =(float*) malloc(DIM*DIM*sizeof(float));

	void *d_in;
	void *d_out;

	cudaMalloc(&d_in,DIM*DIM*sizeof(float));
	cudaMalloc(&d_out,DIM*DIM*sizeof(float));

	int value = 1;
	for (int i = 0; i < DIM; i++) {
		for (int j = 0; j < DIM; j++) {
			h_in[i + j*DIM] = value++;
		}
	}

	startClock("copy in");
	cudaMemcpy(d_in,h_in,DIM*DIM*sizeof(float),cudaMemcpyHostToDevice);
	stopClock("copy in");

	startClock("compute");
	cudaTransposeSerial<<<1,1>>>((float*)d_out,(float*)d_in,DIM);
	cudaThreadSynchronize();
	stopClock("compute");
		
	startClock("copy out");
	cudaMemcpy(h_out,d_out,DIM*DIM*sizeof(float),cudaMemcpyDeviceToHost);
	stopClock("copy out");

	// sanity check

	for (int i = 0; i < DIM; i++) {
		for (int j = 0; j < DIM; j++) {
			if (h_in[i + j*DIM] != h_out[i*DIM + j]) {
				printf("ERROR");
				exit(1);
			}
		}
	}
	free(h_in);
	free(h_out);
	cudaFree(d_in);
	cudaFree(d_out);

	printClock("copy in");
	printClock("compute");
	printClock("copy out");
}

__global__ void cudaTransposeSerial(float* out, float* in, int size) {
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			out[j + i*size] = in[i + j*size];
		}
	}
}	

