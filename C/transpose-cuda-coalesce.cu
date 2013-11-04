/*
 * transpose an array - using a cuda device, use coalescing
 */

#include <stdio.h>
#include <stdlib.h>

__global__ void cudaTransposeCoalesce(float* out, float *in, int size);
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

	int tileSize = 8;
	int tempMem = tileSize*tileSize*sizeof(float);
	dim3 blocks(128,128,1);
	dim3 threads(8,8);

	cudaTransposeCoalesce<<<blocks,threads,tempMem>>>((float*)d_out,(float*)d_in,DIM);
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

__global__ void cudaTransposeCoalesce(float* out, float* in, int size) {

	__shared__ float shared[1024];


	// starting points inside the input data

	int iStart = blockDim.x*blockIdx.x;
	int jStart = blockDim.y*blockIdx.y;

	// let adjacent threads pick up adjacent items from input
	// transpose on the fly

	float data = in[iStart + threadIdx.x + (jStart + threadIdx.y)*size];
	shared[threadIdx.y + threadIdx.x*blockDim.x] = data;

	__syncthreads();

	// ok now put them back to out, but travel with the grain!

	int temp = jStart;
	jStart = iStart;
	iStart = temp;

	out[iStart + threadIdx.x + (jStart + threadIdx.y)*size] = shared[threadIdx.x + threadIdx.y*blockDim.y];
}

