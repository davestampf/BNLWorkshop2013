/*
 * A very simple cuda implementation of reduce. Uses an array of 1024x1024
 * items which are summed into a 1024 array and then summed into a value.
 */

#include <stdio.h>
#include <stdlib.h>

/*
 * this kernel will sum all of the data from in into out - at
 * least as far as the block will carry you
 */

__global__ void reduce(float* out, float* in, int size); 

void startClock(char*);
void stopClock(char*);
void printClock(char*);

int main(int argc, char** argv) {

        int size = 1024*1024;
        printf("size = %d\n",size);
	
	void *d_in;	// device data
	void *d_mid;	// device data - middle results
	void *d_out;	// device data - the answer

	float *h_in;	// host data
	float h_out;

	int numBlocks = 1024;

	cudaMalloc(&d_in,size*sizeof(float));
	cudaMalloc(&d_mid,numBlocks*sizeof(float));
	cudaMalloc(&d_out,sizeof(float));

	h_in = (float*) malloc(size*sizeof(float));

	for (int i = 0; i < size; i++) {
		h_in[i] = 1;
	}

	startClock("copy data to device");	
	cudaMemcpy(d_in,h_in,size*sizeof(float),cudaMemcpyHostToDevice);
	stopClock("copy data to device");	

	startClock("compute");
	
	// use max threads/block and the required # of blocks AND
	// ask for some shared memory

	reduce<<<1024,1024,1024>>>((float*) d_mid,(float*) d_in,size);
	reduce<<<1,1024,1024>>>((float*)d_out,(float*)d_mid,1024);
	cudaThreadSynchronize();

	stopClock("compute");
	
	startClock("copy data to host");
	h_out = -17;
	cudaMemcpy(&h_out,d_out,sizeof(float),cudaMemcpyDeviceToHost);
	stopClock("copy data to host");

	printf("The total is %f\n",h_out);
	free(h_in);
	cudaFree(d_in);
	cudaFree(d_out);

	printClock("copy data to device");
	printClock("compute");
	printClock("copy data to host");
}

