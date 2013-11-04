/*
 * a simple test of the scan kernel.
 */

#include <stdio.h>
#include <stdlib.h>

__global__ void cudaScan(float* out, float *in, int size);
void startClock(char*);
void stopClock(char*);
void printClock(char*);

int main(int argc, char** argv) {

	if (argc < 2) {
		printf("Usage: %s size-of-array\n",argv[0]);
		exit(1);
	}
	int size = atoi(argv[1]);
	cudaDeviceProp props;
	cudaGetDeviceProperties (&props,0);
	if (size > props.maxThreadsPerBlock) {
		fprintf(stderr,"At most %d elements for small version\n",
			props.maxThreadsPerBlock);
		exit(1);
	}
	printf("size = %d\n",size);

	/* host based data */

	float *h_in;
	float *h_out;

	/* device based data */

	float *d_in;
	float *d_out;

	h_in = (float*) malloc(size*sizeof(float));
	h_out =(float*) malloc(size*sizeof(float));

	cudaMalloc(&d_in,size*sizeof(float));
	cudaMalloc(&d_out,size*sizeof(float));

	for (int i = 0; i < size; i++) {
		h_in[i] = 1.0;
	}

	startClock("copy data to device");
	cudaMemcpy(d_in,h_in,size*sizeof(float),cudaMemcpyHostToDevice);
	stopClock("copy data to device");

	startClock("compute");

	cudaScan<<<1,size,2*size*sizeof(float)>>>(d_out,d_in,size);
	cudaThreadSynchronize();

	stopClock("compute");

	startClock("copy data from device");
	cudaMemcpy(h_out,d_out,size*sizeof(float),cudaMemcpyDeviceToHost);
	stopClock("copy data from device");

	float sum = 0.0f;
	for (int i = 0; i < size; i++) {
		sum += h_in[i];
		printf("%d %f -> %f (%f)\n",i,h_in[i],h_out[i],sum);
	}

	free(h_in);
	free(h_out);

	cudaFree(d_in);
	cudaFree(d_out);

	printClock("copy data to device");
	printClock("compute");
	printClock("copy data from device");
}

