/*
 * Map the square function sqeezing all of this data through 1 CPU.
 */

#include <stdio.h>
#include <stdlib.h>

void nonCudaMap(float* out, float *in, int size);
void startClock(char*);
void stopClock(char*);
void printClock(char*);

float square(float x) {
	return x*x;
}

int main(int argc, char** argv) {

	if (argc < 2) {
		printf("Usage: %s #-of-floats\n",argv[0]);
		exit(1);
	}
	int size = atoi(argv[1]);
	printf("size = %d\n",size);

	float *h_in, *h_out;

	h_in = (float*) malloc(size*sizeof(float));
	h_out =(float*) malloc(size*sizeof(float));

	for (int i = 0; i < size; i++) {
		h_in[i] = i;
	}

	startClock("compute");
	nonCudaMap(h_out,h_in,size);
	stopClock("compute");
		
	for (int i = 0; i < size; i++) {
		printf("%f -> %f\n",h_in[i],h_out[i]);
	}

	free(h_in);
	free(h_out);

	printClock("compute");
}

void nonCudaMap(float* out, float* in, int size) {
	for (int i = 0; i < size; i++) {
		out[i] = square(in[i]);
	}
}	

