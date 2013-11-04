/*
 * reduce an array of 1's by the sum
 */

#include <stdio.h>
#include <stdlib.h>

void nonCudaReduce(float* out, float *in, int size);
void startClock(char*);
void stopClock(char*);
void printClock(char*);

int main(int argc, char** argv) {

	if (argc < 2) {
		printf("Usage: %s #-of-floats\n",argv[0]);
		exit(1);
	}
	int size = atoi(argv[1]);
	printf("size = %d\n",size);

	float *h_in;
	float h_out;

	h_in = (float*) malloc(size*sizeof(float));

	for (int i = 0; i < size; i++) {
		h_in[i] = 1;
	}

	startClock("compute");
	nonCudaReduce(&h_out,h_in,size);
	stopClock("compute");
	
	printf("The sum is %f\n",h_out);	

	free(h_in);

	printClock("compute");
}

void nonCudaReduce(float* out, float* in, int size) {
	*out = 0.0;
	for (int i = 0; i < size; i++) {
		*out += in[i];
	}
}	

