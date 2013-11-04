/*
 * a simple serial scan.
 */

#include <stdio.h>
#include <stdlib.h>

void scan(float* out, float* in, int size);
void startClock(char*);
void stopClock(char*);
void printClock(char*);

int main(int argc, char** argv) {

	if (argc < 2) {
		printf("Usage: %s size-of-array\n",argv[0]);
		exit(1);
	}
	int size = atoi(argv[1]);
	printf("size = %d\n",size);

	/* host based data */

	float *h_in;
	float *h_out;

	h_in = (float*) malloc(size*sizeof(float));
	h_out =(float*) malloc(size*sizeof(float));

	for (int i = 0; i < size; i++) {
		h_in[i] = 1.0;
	}

	startClock("compute");
	scan(h_out,h_in,size);
	stopClock("compute");

	for (int i = 0; i < size; i++) {
		printf("%d %f -> %f\n",i,h_in[i],h_out[i]);
	}

	free(h_in);
	free(h_out);

	printClock("compute");
}

void scan(float* out, float* in, int size) {
	out[0] = 0;
	for (int i = 1; i < size; i++) {
		out[i] = out[i-1] + in[i-1];
	}
}

