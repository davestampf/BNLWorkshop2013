/*
 * transpose an array
 */

#include <stdio.h>
#include <stdlib.h>

void nonCudaTranspose(float* out, float *in, int size);
void startClock(char*);
void stopClock(char*);
void printClock(char*);

#define DIM 1024 

int main(int argc, char** argv) {

	float *h_in;
	float *h_out;

	h_in = (float*) malloc(DIM*DIM*sizeof(float));
	h_out =(float*) malloc(DIM*DIM*sizeof(float));

	int value = 1;
	for (int i = 0; i < DIM; i++) {
		for (int j = 0; j < DIM; j++) {
			h_in[i + j*DIM] = value++;
		}
	}

/*	for (int i = 0; i < DIM; i++) {
		for (int j = 0; j < DIM; j++) {
			printf("%10.2f ",h_in[i+j*DIM]);
		}
		printf("\n");
	}
*/
	startClock("compute");
	nonCudaTranspose(h_out,h_in,DIM);
	stopClock("compute");
		
/*	for (int i = 0; i < DIM; i++) {
		for (int j = 0; j < DIM; j++) {
			printf("%10.2f ",h_out[i+j*DIM]);
		}
		printf("\n");
	}
*/
	free(h_in);
	free(h_out);

	printClock("compute");
}

void nonCudaTranspose(float* out, float* in, int size) {
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			out[j + i*size] = in[i + j*size];
		}
	}
}	

