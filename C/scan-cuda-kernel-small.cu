/*
 * a simple scan program. compute the partial sums of
 * the elements of the input array. This version uses
 * only one block
 */

__global__ void cudaScan(float* d_out, float* d_in, int n) {

	// shared array allocated by the launch of the kernel

	extern __shared__ float temp[];
	
	int threadId = threadIdx.x;
	if (threadId >=  n) return;

	int fromBuffer = 1;
	int toBuffer = 0;

	// make a local copy of the data

	temp[threadId] = d_in[threadId];
	__syncthreads();


	int maxOffset =(int)ceil(log2(1.0f*n));
	maxOffset = pow(2.0f,1.0f*maxOffset);
	for (int offset = 1; offset < maxOffset; offset *= 2) {

		fromBuffer = 1-fromBuffer;
		toBuffer = 1-toBuffer;
		if (threadId >= offset) {
			temp[toBuffer*n + threadId] = temp[fromBuffer*n + threadId - offset] +
				temp[fromBuffer*n + threadId]; 
		} else {
			temp[toBuffer*n+ threadId] = temp[fromBuffer*n + threadId];
		}
		__syncthreads();
	}

	d_out[threadId] = temp[toBuffer*n + threadId];
}

