/*
 * reduction kernel. Initially, each thread will copy 1 item of data
 * from global to shared memory. Then will will do the binary tree dance. 
 */

__global__ void reduce(float* out, float* in, int size) {
	__shared__ float temp[1024];

	int index = blockDim.x*blockIdx.x + threadIdx.x;
	int myId = threadIdx.x;

	if (index >= size) return;

	// move data to shared memory for speed

	temp[myId] = in[index];
	__syncthreads();

	int stride = blockDim.x/2;
	while (stride >= 1) {
		if (myId < stride) {
			temp[myId] += temp[myId + stride];
		}
		__syncthreads();
		stride = stride/2;
	}	
	out[blockIdx.x] = temp[0];
}

