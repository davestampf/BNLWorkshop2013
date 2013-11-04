/*
 * squaring map kernel that runs in N linear blocks
 */

/*
 * runs on and callable from the device
 */

__device__ float square(float x) {
	return x*x;
}

__global__ void map(float* out, float* in, int size) {
/*	int index = blockDim.x*blockIdx.x + threadIdx.x;
	if (index >= size) return;
	out[index] = square(in[index]);
*/
}

