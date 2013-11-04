// 2 dimensional powder diffraction OPTIMIZED WITH LINEAR INDEXING

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <complex.h>
#include <cuComplex.h>

// #define COMPARE_SERIAL
#define DO_OUTPUT

/*__global__ void get1(double *d) {
	*d = 1;
}
*/
__global__ void cudaComp(double* Qx,double* Qy,double* Qz,double* U_L,double* V_L,double* W_L,double* Int_2D,int pixel_n,int M) {

	// compute the value for o and p based upon the blockIdx, blockIdy, threadIdx, threadIdy
	int p = threadIdx.x + blockIdx.x*blockDim.x;
	int o = threadIdx.y + blockIdx.y*blockDim.y;
	if (o >= pixel_n || p >= pixel_n) return;
	double qx, qy, qz;
	qx = Qx[o*pixel_n + p];
	qy = Qy[o*pixel_n + p];
	qz = Qz[o*pixel_n + p];

	cuDoubleComplex Amp = make_cuDoubleComplex(0,0);
	double PI = 4.*atan(1.0);
	for (int a = 0; a < M; a++) {
		double exp = 2*PI*(qx*U_L[a] + qy*V_L[a] + qz*W_L[a]);
		double cosValue, sinValue;
		sincos(exp,&sinValue,&cosValue); 
		Amp = cuCadd(Amp,make_cuDoubleComplex(cosValue,sinValue));
	}
	Int_2D[o*pixel_n+p] = Int_2D[o*pixel_n+p] + cuCreal(cuCmul(Amp, cuConj(Amp)));
}

int main()
{
    FILE *fp;  // pointer to the file to write the results
    FILE *coords;   //pointer to the file that contains atomic coordinates
    FILE *angles;
    double lambda = 1.0/10000000000;   // wavelength of illumination
    double PI = 4.*atan(1.0);
    time_t start,stop;
    double diff;

    int i;
    int k;
    int a;
    int b;

	int M = 736;        //no of atoms

    int no_rot = 1000;                  // how many rotations of the crystal will we consider to find the average intensity

    double *x,*y,*z,*alpha_x,*alpha_z,*alpha_y,*detx,*dety;
    x = (double*)calloc(M, sizeof(double));    // x coordinates of the atoms
    y = (double*)calloc(M, sizeof(double));
    z = (double*)calloc(M, sizeof(double));

    double D = 0.2/3;                                    // sample to detector distance

    int pixel_n = 360;                                  // at how many points we will compute the diffracted intensity

    alpha_x = (double*)calloc(no_rot,sizeof(double));  // sample rotation angle around local x axis
    alpha_z = (double*)calloc(no_rot,sizeof(double));
    alpha_y = (double*)calloc(no_rot,sizeof(double));

    double k_i[3];                                     // incoming wave vector
    detx = (double*)calloc(pixel_n,sizeof(double));
    dety = (double*)calloc(pixel_n,sizeof(double));
    int o,p;

 double *Qx,*Qy,*Qz;
 
    Qx = (double*)calloc(pixel_n*pixel_n, sizeof(double));
    Qy = (double*)calloc(pixel_n*pixel_n, sizeof(double));
    Qz = (double*)calloc(pixel_n*pixel_n, sizeof(double));

    
    double SL[9];
    double *U_L,*V_L,*W_L;
    U_L = (double*)calloc(M,sizeof(double));
    V_L = (double*)calloc(M,sizeof(double));
    W_L = (double*)calloc(M,sizeof(double));

    double *Int_2D;
    Int_2D = (double*)calloc(pixel_n*pixel_n, sizeof(double)); 
    double *cuda_Int_2D;
    cuda_Int_2D = (double*)calloc(pixel_n*pixel_n, sizeof(double)); 
    
    double complex Amp = 0;

    time(&start);

   angles = fopen("100randomAngles.txt","r+");
   if(angles == NULL)
{
    printf("Angle File can not be opened\n");
    getchar();
    exit(0);
}
    for (i=0; i<no_rot; i++){
    fscanf(angles,"%lf%lf%lf", &alpha_x[i], &alpha_y[i], &alpha_z[i]);
}
   fclose(angles);

  coords = fopen("736fccSph.txt","r+");

  if(coords == NULL)
 {
    printf("File can not be opened\n");
    getchar();
    exit(0);
 }

  for (i=0; i<M; i++){
    fscanf(coords,"%lf%lf%lf", &x[i], &y[i], &z[i]);
 }
fclose(coords);

    k_i[0]=1/lambda*0;
    k_i[1]=1/lambda*0;
    k_i[2]=1/lambda*1;

        for (o=0; o<pixel_n; o++)
    {

        detx[o] = -0.2 + o*(0.2-(-0.2))/(pixel_n-1);
        dety[o] = -0.2 + o*(0.2-(-0.2))/(pixel_n-1);

    }

 // LINEAR ARRAY ASSIGNING VALUES FOR 3 COMPONENTS OF Q VECTOR

     for (o=0; o<pixel_n; o++)
     {
         for(p=0; p<pixel_n; p++)
         {

             Qx[o*pixel_n+p] = detx[p]/lambda/sqrt(detx[p]*detx[p] + dety[o]*dety[o] + D*D) - k_i[0];
             Qy[o*pixel_n+p] = dety[o]/lambda/sqrt(detx[p]*detx[p] + dety[o]*dety[o] + D*D) - k_i[1];
             Qz[o*pixel_n+p] = D/lambda/sqrt(detx[p]*detx[p] + dety[o]*dety[o] + D*D) - k_i[2];
                  
          }
      }

	double *d_Qx;
	double *d_Qy;
	double *d_Qz;
	double *d_U_L, *d_V_L, *d_W_L;
	double *d_Int_2D;
	cudaError_t err;

	err = cudaMalloc(&d_Qx,pixel_n*pixel_n*sizeof(double));
	if (err != cudaSuccess) {
		fprintf(stderr,"Fail to allocate d_Qx on device %s\n",cudaGetErrorString(err));
		exit(1);
	}
	err = cudaMalloc(&d_Qy,pixel_n*pixel_n*sizeof(double));
	if (err != cudaSuccess) {
		fprintf(stderr,"Fail to allocate d_Qy on device %s\n",cudaGetErrorString(err));
		exit(1);
	}
	err = cudaMalloc(&d_Qz,pixel_n*pixel_n*sizeof(double));
	if (err != cudaSuccess) {
		fprintf(stderr,"Fail to allocate d_Qz on device %s\n",cudaGetErrorString(err));
		exit(1);
	}

	err = cudaMalloc(&d_U_L,M*sizeof(double));
	if (err != cudaSuccess) {
		fprintf(stderr,"Fail to allocate d_U_L on device %s\n",cudaGetErrorString(err));
		exit(1);
	}
	err = cudaMalloc(&d_V_L,M*sizeof(double));
	if (err != cudaSuccess) {
		fprintf(stderr,"Fail to allocate d_V_L on device %s\n",cudaGetErrorString(err));
		exit(1);
	}
	err = cudaMalloc(&d_W_L,M*sizeof(double));
	if (err != cudaSuccess) {
		fprintf(stderr,"Fail to allocate d_W_L on device %s\n",cudaGetErrorString(err));
		exit(1);
	}

	err = cudaMalloc(&d_Int_2D,pixel_n*pixel_n*sizeof(double));
	if (err != cudaSuccess) {
		fprintf(stderr,"Fail to allocate d_Int_2D on device %s\n",cudaGetErrorString(err));
		exit(1);
	}

	//printf("copying static data to device\n");
	err = cudaMemcpy(d_Qx,Qx,pixel_n*pixel_n*sizeof(double),cudaMemcpyHostToDevice);
	if (err != cudaSuccess) {
		fprintf(stderr,"Fail to copy to d_Qx on device %s\n",cudaGetErrorString(err));
		exit(1);
	}
	err = cudaMemcpy(d_Qy,Qy,pixel_n*pixel_n*sizeof(double),cudaMemcpyHostToDevice);
	if (err != cudaSuccess) {
		fprintf(stderr,"Fail to copy d_Qy on device %s\n",cudaGetErrorString(err));
		exit(1);
	}
	err = cudaMemcpy(d_Qz,Qz,pixel_n*pixel_n*sizeof(double),cudaMemcpyHostToDevice);
	if (err != cudaSuccess) {
		fprintf(stderr,"Fail to copy d_Qz on device %s\n",cudaGetErrorString(err));
		exit(1);
	}
	err = cudaMemcpy(d_Int_2D,Int_2D,pixel_n*pixel_n*sizeof(double),cudaMemcpyHostToDevice);
	if (err != cudaSuccess) {
		fprintf(stderr,"Fail to copy d_Int_2D on device %s\n",cudaGetErrorString(err));
		exit(1);
	}


   for (k=0; k<no_rot; k++)
    {

      SL[0] = cos(alpha_z[k])*cos(alpha_y[k]);
      SL[1] = sin(alpha_z[k])*cos(alpha_y[k]);
      SL[2] = sin(alpha_y[k]);
      SL[3] = -cos(alpha_x[k])*sin(alpha_z[k])-sin(alpha_x[k])*sin(alpha_y[k])*cos(alpha_z[k]);
      SL[4] = cos(alpha_z[k])*cos(alpha_x[k])-sin(alpha_x[k])*sin(alpha_y[k])*sin(alpha_z[k]);
      SL[5] = sin(alpha_x[k])*cos(alpha_y[k]);
      SL[6] = sin(alpha_x[k])*sin(alpha_z[k])-sin(alpha_y[k])*cos(alpha_x[k])*cos(alpha_z[k]);
      SL[7] = -sin(alpha_x[k])*cos(alpha_z[k])-sin(alpha_y[k])*cos(alpha_x[k])*sin(alpha_z[k]);
      SL[8] = cos(alpha_x[k])*cos(alpha_y[k]);


      for (a=0; a<M; a++)
      {
          U_L[a] = SL[0]*x[a] + SL[1]*y[a] + SL[2]*z[a];
          V_L[a] = SL[3]*x[a] + SL[4]*y[a] + SL[5]*z[a];
          W_L[a] = SL[6]*x[a] + SL[7]*y[a] + SL[8]*z[a];

      }
/*
 * This is the main time sink. pixel_n is on the order of 360, but M is equal to the number of atoms - probably in the thousands.
 * First, copy over the data to the device
 */

	//printf("copying dynamic data to device\n");
	err = cudaMemcpy(d_U_L,U_L,M*sizeof(double),cudaMemcpyHostToDevice);
	if (err != cudaSuccess) {
		fprintf(stderr,"Fail to copy to d_U_L on device %s\n",cudaGetErrorString(err));
		exit(1);
	}
	err = cudaMemcpy(d_V_L,V_L,M*sizeof(double),cudaMemcpyHostToDevice);
	if (err != cudaSuccess) {
		fprintf(stderr,"Fail to copy to d_V_L on device %s\n",cudaGetErrorString(err));
		exit(1);
	}
	err = cudaMemcpy(d_W_L,W_L,M*sizeof(double),cudaMemcpyHostToDevice);
	if (err != cudaSuccess) {
		fprintf(stderr,"Fail to copy to d_W_L on device %s\n",cudaGetErrorString(err));
		exit(1);
	}
	
	int numThreadsX = 32;
	int numThreadsY = 32;
	dim3 threadBlocks = dim3(numThreadsX,numThreadsY,1);
	dim3 grid = dim3((pixel_n + numThreadsX-1)/numThreadsX,(pixel_n + numThreadsY - 1)/numThreadsY,1);
	//printf("Launching Kernel\n");
	//printf("grid = %d %d %d\n",grid.x,grid.y,grid.z);
	//printf("threadBlocks = %d %d %d\n",threadBlocks.x,threadBlocks.y,threadBlocks.z);

	cudaComp<<<grid,threadBlocks>>>((double*)d_Qx,(double*)d_Qy,(double*)d_Qz,(double*)d_U_L,(double*)d_V_L,(double*)d_W_L,(double*)d_Int_2D,pixel_n,M);

	err = cudaGetLastError();
	if (err != cudaSuccess) {
		printf("launch failure: %s\n",cudaGetErrorString(err));
		exit(1);
	}

#ifdef COMPARE_SERIAL
	printf("Serial Computation\n");
     for (o=0; o<pixel_n; o++)
      {
          for (p=0; p<pixel_n; p++)
          {

              for(a=0; a<M; a++)
              {
                         
                             Amp = Amp + cexp(2*PI*I* (Qx[o*pixel_n+p]*U_L[a] + Qy[o*pixel_n+p]*V_L[a] + Qz[o*pixel_n+p]*W_L[a]));
                         
              }
                          
                         Int_2D[o*pixel_n+p] = Int_2D[o*pixel_n+p] + creal((Amp* conj(Amp)));
                          Amp = 0;
          }
     }
#endif
	
}


	/* copy back the results */

	//printf("Copy back\n");
	err = cudaMemcpy(cuda_Int_2D,d_Int_2D,pixel_n*pixel_n*sizeof(double),cudaMemcpyDeviceToHost);
	if (err != cudaSuccess) {
		fprintf(stderr,"Fail to copy c_Int_2D back to cuda_Int_2D: %s\n",cudaGetErrorString(err));
		exit(1);
	}

#ifdef COMPARE_SERIAL
	printf("Compare answers\n");
	double maxpct = 0.0;
	double maxabs = 0.0;
	for (o=0; o<pixel_n;o++) {
		for (p=0; p<pixel_n;p++) {
			double diff = Int_2D[o*pixel_n+p] - cuda_Int_2D[o*pixel_n+p];
			double pct;
			if (abs(Int_2D[o*pixel_n+p]) < 1.0e-6) {
				pct = abs(cuda_Int_2D[o*pixel_n+p]);
			} else {
				pct = abs(diff/Int_2D[o*pixel_n + p]);
			}
			if (pct > .0001) printf("%d %d %e %e %e %e\n",o,p,cuda_Int_2D[o*pixel_n+p],Int_2D[o*pixel_n+p],diff,pct);
			if (pct > maxpct) maxpct = pct;
			if (abs(diff) > maxabs) maxabs = abs(diff);
		}
	}
	printf("Max pct difference = %e\n",100*maxpct);
	printf("Max abs difference = %e\n",maxabs);
#endif

#ifdef DO_OUTPUT
fp = fopen("736atom10rots_june19_2.txt","w");

if(fp==NULL)
{
   printf("I couldn't open results.txt for writing.\n");
   exit(0);
}

for (a=0; a<pixel_n; a++)
{

    for (b=0; b<pixel_n; b++)
    {
          
      fprintf(fp, "%f\t", cuda_Int_2D[a*pixel_n+b]/100);
      //fprintf(fp, "%d %d %f\t %f\n", a, b, Int_2D[a*pixel_n+b]/10, cuda_Int_2D[a*pixel_n+b]/10);
    }
    //printf("\n");
    fprintf(fp,"\n");

}


fclose(fp);
#endif

time(&stop);
diff = difftime(stop,start);
printf("Intensity is averaged over # of rotations and Time elapsed is %f seconds\n", diff);

free (detx);
free(dety);
free (alpha_x);
free(alpha_y);
free(alpha_z);


    return 0;

}
