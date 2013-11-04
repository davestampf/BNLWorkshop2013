// 2 dimensional powder diffraction OPTIMIZED WITH LINEAR INDEXING

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <complex.h>

int main()
{
    FILE *fp;  // pointer to the file to write the results
    FILE *coords;   //pointer to the file that contains atomic coordinates
    FILE *angles;
    double lambda = 1.0/10000000000;   // wavelength of illumination
    const double PI = 3.14159265;
    time_t start,stop;
    double diff;

    int j;
    int i;
    int k;
    int a;
    int b;
    double complex temp = 0.0;

    int M = 736;        //no of atoms

    int no_rot = 10;                  // how many rotations of the crystal will we consider to find the average intensity

    double *x,*y,*z,*alpha_x,*alpha_z,*alpha_y,*detx,*dety;
    x = (double*)calloc(M, sizeof(double));    // x coordinates of the atoms
    y = (double*)calloc(M, sizeof(double));
    z = (double*)calloc(M, sizeof(double));

    double f = 1.0;                                   // atomic scattering factor--1 for constant
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
    
    double complex Amp = 0;

    time(&start);

   angles = fopen("100randomAngles.txt","r+");
   if(angles == NULL) {
   	printf("Angle File can not be opened\n");
   	getchar();
   	exit(0);
   }
   for (i=0; i<no_rot; i++){
      fscanf(angles,"%lf%lf%lf", &alpha_x[i], &alpha_y[i], &alpha_z[i]);
   }
   fclose(angles);

  coords = fopen("736fccSph.txt","r+");

  if(coords == NULL) {
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

      for (o=0; o<pixel_n; o++)
      {
          for (p=0; p<pixel_n; p++)
          {

              for(a=0; a<M; a++)
              {
                         
                             Amp = Amp + cexp(2*PI*I* (Qx[o*pixel_n+p]*U_L[a] + Qy[o*pixel_n+p]*V_L[a] + Qz[o*pixel_n+p]*W_L[a]));
                         
              }
                          
                         Int_2D[o*pixel_n+p] = Int_2D[o*pixel_n+p] + (Amp* conj(Amp));
                          Amp = 0;
          }
     }

}

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
          
      fprintf(fp, "%f\t", Int_2D[a*pixel_n+b]/10);
    }
    fprintf(fp,"\n");

}


fclose(fp);

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
