/* Naive off reader for point cloud from .off file */
/* Only works for TP5 */
/* Point cloud is stored in points */
/* The size of points is N */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "imageFormationUtils.h"


struct point3d *readOff(char *fileName, int *N){
	struct point3d *points;

	FILE *ifp = fopen(fileName, "r");
	if (ifp == NULL) {
  		fprintf(stderr, "Can't open input file in.list!\n");
  		exit(1);
	}

	char fooChar[16];
	int N_v, fooInt;
	fscanf(ifp, "%s", fooChar);
	//printf("%s \n", fooChar);
	fscanf(ifp, "%d %d %d", &N_v, &fooInt, &fooInt);
	*N = N_v;
	points = (struct point3d *) malloc(N_v * sizeof(struct point3d));
	//printf("%d %d %d\n", N_v, fooInt, fooInt);
	int i = 0;
	for(i = 0; i < N_v; i++){
		fscanf(ifp, "%f %f %f %d %d %d %d", &(points[i].x), &(points[i].y), &(points[i].z), &(points[i].r), &(points[i].g), &(points[i].b), &fooInt);
		//printf("%f %f %f %d %d %d\n", points[i].x, points[i].y, points[i].z, points[i].r, points[i].g, points[i].b);
	}

	return points;
}

void centerThePCL(struct point3d *points, int N){
  int i = 0;
  float xMin = 100000, xMax = -1000000, yMin = 1000000, yMax = -10000000, zMin = 100000000, zMax = -10000000;
  for(i = 0; i < N; i++){
    if(points[i].x > xMax)
      xMax = points[i].x;
    else if (points[i].x < xMin)
      xMin = points[i].x;

    if(points[i].y > yMax)
      yMax = points[i].y;
    else if (points[i].y < yMin)
      yMin = points[i].y;

    if(points[i].z > zMax)
      zMax = points[i].z;
    else if (points[i].z < zMin)
      zMin = points[i].z;
  }

  float xMid = (xMin+xMax)/2, yMid = (yMin+yMax)/2, zMid = (zMin+zMax)/2;
  for(i = 0; i < N; i++){
    points[i].x = points[i].x-xMid;
    points[i].y = points[i].y-yMid;
    points[i].z = points[i].z-zMid;
  }
}


void Rx(double alpha, float * R){

  int i = 0, j = 0;
  for(i = 0; i < 4; i++)
    for(j = 0; j < 4; j++)
      R[i*4+j] = 0;
    // 1      0     0    0
    //   0     cosA -sinA  0
    //   0     sinA  cosA  0
    //   0      0     0    1
      
  R[0] = 1;

  R[5] = cos(alpha);
  R[6] = -sin(alpha);

  R[9] = sin(alpha);
  R[10] = cos(alpha);

  R[15] = 1;

  return;
}

void Ry(float alpha, float *R){

  int i = 0, j = 0;
  for(i = 0; i < 4; i++)
    for(j = 0; j < 4; j++)
      R[i*4+j] = 0;
  /* cosA    0    sinA   0
      0      1     0     0
    -sinA    0    cosA   0
      0      0     0     1
      */
  R[0] = cos(alpha);
  R[2] = sin(alpha);

  R[5] = 1;

  R[8] = -sin(alpha);
  R[10] = cos(alpha);

  R[15] = 1;

  return;
}

void Rz(float alpha, float *R){

  int i = 0, j = 0;
  for(i = 0; i < 4; i++)
    for(j = 0; j < 4; j++)
      R[i*4+j] = 0;
  /* cosA  -sinA   0   0
     sinA   cosA   0   0
      0      0     1   0
      0      0     0   1
      */
  R[0] = cos(alpha);
  R[1] = -sin(alpha);

  R[4] = sin(alpha);
  R[5] = cos(alpha);

  R[10] = 1;
  R[15] = 1;

  return;
}

void matMul(float *A, float *B, float *result){
  float tmpA[16], tmpB[16];
  int i = 0, j = 0;
  for(i = 0; i < 4; i++){
    for(j = 0; j < 4; j++){
      tmpA[i*4+j] = A[i*4+j];
      tmpB[i*4+j] = B[i*4+j];
    }
  }
  for(i = 0; i < 4; i++){
    for(j = 0; j < 4; j++){
      result[i*4+j] = tmpA[i*4+0]*tmpB[0*4+j]+tmpA[i*4+1]*tmpB[1*4+j]+tmpA[i*4+2]*tmpB[2*4+j]+tmpA[i*4+3]*tmpB[3*4+j];
    }
  }
  return;
}

void computeTrans(float gama, float beta, float alpha, float T_x, float T_y, float T_z, float *result){
  float R_x[16], R_y[16], R_z[16];
  Rx(gama, R_x);
  Ry(beta, R_y);
  Rz(alpha, R_z);

  matMul(R_y, R_z, result);
  matMul(R_x, result, result);
  result[3] = T_x;
  result[7] = T_y;
  result[11] = T_z;
  return;
}

struct point2d *pinhole_projection(struct point3d *points, int N, int f)
{
    struct point2d* result = malloc(N*sizeof (struct point2d));

    for(int i=0;i<N;i++)
    {
        result[i].r = points[i].r;
        result[i].g = points[i].g;
        result[i].b = points[i].b;

        result[i].x = points[i].x/(1+points[i].z/f);
        result[i].y = points[i].y/(1+points[i].z/f);
//        printf("%f %f\n", result[i].x, result[i].y);

    }
    return result;
}

image_structure_rgb_t *pixel_coordinate_projection(struct point2d *points, int N, int height, int width, int u0, int v0, int alpha_u, int alpha_v)
{
    image_structure_rgb_t *new_image = malloc(sizeof(image_structure_rgb_t));
    unsigned char *r = (unsigned char*) calloc(width* height, sizeof (unsigned char));
    unsigned char *g = (unsigned char*) calloc(width* height, sizeof (unsigned char));
    unsigned char *b = (unsigned char*) calloc(width* height, sizeof (unsigned char));
    int u, v;
    for(int i=0;i<N;i++)
    {
        u = (int)(points[i].x*alpha_u + u0);
        v = (int)(points[i].y*alpha_v + v0);
//        printf("%d %d", u, v);

        if(!(u<0 || u>=height || v<0 || v>=width))
        {
            r[u*width+v] = points[i].r;
            g[u*width+v] = points[i].g;
            b[u*width+v] = points[i].b;

        }

    }


    new_image->maxval = 255;
    new_image->rows = height;
    new_image->cols = width;
    new_image->red = r;
    new_image->green = g;
    new_image->blue = b;

    return new_image;
}

void write_image_to_file_rgb(image_structure_rgb_t *img, FILE * fp)
{
//    FILE *fp = fopen(file_name, "w");
    fprintf(fp, "P6\n");

    fprintf(fp, "%d %d \n", img->cols, img->rows);
    fprintf(fp, "%d\n",img->maxval);

    for(int i=0; i < img->rows; i++)
        for(int j=0; j < img->cols ; j++)
        {
            fprintf(fp,"%c",img->red[i * img->cols + j]);
            fprintf(fp,"%c",img->green[i * img->cols + j]);
            fprintf(fp,"%c",img->blue[i * img->cols + j]);


        }
    fclose(fp);
}

struct point3d *rigid_transformation(struct point3d *points, int N, float alpha, float beta, float gamma)
{
    float *rx = calloc(16,sizeof(float));
    float *ry = calloc(16,sizeof(float));
    float *rz = calloc(16,sizeof(float));

    Rx(gamma, rx);Ry(beta, ry);Rz(alpha, rz);

    float *result_p = calloc(16,sizeof(float));
    float *result = calloc(16,sizeof(float));

    matMul(rz,ry,result_p);
    matMul(result_p,rx,result);

    struct point3d* final_result = malloc(N*sizeof (struct point3d));
    for(int i=0;i<N;i++)
    {
        final_result[i].r = points[i].r;
        final_result[i].g = points[i].g;
        final_result[i].b = points[i].b;

        final_result[i].x = result[0]*points[i].x + result[1]*points[i].y + result[2]*points[i].z + result[3];
        final_result[i].y = result[4]*points[i].x + result[5]*points[i].y + result[6]*points[i].z + result[7];
        final_result[i].z = result[8]*points[i].x + result[9]*points[i].y + result[10]*points[i].z + result[11];

    }

    return final_result;



}


struct point3d *filter_points(struct point3d *points, int *N)
{
    struct point3d* final_result = malloc(*N*sizeof (struct point3d));
    printf("%d\n", *N);

    int idx = 0, flag = 0;
    for(int i=0;i<*N;i++)
    {
        flag = 0;
        for(int j=0;j<idx;j++)
        {
//            printf("bonjour from above\n");
            if(final_result[j].x == points[i].x && final_result[j].y == points[i].y)
            {
//                printf("bonjourr\n");
                if(final_result[j].z > points[i].z) final_result[j].z = points[i].z;
                flag=1;
                break;
            }
        }

        if(flag == 0)
        {
            final_result[idx].x = points[i].x;
            final_result[idx].y = points[i].y;
            final_result[idx].z = points[i].z;

            final_result[idx].r = points[i].r;
            final_result[idx].g = points[i].g;
            final_result[idx].b = points[i].b;


            idx++;

        }
    }

    *N = idx;
    printf("%d\n", *N);
    return final_result;
}

struct point2d *another_filter_points(struct point3d *points, int N, int f)
{
    struct point2d* final_result = malloc(N*sizeof (struct point2d));
    struct point3d* aux = malloc(N*sizeof (struct point3d));
    int idx = 0, flag = 0;
    float aux_x, aux_y;
//    for(int i=0;i<*N;i++)
//    {
//        flag = 0;
//        aux_x = points[i].x/(1+points[i].z/f);
//        aux_y = points[i].y/(1+points[i].z/f);
//        for(int j=0;j<idx;j++)
//        {
//
//            if(final_result[j].x == aux_x && final_result[j].y == aux_y)
//            {
//                if(aux[j].z > points[i].z)
//                {
//                    final_result[idx].r = points[i].r;
//                    final_result[idx].g = points[i].g;
//                    final_result[idx].b = points[i].b;
//                }
//                flag=1;
//                break;
//            }
//        }
//
//        if(flag == 0)
//        {
//
//            final_result[idx].x = points[i].x/(1+points[i].z/f);
//            final_result[idx].y = points[i].y/(1+points[i].z/f);
//            aux[idx].z = aux[i].z;
//
//            final_result[idx].r = points[i].r;
//            final_result[idx].g = points[i].g;
//            final_result[idx].b = points[i].b;
//
//
//            idx++;
//
//        }
//
//    }
    for(int i=0;i<N;i++)
    {
        aux[i].r = points[i].r;
        aux[i].g = points[i].g;
        aux[i].b = points[i].b;

        aux[i].x = points[i].x/(1+points[i].z/f);
        aux[i].y = points[i].y/(1+points[i].z/f);
        aux[i].z = points[i].z;
//        printf("%f %f\n", result[i].x, result[i].y);

    }
    int r,g,b,mz=0;
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            if(aux[i].x == aux[j].x && aux[i].y == aux[j].y)
            {
                if(aux[j].z >mz)
                {
                    r = aux[j].r;
                    g = aux[j].g;
                    b = aux[j].b;

                    mz = aux[j].z;
                }

            }
        }
        for(int j=0;j<N;j++)
        {
            if(aux[i].x == aux[j].x && aux[i].y == aux[j].y && aux[i].z < mz)
            {

                   final_result[j].x = aux[j].x;
                    final_result[j].y = aux[j].y;
                    final_result[j].r = 0;
                    final_result[j].g = 0;
                    final_result[j].b = 0;

            }
            else
            {
                final_result[j].x = aux[j].x;
                final_result[j].y = aux[j].y;
                final_result[j].r = aux[j].r;
                final_result[j].g = aux[j].g;
                final_result[j].b = aux[j].b;
            }
        }
        }



    return final_result;
}