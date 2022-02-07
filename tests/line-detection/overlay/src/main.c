#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include "routinesCPU.h"

/* Time */
#include <sys/time.h>
#include <sys/resource.h>        


int main(int argc, char **argv)
{
	const int width = 1920, height = 1080;
	uint8_t imtmp[width*height], im[width*height];

	float sin_table[180], cos_table[180];
	int nlines=0; 
	int x1[10], x2[10], y1[10], y2[10];
	int l;
	double t0, t1;


	/* Read images */
	FILE *fp = fopen("img0.bin", "rb");
	size_t read = fread(imtmp, sizeof(uint8_t), width*height, fp);
	fclose(fp);
	fp = fopen("img0_bw.bin", "rb");
	read = fread(im, sizeof(uint8_t), width*height, fp);
	fclose(fp);

	init_cos_sin_table(sin_table, cos_table, 180);	

	// Create temporal buffers 
	uint8_t *imEdge = (uint8_t *)malloc(sizeof(uint8_t) * width * height);
	int *NR = (int *)malloc(sizeof(int) * width * height);
	float *G = (float *)malloc(sizeof(float) * width * height);
	int *phi = (int *)malloc(sizeof(int) * width * height);
	int *Gx = (int *)malloc(sizeof(int) * width * height);
	int *Gy = (int *)malloc(sizeof(int) * width * height);
	uint8_t *pedge = (uint8_t *)malloc(sizeof(uint8_t) * width * height);

	//Create the accumulators
	float hough_h = ((sqrt(2.0) * (float)(height>width?height:width)) / 2.0);
	int accu_height = hough_h * 2.0; // -rho -> +rho
	int accu_width  = 180;
	uint32_t *accum = (uint32_t*)malloc(accu_width*accu_height*sizeof(uint32_t));

	//Execute on CPU
	line_asist_CPU(im, height, width, 
		imEdge, NR, G, phi, Gx, Gy, pedge,
		sin_table, cos_table,
		accum, accu_height, accu_width,
		x1, y1, x2, y2, &nlines);

	for (int l=0; l<nlines; l++)
		printf("(x1,y1)=(%d,%d) (x2,y2)=(%d,%d)\n", x1[l], y1[l], x2[l], y2[l]);

}

