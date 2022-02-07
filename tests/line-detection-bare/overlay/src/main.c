#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include "routinesCPU.h"
#include "img0.h"

/* Time */
#include <sys/time.h>
#include <sys/resource.h>        


int main(int argc, char *argv[])
{
	const int width = 1920, height = 1080;
	uint8_t * imtmp;
	uint8_t * im;

	float sin_table[180], cos_table[180];
	int nlines=0; 
	int x1[10], x2[10], y1[10], y2[10];
	int l;
	double t0, t1;


	/* Read images */
	imtmp = getImtmp();
	im    = getImBW();

	init_cos_sin_table(sin_table, cos_table, 180);	

	// Create temporal buffers 
	uint8_t *imEdge = (uint8_t *)malloc(sizeof(uint8_t) * width * height);

	//Create the accumulators
	const float hough_h = ((sqrt(2.0) * (float)(height>width?height:width)) / 2.0);
	const int accu_height = hough_h * 2.0; // -rho -> +rho
	const int accu_width  = 180;
	uint32_t accum[accu_width*accu_height];

	//Execute on CPU
	line_asist_CPU(im,
		imEdge,
		sin_table, cos_table,
		accum, accu_height, accu_width,
		x1, y1, x2, y2, &nlines);

	for (int l=0; l<nlines; l++)
		printf("(x1,y1)=(%d,%d) (x2,y2)=(%d,%d)\n", x1[l], y1[l], x2[l], y2[l]);
	
	return 0;
}

