#include <stdio.h>
//#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include "routinesCPU.h"
     

#define width 1920
#define height 1080
#define accu_height 2715
#define accu_width 180

int main(int argc, char **argv)
{
	int * imtmp;
	int * im;

	float sin_table[180], cos_table[180];
	int nlines=0; 
	int x1[10], x2[10], y1[10], y2[10];
	int l;
	double t0, t1;


	/* Read images */
	imtmp = getImtmp();
	imtmp = getImBW();

	init_cos_sin_table(sin_table, cos_table, 180);	

	// Create temporal buffers 
	uint8_t imEdge[width * height];

	//Create the accumulators
	uint32_t accum[accu_width*accu_height];

	//Execute on CPU
	line_asist_CPU(im,
		imEdge,
		sin_table, cos_table,
		accum, accu_height, accu_width,
		x1, y1, x2, y2, &nlines);

	for (int l=0; l<nlines; l++)
		printf("(x1,y1)=(%d,%d) (x2,y2)=(%d,%d)\n", x1[l], y1[l], x2[l], y2[l]);

}

