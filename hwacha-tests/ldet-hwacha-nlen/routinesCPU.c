#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "routinesCPU.h"
#include "vec_mul.h"

#define DEG2RAD 0.017453f

#define width 1920
#define height 1080

int NR[width * height];
float G[width * height]; 
int phi[width * height];
int Gx[width * height];
int Gy[width * height];
uint8_t pedge[width * height];

unsigned long read_cycles(void)
{
  unsigned long cycles;
  asm volatile ("rdcycle %0" : "=r" (cycles));
  return cycles;
}

/*unsigned long read_time(void)
{
  unsigned long time;
  asm volatile ("rdtime %0" : "=r" (time));
  return time;
}*/

unsigned long read_instret(void)
{
  unsigned long instret;
  asm volatile ("rdinstret %0" : "=r" (instret));
  return instret;
}

#define VLen_NR 25
#define VLen_G 20
#define nElemsxIt 8 //OJO: nElemsxIt debe ser divisor de width*height (pq no tengo en cuenta strip)

void next_NR(int im[], int * i, int * j, int data[]){
	int iAux = *i, jAux = *j;
	for (int k = 0; k < nElemsxIt; k++){
		for(int ii = iAux - 2; ii <= iAux + 2; ii++)
			for(int jj = jAux - 2; jj <= jAux + 2; jj++)
				data[jj] = im[ii*width + jj];

		jAux += 1;
		if(jAux == width - 2){
			//He acabado con el i
			iAux++;
			if(iAux == height - 2){
				//He acabado con todo
				return;
			}
		}
	}

	*i = iAux;
	*j = jAux;
}

void next_Gx(int im[], int * i, int * j, int data[]){
	int iAux = *i, jAux = *j;
	for (int k = 0; k < nElemsxIt; k++){
		for (int ii = iAux - 2; ii <= iAux + 2; ii++)
			for (int jj = jAux - 2; jj <= jAux + 2; jj++){
				if (jj == jAux) jj++;
				data[jj] = im[ii*width + jj];
			}

		jAux += 1;
		if (jAux == width - 2){
			//He acabado con el i
			iAux++;
			if (iAux == height - 2){
				//He acabado con todo
				return;
			}
		}
	}

	*i = iAux;
	*j = jAux;
}

void next_Gy(int im[], int * i, int * j, int data[]){
	int iAux = *i, jAux = *j;
	for (int k = 0; k < nElemsxIt; k++){
		for (int ii = iAux - 2; ii <= iAux + 2; ii++){
			if (ii == iAux) ii++;
			for (int jj = jAux - 2; jj <= jAux + 2; jj++)
				data[jj] = im[ii*width + jj];
		}

		jAux += 1;
		if (jAux == width - 2){
			//He acabado con el i
			iAux++;
			if (iAux == height - 2){
				//He acabado con todo
				return;
			}
		}
	}

	*i = iAux;
	*j = jAux;
}


void canny(int *im, uint8_t *image_out, float level)
{
	float PI = 3.141593;
	float phi_aux[width * height];

	float lowthres, hithres;
	

	int mask_NR_base[VLen_NR] ={
		2,4,5,4,2,
		4,9,12,9,4,
		5,12,15,12,
		5,4,9,12,9,
		4,2,4,5,4,2};


	int mask_NR[VLen_NR * nElemsxIt];
	for(int i = 0; i < VLen_NR * nElemsxIt; i++){
		mask_NR[i] = mask_NR_base[i%VLen_NR];
	}

	int i=2,j=2;
	int termino = 0;
	int data_NR[VLen_NR * nElemsxIt];
	while(termino == 0){
		int iAux = i, jAux = j;
		next_NR(im, &i, &j, data_NR);
		//vec_mul_asm(N_DATA, in_data, inOut_data)
		vec_mul_asm(VLen_NR * nElemsxIt, mask_NR, data_NR);

		//Suma los resultados
		for (int k = 0; k < nElemsxIt; k++){
			int value = 0;
			for (int i = 0; i < VLen_NR; i++){
				value += data_NR[k*VLen_NR + i];
			}

			NR[iAux*width+jAux] = value/159;
			jAux++;
			if(jAux >= width - 2){
				jAux = 0;
				iAux++;
				if(iAux >= height -2 ){
					termino = 1;
					break;
				} 
			}
		}
	}


	int mask_Gx_base[VLen_G] = {
		1,2,-2,-1,
		4,8,-8,-4,
		6,12,-12,-6,
		4,8,-8,-4,
		1,2,-2,-1};

	int mask_Gy_base[VLen_G] = {
		-1,-4,-6,-4, -1,
		-2,-8,-12,-8,-2,
		2,8,12,8,2,
		1,4,6,4,1};


	int mask_Gx[VLen_G * nElemsxIt];
	for(int i = 0; i < VLen_G * nElemsxIt; i++){
		mask_Gx[i] = mask_Gx_base[i%VLen_G];
	}

	int mask_Gy[VLen_G * nElemsxIt];
	for(int i = 0; i < VLen_G * nElemsxIt; i++){
		mask_Gy[i] = mask_Gy_base[i%VLen_G];
	}


	i=2,j=2;
	termino = 0;
	int data_Gx[VLen_G * nElemsxIt], data_Gy[VLen_G * nElemsxIt];
	while(termino == 0){
		int iAux = i, jAux = j;
		next_Gx(NR, &i, &j, data_Gx);
		//vec_mul_asm(N_DATA, in_data, inOut_data)
		vec_mul_asm(VLen_G * nElemsxIt, mask_Gx, data_Gx);

		next_Gy(NR, &i, &j, data_Gy);	
		vec_mul_asm(VLen_G * nElemsxIt, mask_Gy, data_Gy);

		int value_Gx = 0, value_Gy = 0;
		//Suma los resultados
		for (int k = 0; k < nElemsxIt; k++){
			value_Gx = 0;
			value_Gy = 0;
			for (int i = 0; i < VLen_G; i++){
				value_Gx += data_Gx[k*VLen_G + i];
				value_Gy += data_Gy[k*VLen_G + i];
			}

			Gx[iAux*width+jAux] = value_Gx;
			Gy[iAux*width+jAux] = value_Gy;

			G[iAux*width+jAux]   = sqrtf((Gx[iAux*width+jAux]*Gx[iAux*width+jAux])+(Gy[iAux*width+jAux]*Gy[iAux*width+jAux]));	//G = √Gx²+Gy²
			phi_aux[iAux*width+jAux] = atan2f(fabs(Gy[iAux*width+jAux]),fabs(Gx[iAux*width+jAux]));

			if(fabs(phi_aux[iAux*width+jAux])<=PI/8 )
				phi[iAux*width+jAux] = 0;
			else if (fabs(phi_aux[iAux*width+jAux])<= 3*(PI/8))
				phi[iAux*width+jAux] = 45;
			else if (fabs(phi_aux[iAux*width+jAux]) <= 5*(PI/8))
				phi[iAux*width+jAux] = 90;
			else if (fabs(phi_aux[iAux*width+jAux]) <= 7*(PI/8))
				phi[iAux*width+jAux] = 135;
			else phi[iAux*width+jAux] = 0;

			jAux++;
			if(jAux >= width - 2){
				jAux = 0;
				iAux++;
				if(iAux >= height -2 ){
					termino = 1;
					break;
				} 
			}
		}
		
	}

	// Edge
	for(i=3; i<height-3; i++)
		for(j=3; j<width-3; j++)
		{
			pedge[i*width+j] = 0;
			if(phi[i*width+j] == 0){
				if(G[i*width+j]>G[i*width+j+1] && G[i*width+j]>G[i*width+j-1]) //edge is in N-S
					pedge[i*width+j] = 1;

			} else if(phi[i*width+j] == 45) {
				if(G[i*width+j]>G[(i+1)*width+j+1] && G[i*width+j]>G[(i-1)*width+j-1]) // edge is in NW-SE
					pedge[i*width+j] = 1;

			} else if(phi[i*width+j] == 90) {
				if(G[i*width+j]>G[(i+1)*width+j] && G[i*width+j]>G[(i-1)*width+j]) //edge is in E-W
					pedge[i*width+j] = 1;

			} else if(phi[i*width+j] == 135) {
				if(G[i*width+j]>G[(i+1)*width+j-1] && G[i*width+j]>G[(i-1)*width+j+1]) // edge is in NE-SW
					pedge[i*width+j] = 1;
			}
		}

	// Hysteresis Thresholding
	lowthres = level/2;
	hithres  = 2*(level);
	int ii, jj;
	for(i=3; i<height-3; i++)
		for(j=3; j<width-3; j++)
		{
			image_out[i*width+j] = 0;
			if(G[i*width+j]>hithres && pedge[i*width+j])
				image_out[i*width+j] = 255;
			else if(pedge[i*width+j] && G[i*width+j]>=lowthres && G[i*width+j]<hithres)
				// check neighbours 3x3
				for (ii=-1;ii<=1; ii++)
					for (jj=-1;jj<=1; jj++)
						if (G[(i+ii)*width+j+jj]>hithres)
							image_out[i*width+j] = 255;
		}
}


void getlines(int threshold, uint32_t *accumulators, int accu_width, int accu_height,  
	float *sin_table, float *cos_table,
	int *x1_lines, int *y1_lines, int *x2_lines, int *y2_lines, int *lines)
{
	int rho, theta, ii, jj;
	uint32_t max;

	for(rho=0;rho<accu_height;rho++)
	{
		for(theta=0;theta<accu_width;theta++)  
		{  

			if(accumulators[(rho*accu_width) + theta] >= threshold)  
			{  
				//Is this point a local maxima (9x9)  
				max = accumulators[(rho*accu_width) + theta]; 
				for(int ii=-4;ii<=4;ii++)  
				{  
					for(int jj=-4;jj<=4;jj++)  
					{  
						if( (ii+rho>=0 && ii+rho<accu_height) && (jj+theta>=0 && jj+theta<accu_width) )  
						{  
							if( accumulators[((rho+ii)*accu_width) + (theta+jj)] > max )  
							{
								max = accumulators[((rho+ii)*accu_width) + (theta+jj)];
							}  
						}  
					}  
				}  

				if(max == accumulators[(rho*accu_width) + theta]) //local maxima
				{
					int x1, y1, x2, y2;  
					x1 = y1 = x2 = y2 = 0;  

					if(theta >= 45 && theta <= 135)  
					{
						if (theta>90) {
							//y = (r - x cos(t)) / sin(t)  
							x1 = width/2;  
							y1 = ((float)(rho-(accu_height/2)) - ((x1 - (width/2) ) * cos_table[theta])) / sin_table[theta] + (height / 2);
							x2 = width;  
							y2 = ((float)(rho-(accu_height/2)) - ((x2 - (width/2) ) * cos_table[theta])) / sin_table[theta] + (height / 2);  
						} else {
							//y = (r - x cos(t)) / sin(t)  
							x1 = 0;  
							y1 = ((float)(rho-(accu_height/2)) - ((x1 - (width/2) ) * cos_table[theta])) / sin_table[theta] + (height / 2);
							x2 = width*2/5;  
							y2 = ((float)(rho-(accu_height/2)) - ((x2 - (width/2) ) * cos_table[theta])) / sin_table[theta] + (height / 2); 
						}
					} else {
						//x = (r - y sin(t)) / cos(t);  
						y1 = 0;  
						x1 = ((float)(rho-(accu_height/2)) - ((y1 - (height/2) ) * sin_table[theta])) / cos_table[theta] + (width / 2);  
						y2 = height;  
						x2 = ((float)(rho-(accu_height/2)) - ((y2 - (height/2) ) * sin_table[theta])) / cos_table[theta] + (width / 2);  
					}
					x1_lines[*lines] = x1;
					y1_lines[*lines] = y1;
					x2_lines[*lines] = x2;
					y2_lines[*lines] = y2;
					(*lines)++;
				}
			}
		}
	}
}

void init_cos_sin_table(float *sin_table, float *cos_table, int n)
{
	int i;
	for (i=0; i<n; i++)
	{
		sin_table[i] = sinf(i*DEG2RAD);
		cos_table[i] = cosf(i*DEG2RAD);
	}
}

void houghtransform(uint8_t *im, uint32_t *accumulators, int accu_width, int accu_height, 
	float *sin_table, float *cos_table)
{
	int i, j, theta;

	float hough_h = ((sqrt(2.0) * (float)(height>width?height:width)) / 2.0);

	for(i=0; i<accu_width*accu_height; i++)
		accumulators[i]=0;	

	float center_x = width/2.0; 
	float center_y = height/2.0;
	for(i=0;i<height;i++)  
	{  
		for(j=0;j<width;j++)  
		{  
			if( im[ (i*width) + j] > 250 ) // Pixel is edge  
			{  
				for(theta=0;theta<180;theta++)  
				{  
					float rho = ( ((float)j - center_x) * cos_table[theta]) + (((float)i - center_y) * sin_table[theta]);
					accumulators[ (int)((round(rho + hough_h) * 180.0)) + theta]++;

				} 
			} 
		} 
	}
}

void line_asist_CPU(int *im, 
	uint8_t *imEdge, 
	float *sin_table, float *cos_table, 
	uint32_t *accum, int accu_height, int accu_width,
	int *x1, int *y1, int *x2, int *y2, int *nlines)
{
	int threshold;
	long cycles, instret;

	/* Canny */

   	cycles = - read_cycles();
   	instret = - read_instret();
	
	canny(im, imEdge, 1000.0f); //level
	
	cycles += read_cycles();
   	instret += read_instret();
	printf("Canny\nCycles: %ld \nInstret: %ld \n", cycles, instret);

	/* hough transform */
	cycles = - read_cycles();
    	instret = - read_instret();
	
	houghtransform(imEdge, accum, accu_width, accu_height, sin_table, cos_table);
	
	cycles += read_cycles();
	instret += read_instret();
	printf("Hough transform\nCycles: %ld \nInstret: %ld \n", cycles, instret);

	if (width>height) threshold = width/6;
	else threshold = height/6;

	cycles = - read_cycles();
    	instret = - read_instret();
	
	getlines(threshold, accum, accu_width, accu_height, 
		sin_table, cos_table,
		x1, y1, x2, y2, nlines);
	
	cycles += read_cycles();
	instret += read_instret();
	printf("Getlines\nCycles: %ld \nInstret: %ld \n", cycles, instret);
}

