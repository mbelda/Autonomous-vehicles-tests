#ifndef ROUTINES_H
#define ROUTINES_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" 
#endif
void init_cos_sin_table(float *sin_table, float *cos_table, int n);

#ifdef __cplusplus
extern "C" 
#endif
void line_asist_CPU(int *im,
	uint8_t *imEdge,
	float *sin_table, float *cos_table, 
	uint32_t *accum, int accu_height, int accu_width,
	int *x1, int *y1, int *x2, int *y2, int *nlines);

#endif

unsigned long read_cycles(void);
unsigned long read_instret(void);
