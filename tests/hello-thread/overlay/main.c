#include <stdio.h>
#include <stdlib.h>
#include <unistd.h> //Header file for sleep(). man 3 sleep for details.
#include <pthread.h>
#include <dirent.h>

/* Time */
#include <sys/time.h>
#include <sys/resource.h>


static struct timeval tv0;
double get_time()
{
	double t;
	gettimeofday(&tv0, (struct timezone*)0);
	t = ((tv0.tv_usec) + (tv0.tv_sec)*1000000);

	return (t);
}

unsigned long read_cycles(void)
{
  unsigned long cycles;
  asm volatile ("rdcycle %0" : "=r" (cycles));
  return cycles;
}


unsigned long read_time(void)
{
  unsigned long time;
  asm volatile ("rdtime %0" : "=r" (time));
  return time;
}

unsigned long read_instret(void)
{
  unsigned long instret;
  asm volatile ("rdinstret %0" : "=r" (instret));
  return instret;
}


// CONSTANTES
#define N_threads 2
#define N 10000000
const int vivo = 100000;
const int N_times = 8;
int rand_A_vals[3*N], rand_B_vals[3*N];

void *myThreadFun(void *vargp)
{
	
	int a,b,c;
	int * nThread = (int *) vargp;
	for(int j=0; j < N_times; j++){
		for (int i=0; i < N; i++){
			a = rand_A_vals[i + *nThread * N];
			b = rand_B_vals[i + *nThread * N];
			c = a + b;
			if (i%vivo == 0) {
				printf("En proceso...%d\n", i);
			}
		}
	}

	return NULL;
}


int main()
{
	
	printf("Create random values");
	
	for (int i=0; i<3*N; i++){
		rand_A_vals[i] = rand();
		rand_B_vals[i] = rand();
	}
	printf("Arquitectura\n");
	system("cat /proc/cpuinfo");	
	printf("Create threads\n");

	long cycles, time, instret;
	double t;

	cycles = - read_cycles();
	time = - read_time();
	instret = - read_instret();
	t = - get_time();

	printf("Cycles: %ld \nTime: %ld \nTime library: %d \nInstret: %ld \n", cycles, time, t, instret);

	pthread_t tid[4];
	int nThread[N_threads];
	for (int i = 0; i < N_threads; i++) {
		nThread[i] = i;
    	pthread_create(&tid[i], NULL, myThreadFun, (void *) &nThread[i]);
    }
	for (int i = 0; i < N_threads; i++)
		pthread_join(tid[i], NULL);	

	cycles += read_cycles();
	time += read_time();
	instret += read_instret();
	t += get_time();


	printf("After threads\n");
	printf("Cycles: %ld \nTime: %ld \nTime library: %d \nInstret: %ld \n", cycles, time, t, instret);

	exit(0);
}


