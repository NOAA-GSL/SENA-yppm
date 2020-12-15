#define _GNU_SOURCE
#include <stdio.h>
#include <unistd.h>  // gethostname, getopt
#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef __APPLE__
#include <sched.h>   // sched_getaffinity
extern void runnable (cpu_set_t *, int *, int *);
void print_affinity_ ()
{
  char hnbuf[64];
  int thread = 0;
  int lo;
  int hi;
  cpu_set_t coremask;
  gethostname (hnbuf, sizeof (hnbuf));
#pragma omp parallel private (thread, coremask, lo, hi)
  {
#ifdef _OPENMP
    thread = omp_get_thread_num ();
#endif
    // Passing zero means use the calling process
    sched_getaffinity (0, sizeof (coremask), &coremask);
    runnable (&coremask, &lo, &hi);
#pragma omp critical
    {
      printf ("Thread %d on %s. (Runnable range: lo=%d hi=%d)\n",
	      thread, hnbuf, lo, hi);
      fflush (stdout);
    }
  }
}
#else
void print_affinity_ ()
{
 printf("print_affinity is not supported on Mac OS\n");
}
#endif 
