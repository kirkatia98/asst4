/* Adapted from file Cycletimer.h in 15-418 code repository */

#if defined(__APPLE__)
//  #if defined(__x86_64__)
  #if 0
    #include <sys/sysctl.h>
  #else
    #include <mach/mach.h>
    #include <mach/mach_time.h>
  #endif // __x86_64__ or not


#elif _WIN32
  #include <windows.h>
  #include <time.h>
#else
  #include <string.h>
  #include <sys/time.h>
#endif

#include <stdio.h>  // fprintf
#include <stdlib.h> // exit
#include <stdint.h>
#include <stdbool.h>

#include "cycletimer.h"



  // This uses the cycle counter of the processor.  Different
  // processors in the system will have different values for this.  If
  // you process moves across processors, then the local time you
  // measure will likely be incorrect.  This is mostly for fine
  // grained measurements where the process is likely to be on the
  // same processor.  For more global things you should use the
  // Time interface.

  // Also note that if you processors' speeds change (i.e. processors
  // scaling) or if you are in a heterogenous environment, you will
  // likely get spurious results.

typedef uint64_t SysClock;

    //////////
    // Return the current CPU time, in terms of clock ticks.
    // Time zero is at some arbitrary point in the past.
static SysClock currentTicks() {
    //#if defined(__APPLE__) && !defined(__x86_64__)
#if defined(__APPLE__)
      return mach_absolute_time();
#elif defined(_WIN32)
      LARGE_INTEGER qwTime;
      QueryPerformanceCounter(&qwTime);
      return qwTime.QuadPart;
#elif defined(__x86_64__)
      unsigned int a, d;
      asm volatile("rdtsc" : "=a" (a), "=d" (d));
      return ((uint64_t) d << 32) + a;
#elif defined(__ARM_NEON__) && 0 // mrc requires superuser.
      unsigned int val;
      asm volatile("mrc p15, 0, %0, c9, c13, 0" : "=r"(val));
      return val;
#else
      timespec spec;
      clock_gettime(CLOCK_THREAD_CPUTIME_ID, &spec);
      return (uint64_t) spec.tv_sec * 1000 * 1000 * 1000 + spec.tv_nsec;
#endif
}

//////////
// Return the conversion from ticks to seconds.
static double secondsPerTick() {
    static bool initialized = false;
    static double secondsPerTick_val;
    if (initialized) return secondsPerTick_val;
#if defined(__APPLE__)
    #if 0
    //  #ifdef __x86_64__
    int args[] = {CTL_HW, HW_CPU_FREQ};
    unsigned int Hz;
    size_t len = sizeof(Hz);
    if (sysctl(args, 2, &Hz, &len, NULL, 0) != 0) {
	fprintf(stderr, "Failed to initialize secondsPerTick_val!\n");
	exit(-1);
    }
    secondsPerTick_val = 1.0 / (double) Hz;
#else
    mach_timebase_info_data_t time_info;
    mach_timebase_info(&time_info);

    // Scales to nanoseconds without 1e-9f
    secondsPerTick_val = (1e-9 * (double) time_info.numer)/(double) time_info.denom;
#endif // x86_64 or not
#elif defined(_WIN32)
      LARGE_INTEGER qwTicksPerSec;
      QueryPerformanceFrequency(&qwTicksPerSec);
      secondsPerTick_val = 1.0/(double) qwTicksPerSec.QuadPart;
#else
      FILE *fp = fopen("/proc/cpuinfo","r");
      char input[1024];
      if (!fp) {
         fprintf(stderr, "cycletimer failed: couldn't find /proc/cpuinfo.");
         exit(-1);
      }
      // In case we don't find it, e.g. on the N900
      secondsPerTick_val = 1e-9;
      while (!feof(fp) && fgets(input, 1024, fp)) {
	  // NOTE(boulos): Because reading cpuinfo depends on dynamic
	  // frequency scaling it's better to read the @ sign first
	  float GHz, MHz;
	  if (strstr(input, "model name")) {
	      char* at_sign = strstr(input, "@");
	      if (at_sign) {
		  char* after_at = at_sign + 1;
		  char* GHz_str = strstr(after_at, "GHz");
		  char* MHz_str = strstr(after_at, "MHz");
		  if (GHz_str) {
		      *GHz_str = '\0';
		      if (1 == sscanf(after_at, "%f", &GHz)) {
			  //printf("GHz = %f\n", GHz);
			  secondsPerTick_val = 1e-9f / GHz;
			  break;
		      }
		  } else if (MHz_str) {
		      *MHz_str = '\0';
		      if (1 == sscanf(after_at, "%f", &MHz)) {
			  //printf("MHz = %f\n", MHz);
			  secondsPerTick_val = 1e-6f / GHz;
			  break;
		      }
		  }
	      }
	  } else if (1 == sscanf(input, "cpu MHz : %f", &MHz)) {
	      //printf("MHz = %f\n", MHz);
	      secondsPerTick_val = 1e-6f / MHz;
	      break;
	  }
      }
      fclose(fp);
#endif
      
      initialized = true;
      return secondsPerTick_val;
}

//////////
// Return the current CPU time, in terms of seconds.
// This is slower than currentTicks().  Time zero is at
// some arbitrary point in the past.
double currentSeconds() {
    return currentTicks() * secondsPerTick();
}

