/*+C
 ************************************************************************
 * This program reads in a mass spec and prints out the weights 
 *
 * This program reads in a CSV file (ignoring the first line), assumed
 * to have two columns: the first is a mass, the second a spectrum
 * intensity. It does NOT presume the masses are equispaced. It then
 * identifies peaks in the spectrum by forming a split-window
 * normalizer around each point, identifying sharp tight peaks, and
 * then prints out the peak value.
 *
 * This it does by reading in points 1 at a time and storing them in a
 * circular buffer. It maintains indices to the first and last entries
 * for the window and the peak areas, incrementing them as necessary
 * as more points are read in.
 ************************************************************************
 */

/* Includes */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* File-Scope Constants, Macros, and Enumerations */

/* This is the maximum number of points we can work on at any one time. */
#define BUFFER_SIZE (1024)


/* #define DPRINTF(x) */
#define DPRINTF(x) printf x

/* This macro applies a circular offset to an index */
#define INDEX(x) ((x)%BUFFER_SIZE)


/* This is the usage error */
#define USAGE(pname) \
  {printf("Usage: %s <-window_size #> <-peak_size #> " \
	  "<-threshold #> Filename\n",pname); exit(1);}

/* File-Scope Type Definitions */

/* File-Scope Variables */
static float masses[BUFFER_SIZE];
static float intensities[BUFFER_SIZE];

/* File-Scope Function Prototypes */
static float maxValue(float *input, int start, int end);

static float maxValue(float *input, int start, int end)
{
  float max = input[INDEX(start)];

  while(start <= end) {
    if (input[INDEX(start)] > max)
      max = input[INDEX(start)];
    start++;
  }

  return(max);
  
}
int main(int argc, char **argv)
{
  char *pname;
  char line[128];

  int num_points=0;
  int current_point = 0;

  int peak_start = 0;
  int window_start = 0;

  int peak_end = 0;
  int window_end = 0;

  float threshold=5.0;

  float peak_size=8;
  float window_size = 16.0;

  float peak_half_size, window_half_size;

  float window_sum, peak_sum, normalizer;

  FILE *fp;

  /* Store the program name */
  pname = argv[0]; argc--; argv++;

  /* Parse the arguments, leaving the last one as the file name */
  while (argc > 1) {

    /* Check for window width */
    if (strcmp(argv[0],"-window_size")) {
      argc--; argv++;
      if (argc == 0 || sscanf(argv[0],"%f",&window_size) != 1)
	USAGE(pname);
      argc--; argv++;
      continue;
    }

    /* Check for peak width */
    if (strcmp(argv[0],"-peak_size")) {
      argc--; argv++;
      if (argc == 0 || sscanf(argv[0],"%f",&peak_size) != 1)
	USAGE(pname);
      argc--; argv++;
      continue;
    }
    /* And the threshold */
    if (strcmp(argv[0],"-threshold")) {
      argc--; argv++;
      if (argc == 0 || sscanf(argv[0],"%f",&threshold) != 1)
	USAGE(pname);
      argc--; argv++;
      continue;
    }

    /* IF we got hear, there was a bad command line argument */
    USAGE(pname);
  }

  peak_half_size = peak_size/2;
  window_half_size = window_size/2;
  DPRINTF((" Detection %.1f %.1f %.1f\n",window_size,peak_size,threshold));
  
  /* Now open the file and get the first line */
  if ((fp = fopen(argv[0],"r")) == NULL ||
      fgets(line,sizeof(line),fp) == NULL) {
    DPRINTF(("Unable to open file or read first line\n"));
    exit(1);
  }

  while(1) {
    if (fgets(line,sizeof(line),fp) == NULL ||
	sscanf(line,"%f,%f",
	       masses+INDEX(num_points),
	       intensities+INDEX(num_points)) != 2) break;
    DPRINTF(("%d: %s",num_points,line));
    
    /* 
     * When the window first gets fully realized, initialize all the
     * indices. We can use the peak_end index as an indicator that the
     * window is not set yet
     */
    if (peak_end == 0) {

      /* Set the new window end to the current point */
      window_end = num_points;
      if (window_end - window_start >= BUFFER_SIZE) {
	printf("Out Of Bounds!!\n");
	exit(1);
      }

      /* Check to see if the window is complete */
      if (masses[INDEX(window_end)] -
	  masses[INDEX(window_start)] > window_size) {

	/* If so, set the current point */
	while(masses[INDEX(current_point)] <
	      masses[INDEX(window_start)] + window_half_size)
	  current_point++;

	/* 
	 * Now check to see if we have enough from the other end: this
	 * may not happen the first few times so we can't set the peak
	 * indices if not
	 */
	if (masses[INDEX(window_end)] -
	    masses[INDEX(current_point)] >= window_half_size) {

	  /* We have a fully realized situation: set the peak limits */
	  for (peak_start = current_point;
	       masses[INDEX(current_point)] -
		 masses[INDEX(peak_start)] < peak_half_size;
	       peak_start--)
	    ;

	  for (peak_end = current_point;
	       masses[INDEX(peak_end)] -
		 masses[INDEX(current_point)] < peak_half_size;
	       peak_end++)
	    ;
	}
      }
    } else {

      /* 
       * We have to update for the new point: first move the bottom of
       * the window up 
       */
      window_end = num_points;
      while(masses[INDEX(window_end)] -
	    masses[INDEX(window_start+2)] > window_size)
	window_start++;

      /* Now move the current point up, processing each point as we go */
      while (masses[INDEX(window_end)] -
	     masses[INDEX(current_point+1)] > window_half_size) {

	/* Reset the peak limits */
	while (masses[INDEX(current_point)] -
	       masses[INDEX(peak_start+2)] > peak_half_size)
	  peak_start++;

	while (masses[INDEX(peak_end)] -
	       masses[INDEX(current_point)] < peak_half_size)
	  peak_end++;

	DPRINTF((" Update current[%d,%d,%d,%d,%d] (%.0f %.0f %.0f %.f) .. ",
		 window_start,peak_start,current_point,peak_end,window_end,
		 maxValue(intensities,window_start, peak_start),
		 maxValue(intensities,peak_start, peak_end),
		 maxValue(intensities,peak_end, window_end),
		 intensities[INDEX(current_point)]));
	if (intensities[INDEX(current_point)] >
	    threshold * maxValue(intensities,window_start, peak_start) &&

	    intensities[INDEX(current_point)] >
	    threshold * maxValue(intensities,peak_end,window_end) &&
	    intensities[INDEX(current_point)] ==
	    maxValue(intensities,peak_start,peak_end))
	  printf("%.6f ",masses[INDEX(current_point)]);
	DPRINTF(("\n"));
	current_point++;
      }
    }

    /* Next point */
    num_points++;
  }
  printf("\n");
}
