/***************************************************************************
 *
 *   File        : main.c
 *   Student Id  : <758397>
 *   Name        : <JUSTIN BUGEJA>
 *
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <string.h>
#include "tasks.h"

#define TIMEDIV 1000.0

int main(int argc, char *argv[]) {

	if (argc < 5) {
		printf("NOT ENOUGH ARGUMENTS\n");
		exit(EXIT_FAILURE);
	}

	char* q2_file = argv[1];
	char* q3_file = argv[2];
	char* q5_file = argv[3];
	double xo = atoi(argv[4]);
	char* q6_file = argv[5];

	/* timing for each task and output running time in ms */
	struct timeval start;
	struct timeval stop;

	/* Question 2 */
	gettimeofday(&start, NULL);
	shockwave(q2_file);
	gettimeofday(&stop, NULL);
	double elapsed_ms = (stop.tv_sec - start.tv_sec) * TIMEDIV;
	elapsed_ms += (stop.tv_usec - start.tv_usec) / TIMEDIV;
	printf("QUESTION 2: %.2f ms\n", elapsed_ms);

	/* Question 3 */
	gettimeofday(&start, NULL);
	linalgbsys(q3_file);
	gettimeofday(&stop, NULL);
	elapsed_ms = (stop.tv_sec - start.tv_sec) * TIMEDIV;
	elapsed_ms += (stop.tv_usec - start.tv_usec) / TIMEDIV;
	printf("QUESTION 3: %.2f ms\n", elapsed_ms);

	/* Question 5 */
	gettimeofday(&start, NULL);
	interp(q5_file, xo);
	gettimeofday(&stop, NULL);
	elapsed_ms = (stop.tv_sec - start.tv_sec) * TIMEDIV;
	elapsed_ms += (stop.tv_usec - start.tv_usec) / TIMEDIV;
	printf("QUESTION 5: %.2f ms\n", elapsed_ms);

	/* Question 6 */
	gettimeofday(&start, NULL);
	heateqn(q6_file);
	gettimeofday(&stop, NULL);
	elapsed_ms = (stop.tv_sec - start.tv_sec) * TIMEDIV;
	elapsed_ms += (stop.tv_usec - start.tv_usec) / TIMEDIV;
	printf("QUESTION 6: %.2f ms\n", elapsed_ms);

	return (EXIT_SUCCESS);
}
