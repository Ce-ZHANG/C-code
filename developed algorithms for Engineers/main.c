/***************************************************************************
 *
 *   File        : main.c
 *   Student Id  : 906554
 *   Name        : SHUOYANG QIN
 *
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <string.h>
#include "tasks.h"

int main(int argc, char *argv[]) {
	
	/* TODO: Parse Command Line Arguments
	DONOT explicitly set arguments to filenames */
	//char* q2_file = argv[1];
	//char* q4_file = argv[1];
	//char* q5_file = argv[0];
	//double xo;
	char* q6_file = argv[1];

	/* TODO: Add timing for each task and output running time in ms */
    
	/* Question 2 */
	//shockwave(q2_file);
	
	/* Question 4 */
	//linalgbsys(q4_file);
	
	/* Question 5 */
	//interp(q5_file,xo);
	
	/* Question 6 */
	heateqn(q6_file);
    
	return (0);
}
