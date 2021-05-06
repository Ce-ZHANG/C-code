/***************************************************************************
 *
 *   File        : task3.c
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
double* ThomasAlgorithm(double *a, double *b, double *c, double *q, int n){
	int i;
	double *aa, *qq;
	double *x;
	aa = (double *)malloc(n*sizeof(double));
	qq = (double *)malloc(n*sizeof(double));
	x = (double *)malloc(n*sizeof(double));

	aa[0] = a[0];
	qq[0] = q[0];
	for(i = 1; i < n; i++){
		aa[i]=a[i]-c[i]/aa[i-1]*b[i-1];
		qq[i]=q[i]-c[i]/aa[i-1]*qq[i-1];
	}
	x[n-1] = qq[n-1]/aa[n-1];
	for(i = n-2; i >= 0; i--){
		x[i] = (qq[i]-b[i]*x[i+1])/aa[i];
	}
	free(aa);
	free(qq);
	return x;
}
void linalgbsys(const char* q4_file)
{
 	double *a,*b,*c,*q,*x;
	/*
	 * Read the data from the q4_file
	 */
    FILE *fp;
    fp = fopen(q4_file,"r");
    fseek(fp, 8L, SEEK_SET);
    for(int i=0; i<5; i++){
    	fscanf(fp,"%lf,%lf,%lf,%lf",&a[i],&b[i],&c[i],&q[i]);
    }
    x = ThomasAlgorithm(a,b,c,q,5);
	/*
	 * Output the result x into out_linalsys.csv
	 */
	FILE *fw;
	fw = fopen("out_linalsys.csv","w");
	fprintf(fw,"x\n");
	for(int i = 0; i < 5; i++){
		fprintf(fw,"%.4f\n",x[i]);
	}
	free(x);
	fclose(fp);
	fclose(fw);

}
