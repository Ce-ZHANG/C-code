/***************************************************************************
 *
 *   File        : task4.c
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
double* Thomas(double *a, double *b, double *c, double *q, int n){
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
double Lagrange(double *x, double *y, double xo){
	double f = 0;
	double xx[3] = {x[1],x[2],x[3]};
	double yy[3] = {y[1],y[2],y[3]};
	f=((xo-xx[1])/(xx[0]-xx[1]))*((xo-xx[2])/(xx[0]-xx[2]))*yy[0]+
	     ((xo-xx[0])/(xx[1]-xx[0]))*((xo-xx[2])/(xx[1]-xx[2]))*yy[1]+
	       ((xo-xx[0])/(xx[2]-xx[0]))*((xo-xx[1])/(xx[2]-xx[1]))*yy[2];
   	return f;
}
double Cubicsplines(double *x, double *y, double xo, int n){
	int i;
	double s;
	double *a, *b, *c, *d, *h;
	a = (double *)malloc(n*sizeof(double));
	c = (double *)malloc(n*sizeof(double));
	h = (double *)malloc((n-1)*sizeof(double));
	b = (double *)malloc((n-1)*sizeof(double));
	d = (double *)malloc((n-1)*sizeof(double));
	double *aa, *bb, *cc, *qq;
	aa = (double *)malloc(n*sizeof(double));
	bb = (double *)malloc(n*sizeof(double));
	cc = (double *)malloc(n*sizeof(double));
	qq = (double *)malloc(n*sizeof(double));
	for (i = 0; i < n; i++){
		a[i] = y[i];
	}
	for (i = 0; i <n-1; i++){
		h[i] = x[i+1]-x[i];
	}
	aa[0] = 1;
	aa[n-1] = 1;
	bb[0] = 0;
	bb[n-1] = 0;
	cc[0] = 0;
	cc[n-1] = 0;
	qq[0] = 0;
	qq[n-1] = 0;
	for (i = 1; i < n-1; i++){
		aa[i] = 2*(h[i-1]+h[i]);
		bb[i] = h[i];
		cc[i] = h[i-1];
		qq[i] = 3*((a[i+1]-a[i])/h[i]+(a[i-1]-a[i])/h[i-1]);
	}
	c = Thomas(aa,bb,cc,qq,4);
	for (i = 0; i < 3; i++){
		d[i] = (c[i+1]-c[i])/(3*h[i]);
		b[i] = (a[i+1]-a[i])/h[i]-(2*c[i]+c[i+1])*h[i]/3;
	}
	s = a[2]+b[2]*(xo-x[2])+c[2]*pow((xo-x[2]),2)+d[2]*pow((xo-x[2]),3);
	free(a);
	free(b);
	free(c);
	free(d);
	free(h);
	free(aa);
	free(bb);
	free(cc);
	free(qq);
	return s;
}
void interp(const char* q5_file, const double xo)
{
	double x[4];
	double y[4];
	double outputy = 0;
    FILE *fp;
    fp = fopen(q5_file,"r");
    fseek(fp, 7L, SEEK_SET);
    for (int i = 0; i < 4; i++){
    	fscanf(fp,"%lf,%lf",&x[i],&y[i]);
    }
    FILE *fw;
    fw = fopen("out_interp.csv","w");
    outputy = Lagrange(x,y,xo);
    fprintf(fw,"lagrange\n");
    fprintf(fw,"%.4f\n",outputy);
    outputy = Cubicsplines(x,y,xo,4);
    fprintf(fw,"cubic\n");
    fprintf(fw,"%.4f\n",outputy);
    fclose(fp);
    fclose(fw);
}
