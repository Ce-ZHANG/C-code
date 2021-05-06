/***************************************************************************
 *
 *   File        : task6.c
 *   Student Id  : 906554
 *   Name        : SHUOYANG QIN
 *
 ***************************************************************************/
#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <sys/time.h>
#include <string.h>
#include "tasks.h"

#define PI 3.1415926

double* TA(double *a, double *b, double *c, double *q, int n){
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


double initial(double x){
	double f;
	if(x >= 0 && x < 0.125){
		f = 0;
	}
	if(x >= 0.125 && x <= 0.375){
		f = 0.5*(1-cos(8*PI*(x-0.125)));
	}
	if(x > 0.375 && x <= 1){
		f = 0;
	}
	return f;
}
/*
 * Fixed ends for explicit routine
 */
double* RHS_Explicit_FixedEnd(int Nt, int Nx, double mu){

	int n, i;
	double x;
	double dt = 2 / (double)Nt;
	double dx = 1 / (double)Nx;
	double **fx = (double **)malloc((Nt+1) * sizeof(double*));
	for(n = 0; n < Nt+1; n++)
	{
		fx[n] = (double *)malloc((Nx+1) * sizeof(double));
	}
	for(n = 0; n < Nt+1; n++){
		fx[n][0] = 0;
		fx[n][Nx] = 0;
   }
   for(i = 0; i < Nx+1; i++){
 	  x = dx * i;
      fx[0][i] = initial(x);
  }
  for(n = 0; n < Nt; n++){
    for(i = 1; i < Nx; i++){
      fx[n+1][i]=fx[n][i]+dt*mu*(fx[n][i+1]-2*fx[n][i]+fx[n][i-1])/pow(dx,2);
      }
	}
  return fx;
}
/*
 * Variables ends for explicit routine
 */
double** RHS_Explicit_VariablesEnd(int Nt, int Nx, double mu){
	int n, i;
	double x;
	double dt = 2 / (double)Nt;
	double dx = 1 / (double)Nx;
	double **fx = (double **)malloc((Nt+1) * sizeof(double*));
	for(n = 0; n < Nt+1; n++){
		fx[n] = (double *)malloc((Nx+1) * sizeof(double));
	}
	for(i = 0; i < Nx+1; i++){
			x = dx * i;
			fx[0][i] = initial(x);
		}
	for(n = 0; n < Nt; n++){
		for(i = 1; i < Nx; i++){
			fx[n+1][i]=fx[n][i]+dt*mu*(fx[n][i+1]-2*fx[n][i]+fx[n][i-1])/pow(dx,2);
		}
		fx[n+1][0] = fx[n][0]+dt*mu*(fx[n][0]-2*fx[n][1]+fx[n][2])/pow(dx,2);
		fx[n+1][Nx]= fx[n][Nx]+dt*mu*(fx[n][Nx]-2*fx[n][Nx-1]+fx[n][Nx-2])/pow(dx,2);
	}
	return fx;
}
/*
 * Fixed ends for Implicit routine
 */
double* RHS_Implicit_FixedEnd(int Nt, int Nx, double mu){
	int n, i;
	double x;
	double dt = 2 / (double)Nt;
	double dx = 1 / (double)Nx;

	double *fn = (double *)malloc((Nx+1) * sizeof(double));
	double *fn1 = (double *)malloc((Nx+1) * sizeof(double));
	double *a = (double *)malloc((Nx+1) * sizeof(double));
	double *b = (double *)malloc((Nx+1) * sizeof(double));
	double *c = (double *)malloc((Nx+1) * sizeof(double));
	for(i = 1; i < Nx; i++){
		a[i] = 2*(dt*mu/pow(dx,2))+1;
		b[i] = -(dt*mu/pow(dx,2));
		c[i] = -(dt*mu/pow(dx,2));
	}
		a[0] = 1;
		a[Nx] = 1;
		b[0] = 0;
		b[Nx] = 0;
		c[0] = 0;
		c[Nx] = 0;
	for(i = 0; i < Nx+1; i++){
		x = dx * i;
		fn[i] = initial(x);
	}
	for(n = 0; n < 100; n++){
		fn[0] = 0 ;
		fn[Nx] = 0;
		for(i = 0; i < Nx+1; i++){
			fn1 = TA(a,b,c,fn,Nx+1);
		}
		fn =fn1;
	}
	free(a);
	free(b);
	free(c);
	free(fn1);
	return fn;
}
void waveeqn(const char* q6_file)
{
	double mu=0;
	int Nt=0;
    int Nx=0;
	int i,n;
    double x[Nx+1];
	FILE *fp;
	fp = fopen(q6_file,"r");
	fseek(fp, 9L, SEEK_SET);
	fscanf(fp,"%lf,%d,%d", &mu, &Nx, &Nt);
	fclose(fp);

    for(i = 0; i < Nx+1; i++){
    	x[i]= i /(double)Nx;
    }

  	FILE *fw1;
	fw1 = fopen("out_waveeqn_EU.csv","w");
	fprintf(fw1,"x,F(x)\n");
    double **f1 = RHS_Explicit_FixedEnd(Nt,Nx,mu);
	for(i = 0; i < Nx+1; i++){
  		fprintf(fw1,"%.6f,%.6f\n",x[i],f1[99][i]);
	}
	fclose(fw1);

	FILE *fw2;
	fw2 = fopen("out_waveeqn_LW.csv","w");
	fprintf(fw2,"x,F(x)\n");
	f1 = RHS_Explicit_VariablesEnd(Nt,Nx,mu);
	for(i = 0; i < Nx+1; i++){
		fprintf(fw2,"%.6f,%.6f\n",x[i],f1[99][i]);
	}
	fclose(fw2);

	FILE *fw3;
	fw3 = fopen("out_heateqn_implicit_fe.csv","w");
	fprintf(fw3,"x,F(x)\n");
	double *f2 = RHS_Implicit_FixedEnd(Nt,Nx,mu);
	for(i = 0; i < Nx+1; i++){
		fprintf(fw3,"%.4f,%.4f\n",x[i],f2[i]);
	}
	fclose(fw3);

	for(n = 0; n < Nt+1; n++)
	{
		free(f1[n]);
	}
	free(f1);
	free(f2);

	//Nx = 80.0;
	//CFL = 1.002;
	//t6.dx = spatial_interval / Nx;  // calculate delta x
	//t6.x = (double*)calloc(Nx + 1, sizeof(double));

	//t6.dt = t6.dx * CFL / c;  // calculate delta t
	//double Nt = time_interval / t6.dt;
	//t6.t = (double*)calloc(Nt + 1, sizeof(double));


	//for (int i = 0; i <= Nx; i++)
	//{
	//	t6.x[i] = t6.dx * i;   //construct array x
	//						   /*printf("%lf\n", t6.x[i]);*/
	//}

	//for (int i = 0; i <= Nt; i++)
	//{
	//	t6.t[i] = t6.dt * i;   //construct array t
	//}


	//initalize(Nx, Nt);  // set the initial conditions of both methods
	//Euler_explicit(Nx, Nt, c);
	//if (Nx == 80.0&&CFL == 1.002) 
	//{
	//	FILE *fw6_2_1;
	//	fw6_2_1 = fopen("Euler_Nx80_CFL1.5.csv", "w");

	//	FILE *fw6_2_2;
	//	fw6_2_2 = fopen("Lax_Nx80_CFL1.5.csv", "w");
	//	for (int j = 2; j <= 10; j = j + 2)
	//	{
	//		int x = j / t6.dt;
	//		fprintf(fw6_2_1, "x, t = %d\n", j);
	//		for (int i = 0; i < Nx + 1; i++)
	//		{
	//			fprintf(fw6_2_1, "%.6f,%.6f\n", t6.x[i], t6.fE[i][x]);
	//		}
	//		/*for (int i = 0; i < Nx + 1; i++)
	//		{
	//		printf("%.6f,%.6f\n", t6.x[i], t6.fE[i][2]);
	//		}*/

	//		Lax(Nx, Nt, c);
	//		fprintf(fw6_2_2, "x, t = %d\n", j);
	//		for (int i = 0; i < Nx + 1; i++)
	//		{
	//			fprintf(fw6_2_2, "%.6f,%.6f\n", t6.x[i], t6.fL[i][x]);
	//		}

	//	}
	//	fclose(fw6_2_1);
	//	fclose(fw6_2_2);
	//}
	//if (Nx == 80.0&&CFL == 1.0) 
	//{
	//	FILE *fw6_2_1;
	//	fw6_2_1 = fopen("Euler_Nx80_CFL1.csv", "w");

	//	FILE *fw6_2_2;
	//	fw6_2_2 = fopen("Lax_Nx80_CFL1.csv", "w");
	//	for (int j = 2; j <= 10; j = j + 2)
	//	{
	//		int x = j / t6.dt;
	//		/*fprintf(fw6_2_1, "x,F(x)\n");*/
	//		for (int i = 0; i < Nx + 1; i++)
	//		{
	//			fprintf(fw6_2_1, "%.6f,%.6f\n", t6.x[i], t6.fE[i][x]);
	//		}
	//		/*for (int i = 0; i < Nx + 1; i++)
	//		{
	//		printf("%.6f,%.6f\n", t6.x[i], t6.fE[i][2]);
	//		}*/

	//		Lax(Nx, Nt, c);
	//		/*fprintf(fw6_2_2, "x,F(x)\n");*/
	//		for (int i = 0; i < Nx + 1; i++)
	//		{
	//			fprintf(fw6_2_2, "%.6f,%.6f\n", t6.x[i], t6.fL[i][x]);
	//		}

	//	}
	//	fclose(fw6_2_1);
	//	fclose(fw6_2_2);
	//}
	//if (Nx == 80.0&&CFL == 0.75)
	//{
	//	FILE *fw6_2_1;
	//	fw6_2_1 = fopen("Euler_Nx80_CFL0.75.csv", "w");

	//	FILE *fw6_2_2;
	//	fw6_2_2 = fopen("Lax_Nx80_CFL0.75.csv", "w");
	//	for (int j = 2; j <= 10; j = j + 2)
	//	{
	//		int x = j / t6.dt;
	//		fprintf(fw6_2_1, "x,F(x)\n");
	//		for (int i = 0; i < Nx + 1; i++)
	//		{
	//			fprintf(fw6_2_1, "%.6f,%.6f\n", t6.x[i], t6.fE[i][x]);
	//		}
	//		/*for (int i = 0; i < Nx + 1; i++)
	//		{
	//		printf("%.6f,%.6f\n", t6.x[i], t6.fE[i][2]);
	//		}*/

	//		Lax(Nx, Nt, c);
	//		fprintf(fw6_2_2, "x,F(x)\n");
	//		for (int i = 0; i < Nx + 1; i++)
	//		{
	//			fprintf(fw6_2_2, "%.6f,%.6f\n", t6.x[i], t6.fL[i][x]);
	//		}

	//	}
	//	fclose(fw6_2_1);
	//	fclose(fw6_2_2);
	//}
	//if (Nx == 80.0&&CFL == 0.25)
	//{
	//	FILE *fw6_2_1;
	//	fw6_2_1 = fopen("Euler_Nx80_CFL0.25.csv", "w");

	//	FILE *fw6_2_2;
	//	fw6_2_2 = fopen("Lax_Nx80_CFL0.25.csv", "w");
	//	for (int j = 2; j <= 10; j = j + 2)
	//	{
	//		int x = j / t6.dt;
	//		fprintf(fw6_2_1, "x,F(x)\n");
	//		for (int i = 0; i < Nx + 1; i++)
	//		{
	//			fprintf(fw6_2_1, "%.6f,%.6f\n", t6.x[i], t6.fE[i][x]);
	//		}
	//		/*for (int i = 0; i < Nx + 1; i++)
	//		{
	//		printf("%.6f,%.6f\n", t6.x[i], t6.fE[i][2]);
	//		}*/

	//		Lax(Nx, Nt, c);
	//		fprintf(fw6_2_2, "x,F(x)\n");
	//		for (int i = 0; i < Nx + 1; i++)
	//		{
	//			fprintf(fw6_2_2, "%.6f,%.6f\n", t6.x[i], t6.fL[i][x]);
	//		}

	//	}
	//	fclose(fw6_2_1);
	//	fclose(fw6_2_2);
	//}
	//if (Nx == 200.0&&CFL == 1.0)
	//{
	//	FILE *fw6_2_1;
	//	fw6_2_1 = fopen("Euler_Nx200_CFL1.csv", "w");

	//	FILE *fw6_2_2;
	//	fw6_2_2 = fopen("Lax_Nx200_CFL1.csv", "w");
	//	for (int j = 2; j <= 10; j = j + 2)
	//	{
	//		int x = j / t6.dt;
	//		fprintf(fw6_2_1, "x,t = %d\n", j);
	//		for (int i = 0; i < Nx + 1; i++)
	//		{
	//			fprintf(fw6_2_1, "%.6f,%.6f\n", t6.x[i], t6.fE[i][x]);
	//		}
	//		/*for (int i = 0; i < Nx + 1; i++)
	//		{
	//		printf("%.6f,%.6f\n", t6.x[i], t6.fE[i][2]);
	//		}*/

	//		Lax(Nx, Nt, c);
	//		fprintf(fw6_2_2, "x,t = %d\n", j);
	//		for (int i = 0; i < Nx + 1; i++)
	//		{
	//			fprintf(fw6_2_2, "%.6f,%.6f\n", t6.x[i], t6.fL[i][x]);
	//		}

	//	}
	//	fclose(fw6_2_1);
	//	fclose(fw6_2_2);
	//}
	//if (Nx == 200.0&&CFL == 0.75)
	//{
	//	FILE *fw6_2_1;
	//	fw6_2_1 = fopen("Euler_Nx200_CFL0.75.csv", "w");

	//	FILE *fw6_2_2;
	//	fw6_2_2 = fopen("Lax_Nx200_CFL0.75.csv", "w");
	//	for (int j = 2; j <= 10; j = j + 2)
	//	{
	//		int x = j / t6.dt;
	//		fprintf(fw6_2_1, "x,F(x)\n");
	//		for (int i = 0; i < Nx + 1; i++)
	//		{
	//			fprintf(fw6_2_1, "%.6f,%.6f\n", t6.x[i], t6.fE[i][x]);
	//		}
	//		/*for (int i = 0; i < Nx + 1; i++)
	//		{
	//		printf("%.6f,%.6f\n", t6.x[i], t6.fE[i][2]);
	//		}*/

	//		Lax(Nx, Nt, c);
	//		fprintf(fw6_2_2, "x,F(x)\n");
	//		for (int i = 0; i < Nx + 1; i++)
	//		{
	//			fprintf(fw6_2_2, "%.6f,%.6f\n", t6.x[i], t6.fL[i][x]);
	//		}

	//	}
	//	fclose(fw6_2_1);
	//	fclose(fw6_2_2);
	//}
	//if (Nx == 200.0&&CFL == 0.25)
	//{
	//	FILE *fw6_2_1;
	//	fw6_2_1 = fopen("Euler_Nx200_CFL0.25.csv", "w");

	//	FILE *fw6_2_2;
	//	fw6_2_2 = fopen("Lax_Nx200_CFL0.25.csv", "w");
	//	for (int j = 2; j <= 10; j = j + 2)
	//	{
	//		int x = j / t6.dt;
	//		fprintf(fw6_2_1, "x, t = %d\n", j);
	//		for (int i = 0; i < Nx + 1; i++)
	//		{
	//			fprintf(fw6_2_1, "%.6f,%.6f\n", t6.x[i], t6.fE[i][x]);
	//		}
	//		/*for (int i = 0; i < Nx + 1; i++)
	//		{
	//		printf("%.6f,%.6f\n", t6.x[i], t6.fE[i][2]);
	//		}*/

	//		Lax(Nx, Nt, c);
	//		fprintf(fw6_2_2, "x,t = %d\n", j);
	//		for (int i = 0; i < Nx + 1; i++)
	//		{
	//			fprintf(fw6_2_2, "%.6f,%.6f\n", t6.x[i], t6.fL[i][x]);
	//		}

	//	}
	//	fclose(fw6_2_1);
	//	fclose(fw6_2_2);
	//}

}
