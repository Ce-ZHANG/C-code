/***************************************************************************
 *
 *   File        : task2.c
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

#define PI      3.141592635
#define diff    0.000000001
#define times   50
#define rad     180

double fun(double M, double gamma, double theta, double beta){
	double radtheta;
	double radbeta;
	radtheta = theta * PI /rad;
	radbeta = beta * PI /rad;
	double fun =  2*(1/tan(radbeta))*(pow(M,2)*pow(sin(radbeta),2)-1)/(pow(M,2)*(gamma+cos(2*radbeta))+2)-tan(radtheta);
	return fun;
}
double dfun(double M, double gamma, double beta){
	double radbeta;
	radbeta = beta * PI /rad;
	double dfun= (4*pow(M,2)*cos(radbeta)*(1/tan(radbeta))*sin(radbeta))/((gamma + cos(2*radbeta))*pow(M,2) + 2) - (2*(pow((1/tan(radbeta)),2) + 1)*(pow(M,2)*pow(sin(radbeta),2) - 1))/((gamma + cos(2*radbeta))*pow(M,2) + 2) + (4*pow(M,2)*sin(2*radbeta)*(1/tan(radbeta))*(pow(M,2)*pow(sin(radbeta),2) - 1))/pow(((gamma + cos(2*radbeta))*pow(M,2) + 2),2);
	return dfun;
}
double Newton_a(double M, double gamma, double theta, double radbeta){

	double radbeta0 = 0;
	double f,df;
	do{
		for (int i = 0; i < times; i++){
			radbeta0 = radbeta;
			f = fun(M,gamma,theta,radbeta0);
			df = dfun(M,gamma,radbeta0);
			radbeta = radbeta0 - f/df;
		}
	}while(fabs(radbeta-radbeta0) >= diff);
	return radbeta;
}

void shockwave(const char* q2_file)
{
	double M;
	double theta;
	double beta_l;
	double beta_u;
	double gamma;
	double outputl;
	double outputu;

	/*
	 * task2.3_a
	 */
	FILE *fp;
	fp = fopen(q2_file,"r");
	fseek(fp, 28L, SEEK_SET);
	fscanf(fp,"%lf,%lf,%lf,%lf,%lf",&M,&theta,&beta_l,&beta_u,&gamma);
	/* Find and output beta_l for task2.3_a */
	printf("%.4f\n",Newton_a(M,gamma,theta,beta_l));
	/* Find and output beta_u for task2.3_a*/
	printf("%.4f\n",Newton_a(M,gamma,theta,beta_u));

	/*
	 * task2.3_b,c
	 */
	FILE *fw;
	fw = fopen("out_shock.csv","w");
	fprintf(fw,"M,theta,beta_lower,beta_upper\n");
	fseek(fp, 3L, 1);
	/* Find and output beta_l and beta_u*/
	while(fscanf(fp,"%lf",&M) == 1){
		for(int theta = 0; theta < 90; theta++){
			outputl = Newton_a(M,gamma,theta,beta_l);
			outputu = Newton_a(M,gamma,theta,beta_u);
			if(outputl<90&&outputl>0&&outputu<90&&outputu>0)
			fprintf(fw,"%.4f,%d,%.4f,%.4f\n",M,theta,outputl,outputu);
		}
	}
	fclose(fp);
	fclose(fw);

}
