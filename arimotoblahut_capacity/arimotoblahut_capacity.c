#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

//#define ERRCHANNEL

int X,Y;

#define ERRT 1e-4

void malloc2d(double ***grid, int nrows, int ncols){
    int i;
    *grid = (double **)malloc( sizeof(double *) * nrows);
	
    if (*grid == NULL){
        printf("ERROR: out of memory\n");
        exit(1);
    }

    for (i=0;i<nrows;i++){
        (*grid)[i] = (double *)malloc( sizeof(double) * ncols);
        if ((*grid)[i] == NULL){
            printf("ERROR: out of memory\n");
            exit(1);
        }
    }
}

void free2d(double ***grid, int nrows){
	int i;
	for(i=0;i<nrows;i++)
		free((*grid)[i]);
	
	free(*grid);
}

void init_memory(double ***p,double **r,double ***q){
	malloc2d(p,Y,X);
	*r = (double *)malloc( sizeof(double)*X);
	malloc2d(q,X,Y);
}

void free_memory(double ***p,double **r,double ***q){
	free(*r);
	free2d(q,Y);
	free2d(p,Y);
}

void init_data(double ***p,double **r,double ***q,int type){
	int i,j;

	if(type==0){
		X=Y=2;
		init_memory(p,r,q);
		printf("\nZ Channel\n");
		(*p)[0][0] = 1;		(*p)[0][1] = 0.5;
		(*p)[1][0] = 0;		(*p)[1][1] = 0.5;
	}else if(type==1){
		X=Y=3;
		printf("\n3x3 Channel\n");
		init_memory(p,r,q);
		(*p)[0][0] = 2.0/3.0;	(*p)[0][1] = 1.0/3.0;	(*p)[0][2] = 0.0;
		(*p)[1][0] = 1.0/3.0;	(*p)[1][1] = 1.0/3.0;	(*p)[1][2] = 1.0/3.0;
		(*p)[2][0] = 0.0;		(*p)[2][1] = 1.0/3.0;	(*p)[2][2] = 2.0/3.0;
	}

	for(j=0;j<X;j++)
		(*r)[j] = 1.0/X;

	for(i=0;i<X;i++)
		for(j=0;j<Y;j++)
			(*q)[i][j] = 0;

}

double arimotoblahutCapacity(double **p,double *r,double **q){
	int i,j,t;
	double sum;
	double *tempY=(double *)malloc(sizeof(double)*Y);
	double *tempX=(double *)malloc(sizeof(double)*X);

	for(j=0;j<Y;j++)
		for(i=tempY[j]=0;i<X;i++)
			tempY[j]+=r[i]*p[j][i];
	
	for(i=0;i<X;i++)
		for(j=0;j<Y;j++)
			q[i][j] = r[i]*p[j][i]/tempY[j];

	for(i=sum=0;i<X;i++){
		for(j=0,tempX[i]=1;j<Y;j++){
			tempX[i]*= pow(q[i][j],p[j][i]);
		}
		sum+=tempX[i];
	}

	for(i=0;i<X;i++)		
		tempX[i] = tempX[i]/sum;

	for(i=0,sum=0;i<X;i++){
		sum+=pow(tempX[i]-r[i],2.0);
	}

	memcpy(r,tempX,sizeof(double)*X);
	free(tempX);free(tempY);

	return sqrt(sum);
}

double calcCapacity(double **p,double *r,double **q){
	int i,j;
	double C = 0;
	for(i=0;i<X;i++){
		for(j=0;j<Y;j++){
			if(r[i]>0 && q[i][j]>0)
				C+=r[i]*p[j][i]*log(q[i][j]/r[i]);
		}
	}
	return C/log(2.0);
}

int main(int argc,char *argv[]){
	int i,type;
	double **p,*r,**q,error;

	for(type=0;type<2;type++){
		init_data(&p,&r,&q,type);
		for(i=0;i<1000;i++){
			error = arimotoblahutCapacity(p,r,q);
			if(error<ERRT){
				printf("Iter: %3d, Channel capacity: %.5lf\n",i,calcCapacity(p,r,q));
				break;
			}
		}
		free_memory(&p,&r,&q);
	}

	return 1;
}

