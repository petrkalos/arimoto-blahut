#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int X;
int XI;

#define ERRT 1e-6

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

void init_memory(double **p,double **r,double ***q,double ***d){
	*p = (double *)malloc( sizeof(double)*X);
	*r = (double *)malloc( sizeof(double)*XI);
	malloc2d(q,XI,X);
	malloc2d(d,X,XI);
}

void free_memory(double **p,double **r,double ***q,double ***d){
	free(*p);free(*r);
	free2d(q,XI);
	free2d(d,X);
}

void init_data(double **p,double **r,double ***q,double ***d,int type){
	int i,j;

	if(type==0){
		X=XI=2;
		init_memory(p,r,q,d);
		(*p)[0] = 0.25;	(*p)[1]=0.75;
		(*d)[0][0] = 0;	(*d)[0][1]=1;
		(*d)[1][0] = 1;	(*d)[1][1]=0;
	}else if(type==1){
		double e=0.25;
		X=2;XI=3;		
		init_memory(p,r,q,d);
		(*p)[0] = 0.5;	(*p)[1] = 0.5;
		(*d)[0][0] = 0;	(*d)[0][1]=1;	(*d)[0][2]=e;
		(*d)[1][0] = 1;	(*d)[1][1]=0;	(*d)[1][2]=e;
	}else if(type==2){
		X=XI=2;
		init_memory(p,r,q,d);
		(*p)[0] = 0.25;	(*p)[1]=0.75;
		(*d)[0][0] = 0;	(*d)[0][1]=2;
		(*d)[1][0] = 1;	(*d)[1][1]=0;
	}

	for(i=0;i<XI;i++)
		(*r)[i] = 1.0/XI;
	
	for(i=0;i<XI;i++)
		for(j=0;j<X;j++)
			(*q)[i][j] = 0;

}

double calcRate(double *p,double *r,double **q){
	int i,j;
	double R = 0;
	for(i=0;i<XI;i++)
		for(j=0;j<X;j++)
			R+=p[j]*q[i][j]*log(q[i][j]/r[i]);
		
	return R/log(2.0);
}

double calcDist(double *p,double **d,double **q){
	int i,j;
	double C = 0;
	for(i=0;i<XI;i++)
		for(j=0;j<X;j++)
			C+=p[j]*q[i][j]*d[j][i];
		
	return C;
}

double arimotoblahutDistortion(double *p,double *r,double **q,double **d,double lamda){
	int i,j;
	double *temp,*c;
	double err;
	
	temp = (double *)malloc(sizeof(double)*X);
	c = (double *)malloc(sizeof(double)*XI);

	for(j=0;j<X;j++)
		for(i=temp[j]=0;i<XI;i++)
			temp[j] += r[i]*exp(-lamda*d[j][i]);
		
	for(i=0;i<XI;i++)
		for(j=0;j<X;j++)
			q[i][j] = r[i]*exp(-lamda*d[j][i])/temp[j];

	for(i=0;i<XI;i++)
		for(j=r[i]=0;j<X;j++)
			r[i]+=p[j]*q[i][j];

	for(i=0;i<XI;i++)
		for(j=c[i]=0;j<X;j++)
			c[i]+=p[j]*exp(-lamda*d[j][i])/temp[j];

	for(i=err=0;i<XI;i++)
		err+=r[i]*fabs(1-c[i]);

	free(c);free(temp);

	return err;
}

int main(int argc,char *argv[]){
	int i,j,type,count;
	double l;
	double lamda[3]={5.0,2.5,1.25};
	double *p,*r,**q,**d,error;
	
	printf("\nType %d\n",0);
	for(j=0;j<3;j++){
		init_data(&p,&r,&q,&d,0);
		l = lamda[j];

		for(i=0;i<1000;i++){
			error = arimotoblahutDistortion(p,r,q,d,l);
			if(error<=ERRT){
				printf("Iter: %3d Lamda %5.2lf,Rate = %.5lf,Distortion = %.5lf\n",i,l,calcRate(p,r,q),calcDist(p,d,q));
				break;
			}
			
		}
		free_memory(&p,&r,&q,&d);
	}


	for(type=1;type<3;type++){
		printf("\nType %d\n",type);
		count=0;
		for(l=0;l<5;l+=0.001){
			init_data(&p,&r,&q,&d,type);
			for(i=0;i<1000;i++){
				error = arimotoblahutDistortion(p,r,q,d,l);
				if(error<=ERRT){
					break;
				}
			}
			count+=i;
			printf("Iter: %3d Lamda %5.2lf,Rate = %.5lf,Distortion = %.5lf\n",i,l,calcRate(p,r,q),calcDist(p,d,q));
			//printf("%5.2lf\t%.5lf\t%.5lf\n",l,calcRate(p,r,q),calcDist(p,d,q));
			free_memory(&p,&r,&q,&d);
		}
		printf("Average iter %.1lf\n",count/(l/0.001));
	}

	return 1;
}

