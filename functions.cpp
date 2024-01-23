#include<iostream>
#include<cmath>
#include<ctime>
#include"matrix.h"
#include"set.h"
using namespace std;


double corr_function(double noise_level, double gamma, double lambda, int i, int j, int n_x, int n_y, double delta, double x_min, double y_min)
{
       //if((i>n_x*n_y&&j<=n_x*n_y)||(j>n_x*n_y&&i<=n_x*n_y)) return 0;
       double x1,y1,x2,y2;
       i%=(n_x*n_y);
       j%=(n_x*n_y);
       int k(i),l(j);
       k%=n_x;
       l%=n_x;
       x1=k*delta+x_min;
       x2=l*delta+x_min;
       y1=(i-k)/(double)n_x*delta+y_min;
       y2=(j-l)/(double)n_x*delta+y_min;
       
       double r=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
       double corr=noise_level*fabs(gamma)*exp(-r/lambda);
       //cout<<"i: "<<i<<" x1:"<<x1<<" y1:"<<y1<<"; j:"<<j<<" x2:"<<x2<<" y2:"<<y2<<"; corr: "<<corr<<endl;
       if(corr<0.01) return 0;
       //else
       return corr;
}

double identity(double noise_level, double gamma,double lambda,int i, int j, int n_x, int n_y, double delta, double x_min, double y_min)
{
     if (i==j) return fabs(gamma)*noise_level;
     else return 0;
}

double zapis(int i, int j, int n_x,int n_y,double delta, double x_min,double y_min)
{
       double x(0);
       if(i>j) return 0;
       else 
       if(rand()%11>=6) return 0.363409875675;
       else return 1;
       /*
       if (i==j) 
       {while(x==0)
       x=rand()%9;
       return x/(double)10;}
       else x=rand()%9;
       //if(x>5) return 0;
       //else
       return x/(double)10;*/
}

double zapis2(int i, int j, int n)
{
       if(i>j) return 0;
       else return((i-1)*n+(j-1)+1);
}

void compare(matrix &A, matrix &B)
{
     int c=1,i,j, n=A.l_elementow;
     if(A.l_elementow!=B.l_elementow) {cout<<"blad l. elementow"<< endl; return;}
     for(i=1;i<=n;++i)
     for(j=1;j<=n;++j)
     {if(A.element(i,j)==B.element(i,j)) continue;
     else {cout<<"BLAD i="<<i<<" j="<<j<<endl; c=0; break;}}
     if(c!=0) cout<<"OK"<<endl;
}

void czasomierz()
{
     static int i(0);
     static double last;
     double x;
     if(i==0)
     {++i;
     last=clock();
     return;}
     x=last;
     last=clock();
     cout<<"czas trwania:"<<(last-x)/1000<<"s\n"<<endl;
}

double multiply(int i, int j, matrix &A, matrix &B)
{
       double suma(0);
       int k;
       for(k=1;k<=A.l_elementow;++k)
       suma+=A.element(i,k)*B.element(k,j);
       return suma;
}

double multiply_T(int i, int j, matrix &A, matrix &B)
{
       double suma(0), x1, x2;
       int k;
       for(k=1;k<=A.l_elementow;++k)
       {x1=A.element(k,i);
       if(x1==0) continue;
       x2=B.element(k,j);
       if(x2==0) continue;
       suma+=x1*x2;}
       return suma;
}

matrix operator*(matrix &A, matrix &B)
{
       matrix C(multiply_T,A,B);
       return C;
}

double transpose_el(int i, int j, matrix &A,matrix &B)
{
       return A.element(j,i);
}
