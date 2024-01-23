#include<iostream>
#include<cstdlib>
#include<cmath>
#include "matrix.h"

using namespace std;

//destruktor
matrix:: ~matrix()
{
         //cout<<"matrix d-tor"<<endl;
         if (blok!=NULL) usuwanie(blok);
}

//konstruktor kopiujacy
matrix:: matrix(const matrix &A): l_elementow(A.l_elementow), blok(NULL)
{
         if(this==&A) return;
         blok=A.copy_struct();
         l_elementow=A.l_elementow;
}

//konstruktor macierzy fcj zwraca element (i,j)
matrix:: matrix(double (*fcja)(double,double,double,int,int,int,int,double,double,double),double noise_level,double gamma,double lambda,int n_x,int n_y,double delta,double x_min, double y_min,int n): l_elementow(n), blok(NULL)
{        
         //cout<<"pocz. c-tor"<<endl;
         int i,j;
         int k(0), N(260000);
         double bufor[260000];
         double element;
         list *wsk;
         
         wsk=new list();
         blok=wsk;
         for(i=1;i<=n;++i)
         for(j=1;j<=n;++j)
         {
                         if(k==N-1)
                         {
                                   wsk=oproznij_bufor(bufor,k,wsk);
                                   wsk->next=new list(); //alokacja next
                                   wsk=wsk->next; //wsk wskazuje na next
                                   k=0;         
                         }
                         element=fcja(noise_level,gamma,lambda,i,j,n_x,n_y,delta,x_min, y_min);
                         if(element==0&&k==0) {++wsk->l_zer; continue;}//a co jak N-1?
                         if(element==0&&k!=0)
                         {
                                       if(bufor[k-1]!=0) {bufor[k]=element; ++k; continue;}
                                       else 
                                       {
                                            wsk=oproznij_bufor(bufor,k,wsk);
                                            wsk->next=new list(); //alokacja next
                                            wsk=wsk->next; //wsk wskazuje na next
                                            ++wsk->l_zer;
                                            k=0;
                                       };
                                       continue;
                         }
                         if(element!=0) {bufor[k]=element; ++k; continue;}
         }
         if(i==n+1||j==n+1) wsk=oproznij_bufor(bufor,k,wsk);//inkrementacja ostatnia:p +1
}

//konstruktor budujacy z dwoch macierzy A i B (do iloczynu)
matrix:: matrix(double (*fcja)(int, int, matrix&, matrix&),matrix &A, matrix &B): blok(NULL)
{        
         //cout<<"pocz. AB-ctor"<<endl;
         if (A.l_elementow!=B.l_elementow) {cout<<"inne rozmiary"<<endl; return;}
         else l_elementow=A.l_elementow;
         int i,j,n=l_elementow;
         int k(0), N(260000);
         double bufor[260000];
         double element;
         list *wsk;
         
         wsk=new list();
         blok=wsk;
         for(i=1;i<=n;++i)
         for(j=1;j<=n;++j)
         {
                         if(k==N-1)
                         {
                                   wsk=oproznij_bufor(bufor,k,wsk);
                                   wsk->next=new list(); //alokacja next
                                   wsk=wsk->next; //wsk wskazuje na next
                                   k=0;         
                         }
                         element=fcja(i,j,A,B);
                         if(element==0&&k==0) {++wsk->l_zer; continue;}//a co jak N-1?
                         if(element==0&&k!=0)
                         {
                                       if(bufor[k-1]!=0) {bufor[k]=element; ++k; continue;}
                                       else 
                                       {
                                            wsk=oproznij_bufor(bufor,k,wsk);
                                            wsk->next=new list(); //alokacja next
                                            wsk=wsk->next; //wsk wskazuje na next
                                            ++wsk->l_zer;
                                            k=0;
                                       };
                                       continue;
                         }
                         if(element!=0) {bufor[k]=element; ++k; continue;}
         }
         if(i==n+1||j==n+1) wsk=oproznij_bufor(bufor,k,wsk);//hehe inkrementacja ostatnia:p +1
}

//kostruktor do rozk³adu cholesky'ego
matrix:: matrix(matrix &A,int a): l_elementow(A.l_elementow), blok(NULL)
{
         //cout<<"pocz. L-ctor"<<endl;
         int i,j;
         int k(0), N(260000),n(l_elementow);
         double bufor[260000];
         double element;
         list *wsk;
         
         wsk=new list();
         blok=wsk;
         for(i=1;i<=n;++i)
         for(j=1;j<=n;++j)
         {
                         if(k==N-1)
                         {
                                   wsk=oproznij_bufor(bufor,k,wsk);
                                   wsk->next=new list(); //alokacja next
                                   wsk=wsk->next; //wsk wskazuje na next
                                   k=0;         
                         }
                         if (i>j) element=0;  //poniewaz macierz jest trojkatna
                         else element=LT(i,j,A,blok,bufor);
                         
                         //cout<<"element: "<<element<<endl;
                         
                         if(element==0&&k==0)
                         {++wsk->l_zer;
                         continue;}//a co jak N-1?
                         if(element==0&&k!=0)
                         {
                                       if(bufor[k-1]!=0) {bufor[k]=element; ++k; continue;}
                                       else 
                                       {
                                            wsk=oproznij_bufor(bufor,k,wsk);
                                            wsk->next=new list(); //alokacja next
                                            wsk=wsk->next; //wsk wskazuje na next
                                            ++wsk->l_zer;
                                            k=0;
                                       };
                                       continue;
                         }
                         if(element!=0) {bufor[k]=element; ++k; continue;}
         }
         if(i==n+1||j==n+1) wsk=oproznij_bufor(bufor,k,wsk);//inkrementacja ostatnia:p +1
}
