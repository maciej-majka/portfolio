#ifndef MATRIX_H
#define MATRIX_H
#include<iostream>
using namespace std;

struct list
{
       list(): next(NULL), l_zer(0),dl_bloku(0), tab(NULL) {}
       double *tab;
       list *next;
       int l_zer,dl_bloku;
};

class matrix
{
      public:
             matrix(): blok(NULL), l_elementow(0) {cout<<"c-tor"<<endl;}//c-tor
             matrix(const matrix &A);
             matrix(double (*fcja)(double,double,double,int,int,int,int,double,double,double),double,double,double,int,int,double,double,double,int);
             matrix(matrix &A,int a);
             matrix(double (*fcja)(int, int, matrix &, matrix &), matrix &A, matrix &B);
             ~matrix();
             void wypisz();
             void wypisz_element(int i, int j);
             double element(int i, int j);
             double LT(int i, int j, matrix &A, list *blok, double *bufor);
             void wypisz_struct();
             //double multiply(int i, int j, matrix &A, matrix &B);
             void solve_vec(double *,double *);
             
             void transpose();
             
             matrix& operator=(matrix const &A);
             //matrix& operator=(matrix A);
             
             list* const copy_struct(list const *wsk) const;
             list* const copy_struct() const;
             
             friend matrix operator*(matrix &A, matrix &B);
             
             int l_elementow;
      private:
              list *blok;
              
              void usuwanie(list *ptr);
              list* oproznij_bufor(double *bufor,int k,list *wsk);
};

#endif
