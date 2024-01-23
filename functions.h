#ifndef FUNCTIONS_H
#define FUNCTIONS_H

//zapisy
double zapis(int , int , int, int, double, double, double);

double zapis2(int i, int j, int n);

//fcja sprawdzajaca
void compare(matrix &A, matrix &B);

//czas
void czasomierz();

//mnozenie
double multiply(int i, int j, matrix &A, matrix &B);

double multiply_T(int i, int j, matrix &A, matrix &B);

//globalny operator
matrix operator*(matrix &A, matrix &B);

//WAZNE: fcje korelacji
double corr_function(double,double,double,int i, int j, int, int, double, double, double);

double identity(double,double,double,int i, int j, int, int, double, double, double);

double transpose_el(int, int, matrix&,matrix&);

void simulation();
#endif
