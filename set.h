#ifndef SET_H
#define SET_H

#include<iostream>
#include<fstream>
#include<string>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include<gsl/gsl_histogram.h>
#include"plot_anim.h"
#include"matrix.h"

using namespace std;

struct bead
{
       bead(): x(0),y(0),fi(0),k(1), k2(0.25),R0(2),L0(3),sigma(1), epsilon(1), m(1), gamma(0),noise_level(0){}
       double m; //masa wezla
       double x,y; // x,y -polozenie, fi -kat nachylenia do wzgledem poprzedniego
       double v_x, v_y;

       double x_o,y_o;
       double fi; // zmienne wewnêtrzne
       double psi;
       double d, d_x, d_y;
       double suma_fx, suma_fy; //sumy sk³adowych x i y si³ dzia³ajacych na wêze³, bez stochastycznych
       double sumaLJ_fx_z, sumaLJ_fy_z, S_x, S_y; //sily zewnetrzne S-stochastyczna, LJ z innymi lancuchami
       double sumaR2_fx, sumaR2_fy, sumaR_fx, sumaR_fy, sumaLJ_fx_w, sumaLJ_fy_w, sumaPSI_fx, sumaPSI_fy;//sily KAT -oddz. katowe, R -rozciaganie, LJ -w samym lancuchu
       double sumaT_fx, sumaT_fy; //sily tarcia
       double gamma;
       double noise_level;
       double epsilon; //³adunek w oddzia³ywaniu LJ
       double k; //stala sprezystosci dla segmentu, opisuje odcinek i - i-1
       double k2; //stala sprezystosci oddzialywan II rzedu miedzy i-tym oraz i-2
       double PSI0;
       double L0; //rownowagowa odleglosc sasiadow II rzedu miedzy i-tym oraz i-2
       double R0; //rownowagowa odleglosc wezlow miedzy i - i-1
       double sigma; //stala do potencjalu LJ

       int kink;
       int rad;

};

struct polymer
{
       polymer(): wezel(NULL), beads_num(0),suma_Fx(0),suma_Fy(0),sumaS_x(0), sumaS_y(0) {}
       double Xcm,Ycm; //wsp. srodka masy;
       double mass; //calkowita masa
       int beads_num; //liczba monomerow
       bead *wezel; //tablica opisujaca wlasnosci wezla
       double suma_Fx,suma_Fy,sumaS_x, sumaS_y;

       polymer& operator=(const polymer&);
};

struct check_data
{
       bool check;
       double h; //krok calkowania
       double new_h; //zmieniony krok calkowania
       double t0; //czas t0 dla ktorego wykonano ostatni backup
       double t_end; //czas po ktorym nastapi powrot do wlasciwej dlugosci kroku h
       double jumpE; //maksymalny skok energii
       double last_energy; //energia ukladu w chwili t0
};

class set
{
      public:
             int NUMBER; //nr porzadkowy symulacji
             string output_dir;
             set(int N, int (*f_dlugosc)(int), bead (*f_wezel)(polymer*,int,int));
             ~set();
             set(int,int,double,double,double,double,double,double,double,double,double,double,double,char,string,int,int,bool);
             void solution(double,double,char);
             ofstream alldata;
             ofstream current_file;
             void file_operation();
             void file_operation2();
             void print_data();
             bool grafika;
             //void file_tmp(string);

             void velocity_dependence(ofstream&);
             void reduced_dynamics(ofstream&);

             double corr_fun(int,int);

             void delete_matrix(matrix*);
             double skala;

             void print_kink_correlations();
             void print_rad_correlations();
             void print_cross_correlations();
             void corr_trans_file();
             double get_corr_element(int,int,int);
             double get_corr_range(char*,char*);
             //sprawdzanie korelacji si³
              void check_force_correlation(int N);
              void force_dynamic_file();

              void print_noise();

              void set_send(FILE*);
              void anim_start(FILE*);


      private:
              double length;
              double level_kT;
              double gamma_coef;

              int iteration_counter;
              double time; //czas
              int polymer_num; //liczba polimerow
              polymer *polymers; //dane o polimerach, na ktorych przeprowadzane sa aktualne obliczenia
              polymer *Runge_k[4]; //pomocnicze wektory do metody Rungego
              polymer *back_up;  //kopia ukladu dla funkcji back_up
              int max_length; //najwiêksza d³ugoœæ ³añcucha

              check_data control; //struktura z parametrami kotrolnymi dla fcji back_up

              char typ_korelacji;
              double sqrt_detS; //pierwiastek z wyznacznika macierzy S
              double matrixL_det(); //fcja oblicza wyznacznik macierzy L
              matrix *correlations; //macierz korelacji
              matrix *corrLT; //macierz trojkatna, pochodzaca z rozkaldu cholesky'ego macierzy korelacji
              double x_max, x_min, y_max, y_min, delta_r; //dane o wielkosci obszaru, w ktorym wprowadzone sa korelacje
              int n_x,n_y;  //podzial dyskretny obszaru
              double *corr_F, *uncorr_F; //tablica przechowujaca wektor skladowych sil stochastycznych skorelowanych i nieskorelowanych
              void find_area(double,double); //fcja okresla obszar siatki na ktorej generowane sa sily stochastyczne
              int matrix_rank; //wymiar macierzy korelcji

              //generator
              gsl_rng *generator;

              double ij_to_xy(int,char); //fcja przelicza indeksy punktu ij na jego polozenie

              //energia
              double U(int,int);
              double U_LJ();
              double U_k();
              double U_psi();
              double Ek();
              double energy();

              //skladowe sil
              void F_LJ(); //fcja oblicza przyczynek od oddzia³ywania LJ
              void F_range(); //fcja oblicza przyczynek od sil sprezystosci
              void F_second(); //fcja oblicza przyczynek od oddzia³ywan z sasiadami II rzedu
              void F_stochastic(); //NIE DAJE EFEKTU! fcja oblicza przyczynek od sil stochastycznych, generujac losowy kat i skorelwoane wartosci sil
              void F_stochastic_v2();// fcja oblicza skorelowane skladowe x i y sil stochastycznych
              void F_stochastic_dynamic(); //dziala TYLKO dla 1 lancucha! dynamiczna zmiana macierzy korelacji przy kazdym kroku
              //void F_psi();//NIEU¯YWANE fcja oblicza sily od potencjalu katowego opartego o kat
              void F_T(); //Przyczynek od sil tarcia
              //funkcje do ca³kowania
              void forces(); //fcja wykonuje wszystkie operacje na silach, wymagane w kroku calkowania
              void suma_F(); //fcja przeprowadzajaca sumowanie odpowiednich skladowych sil we wszystkich wezlach
              void mass_sum();
              void d_innear(); //fcja oblicza odleglosci miedzy wezlami
              void psi_innear(); //fcja oblicza katy psi miedzy wezlami (nie rozroznia kierunku zgiecia)

              void CM();
              void relaxation(double multiplayer,int iterations,double h);
              void integration_Verlet(double h);
              void integration_Runge(double h,double);
              void back_up_control(double,double*,double);
              void fi_to_psi();
              void psi_to_fi();

              //histogramy obsługa
              void hist_normalize(gsl_histogram*);
              void hist_to_file(gsl_histogram*, int, ofstream&);
              //histogram odległości
              gsl_histogram *hist_d;
              int hist_d_bin_num;
              void hist_d_update();
              //histogram kątów
              gsl_histogram *hist_psi;
              int hist_psi_bin_num;
              void hist_psi_update();
              //analiza kątów
              void angles_dat(ofstream&);
              //analiza korelacji
              int delta_corr;
              void corr(gsl_histogram*,gsl_histogram*,double**,char);
              //korelacja katow
              double **ang;
              gsl_histogram *hist_ang;

              //korelacja dlugosci segmentow
              double **rad;
              gsl_histogram *hist_rad;

              //analiza liniowości
              gsl_histogram *hist_lin;
              void def_lin(double);

              //funkcje inne
              void f_max_length();
              void cartesian_to_innear();

             //badanie korelacji przejsc miedzy minimami katowymi
              void def_kink(); //fcja rozpoznaje w ktora strone jest zgiety polimer
              int *last_kink;
              double *kink_corr;
              int *kink_change;
              int kink_t_range;
              void kink_correlation(double);
              //badanie korelacji przejsc miedzy minimami potencja³u kwadratowy-LJ
              double rad_max;
              double find_rad_max();
              void def_rad();
              int *last_rad;
              double *rad_corr;
              int *rad_change;
              int rad_t_range;
              void rad_correlation(double);

              //badanie korelacji krzy¿owych miedzy przejsciami k¹towymi i promieniowymi
              int cross_t_range;
              double *cross_corr;
              void cross_correlation(double);

              int corr_index(int,int,int);

              double **sily;
              int **czestosc;
              void check_force_dynamics();
};

#endif
