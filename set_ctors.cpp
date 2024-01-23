#include<iostream>
#include<string>
#include<cmath>
#include<ctime>
#include<cstdlib>
#include<gsl/gsl_rng.h>
#include"set.h"
#include"plot_anim.h"
#include"matrix.h"
#include"functions.h"
using namespace std;

set:: ~set()
{
      int i,j,k;

      delete [] kink_corr;
      delete [] kink_change;
      delete [] last_kink;

      delete [] rad_corr;
      delete [] rad_change;
      delete [] last_rad;

      delete [] cross_corr;
      /*
      for(i=0;i<100;++i)
      {delete [] sily[i];
      delete [] czestosc[i];}
      delete [] sily;
      delete [] czestosc;
      */

      gsl_rng_free(generator);
      gsl_histogram_free(hist_d);
      gsl_histogram_free(hist_psi);
      gsl_histogram_free(hist_lin);
      gsl_histogram_free(hist_ang);
      gsl_histogram_free(hist_rad);

      for(i=0;i<polymers[0].beads_num;++i)
      {
          delete [] ang[i];
          delete [] rad[i];
      }
      delete [] ang;
      delete [] rad;

      if(polymers==NULL) return;

      for(i=0;i<polymer_num;++i)
      {
                                if(polymers[i].wezel!=NULL) delete[](polymers[i].wezel);
                                //if(back_up[i].wezel!=NULL) delete[](back_up[i].wezel);
      }
      delete[](polymers);
      //delete[] (back_up);

      //if(corr_F!=NULL) delete[] (corr_F);
      //if(uncorr_F!=NULL) delete[] (uncorr_F);

      for(k=0;k<4;++k)
      {
                      for(i=0;i<polymer_num;++i)
                      delete [] Runge_k[k][i].wezel;
                      delete [] Runge_k[k];
      }

      if(correlations!=NULL) delete correlations;
      if(corrLT!=NULL) delete corrLT;

      //cout<<"SET d-tor: koniec"<<endl;

}

set:: set(int N, int (*f_dlugosc)(int), bead (*f_wezel)(polymer*,int,int)): polymer_num(N), polymers(NULL), grafika(false)
{
      int i,j;
      polymers=new polymer [N];
      for(i=0;i<N;++i)
      {
                      polymers[i].beads_num=f_dlugosc(i);
                      if(polymers[i].beads_num!=0) polymers[i].wezel=new bead [polymers[i].beads_num];
      }
      for(i=0;i<N;++i)
      for(j=0;j<polymers[i].beads_num;++j)
                                                 polymers[i].wezel[j]=f_wezel(polymers,i,j); //fcja okreslajaca polozenie poczatkowe
}

set:: set(int liczba, int dlugosc, double R0, double k1, double psi0, double k2, double epsilon, double sigma, double lambda,double gamma, double level, double margin, double delta,char c,string katalog,int seed,int numer,bool graf): NUMBER(numer), time(0), iteration_counter(0), grafika(graf)
{
      //cout<<"SET-ctor"<<endl;
      polymer_num=liczba;
      polymers=new polymer[liczba];
      int i,j;
      katalog.erase(0,3);
      output_dir=katalog;
      double A;

      srand(clock());

      for(i=0;i<liczba;++i)
      {
                          polymers[i].wezel=new bead[dlugosc];
                          polymers[i].beads_num=dlugosc;

                          for(j=0;j<dlugosc;++j)
                          {

                                                /*
                                                double dk;
                                                dk=(double)2/(double)dlugosc;

                                                polymers[i].wezel[j].x=(j-dlugosc/(double)2)*R0;
                                                if(a!=0) A=R0*dlugosc/a*0.5;
                                                else A=1;
                                                polymers[i].wezel[j].y=A*cos(a*3.1415/((double)dlugosc*R0)*polymers[i].wezel[j].x)-A/2+A*i*1.3;
                                                */
                                                polymers[i].wezel[j].k=k1;
                                                polymers[i].wezel[j].k2=k2;
                                                polymers[i].wezel[j].PSI0=psi0;
                                                polymers[i].wezel[j].m=1;
                                                polymers[i].wezel[j].epsilon=epsilon;//145;
                                                polymers[i].wezel[j].R0=R0;
                                                polymers[i].wezel[j].sigma=sigma;
                                                polymers[i].wezel[j].v_x=0;
                                                polymers[i].wezel[j].v_y=0;
                                                polymers[i].wezel[j].L0=polymers[i].wezel[j].R0*sqrt(2*(1-cos(polymers[i].wezel[j].PSI0)));
                                                polymers[i].wezel[j].gamma=gamma;
                                                polymers[i].wezel[j].noise_level=level;
                                                polymers[i].wezel[j].S_x=0;
                                                polymers[i].wezel[j].S_y=0;


                                                if(j==0)
                                                {
                                                    polymers[i].wezel[j].x=0;
                                                    polymers[i].wezel[j].y=0;
                                                }
                                                if(j==1)
                                                {
                                                    double a=(rand()%360)*6.28/360.0;
                                                    polymers[i].wezel[j].x=R0*cos(a);
                                                    polymers[i].wezel[j].y=R0*sin(a);
                                                }

                                                if(j>1)
                                                {
                                                    double a1,b1,a2,b2,c;
                                                    a1=polymers[i].wezel[j-1].x-polymers[i].wezel[j-2].x;
                                                    b1=polymers[i].wezel[j-1].y-polymers[i].wezel[j-2].y;
                                                    a2=((rand()%100)/100.0*2.0-1)*fabs(b1);
                                                    //cout<<a1<<" "<<b1<<" "<<a2<<endl;
                                                    c=sqrt(R0*R0-a2*a2);
                                                    //b2=c;
                                                    if((a1*a2+b1*c)>(a1*a2-b1*c)) b2=c;
                                                    else b2=-c;
                                                    polymers[i].wezel[j].x=polymers[i].wezel[j-1].x+a2;
                                                    polymers[i].wezel[j].y=polymers[i].wezel[j-1].y+b2;
                                                }
                                                //cout<<j<<":"<<polymers[i].wezel[j].x<<" "<<polymers[i].wezel[j].y<<endl;

                                                //po³ozenia random
                                                /*
                                                {
                                                    polymers[i].wezel[j].x=rand()%100-50;
                                                    polymers[i].wezel[j].y=rand()%100-50;
                                                }
                                                */




                          }
      }
      //cout<<"alokowano"<<endl;

      mass_sum();
      CM();
      d_innear();
      psi_innear();
      suma_F();
      f_max_length();

      //stoch_force=new  double [matrix_rank];

      //cout<<"tutaj:"<<endl;
      for(i=0;i<4;++i)
      {
                      Runge_k[i]=new polymer[liczba];
                      for(j=0;j<liczba;++j)
                      Runge_k[i][j]=polymers[j];
      }
      //alokacja pamieci pod wektory backup'u
      //back_up=new polymer[liczba];

      //print_data();
      /*
      find_area(margin,delta);
      cout<<"xmin: "<<x_min<<" xmax:"<<x_max<<endl;
      cout<<"ymin: "<<y_min<<" ymax:"<<y_max<<endl;
      cout<<"matrix rank: "<<matrix_rank<<endl;
      cout<<"margin: "<<margin<<" delta: "<<delta<<endl;

      //skala=a

      uncorr_F=new double[matrix_rank];
      corr_F=new double[matrix_rank];
      for(i=0;i<matrix_rank;++i)
      {
                                uncorr_F[i]=0;
                                corr_F[i]=0;
      }

      for(i=0;i<matrix_rank/2;++i)
      for(j=0;j<matrix_rank/2;++j)
      corr_function(i,j,n_x,n_y,delta,x_min,y_min);
      */

      if(level>0)
      {
                 /*
                 if(c=='c') correlations=new matrix(corr_function,level,gamma,lambda,n_x,n_y,delta_r,x_min,y_min,matrix_rank);
                 if(c=='i') correlations=new matrix(identity,level,gamma,lambda,n_x,n_y,delta_r,x_min,y_min,matrix_rank);
                 if(c=='c') cout<<"korelacje eksponecjalne"<<endl;
                 if(c=='i') cout<<"korelacje jednostkowe"<<endl;
                 */
                 correlations=NULL;
                 corrLT=NULL;
                 //correlations->wypisz();
                 if(correlations!=NULL)
                 {
                                       corrLT=new matrix(*correlations,1);
                                       corrLT->transpose();
                                       //corrLT->wypisz();
                                       delete_matrix(correlations);
                                       correlations=NULL;
                                       //sqrt_detS=matrixL_det();
                                       //cout<<"sqrt_det: "<<sqrt_detS<<endl;
                 }

                 typ_korelacji=c;
      }
      else
      {
          correlations=NULL;
          corrLT=NULL;
          cout<<"szum wylaczony"<<endl;
      }

      length=lambda;
      level_kT=level;
      gamma_coef=gamma;

      generator=gsl_rng_alloc(gsl_rng_mt19937);
      gsl_rng_set(generator,seed);

      //ten kod sluzy do dynamicznego sprawdzania korelacji sil
      /*
      sily=new double* [100];
      czestosc=new int *[100];
      for(i=0;i<100;++i)
      {sily[i]=new double[100];
      czestosc[i]=new int[100];}

      for(i=0;i<100;++i)
      for(j=0;j<100;++j)
      {sily[i][j]=0;
      czestosc[i][j]=0;}
      */


      forces();
      /*
      for(i=0;i<matrix_rank;++i)
      cout<<"i: "<<i<<" unc:"<<uncorr_F[i]<<" c:"<<corr_F[i]<<endl;
      */

      //operacje dla kink
      kink_t_range=20;
      last_kink=new int[polymers[0].beads_num-2];
      int wym=(polymers[0].beads_num-2)*kink_t_range;
      cout<<"wym:"<<wym<<endl;
      kink_change=new int[wym];
      kink_corr=new double[(polymers[0].beads_num-2)*kink_t_range];

      //zerowanie macierzy korelacji kink i macierzy zmiany kink
      for(i=0;i<(polymers[0].beads_num-2)*kink_t_range;++i)
      {
                                            kink_corr[i]=0;
                                            kink_change[i]=0;
      }

      //operacje dla rad
      rad_max=find_rad_max();

      rad_t_range=20;
      last_rad=new int[polymers[0].beads_num-1];

      rad_change=new int[(polymers[0].beads_num-1)*rad_t_range];
      rad_corr=new double[(polymers[0].beads_num-1)*rad_t_range];

      //zerowanie macierzy korelacji kink i macierzy zmiany kink
      for(i=0;i<(polymers[0].beads_num-1)*rad_t_range;++i)
      {
                                                    rad_corr[i]=0;
                                                    rad_change[i]=0;
      }
      //operacje dla cross
      if(kink_t_range>rad_t_range) cross_t_range=kink_t_range;
      else cross_t_range=rad_t_range;

      cross_corr=new double[(polymers[0].beads_num-2)*cross_t_range];

      for(i=0;i<(polymers[0].beads_num-2)*cross_t_range;++i)
      {
                                            cross_corr[i]=0;
      }

      hist_d_bin_num=100;
      hist_d=gsl_histogram_alloc(hist_d_bin_num);
      gsl_histogram_set_ranges_uniform(hist_d,0,3*R0);
      hist_psi_bin_num=360;
      hist_psi=gsl_histogram_alloc(hist_psi_bin_num);
      gsl_histogram_set_ranges_uniform(hist_psi,0,360);

      hist_lin=gsl_histogram_alloc(50);
      gsl_histogram_set_ranges_uniform(hist_lin,1,51);

      delta_corr=100;
      ang=new double*[dlugosc];
      rad=new double*[dlugosc];

      for(i=0;i<dlugosc;++i)
      {
          ang[i]=new double[delta_corr];
          rad[i]=new double[delta_corr];
      }

      for(i=0;i<delta_corr;++i)
      for(j=0;j<dlugosc;++j)
      {
          ang[j][i]=0;
          rad[j][i]=0;
      }

      hist_ang=gsl_histogram_alloc(delta_corr);
      gsl_histogram_set_ranges_uniform(hist_ang,0,delta_corr);

      hist_rad=gsl_histogram_alloc(delta_corr);
      gsl_histogram_set_ranges_uniform(hist_rad,0,delta_corr);

      //cout<<"SET c-tor koniec"<<endl;
}
