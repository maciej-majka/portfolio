#include<iostream>
#include<cmath>
#include<string>
#include<fstream>
#include<gsl/gsl_histogram.h>
#include"set.h"
using namespace std;

void set:: corr(gsl_histogram *val,gsl_histogram *corr, double **tab,char c)
{
    int i,j;
    double m=0,l,u;
    double s=gsl_histogram_sum(val);

    for(i=0;i<gsl_histogram_bins(val);++i)
    {
        gsl_histogram_get_range(val,i,&l,&u);
        m+=gsl_histogram_get(val,i)*l/((double)s);
    }

    for(i=delta_corr-1;i>0;--i)
    for(j=0;j<polymers[0].beads_num;++j)
    {
        tab[j][i]=tab[j][i-1];
    }
    if(c=='a')
    {
        for(j=1;j<polymers[0].beads_num-1;++j)
        tab[j][0]=polymers[0].wezel[j].psi;
    }
    if(c=='r')
    {
        for(j=1;j<polymers[0].beads_num-1;++j)
        tab[j][0]=polymers[0].wezel[j].d;
    }

    for(i=0;i<delta_corr;++i)
    for(j=1;j<polymers[0].beads_num-1;++j)
    {
        gsl_histogram_accumulate(corr,i,(tab[j][0]-m)*(tab[j][i]-m));
    }

}

void set:: def_lin(double delta_psi)
{
    int lin_len(1),i;
    for(i=1;i<polymers[0].beads_num-1;++i)
    {
        if(fabs(polymers[0].wezel[i].psi-180)<delta_psi)
        lin_len++;
        else
        {
            gsl_histogram_increment(hist_lin,lin_len);
            lin_len=1;
        }
    }
}

void set:: angles_dat(ofstream &plik)
{
    int j;
    for(j=2;j<polymers[0].beads_num-2;++j)
    plik<<polymers[0].wezel[j].psi<<" ";
    plik<<endl;
}

void set:: hist_d_update()
{
    int j;
    for(j=2;j<polymers[0].beads_num-2;++j)
    {
        gsl_histogram_increment(hist_d,polymers[0].wezel[j].d);
    }
}

void set:: hist_normalize(gsl_histogram *hist)
{
    double norm=gsl_histogram_sum(hist);
    gsl_histogram_scale(hist,1/(double)norm);
}

void set:: hist_to_file(gsl_histogram *hist,int hist_bin_num, ofstream& plik)
{
    int i;
    int N=hist_bin_num;
    double x,l;
    for(i=0;i<N;++i)
    {
        gsl_histogram_get_range(hist,i,&l,&x);
        plik<<x<<" "<<gsl_histogram_get(hist,i)<<endl;
    }
}

void set:: hist_psi_update()
{
    int j;
    for(j=2;j<polymers[0].beads_num-2;++j)
    {
        gsl_histogram_increment(hist_psi,polymers[0].wezel[j].psi);
    }
}

int set:: corr_index(int typ,int i, int j)
{
    //cout<<"c_i "<<i<<" "<<j<<endl;
    if(typ==1) return i*(polymers[0].beads_num-2)+j;
    if(typ==2) return i*(polymers[0].beads_num-1)+j;
    if(typ==3) return i*(polymers[0].beads_num-2)+j;
}

double set:: get_corr_element(int typ, int i, int t)
{
       double *tmp;
       char *nazwa;
       if(typ==1) tmp=kink_corr;
       if(typ==2) tmp=rad_corr;
       if(typ==3) tmp=cross_corr;
       return tmp[corr_index(typ,t,i)];
}


void set:: def_kink()
{
     int i,j;
     double theta,new_y;
     //cout<<"kink: ";
     for(i=0;i<polymer_num;++i)
     for(j=1;j<polymers[i].beads_num-1;++j)
     {
                                           if(time!=0&&i==0) last_kink[j-1]=polymers[i].wezel[j].kink;
                  theta=polymers[i].wezel[j].fi; //atan2(polymers[i].wezel[j].d_y,polymers[i].wezel[j].d_x);
                  //cout<<"\ntheta: "<<theta<<endl;
                  new_y=-sin(theta)*(polymers[i].wezel[j+1].d_x+polymers[i].wezel[j].d_x)+cos(theta)*(polymers[i].wezel[j+1].d_y+polymers[i].wezel[j].d_y); //shit, by³o x i y - to duuzy blad
                  if(new_y>0) polymers[i].wezel[j].kink=1;
                  else polymers[i].wezel[j].kink=-1;
                  //cout<<polymers[i].wezel[j].kink<<" ";
     }
     //cout<<endl;
}

double set::find_rad_max()
{
       int i;
       double k1,sigma, epsilon,R0;
       k1=polymers[0].wezel[1].k;
       sigma=polymers[0].wezel[1].sigma;
       epsilon=fabs(polymers[0].wezel[1].epsilon);
       R0=polymers[0].wezel[1].R0;

       double r=pow(2,1/6.0)*sigma,r_max(0),U,Umax(0);
       for(i=0;i<1000;++i)
       {
                          U=epsilon*(pow(sigma/r,12)-pow(sigma/r,6))+k1*(r-R0)*(r-R0)/2.0;
                          //if(i%100==0) cout<<U<<" "<<r<<endl;
                          if(U>Umax)
                          {
                                    r_max=r;
                                    Umax=U;
                          }
                          r+=(R0-sigma)/1000.0;
       }
       cout<<"r_max="<<r_max<<" Umax="<<Umax<<endl;
       return r_max;
}

void set:: def_rad()
{
     int i,j;
     for(i=0;i<polymer_num;++i)
     for(j=1;j<polymers[i].beads_num;++j)
     {
                                         if(time!=0&&i==0) last_rad[j-1]=polymers[i].wezel[j].rad;
                                         if(polymers[i].wezel[j].d>rad_max) polymers[i].wezel[j].rad=1;
                                         else polymers[i].wezel[j].rad=-1;
     }
}

void set:: kink_correlation(double T_END)
{
     int i,j;
     //przesuwam i+1 chwile na miejsce i-tej
     for(j=0;j<polymers[0].beads_num-2;++j)
     for(i=kink_t_range-2;i>=0;--i)
     kink_change[corr_index(1,i+1,j)]=kink_change[corr_index(1,i,j)];

     //wyznaczam nowa zmiane
     //cout<<endl;
     for(j=1;j<polymers[0].beads_num-1;++j)
     kink_change[corr_index(1,0,j-1)]=fabs(polymers[0].wezel[j].kink-last_kink[j-1])/2.0;
     //cout<<endl;

     //wyznaczenie przyczynku korelacji t+h*delta_t, j+k-ty wezel
     int h,k;
     double delta,norm;
     for(k=0;k<polymers[0].beads_num-2;++k)
     for(h=0;h<kink_t_range;++h)
     {
                      delta=0;
                      norm=0;
                      for(j=0;j+k<polymers[0].beads_num-2;++j)
                      {
                                                              delta+=kink_change[corr_index(1,0,j)]*kink_change[corr_index(1,h,j+k)];
                                                              ++norm;
                      }
                      if(norm>0) delta/=norm*T_END;
                      if(h>0&&k>0) delta/=2;
                      kink_corr[corr_index(1,h,k)]+=delta;
     }


     for(k=1;k<polymers[0].beads_num-2;++k)
     for(h=1;h<kink_t_range;++h)
     {
                      delta=0;
                      norm=0;
                      for(j=polymers[0].beads_num-3;j-k>=0;--j)
                      {
                                                              delta+=kink_change[corr_index(1,0,j)]*kink_change[corr_index(1,h,j-k)];
                                                              ++norm;
                      }
                      if(norm>0) delta/=norm*T_END*2;
                      kink_corr[corr_index(1,h,k)]+=delta;
     }

}

void set:: rad_correlation(double T_END)
{
     int i,j;
     //przesuwam i+1 chwile na miejsce i-tej
     for(j=0;j<polymers[0].beads_num-1;++j)
     for(i=rad_t_range-2;i>=0;--i)
     rad_change[corr_index(2,i+1,j)]=rad_change[corr_index(2,i,j)];

     //wyznaczam nowa zmiane
     //cout<<endl;
     for(j=1;j<polymers[0].beads_num;++j)
     rad_change[corr_index(2,0,j-1)]=fabs(polymers[0].wezel[j].rad-last_rad[j-1])/2.0;
     //cout<<endl;

     //wyznaczenie przyczynku korelacji t+h*delta_t, j+k-ty wezel
     int h,k;
     double delta,norm;
     for(k=0;k<polymers[0].beads_num-1;++k)
     for(h=0;h<rad_t_range;++h)
     {
                      delta=0;
                      norm=0;
                      for(j=0;j+k<polymers[0].beads_num-1;++j)
                      {
                                                              delta+=rad_change[corr_index(2,0,j)]*rad_change[corr_index(2,h,j+k)];
                                                              ++norm;
                      }
                      if(norm>0) delta/=norm*T_END;
                      rad_corr[corr_index(2,h,k)]+=delta;
     }

     for(k=1;k<polymers[0].beads_num-1;++k)
     for(h=1;h<rad_t_range;++h)
     {
                      delta=0;
                      norm=0;
                      for(j=polymers[0].beads_num-2;j-k>=0;--j)
                      {
                                                              delta+=rad_change[corr_index(2,0,j)]*rad_change[corr_index(2,h,j-k)];
                                                              ++norm;
                      }
                      if(norm>0) delta/=norm*T_END;
                      rad_corr[corr_index(2,h,k)]+=delta;
     }
     //cout<<"\ntutaj"<<endl;
}

void set:: cross_correlation(double T_END)
{
     int h,k,j;
     double delta,norm;

     for(k=0;k<polymers[0].beads_num-2;++k)
     for(h=0;h<cross_t_range;++h)
     {
                      delta=0;
                      norm=0;
                      for(j=0;j+k<polymers[0].beads_num-2;++j)
                      {
                                                              delta+=kink_change[corr_index(3,0,j)]*rad_change[corr_index(3,h,j+k)];
                                                              ++norm;
                      }
                      if(norm>0) delta/=norm*T_END;
                      cross_corr[corr_index(3,h,k)]+=delta;
     }


     for(k=1;k<polymers[0].beads_num-2;++k)
     for(h=1;h<kink_t_range;++h)
     {
                      delta=0;
                      norm=0;
                      for(j=polymers[0].beads_num-3;j-k>=0;--j)
                      {
                                                              delta+=kink_change[corr_index(3,0,j)]*rad_change[corr_index(3,h,j-k)];
                                                              ++norm;
                      }
                      if(norm>0) delta/=norm*T_END;
                      kink_corr[corr_index(3,h,k)]+=delta;
     }
}


void set::print_kink_correlations()
{
     int h,k;
     cout<<"\nmacierz korelacji przejsc katowych:\n"<<endl;
     for(h=0;h<kink_t_range;++h)
     {
                      cout<<" ";
                      for(k=0;k<polymers[0].beads_num-2;++k)
                      {
                                                            cout<<kink_corr[corr_index(1,h,k)]<<" ";
                      }
                      cout<<endl;
     }
     cout<<endl;
}

void set::print_rad_correlations()
{
     int h,k;
     cout<<"\nmacierz korelacji przejsc promieniowych:\n"<<endl;
     for(h=0;h<rad_t_range;++h)
     {
                      cout<<" ";
                      for(k=0;k<polymers[0].beads_num-1;++k)
                      {
                                                            cout<<rad_corr[corr_index(2,h,k)]<<" ";
                      }
                      cout<<endl;
     }
     cout<<endl;
}

void set::print_cross_correlations()
{
     int h,k;
     cout<<"\nmacierz korelacji przejsc krzyzowych:\n"<<endl;
     for(h=0;h<cross_t_range;++h)
     {
                      cout<<" ";
                      for(k=0;k<polymers[0].beads_num-2;++k)
                      {
                                                            cout<<cross_corr[corr_index(3,h,k)]<<" ";
                      }
                      cout<<endl;
     }
     cout<<endl;
}

void set:: corr_trans_file()
{
     char nazwa[30];
     sprintf(nazwa,"wyniki/korelacje/c%i.dat",NUMBER);
     ofstream plik;

     plik.open(nazwa);

     int i,j;
     for(i=0;i<kink_t_range;++i)
     {
                                for(j=0;j<polymers[0].beads_num-2;++j)
                                plik<<kink_corr[corr_index(1,i,j)]<<" ";
                                plik<<endl;
     }
     plik<<endl;
     /*
     for(i=0;i<rad_t_range;++i)
     {
                                for(j=0;j<polymers[0].beads_num-1;++j)
                                plik<<rad_corr[corr_index(2,i,j)]<<" ";
                                plik<<endl;
     }
     plik<<endl;
     for(i=0;i<cross_t_range;++i)
     {
                                for(j=0;j<polymers[0].beads_num-2;++j)
                                plik<<cross_corr[corr_index(3,i,j)]<<" ";
                                plik<<endl;
     }
     */
     plik.close();
}

double set:: get_corr_range(char *typ,char *wym)
{
       cout<<"range"<<endl;
       if(wym=="czas")
       {
                      cout<<"czas"<<endl;
                      if(typ=="kink") return kink_t_range;
                      if(typ=="rad") return rad_t_range;
                      if(typ=="cross") return cross_t_range;
       }
       if(wym=="wezel")
       {
                       if(typ=="kink") return polymers[0].beads_num-2;
                       if(typ=="rad") return polymers[0].beads_num-1;
                       if(typ=="cross") return polymers[0].beads_num-2;
       }
}

void set:: check_force_correlation(int N)
{
     int n,i,j,k,l,index1,index2;
     double **macierz;
     macierz=new double*[n_x];
     for(i=0;i<n_x;++i)
     macierz[i]=new double[n_y];

     for(i=0;i<n_x;++i)
     for(j=0;j<n_y;++j)
     macierz[i][j]=0;

     double norm(0),delta(0);
     for(n=0;n<N;++n)
     {
                     F_stochastic_v2();
                     for(k=0;k<n_x;++k)
                     for(l=0;l<n_y;++l)
                     {
                                       norm=0;
                                       delta=0;
                                       for(i=0;i+k<n_x;++i)
                                       for(j=0;j+l<n_y;++j)
                                       {
                                                           index1=j*n_x+i;
                                                           index2=(j+l)*n_x+i+k;
                                                           delta+=corr_F[index1]*corr_F[index2];
                                                           ++norm;
                                       }
                                       macierz[k][l]+=delta/(norm*N);
                                       norm=0;
                                       delta=0;
                                       //if(k!=0&&l!=0)
                                       {
                                                     for(i=n_x-1;i-k>=0;--i)
                                                     for(j=0;j+l<n_y;++j)
                                                     {
                                                           index1=j*n_x+i;
                                                           index2=(j+l)*n_x+i-k;
                                                           delta+=corr_F[index1]*corr_F[index2];
                                                           ++norm;
                                                     }
                                                     macierz[k][l]+=delta/(norm*N*set::delta_r);
                                       }
                     }
     }
     ofstream plik;
     plik.open("sily.dat");
     for(i=0;i<n_x;++i)
     {
                       for(j=0;j<n_y;++j)
                       plik<<macierz[i][j]<<" ";
                       plik<<endl;
     }
     plik.close();
     for(i=0;i<n_x;++i)
     delete [] macierz[i];
     delete [] macierz;
}

void set::check_force_dynamics()
{
     int i,j,n_x,n_y;
     double x_i,y_i,x_j,y_j;
     for(i=0;i<polymers[0].beads_num;++i)
     for(j=i;j<polymers[0].beads_num;++j)
     {
                                         x_i=polymers[0].wezel[i].x;
                                         y_i=polymers[0].wezel[i].y;
                                         x_j=polymers[0].wezel[j].x;
                                         y_j=polymers[0].wezel[j].y;
                                         n_x=fabs(x_i-x_j);
                                         n_y=fabs(y_i-y_j);
                                         if(n_x>99||n_y>99) continue;
                                         sily[n_y][n_x]+=polymers[0].wezel[i].S_x*polymers[0].wezel[j].S_x;
                                         ++czestosc[n_y][n_x];
     }
}

void set:: force_dynamic_file()
{
     int i,j;
     ofstream plik;
     plik.open("sily_dyn.dat");
     double norma(0);

     for(i=0;i<100;++i)
     for(j=0;j<100;++j)
     {
                       if(czestosc[i][j]!=0) sily[i][j]/=(double)czestosc[i][j];
     }

     for(i=0;i<100;++i)
     for(j=0;j<100;++j)
     norma+=sily[i][j];
     cout<<"norma="<<norma<<endl;
     for(i=0;i<100;++i)
     for(j=0;j<100;++j)
     sily[i][j]/=norma;

     for(i=0;i<100;++i)
     {
                       for(j=0;j<100;++j)
                       plik<<sily[i][j]<<" ";
                       plik<<endl;
     }
     plik.close();
}

void set::print_noise()
{
      ofstream plik;
      plik.open("szum.dat");
      plik<<"x y Fx Fy"<<endl;

      double Fx[n_x*n_y],Fy[n_x*n_y];
      double alfa=corrLT->element(1,1);

      int i;

      F_stochastic_v2();
      cout<<"tutaj"<<endl;
      for(i=0;i<n_x*n_y;++i)
      {
                            if(typ_korelacji=='i') Fx[i]=alfa*gsl_ran_gaussian(generator,1);
                            else Fx[i]=corr_F[i];

      }

      F_stochastic_v2();
      for(i=0;i<n_x*n_y;++i)
      {
                            if(typ_korelacji=='i') Fy[i]=alfa*gsl_ran_gaussian(generator,1);
                            else Fy[i]=corr_F[i];
      }

      double x,y,x1,y1;
      for(i=0;i<n_x*n_y;++i)
      {
                            x=i%n_x*delta_r;
                            y=(i-i%n_x)/(double)n_y*delta_r;
                            x1=x+Fx[i];
                            y1=y+Fy[i];
                            plik<<x<<" "<<y<<" "<<x1<<" "<<y1<<endl;
      }

      plik.close();
}
