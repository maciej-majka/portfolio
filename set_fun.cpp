#include<iostream>
#include<cmath>
#include<string>
#include<fstream>
#include<gsl/gsl_math.h>
#include"set.h"
using namespace std;

void set:: mass_sum()
{
     //cout<<"mass sum"<<endl;
      int i,j;
      for(i=0;i<polymer_num;++i)
      {
                                polymers[i].mass=0;
      for(j=0;j<polymers[i].beads_num;++j)
      {
                                          polymers[i].mass+=polymers[i].wezel[j].m;
      }
      }
}

void set:: CM()
{
     //cout<<"CM"<<endl;
      int j,i;
      double suma_x(0), suma_y(0),suma_m(0);
      for(i=0;i<polymer_num;++i)
      {
                                suma_x=0;
                                suma_y=0;
                                suma_m=0;
                                for(j=0;j<polymers[i].beads_num;++j)
                                {
                                                                    suma_x+=polymers[i].wezel[j].m*polymers[i].wezel[j].x;
                                                                    suma_y+=polymers[i].wezel[j].m*polymers[i].wezel[j].y;
                                                                    suma_m+=polymers[i].wezel[j].m;
                                }
                                polymers[i].Xcm=suma_x/suma_m;
                                polymers[i].Ycm=suma_y/suma_m;
      }
}

void set:: integration_Verlet(double h)
{
     int i,j;

     forces(); //obliczenie sil w chwili iteracja*h

     for(i=0;i<polymer_num;++i)
     for(j=0;j<polymers[i].beads_num;++j)
     {
                                         polymers[i].wezel[j].v_x+=(h*polymers[i].wezel[j].suma_fx+sqrt(h)*polymers[i].wezel[j].S_x)/polymers[i].wezel[j].m;
                                         polymers[i].wezel[j].v_y+=(h*polymers[i].wezel[j].suma_fy+sqrt(h)*polymers[i].wezel[j].S_y)/polymers[i].wezel[j].m;

                                         polymers[i].wezel[j].x+=h*polymers[i].wezel[j].v_x;
                                         polymers[i].wezel[j].y+=h*polymers[i].wezel[j].v_y;

     }

     CM();
     d_innear();
     psi_innear();
}

void set:: integration_Runge(double h,double gamma)
{
     //cout<<"Runge poczatek"<<endl;
     int i,j,k;
     double a,s;

     forces(); //obliczenie sil w chwili iteracja*h

     for(i=0;i<polymer_num;++i) //kopiowanie do wektora posredniego k0
     Runge_k[0][i]=polymers[i];


     for(k=1;k<4;++k)
     {
                     //cout<<"k"<<k<<endl;
                     if(k==1||k==2) a=0.5;
                     else a=1;
                     //cout<<"a="<<a<<endl;
                     for(i=0;i<polymer_num;++i)
                     for(j=0;j<polymers[i].beads_num;++j)
                     {
                                                         polymers[i].wezel[j].v_x+=a*(h*Runge_k[k-1][i].wezel[j].suma_fx+sqrt(h)*polymers[i].wezel[j].S_x)/polymers[i].wezel[j].m;
                                                         polymers[i].wezel[j].v_y+=a*(h*Runge_k[k-1][i].wezel[j].suma_fy+sqrt(h)*polymers[i].wezel[j].S_y)/polymers[i].wezel[j].m;
                                                         polymers[i].wezel[j].x+=a*h*Runge_k[k-1][i].wezel[j].v_x;
                                                         polymers[i].wezel[j].y+=a*h*Runge_k[k-1][i].wezel[j].v_y;
                     }
                     d_innear();
                     //psi_innear();
                     forces();
                     //print_data();
                     for(i=0;i<polymer_num;++i)
                     {
                                               Runge_k[k][i]=polymers[i];
                                               polymers[i]=Runge_k[0][i];
                     }

     }

     for(i=0;i<polymer_num;++i)
     for(j=0;j<polymers[i].beads_num;++j)
     for(k=0;k<4;++k)
     {
                     //cout<<"koncowa suma k="<<k<<"suma_fx:"<<Runge_k[k][i].wezel[j].suma_fx<<endl;
                     if(k==0||k==3) s=1;
                     else s=2;
                     polymers[i].wezel[j].v_x+=((double)1/(double)6)*s*(h*Runge_k[k][i].wezel[j].suma_fx+sqrt(h)*Runge_k[k][i].wezel[j].S_x);
                     polymers[i].wezel[j].v_y+=((double)1/(double)6)*s*(h*Runge_k[k][i].wezel[j].suma_fy+sqrt(h)*Runge_k[k][i].wezel[j].S_y);
                     polymers[i].wezel[j].x+=((double)1/(double)6)*h*s*Runge_k[k][i].wezel[j].v_x;
                     polymers[i].wezel[j].y+=((double)1/(double)6)*h*s*Runge_k[k][i].wezel[j].v_y;
     }

     CM();
     d_innear();
     //psi_innear();
}


void set:: solution(double h,double T_END, char typ)
{
     bool running=true;
     cout<<"rozwiazanie metoda";
     if(typ=='V') cout<<" Verleta"<<endl;
     if(typ=='R') cout<<" Rungego-Kutty 4 rzedu"<<endl;

     int i;
     double energia(1), Ekin,ULJ,UR,Upsi,suma,E1=1,E2=2,dE=1;

     char nazwa[40];
     ofstream energia_plik,zaleznosc_plik,red_dyn,plik_CM,plik_hist_d,plik_hist_psi,plik_ang,plik_hist_lin,plik_ang_corr,plik_rad_corr;
     sprintf(nazwa,"wyniki/energia/e%i.dat",NUMBER);
     energia_plik.open(nazwa);
     energia_plik<<"time "<<"E "<<"E_kin "<<"E_p "<<"ULJ "<<"UR "<<"Upsi "<<endl;

     sprintf(nazwa,"wyniki/synch/s%i.dat",NUMBER);
     zaleznosc_plik.open(nazwa);
     zaleznosc_plik<<"time ";
     for(i=1;i<10;++i)
     zaleznosc_plik<<"v"<<i<<" ";
     zaleznosc_plik<<endl;

     sprintf(nazwa,"wyniki/red_dyn/r%i.dat",NUMBER);
     red_dyn.open(nazwa);
     red_dyn<<"time "<<"av_kink "<<"d_kink "<<"av_rad "<<"d_rad"<<endl;

     sprintf(nazwa,"wyniki/CM/CM%i.dat",NUMBER);
     plik_CM.open(nazwa);
     plik_CM<<"time "<<"X_CM "<<"Y_CM "<<"R2_CM"<<endl;

     sprintf(nazwa,"wyniki/hist_d/hd%i.dat",NUMBER);
     plik_hist_d.open(nazwa);

     sprintf(nazwa,"wyniki/hist_psi/hp%i.dat",NUMBER);
     plik_hist_psi.open(nazwa);

     sprintf(nazwa,"wyniki/angles/ang%i.dat",NUMBER);
     plik_ang.open(nazwa);

     sprintf(nazwa,"wyniki/lines/lin%i.dat",NUMBER);
     plik_hist_lin.open(nazwa);

     sprintf(nazwa,"wyniki/ang_corr/ac%i.dat",NUMBER);
     plik_ang_corr.open(nazwa);

     sprintf(nazwa,"wyniki/rad_corr/rc%i.dat",NUMBER);
     plik_rad_corr.open(nazwa);

     mass_sum();
     CM();
     forces();
     energia=energy();

     control.h=h;
     control.t0=0;
     control.last_energy=energia;
     control.check=false;
     control.jumpE=20;

     //grafika.start(polymer_num);
     //file_tmp(grafika.file_path());
     FILE *anim;

     if(grafika==true)
     {
        anim=popen("f:////doktorat////magisterska////animation////animation.exe","w");
        if(anim==NULL) cout<<"nie znaleziono animation.exe"<<endl;
     }

     if(grafika) anim_start(anim);

     cout<<"T_END="<<T_END<<endl;
     while(running)
     {
                   //relaksacja, stan poczatkowy
                   relaxation(10,640,h);
                   if(time>T_END) running=false;
                   /*
                   if(energia<1E-9) running=false;
                   if(time>T_END&&energia<1) running=false;
                   if(time>T_END&&energia>=1) if(dE<1E-9) running=false;
                   */
                     if(iteration_counter%128==0) cout<<"\rtime = "<<time<<" h ="<<h<<" E="<<energia<<" dE="<<dE<<"     ";

                     ++iteration_counter;

                     if(typ=='V') integration_Verlet(h);
                     if(typ=='R') integration_Runge(h,0);

                     if(iteration_counter%128==0)
                     {
                        E1=E2;
                        energia=energy();
                        E2=energia;
                        dE=fabs(E2-E1);
                        energia_plik<<time<<" "<<energia<<endl;
                     }


                     //back_up_control(energia,&h,control.h*4);
                     //if(h<1E-5) break;

                     time+=h;

                     /*
                     if(iteration_counter%128==0)
                     {
                                                 //cout<<"zapis energii"<<endl;
                          Ekin=Ek();
                          ULJ=U_LJ();
                          UR=U_k();
                          Upsi=U_psi();
                          suma=Ekin+ULJ+UR+Upsi;
                          energia_plik<<time<<" "<<suma<<" "<<Ekin<<" "<<suma-Ekin<<" "<<ULJ<<" "<<UR<<" "<<Upsi<<endl;
                     }*/


                     if(iteration_counter%128==0)
                     {
                                                 double R2=polymers[0].Xcm*polymers[0].Xcm+polymers[0].Ycm*polymers[0].Ycm;
                                                 plik_CM<<time<<" "<<polymers[0].Xcm<<" "<<polymers[0].Ycm<<" "<<R2<<endl;;
                     }

                     if(iteration_counter%32==0) //katy potrzebne jedynie do statystyk, sily od nich nie zaleza
                     psi_innear();

                     /*
                     if(iteration_counter%128==0)
                     {
                                                 //def_kink();
                                                 //def_rad();
                     }
                     */
                     if(time>99&&iteration_counter%128==0)
                     {
                                                         corr(hist_psi,hist_ang,ang,'a');
                                                         corr(hist_d,hist_rad,rad,'r');
                                                         //kink_correlation(T_END);
                                                         //rad_correlation(T_END);
                                                         //cross_correlation(T_END);
                     }
                     if(iteration_counter%128==0)
                     {
                                                 //cout<<"zapis synch i red"<<endl;
                                                velocity_dependence(zaleznosc_plik);
                                                //reduced_dynamics(red_dyn);
                                                //angles_dat(plik_ang);
                     }
                     if(iteration_counter%32==0&&(time>99))
                     {
                                                hist_d_update();
                                                hist_psi_update();
                                                def_lin(20);
                     }
                     /*
                     file_tmp(grafika.file_path());
                     if(iteration_counter%32==0)
                     {
                     grafika.waiting(0.0001);
                     grafika.refresh();
                     grafika.waiting(0.0001);
                     }
                     */
                     if(iteration_counter%32==0&&grafika==true)
                     {
                                                set_send(anim);
                     }
     }
     cout<<endl;
     file_operation2();
     //grafika.end();
     energia_plik.close();
     zaleznosc_plik.close();
     plik_CM.close();

     hist_normalize(hist_d);
     hist_normalize(hist_psi);
     hist_normalize(hist_lin);
     hist_normalize(hist_ang);

     double norm=1/gsl_histogram_get(hist_ang,0);
     gsl_histogram_scale(hist_ang,norm);

     norm=1/gsl_histogram_get(hist_rad,0);
     gsl_histogram_scale(hist_rad,norm);

     hist_to_file(hist_d,hist_d_bin_num,plik_hist_d);
     hist_to_file(hist_psi,hist_psi_bin_num,plik_hist_psi);
     hist_to_file(hist_lin,50,plik_hist_lin);
     hist_to_file(hist_ang,delta_corr,plik_ang_corr);
     hist_to_file(hist_rad,delta_corr,plik_rad_corr);

     plik_hist_d.close();
     plik_hist_psi.close();

     plik_ang.close();
     plik_ang_corr.close();
     plik_rad_corr.close();

     if(grafika==true) pclose(anim);
}

void set:: back_up_control(double energia,double *h,double delta_t)
{
     int i;
     if(fabs(energia-control.last_energy)>control.jumpE)
     {
                                                         if(control.check==false)
                                                         {
                                                                                 control.t_end=time+delta_t;
                                                                                 control.check=true;
                                                         }
                                                         for(i=0;i<polymer_num;++i)
                                                         polymers[i]=back_up[i];
                                                         time=control.t0;
                                                         control.new_h=(*h)/((double)8);
                                                         *h=control.new_h;
                                                         cout<<"\n !utrata stabilnosc! zmniejszenie kroku t0="<<time<<" h="<<*h<<endl;
                                                         return;
     }
     if(control.check==true)
     {
                            if(time>=control.t_end)
                            {
                            *h=control.h;
                            control.check=false;
                            cout<<"\nodzyskano stabilnosc, h="<<control.h<<endl;
                            }
                            return;
     }
     if(control.check==false)
     {
                                                   if((control.t0+delta_t)<=time)
                                                   {
                                                                              for(i=0;i<polymer_num;++i)
                                                                              back_up[i]=polymers[i];
                                                                              control.t0+=delta_t;
                                                                              control.last_energy=energia;
                                                                              //cout<<"\nzapisano back up"<<endl;
                                                                              return;
                                                   }
                                                   else return;
     }
}


void set:: d_innear()
{
     int i,j;
     for(i=0;i<polymer_num;++i)
     {
                               polymers[i].wezel[0].d_x=polymers[i].wezel[0].x-polymers[i].Xcm;
                               polymers[i].wezel[0].d_y=polymers[i].wezel[0].y-polymers[i].Ycm;
                               polymers[i].wezel[0].d=sqrt(polymers[i].wezel[0].d_x*polymers[i].wezel[0].d_x+polymers[i].wezel[0].d_y*polymers[i].wezel[0].d_y);
                               for(j=1;j<polymers[i].beads_num;++j)
                               {
                                                                   polymers[i].wezel[j].d_x=polymers[i].wezel[j].x-polymers[i].wezel[j-1].x;
                                                                   polymers[i].wezel[j].d_y=polymers[i].wezel[j].y-polymers[i].wezel[j-1].y;
                                                                   polymers[i].wezel[j].d=sqrt(polymers[i].wezel[j].d_x*polymers[i].wezel[j].d_x+polymers[i].wezel[j].d_y*polymers[i].wezel[j].d_y);
                               }
     }
}

void set:: psi_innear()
{
     int i,j;
     double deg=360.0/(2.0*M_PI);
     for(i=0;i<polymer_num;++i)
     {
            polymers[i].wezel[1].fi=atan2(polymers[i].wezel[1].d_x,polymers[i].wezel[1].d_y);
            for(j=2;j<polymers[i].beads_num;++j)
            {
                      //polymers[i].wezel[j].fi=atan2(polymers[i].wezel[j].d_x,polymers[i].wezel[j].d_y);
                      polymers[i].wezel[j].fi=atan2(polymers[i].wezel[j].d_x,polymers[i].wezel[j].d_y);
                      polymers[i].wezel[j-1].psi=(M_PI-(polymers[i].wezel[j-1].fi-polymers[i].wezel[j].fi))*deg;
            }
     }
}

void set:: find_area(double margin, double delta)
{
     int i,j;
     x_max=polymers[0].wezel[0].x;
     y_max=polymers[0].wezel[0].y;
     x_min=polymers[0].wezel[0].x;
     y_min=polymers[0].wezel[0].y;
     delta_r=fabs(delta);

     for(i=0;i<polymer_num;++i)
     for(j=0;j<polymers[i].beads_num;++j)
     {
                                         if(polymers[i].wezel[j].x>x_max) x_max=polymers[i].wezel[j].x;
                                         if(polymers[i].wezel[j].x<x_min) x_min=polymers[i].wezel[j].x;
                                         if(polymers[i].wezel[j].y>y_max) y_max=polymers[i].wezel[j].y;
                                         if(polymers[i].wezel[j].y<y_min) y_min=polymers[i].wezel[j].y;
     }

     x_max+=fabs(margin);
     y_max+=fabs(margin);
     x_min-=fabs(margin);
     y_min-=fabs(margin);

     n_x=1+static_cast<int>(x_max-x_min)/delta;
     n_y=1+static_cast<int>(y_max-y_min)/delta;
     cout<<"nx: "<<n_x<<" n_y: "<<n_y<<endl;
     matrix_rank=n_x*n_y;
}

double set:: corr_fun(int i, int j)
{
       double x1,y1,x2,y2;
       x1=ij_to_xy(i,'x');
       y1=ij_to_xy(i,'y');

       x2=ij_to_xy(j,'x');
       y2=ij_to_xy(j,'y');

       //cout<<i<<" x1:"<<x1<<" y1:"<<y1<<endl;


       double r=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
       double corr;
       corr=exp(-r);
       if(corr<0.2) return 0;
       else return corr;
}

double set::ij_to_xy(int i, char c)
{
       i%=(n_x*n_y);
       int k(i);
       k%=n_x;
       if(c=='x') return k*delta_r+x_min;
       else return (i-k)/(double)n_x*delta_r+y_min;
}

double set::matrixL_det()
{
       if(corrLT==NULL)
       {
                       cout<<"det: brak macierzy"<<endl;
                       return 0;
       }

       double det(1);
       int i;
       for(i=1;i<=matrix_rank;++i)
       {
                                  //cout<<corrLT->element(i,i)<<endl;
                                det*=corrLT->element(i,i);
       }
       return det;
}

void set::relaxation(double multiplyier, int iterations, double h)
{
      if(time==0)
      {
              int i,j;
              for(i=0;i<polymer_num;++i)
              for(j=0;j<polymers[i].beads_num;++j)
              {
                                                  polymers[i].wezel[j].gamma*=multiplyier;
              }
              cout<<"\nannealing on"<<endl;
      }
      if(time==iterations*h)
      {
              int i,j;
              for(i=0;i<polymer_num;++i)
              for(j=0;j<polymers[i].beads_num;++j)
              {
                                                  polymers[i].wezel[j].gamma/=multiplyier;
              }
              cout<<"\nannealing off"<<endl;
      }
}

/*
void set:: back_up_write()
{
     int i,j;
     ofstream plik;
     plik.open("backup.dat","wb");
     plik.write(
     plik.write(time, sizeof(double));
     plik.write(polymer_num,sizeof(int));
     for(i=0;i<polymer_num;++i)
     {
                               plik.write(polymers[i],sizeof(polymer));
                               for(j=0;j<polymers[i].beads_num;++j)
                               {
                                                                   plik.write(polymers[i].wezel[j],sizeof(bead));
                               }
     }
     plik.close();
}

void set:: back_up_read()
{
     int i,j;
     bead tmp;
     polymer tmp_pol;

     ofstream plik;
     plik.open("back.up","rb");
     plik.read(time,sizeof(double));
     plik.read(polymer_num,sizeof(int));
     for(i=0;i<polymer_num;++i)
     {
                               plik.read(tmp_pol,sizeof(polymer));
                               for(j=0;j<polymers[i].beads_num;++j)
                               {
                                                                   plik.read(polymers[i].wezel[j],sizeof(bead));
                               }
     }
}
*/
