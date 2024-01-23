#include<iostream>
#include<cstdlib>
#include<string>
#include<ctime>
#include<cmath>
#include"set.h"
#include"functions.h"
#include<gsl/gsl_rng.h>
using namespace std;

double new_lambda(double lambda)
{
       return lambda+5;
}

double new_kT(double gamma, double kT)
{
       return kT+10;
}

double *alloc_corr(int x, int y)
{
       double *tmp=new double [x*y];
       int i;
       for(i=0;i<x*y;++i)
       tmp[i]=0;
       return tmp;
}

void add_corr_result(double *corr1,double *corr2, double *corr3,set *system,int N, int num_beads, int t_range)
{
     int i,j;
     for(i=0;i<t_range;++i)
     for(j=0;j<num_beads-2;++j)
     corr1[i*(num_beads-2)+j]+=system->get_corr_element(1,j,i)/(double)N;

     for(i=0;i<t_range;++i)
     for(j=0;j<num_beads-1;++j)
     corr2[i*(num_beads-2)+j]+=system->get_corr_element(2,j,i)/(double)N;

     for(i=0;i<t_range;++i)
     for(j=0;j<num_beads-2;++j)
     corr3[i*(num_beads-2)+j]+=system->get_corr_element(3,j,i)/(double)N;
     cout<<"usredniono korealcje"<<endl;
}

void corr_to_file(int x, int y,double *corr1,double *corr2,double *corr3,int nr)
{
     cout<<"zapis do pliku"<<endl;
     int i;
     ofstream plik;
     char str[20];
     sprintf(str,"wyniki/sr_kor/avc%i.dat",nr);
     plik.open(str);
     for(i=0;i<(x-2)*y;++i)
     {
                           if(i%(x-2)==0&&i!=0) plik<<endl;
                           plik<<corr1[i]<<" ";

     }
     plik<<endl;
     plik<<endl;
     for(i=0;i<(x-1)*y;++i)
     {
                           if(i%(x-1)==0&&i!=0) plik<<endl;
                           plik<<corr2[i]<<" ";
     }
     plik<<endl;
     plik<<endl;
     for(i=0;i<(x-2)*y;++i)
     {
                           if(i%(x-2)==0&&i!=0) plik<<endl;
                           plik<<corr3[i]<<" ";
     }
     plik<<endl;
     plik.close();
}

void delete_av_correlations(double *corr1,double *corr2,double *corr3)
{
     delete [] corr1;
     delete [] corr2;
     delete [] corr3;
}

void simulation()
{
     double T_END(100),dt=1/128.0;
     double k1(1),R0(5),psi0(2),k2(2),sigma(4.85),epsilon(-15),lambda,kT,gamma(-5),margin=1,delta=5;
     int N_chain(1), num_beads(40), seed, number(1),av_num(1),i,Rep_num(2);
     char corr;
     string katalog;

     double *av_kink_corr,*av_rad_corr,*av_cross_corr;

     gsl_rng *gen_seed=gsl_rng_alloc(gsl_rng_mt19937);
     gsl_rng_set(gen_seed,CLOCKS_PER_SEC*clock()*1000);

     ofstream log;
     log.open("wyniki/log.dat");
     log<<"nr "<<"av_num "<<"N "<<"n "<<"k1 "<<"R0 "<<"k2 "<<"psi0 "<<"epsilon "<<"sigma "<<"lambda "<<"gamma "<<"2kT "<<"margin "<<"delta "<<"corr "<<"seed"<<endl;
     set *uklad;
     //uklad(l.lancuchow,l.wezlow,R0,k1,psi0,k2,epsilon,sigma,lambda,gamma,noise,a,margin,delta,typ korelacji,katalog,seed)

     //for(epsilon=-10;epsilon>-30;epsilon+=10)
     //for(k2=0.1;k2<3;k2+=0.5)
     {
                             for(lambda=5;lambda<=10;lambda=new_lambda(lambda))
                             for(kT=5;kT<=45;kT=new_kT(gamma,kT))
                             {
                                                                av_kink_corr=alloc_corr(num_beads-2,20);
                                                                av_rad_corr=alloc_corr(num_beads-1,20);
                                                                av_cross_corr=alloc_corr(num_beads-2,20);

                                                                    for(i=0;i<Rep_num;++i)
                                                                    {
                                                                        seed=(int)30000*gsl_rng_uniform(gen_seed);
                                                                        if(lambda==0) corr='i';
                                                                        else corr='c';
                                                                        cout<<"Nr: "<<number<<endl;
                                                                        cout<<"l. lancuchow="<<N_chain<<" l. wezlow="<<num_beads<<endl;
                                                                        cout<<"potencjaly: k1="<<k1<<" R0="<<R0<<" k2="<<k2<<" epsilon="<<epsilon<<" sigma="<<sigma<<endl;
                                                                        cout<<"szum: 2kT="<<kT<<" gamma="<<gamma<<" lambda="<<lambda<<" seed="<<seed<<endl;

                                                                        log<<number<<" "<<av_num<<" "<<N_chain<<" "<<num_beads<<" "<<R0<<" "<<k1<<" "<<psi0<<" "<<k2<<" "<<epsilon<<" "<<sigma<<" "<<lambda<<" "<<gamma<<" "<<kT<<" "<<margin<<" "<<delta<<" "<<corr<<" "<<seed<<endl;
                                                                        uklad=new set(N_chain,num_beads,R0,k1,psi0,k2,epsilon,sigma,lambda,gamma,kT,margin,delta,corr,katalog,seed,number,false);
                                                                        uklad->solution(dt,T_END,'R');
                                                                        uklad->corr_trans_file();
                                                                        add_corr_result(av_kink_corr,av_rad_corr,av_cross_corr,uklad,Rep_num,num_beads,20);
                                                                        delete uklad;
                                                                        //cout<<"tutaj"<<endl;
                                                                        uklad=NULL;
                                                                        ++number;
                                                                        cout<<"****************"<<endl;
                                                                    }
                                                                    corr_to_file(num_beads,20,av_kink_corr,av_rad_corr,av_cross_corr,av_num);
                                                                    delete_av_correlations(av_kink_corr,av_rad_corr,av_cross_corr);
                                                                    ++av_num;
                             }
     }
     log.close();
     gsl_rng_free(gen_seed);
}

