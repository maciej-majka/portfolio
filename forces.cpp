#include<iostream>
#include<cmath>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include"set.h"
using namespace std;

void set:: F_LJ()
{
     //cout<<"F_LJ"<<endl;
      int i,j,k,l;
      double r2, x1,x2,y1,y2, val1,val2, sigma, epsilon;

      for(i=0;i<polymer_num;++i)
      for(j=0;j<polymers[i].beads_num;++j)
      {
                                          polymers[i].wezel[j].sumaLJ_fx_z=0;
                                          polymers[i].wezel[j].sumaLJ_fy_z=0;
                                          polymers[i].wezel[j].sumaLJ_fx_w=0;
                                          polymers[i].wezel[j].sumaLJ_fy_w=0;
      }

      for(i=0;i<polymer_num;++i)
      for(j=0;j<polymers[i].beads_num;++j)
      {
                                           x1=polymers[i].wezel[j].x;
                                           y1=polymers[i].wezel[j].y;
                                           sigma=polymers[i].wezel[j].sigma;
                                           epsilon=polymers[i].wezel[j].epsilon;
                                           for(k=i;k<polymer_num;++k)
                                           for(l=j+1;l<polymers[i].beads_num;++l)
                                           {
                                                                                  x2=polymers[k].wezel[l].x;
                                                                                  y2=polymers[k].wezel[l].y;
                                                                                  r2=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);
                                                                                  val1=pow(1/r2,7);
                                                                                  val2=pow(1/r2,4);

                                                                                  if(i!=k)
                                                                                  {
                                                                                  polymers[i].wezel[j].sumaLJ_fx_z+=epsilon*(x2-x1)*(12*pow(sigma,12)*val1-6*pow(sigma,6)*val2);
                                                                                  polymers[i].wezel[j].sumaLJ_fy_z+=epsilon*(y2-y1)*(12*pow(sigma,12)*val1-6*pow(sigma,6)*val2);

                                                                                  polymers[k].wezel[l].sumaLJ_fx_z+=epsilon*(x1-x2)*(12*pow(sigma,12)*val1-6*pow(sigma,6)*val2);
                                                                                  polymers[k].wezel[l].sumaLJ_fy_z+=epsilon*(y1-y2)*(12*pow(sigma,12)*val1-6*pow(sigma,6)*val2);
                                                                                  }
                                                                                  else
                                                                                  {
                                                                                  polymers[i].wezel[j].sumaLJ_fx_w+=epsilon*(x2-x1)*(12*pow(sigma,12)*val1-6*pow(sigma,6)*val2);
                                                                                  polymers[i].wezel[j].sumaLJ_fy_w+=epsilon*(y2-y1)*(12*pow(sigma,12)*val1-6*pow(sigma,6)*val2);

                                                                                  polymers[k].wezel[l].sumaLJ_fx_w+=epsilon*(x1-x2)*(12*pow(sigma,12)*val1-6*pow(sigma,6)*val2);
                                                                                  polymers[k].wezel[l].sumaLJ_fy_w+=epsilon*(y1-y2)*(12*pow(sigma,12)*val1-6*pow(sigma,6)*val2);
                                                                                  }
                                           }
      }
}

void set:: F_range() //sily sprezystosci segmentu
{
     //cout<<"F_range"<<endl;
     int i,j;
     double Fr, fx_ij, fy_ij,fx_ij1,fy_ij1;
     //zerowanie skladowych
     for(i=0;i<polymer_num;++i)
     for(j=0;j<polymers[i].beads_num;++j)
     {
                                           polymers[i].wezel[j].sumaR_fx=0;
                                           polymers[i].wezel[j].sumaR_fy=0;
     }

     for(i=0;i<polymer_num;++i)
     for(j=1;j<polymers[i].beads_num;++j)
     {
                                           fx_ij=polymers[i].wezel[j].k*(polymers[i].wezel[j].d-polymers[i].wezel[j].R0)*polymers[i].wezel[j].d_x/polymers[i].wezel[j].d;
                                           fy_ij=polymers[i].wezel[j].k*(polymers[i].wezel[j].d-polymers[i].wezel[j].R0)*polymers[i].wezel[j].d_y/polymers[i].wezel[j].d;
                                           /*
                                           if(j<polymers[i].beads_num-1)
                                           {
                                                  fx_ij1=polymers[i].wezel[j+1].k*(polymers[i].wezel[j+1].d-polymers[i].wezel[j+1].R0)*polymers[i].wezel[j+1].d_x/polymers[i].wezel[j+1].d;
                                                  fy_ij1=polymers[i].wezel[j+1].k*(polymers[i].wezel[j+1].d-polymers[i].wezel[j+1].R0)*polymers[i].wezel[j+1].d_y/polymers[i].wezel[j+1].d;
                                           }
                                           else
                                           {
                                               fx_ij1=0;
                                               fy_ij1=0;
                                           }
                                           */
                                           polymers[i].wezel[j-1].sumaR_fx+=fx_ij;
                                           polymers[i].wezel[j-1].sumaR_fy+=fy_ij;
                                           polymers[i].wezel[j].sumaR_fx-=fx_ij;
                                           polymers[i].wezel[j].sumaR_fy-=fy_ij;
     }
}

void set:: F_T()
{
     int i,j;
     for(i=0;i<polymer_num;++i)
     for(j=0;j<polymers[i].beads_num;++j)
     {
                                         polymers[i].wezel[j].sumaT_fx=polymers[i].wezel[j].gamma*polymers[i].wezel[j].v_x;
                                         polymers[i].wezel[j].sumaT_fy=polymers[i].wezel[j].gamma*polymers[i].wezel[j].v_y;
     }
}

double set:: U(int i, int j)
{
      int k,l;
      double sumaLJ(0), sumaPSI(0), sumaR(0);
      double s, r, r_x, r_y;
      for(k=0;k<polymer_num;++k)
      for(l=0;l<polymers[k].beads_num;++l)
      {
                                          if((k==i)&&(l==j)) continue;
                                          s=polymers[i].wezel[j].sigma;
                                          r_x=polymers[i].wezel[j].x-polymers[k].wezel[l].x;
                                          r_y=polymers[i].wezel[j].y-polymers[k].wezel[l].y;
                                          r=r_x*r_x+r_y*r_y;
                                          sumaLJ+=0.5*polymers[i].wezel[j].epsilon*(pow(s,12)/pow(r,6)-pow(s,6)/pow(fabs(r),3));
                                          //cout<<l<<","<<k<<" s:"<<s<<" r:"<<r<<"energia sumaLJ "<<sumaLJ<<endl;
                                          //getchar();
      }
      if(j>0) sumaR=0.5*polymers[i].wezel[j].k*(polymers[i].wezel[j].d-polymers[i].wezel[j].R0)*(polymers[i].wezel[j].d-polymers[i].wezel[j].R0);
      //if(j>1) sumaPSI=polymers[i].wezel[j].k2*cos(polymers[i].wezel[j].psi-polymers[i].wezel[j].PSI0);//*(polymers[i].wezel[j].psi-polymers[i].wezel[j].PSI0);
      if(j>1) sumaPSI=polymers[i].wezel[j].k2*pow(sqrt((polymers[i].wezel[j].x-polymers[i].wezel[j-2].x)*(polymers[i].wezel[j].x-polymers[i].wezel[j-2].x)+(polymers[i].wezel[j].y-polymers[i].wezel[j-2].y)*(polymers[i].wezel[j].y-polymers[i].wezel[j-2].y))-polymers[i].wezel[j].L0,2);
      return (sumaLJ+sumaR+sumaPSI);
}

double set::U_LJ()
{
      int k,l,i,j;
      double sumaLJ(0);
      double s, r, r_x, r_y;
      for(i=0;i<polymer_num;++i)
      for(j=0;j<polymers[i].beads_num;++j)
      for(k=0;k<polymer_num;++k)
      for(l=0;l<polymers[k].beads_num;++l)
      {
                                          if((k==i)&&(l==j)) continue;
                                          s=polymers[i].wezel[j].sigma;
                                          r_x=polymers[i].wezel[j].x-polymers[k].wezel[l].x;
                                          r_y=polymers[i].wezel[j].y-polymers[k].wezel[l].y;
                                          r=r_x*r_x+r_y*r_y;
                                          sumaLJ+=0.5*(pow(s,12)/pow(r,6)-pow(s,6)/pow(r,3));
                                          //cout<<l<<","<<k<<" s:"<<s<<" r:"<<r<<"energia sumaLJ "<<sumaLJ<<endl;
                                          //getchar();
      }
      return sumaLJ;
}

double set::U_k()
{
       double sumaR;
       int i,j;
       for(i=0;i<polymer_num;++i)
       for(j=0;j<polymers[i].beads_num;++j)
       {
                                           if(j>0) sumaR=0.5*polymers[i].wezel[j].k*(polymers[i].wezel[j].d-polymers[i].wezel[j].R0)*(polymers[i].wezel[j].d-polymers[i].wezel[j].R0);
       }
       return sumaR;
}

double set:: U_psi()
{
       int i,j;
       double sumaPSI;
       for(i=0;i<polymer_num;++i)
       for(j=0;j<polymers[i].beads_num;++j)
       {
                                           if(j>1) sumaPSI=polymers[i].wezel[j].k2*cos(polymers[i].wezel[j].psi-polymers[i].wezel[j].PSI0);//*(polymers[i].wezel[j].psi-polymers[i].wezel[j].PSI0);
       }
       return sumaPSI;
}

double set::Ek()
{
       int i,j;
       double sumaEk;
       for(i=0;i<polymer_num;++i)
       for(j=0;j<polymers[i].beads_num;++j)
       {
                                           sumaEk+=0.5*polymers[i].wezel[j].m*(polymers[i].wezel[j].v_x*polymers[i].wezel[j].v_x+polymers[i].wezel[j].v_y*polymers[i].wezel[j].v_y);
       }
       return sumaEk;
}

double set:: energy()
{
     int i,j;
     double sumaEk(0), sumaU(0);

     for(i=0;i<polymer_num;++i)
     for(j=0;j<polymers[i].beads_num;++j)
     {
                                       sumaEk+=0.5*polymers[i].wezel[j].m*(polymers[i].wezel[j].v_x*polymers[i].wezel[j].v_x+polymers[i].wezel[j].v_y*polymers[i].wezel[j].v_y);
                                       sumaU+=U(i,j);
     }
     return (sumaU+sumaEk);
}

void set:: F_stochastic()
{
     if(polymers[0].wezel[0].noise_level==0) return;
     //cout<<"F stochastic"<<endl;
     int i,j,nx,ny,index;
     double x,y,fi;
     //cout<<"tutaj1"<<endl;

     for(i=0;i<matrix_rank;++i)
     uncorr_F[i]=fabs(gsl_ran_gaussian(generator,1));

     //cout<<"tutaj"<<endl;
     corrLT->solve_vec(uncorr_F,corr_F);

     for(i=0;i<polymer_num;++i)
     for(j=0;j<polymers[i].beads_num;++j)
     {
                                         nx=(polymers[i].wezel[j].x-x_min)/delta_r;
                                         ny=(polymers[i].wezel[j].y-y_min)/delta_r;
                                         nx%=n_x;
                                         ny%=n_y;
                                         index=ny*n_y+nx;
                                         fi=6.283185*gsl_rng_uniform(generator);
                                         //cout<<"fi="<<fi<<" nx: "<<nx<<" ny: "<<ny<<endl;
                                         polymers[i].wezel[j].S_x=corr_F[index]*cos(fi);
                                         polymers[i].wezel[j].S_y=corr_F[index]*sin(fi);
     }


}

void set:: F_stochastic_v2()
{
     if(polymers[0].wezel[0].noise_level==0) return;
     //cout<<"F stochastic"<<endl;
     int i,j,nx,ny,index00,index01,index11,index10;
     double x,y;
     double w1,w2,w3,w4,xi;
     //cout<<"tutaj1"<<endl;


     if(typ_korelacji=='i')
     {
                           double alfa=corrLT->element(1,1);
                           for(i=0;i<polymer_num;++i)
                           for(j=0;j<polymers[i].beads_num;++j)
                           {
                                                               //cout<<alfa<<endl;
                                                               polymers[i].wezel[j].S_x=alfa*gsl_ran_gaussian(generator,1);
                                                               polymers[i].wezel[j].S_y=alfa*gsl_ran_gaussian(generator,1);
                           }
     }

     if(typ_korelacji=='c')
     {
                           for(i=0;i<matrix_rank;++i)
                           uncorr_F[i]=gsl_ran_gaussian(generator,1);

                           //cout<<"tutaj"<<endl;
                           corrLT->solve_vec(uncorr_F,corr_F);

                           for(i=0;i<polymer_num;++i)
                           for(j=0;j<polymers[i].beads_num;++j)
                           {
                                         nx=(polymers[i].wezel[j].x-x_min)/delta_r;
                                         ny=(polymers[i].wezel[j].y-y_min)/delta_r;
                                         nx%=n_x;
                                         ny%=n_y;
                                         index00=ny*n_x+nx;
                                         index01=ny*n_x+(nx+1)%n_x;
                                         index11=((ny+1)%n_y)*n_x+(nx+1)%n_x;
                                         index10=((ny+1)%n_y)*n_x+nx;
                                         w1=(polymers[i].wezel[j].x-(x_min+nx*delta_r))/delta_r;
                                         w2=1.0-w1;
                                         w3=(polymers[i].wezel[j].y-(y_min+ny*delta_r))/delta_r;
                                         w4=1.0-w3;
                                         //cout<<"i="<<i<<" j="<<j<<" nx: "<<nx<<" ny: "<<ny<<" index="<<index<<" S="<<corr_F[index]<<endl;
                                         polymers[i].wezel[j].S_x=w1*w3*corr_F[index00]+w2*w3*corr_F[index01]+w2*w4*corr_F[index11]+w1*w4*corr_F[index10];
                           }

                           for(i=0;i<matrix_rank;++i)
                           uncorr_F[i]=gsl_ran_gaussian(generator,1);

                           //cout<<"tutaj"<<endl;
                           corrLT->solve_vec(uncorr_F,corr_F);

                           for(i=0;i<polymer_num;++i)
                           for(j=0;j<polymers[i].beads_num;++j)
                           {
                                         nx=(polymers[i].wezel[j].x-x_min)/delta_r;
                                         ny=(polymers[i].wezel[j].y-y_min)/delta_r;
                                         nx%=n_x;
                                         ny%=n_y;
                                         index00=ny*n_x+nx;
                                         index01=ny*n_x+(nx+1)%n_x;
                                         index11=((ny+1)%n_y)*n_x+(nx+1)%n_x;
                                         index10=((ny+1)%n_y)*n_x+nx;
                                         w1=(polymers[i].wezel[j].x-(x_min+nx*delta_r))/delta_r;
                                         w2=1.0-w1;
                                         w3=(polymers[i].wezel[j].y-(y_min+ny*delta_r))/delta_r;
                                         w4=1.0-w3;
                                         polymers[i].wezel[j].S_y=w2*w4*corr_F[index00]+w1*w4*corr_F[index01]+w1*w3*corr_F[index11]+w2*w3*corr_F[index10];
                                         /*
                                         cout<<"\nx="<<polymers[i].wezel[j].x<<" y="<<polymers[i].wezel[j].y<<" S="<<polymers[i].wezel[j].S_y<<endl;
                                         cout<<"x00="<<x_min+nx*delta_r<<" y00="<<y_min+ny*delta_r<<" w1="<<w1<<" S00="<<corr_F[index00]<<endl;
                                         cout<<"x01="<<x_min+((nx+1)%n_x)*delta_r<<" y01="<<y_min+ny*delta_r<<" w2="<<w2<<" S01="<<corr_F[index01]<<endl;
                                         cout<<"x01="<<x_min+((nx+1)%n_x)*delta_r<<" y01="<<y_min+((ny+1)%n_y)*delta_r<<" w3="<<w3<<" S11="<<corr_F[index11]<<endl;
                                         cout<<"x01="<<x_min+nx*delta_r<<" y01="<<y_min+((ny+1)%n_y)*delta_r<<" w4="<<w4<<" S10="<<corr_F[index10]<<endl;
                                         */

                           }
     }

}

void set:: F_second()
{
     //cout<<"F_second"<<endl;
     int i,j;
     double fx, fy, L;
     //zerowanie skladowych
     for(i=0;i<polymer_num;++i)
     for(j=0;j<polymers[i].beads_num;++j)
     {
                                           polymers[i].wezel[j].sumaPSI_fx=0;
                                           polymers[i].wezel[j].sumaPSI_fy=0;
     }

     for(i=0;i<polymer_num;++i)
     for(j=2;j<polymers[i].beads_num;++j)
     {
                                           //obliczam skladowe sily sprezystosci
                                           L=sqrt((polymers[i].wezel[j].x-polymers[i].wezel[j-2].x)*(polymers[i].wezel[j].x-polymers[i].wezel[j-2].x)+(polymers[i].wezel[j].y-polymers[i].wezel[j-2].y)*(polymers[i].wezel[j].y-polymers[i].wezel[j-2].y));
                                           fx=(L-polymers[i].wezel[j].L0)/L*(polymers[i].wezel[j].x-polymers[i].wezel[j-2].x)*polymers[i].wezel[j].k2;
                                           fy=(L-polymers[i].wezel[j].L0)/L*(polymers[i].wezel[j].y-polymers[i].wezel[j-2].y)*polymers[i].wezel[j].k2;

                                           //obliczam sile dla i-tego wezla
                                           polymers[i].wezel[j-2].sumaPSI_fx+=fx;
                                           polymers[i].wezel[j-2].sumaPSI_fy+=fy;
                                           //obliczam sile dla i+2 -ego wezla
                                           polymers[i].wezel[j].sumaPSI_fx-=fx;
                                           polymers[i].wezel[j].sumaPSI_fy-=fy;
     }

}


void set:: suma_F()  //oblicza sumy sil na wezel i ca³y polimer
{
      int i,j;
      for(i=0;i<polymer_num;++i)
      for(j=0;j<polymers[i].beads_num;++j)
      {
                                          //sumowanie sil zewnetrzych (LJ) i wewnetrznych (LJ, R, PSI)
                                          polymers[i].wezel[j].suma_fx=polymers[i].wezel[j].sumaLJ_fx_z+polymers[i].wezel[j].sumaLJ_fx_w+polymers[i].wezel[j].sumaR_fx+polymers[i].wezel[j].sumaPSI_fx+polymers[i].wezel[j].sumaT_fx;
                                          polymers[i].wezel[j].suma_fy=polymers[i].wezel[j].sumaLJ_fy_z+polymers[i].wezel[j].sumaLJ_fy_w+polymers[i].wezel[j].sumaR_fy+polymers[i].wezel[j].sumaPSI_fy+polymers[i].wezel[j].sumaT_fy;
      }

      //obliczanie sil zewnetrznych dzialajacych na polimer (LJ i stochastyczne)
      for(i=0;i<polymer_num;++i)
      {
                                polymers[i].suma_Fx=0;
                                polymers[i].suma_Fy=0;

                                for(j=0;j<polymers[i].beads_num;++j)
                                {
                                          polymers[i].suma_Fx+=polymers[i].wezel[j].sumaLJ_fx_z;
                                          polymers[i].suma_Fy+=polymers[i].wezel[j].sumaLJ_fy_z;
                                          polymers[i].sumaS_x+=polymers[i].wezel[j].S_x;
                                          polymers[i].sumaS_y+=polymers[i].wezel[j].S_y;
                                }
      }
}

void set:: forces()
{
     F_LJ();
     F_range();
     F_second();
     //F_stochastic_v2();
     F_stochastic_dynamic();
     F_T();
     suma_F();
}

/*
void set:: F_psi()
{
     //cout<<"F_psi"<<endl;
     int i,j;
     double dpsi_ij_x(0),dpsi_ij1_x(0),dpsi_ij2_x(0),dpsi_ij_y(0),dpsi_ij1_y(0),dpsi_ij2_y(0),alfa_ij(0),alfa_ij1(0),alfa_ij2(0);
     for(i=0;i<polymer_num;++i)
     for(j=0;j<polymers[i].beads_num;++j)
     {
                                           polymers[i].wezel[j].sumaPSI_fx=0;
                                           polymers[i].wezel[j].sumaPSI_fy=0;
     }

     for(i=0;i<polymer_num;++i)
     for(j=0;j<polymers[i].beads_num;++j)
     {
                                         //cout<<"j: "<<j<<endl;
                                         if(j>1)
                                         {
                                                alfa_ij=-(polymers[i].wezel[j-1].d_x*polymers[i].wezel[j].d_x+polymers[i].wezel[j-1].d_y*polymers[i].wezel[j].d_y)/(polymers[i].wezel[j-1].d*polymers[i].wezel[j].d);
                                                if(alfa_ij<=-1)
                                                {
                                                               alfa_ij+=1E-14;
                                                               //cout<<"tu ij alfa="<<alfa_ij<<endl;
                                                }
                                         }

                                         if((j<polymers[i].beads_num-1)&&(j>0))
                                         {
                                                alfa_ij1=-(polymers[i].wezel[j].d_x*polymers[i].wezel[j+1].d_x+polymers[i].wezel[j].d_y*polymers[i].wezel[j+1].d_y)/(polymers[i].wezel[j].d*polymers[i].wezel[j+1].d);
                                                if(alfa_ij1<=-1)
                                                {
                                                                alfa_ij1+=1E-14;
                                                                //cout<<"tu ij1 alfa="<<alfa_ij1<<endl;
                                                }
                                         }

                                         if(j<polymers[i].beads_num-2)
                                         {
                                                alfa_ij2=-(polymers[i].wezel[j+1].d_x*polymers[i].wezel[j+2].d_x+polymers[i].wezel[j+1].d_y*polymers[i].wezel[j+2].d_y)/(polymers[i].wezel[j+1].d*polymers[i].wezel[j+2].d);
                                                if(alfa_ij2<=-1)
                                                 {
                                                                alfa_ij2+=1E-14;
                                                                //cout<<"tu ij2 alfa="<<alfa_ij2<<endl;
                                                 }
                                         }

                                         if(j>1)
                                         {
                                                dpsi_ij_x=(polymers[i].wezel[j-1].d_x/(polymers[i].wezel[j-1].d*polymers[i].wezel[j].d)-alfa_ij*polymers[i].wezel[j].d_x/(polymers[i].wezel[j].d*polymers[i].wezel[j].d))/sqrt(1-alfa_ij*alfa_ij);
                                                dpsi_ij_y=(polymers[i].wezel[j-1].d_y/(polymers[i].wezel[j-1].d*polymers[i].wezel[j].d)-alfa_ij*polymers[i].wezel[j].d_y/(polymers[i].wezel[j].d*polymers[i].wezel[j].d))/sqrt(1-alfa_ij*alfa_ij);
                                                /*
                                                cout<<"  dpsi_ij_x: "<<dpsi_ij_x<<endl;
                                                cout<<"j-1: d="<<polymers[i].wezel[j-1].d<<" d_x="<<polymers[i].wezel[j-1].d_x<<"d_y="<<polymers[i].wezel[j-1].d_y<<endl;
                                                cout<<"j: d="<<polymers[i].wezel[j].d<<" d_x="<<polymers[i].wezel[j].d_x<<"d_y="<<polymers[i].wezel[j].d_y<<endl;
                                                cout<<"alfa ij:"<<alfa_ij<<endl;

                                         }
                                         else
                                         {
                                             dpsi_ij_x=0;
                                             dpsi_ij_y=0;
                                         }

                                         if((j<polymers[i].beads_num-1)&&(j>0))
                                         {
                                                                      dpsi_ij1_x=((polymers[i].wezel[j+1].d_x-polymers[i].wezel[j].d_x)/(polymers[i].wezel[j+1].d*polymers[i].wezel[j].d)-alfa_ij1*(polymers[i].wezel[j].d_x/(polymers[i].wezel[j].d*polymers[i].wezel[j].d)-polymers[i].wezel[j+1].d_x/(polymers[i].wezel[j+1].d*polymers[i].wezel[j+1].d)))/sqrt(1-alfa_ij1*alfa_ij1);
                                                                      dpsi_ij1_y=((polymers[i].wezel[j+1].d_y-polymers[i].wezel[j].d_y)/(polymers[i].wezel[j+1].d*polymers[i].wezel[j].d)-alfa_ij1*(polymers[i].wezel[j].d_y/(polymers[i].wezel[j].d*polymers[i].wezel[j].d)-polymers[i].wezel[j+1].d_y/(polymers[i].wezel[j+1].d*polymers[i].wezel[j+1].d)))/sqrt(1-alfa_ij1*alfa_ij1);

                                                                      //cout<<"  dpsi_ij1_x: "<<dpsi_ij1_x<<endl;
                                                                      //cout<<"alfa ij1:"<<alfa_ij1<<endl;

                                         }
                                         else
                                         {
                                             dpsi_ij1_x=0;
                                             dpsi_ij1_y=0;
                                         }

                                         if(j<polymers[i].beads_num-2)
                                         {
                                                                      dpsi_ij2_x=-(polymers[i].wezel[j+2].d_x/(polymers[i].wezel[j+1].d*polymers[i].wezel[j+2].d)-alfa_ij2*polymers[i].wezel[j+1].d_x/(polymers[i].wezel[j+1].d*polymers[i].wezel[j+1].d))/sqrt(1-alfa_ij2*alfa_ij2);
                                                                      dpsi_ij2_y=-(polymers[i].wezel[j+2].d_y/(polymers[i].wezel[j+1].d*polymers[i].wezel[j+2].d)-alfa_ij2*polymers[i].wezel[j+1].d_y/(polymers[i].wezel[j+1].d*polymers[i].wezel[j+1].d))/sqrt(1-alfa_ij2*alfa_ij2);
                                                                      /*
                                                                      cout<<"  dpsi_ij2_x: "<<dpsi_ij2_x<<endl;
                                                                      cout<<"j+1: d="<<polymers[i].wezel[j+1].d<<" d_x="<<polymers[i].wezel[j+1].d_x<<"d_y="<<polymers[i].wezel[j+1].d_y<<endl;
                                                                      cout<<"j+2: d="<<polymers[i].wezel[j+2].d<<" d_x="<<polymers[i].wezel[j+2].d_x<<"d_y="<<polymers[i].wezel[j+2].d_y<<endl;
                                                                      cout<<"alfa ij2:"<<alfa_ij2<<endl;

                                         }
                                         else
                                         {
                                             dpsi_ij2_x=0;
                                             dpsi_ij2_y=0;
                                         }

                                         if(j>1)
                                         {
                                                polymers[i].wezel[j].sumaPSI_fx-=polymers[i].wezel[j].k2*sin(acos(alfa_ij)-polymers[i].wezel[j].PSI0)*dpsi_ij_x;
                                                polymers[i].wezel[j].sumaPSI_fy-=polymers[i].wezel[j].k2*sin(acos(alfa_ij)-polymers[i].wezel[j].PSI0)*dpsi_ij_y;
                                                //cout<<"f_psi "<<j<<" "<<polymers[i].wezel[j].sumaPSI_fx<<endl;
                                         }

                                         if((j<polymers[i].beads_num-1)&&(j>0))
                                         {
                                                                      polymers[i].wezel[j].sumaPSI_fx-=polymers[i].wezel[j+1].k2*sin(acos(alfa_ij1)-polymers[i].wezel[j+1].PSI0)*dpsi_ij1_x;
                                                                      polymers[i].wezel[j].sumaPSI_fy-=polymers[i].wezel[j+1].k2*sin(acos(alfa_ij1)-polymers[i].wezel[j+1].PSI0)*dpsi_ij1_y;
                                                                      //cout<<"f_psi "<<j<<" "<<polymers[i].wezel[j].sumaPSI_fx<<endl;
                                         }
                                         if(j<polymers[i].beads_num-2)
                                         {
                                                                      polymers[i].wezel[j].sumaPSI_fx-=polymers[i].wezel[j+2].k2*sin(acos(alfa_ij2)-polymers[i].wezel[j+2].PSI0)*dpsi_ij2_x;
                                                                      polymers[i].wezel[j].sumaPSI_fy-=polymers[i].wezel[j+2].k2*sin(acos(alfa_ij2)-polymers[i].wezel[j+2].PSI0)*dpsi_ij2_y;
                                                                      //cout<<"f_psi "<<j<<" "<<polymers[i].wezel[j].sumaPSI_fx<<endl;
                                         }
     }
}
*/

void set::F_stochastic_dynamic()
{
     int w=polymers[0].beads_num;
     double S[w][w],L[w][w],fx[w],fy[w];
     int i,j;
     double r,x_i,y_i,x_j,y_j;

     if(typ_korelacji=='i')
     {
                           double alfa=sqrt(level_kT*fabs(gamma_coef));
                           for(i=0;i<polymer_num;++i)
                           for(j=0;j<polymers[i].beads_num;++j)
                           {
                                                               //cout<<alfa<<endl;
                                                               polymers[i].wezel[j].S_x=alfa*gsl_ran_gaussian(generator,1);
                                                               polymers[i].wezel[j].S_y=alfa*gsl_ran_gaussian(generator,1);
                           }
     }

     if(typ_korelacji!='c') return;

     int strip[polymers[0].beads_num];
     double lev(1);
     double coeff=level_kT*fabs(gamma_coef);

     //for(i=0;i<w;++i)
     //strip[i]=-1;

     for(i=0;i<polymers[0].beads_num;++i)
     {
            S[i][i]=coeff;

            for(j=i+1;j<polymers[0].beads_num;++j)
            {
                                         x_i=polymers[0].wezel[i].x;
                                         y_i=polymers[0].wezel[i].y;
                                         x_j=polymers[0].wezel[j].x;
                                         y_j=polymers[0].wezel[j].y;
                                         r=sqrt((x_i-x_j)*(x_i-x_j)+(y_i-y_j)*(y_i-y_j));
                                         lev=exp(-r/length);
                                         //if(lev<0.1&&(strip[j]==i-1)) strip[j]=i; //okresla od ktorego miejsca pasmo jest istotnie zapelnione
                                         S[i][j]=coeff*lev;
                                         S[j][i]=S[i][j];
            }
     }

     //for(i=0;i<w;++i)
     //if(strip[i]==-1) strip[i]=0;

    /*
     if(iteration_counter%128==64)
     {
     cout<<"S"<<endl;
     for(j=0;j<w;++j)
     cout<<strip[j]<<" ";
     cout<<endl;

     for(i=0;i<w;++i)
     {
                     for(j=0;j<w;++j)
                     cout<<S[i][j]<<" ";
                     cout<<endl;
     }
     }
     */

     int k,m;
     double suma(0);
     for(i=0;i<polymers[0].beads_num;++i)
     for(j=0;j<=i;++j)
     {
                                         suma=0;
                                         if(i==j)
                                         {
                                                 for(k=0/*strip[i]*/;k<i;++k)
                                                 suma+=L[i][k]*L[i][k];
                                                 //cout<<suma<<endl;
                                                 L[i][i]=sqrt(S[i][i]-suma);
                                         }
                                         else
                                         {
                                             //m=max(strip[i],strip[j]);
                                             for(k=0;k<j;++k)
                                             suma+=L[i][k]*L[j][k];
                                             //cout<<suma<<" "<<L[i][i]<<endl;
                                             L[i][j]=(S[i][j]-suma)/L[j][j];
                                         }
     }

     for(i=0;i<w;++i)
     {
                     fx[i]=gsl_ran_gaussian(generator,1);
                     fy[i]=gsl_ran_gaussian(generator,1);
                     polymers[0].wezel[i].S_x=0;
                     polymers[0].wezel[i].S_y=0;
                     //cout<<fx[i]<<endl;
     }

     for(i=0;i<w;++i)
     for(j=0;j<=i;++j)
     {
                    polymers[0].wezel[i].S_x+=L[i][j]*fx[j];
                    polymers[0].wezel[i].S_y+=L[i][j]*fy[j];
     }

     //check_force_dynamics();


     //for(i=0;i<w;++i)
     //cout<<fx[i]<<" "<<polymers[0].wezel[i].S_x<<" "<<fy[i]<<" "<<polymers[0].wezel[i].S_y<<endl;
    /*
     if(iteration_counter%128==16)
     {
     cout<<"\nL"<<endl;
     for(i=0;i<w;++i)
     {
                     for(j=0;j<w;++j)
                     cout<<L[i][j]<<" ";
                     cout<<endl;
     }
     }
    */
}
