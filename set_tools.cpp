#include<iostream>
#include<cmath>
#include<string>
#include<fstream>
#include"set.h"
#include"plot_anim.h"
using namespace std;

void set:: print_data()
{
     cout<<"liczba lancuchow:"<<polymer_num<<endl;
     int i,j;
     for(i=0;i<polymer_num;++i)
     {
                               cout<<"\nPOLYMER nr "<<i<<": liczba wezlow: "<<polymers[i].beads_num<<endl;
                               //cout<<"    Theta: "<<polymers[i].THETA<<endl;
                               //cout<<"    Omega: "<<polymers[i].OMEGA<<endl;
                               cout<<"    Stoch. :"<<polymers[i].sumaS_x<<" , "<<polymers[i].sumaS_y<<endl;
                               cout<<"    Sily_zew: "<<polymers[i].suma_Fx<<" , "<<polymers[i].suma_Fy<<endl;
                               cout<<"    Mass centre: Xcm: "<<polymers[i].Xcm<<" , Ycm: "<<polymers[i].Ycm<<endl;
                               for(j=0;j<polymers[i].beads_num;++j)
                               {
                                         cout<<"\nWEZEL "<<j<<":"<<endl;
                                         cout<<"  m: "<<polymers[i].wezel[j].m<<" k: "<<polymers[i].wezel[j].k<<" k2: "<<polymers[i].wezel[j].k2<<" sigma: "<<polymers[i].wezel[j].sigma<<" epsilon: "<<polymers[i].wezel[j].epsilon<<" L0: "<<polymers[i].wezel[j].L0<<endl;
                                         cout<<"  ZMIENNE:"<<endl;
                                         cout<<"    x: "<<polymers[i].wezel[j].x<<" y: "<<polymers[i].wezel[j].y<<endl;
                                         cout<<"    d: "<<polymers[i].wezel[j].d<<endl;
                                         cout<<"    d_x: "<<polymers[i].wezel[j].d_x<<" d_y: "<<polymers[i].wezel[j].d_y<<endl;
                                         cout<<"    psi: "<<polymers[i].wezel[j].psi<<endl;
                                         cout<<"  SILY:"<<endl;
                                         cout<<"    stoch.: "<<polymers[i].wezel[j].S_x<<" , "<<polymers[i].wezel[j].S_y<<endl;
                                         cout<<"    LJ_zew: "<<polymers[i].wezel[j].sumaLJ_fx_z<<" , "<<polymers[i].wezel[j].sumaLJ_fy_z<<endl;
                                         cout<<"    LJ_wew: "<<polymers[i].wezel[j].sumaLJ_fx_w<<" , "<<polymers[i].wezel[j].sumaLJ_fy_w<<endl;
                                         cout<<"    Range1: "<<polymers[i].wezel[j].sumaR_fx<<" , "<<polymers[i].wezel[j].sumaR_fy<<endl;
                                         cout<<"    F_psi: "<<polymers[i].wezel[j].sumaPSI_fx<<" , "<<polymers[i].wezel[j].sumaPSI_fy<<endl;
                                         cout<<"  PREDKOSCI:"<<endl;
                                         cout<<"    v_x: "<<polymers[i].wezel[j].v_x<<" , v_y: "<<polymers[i].wezel[j].v_y<<endl;
                                         /*
                                         cout<<"  ZMIENNE WEWNETRZNE:"<<endl;
                                         cout<<"    psi: "<<polymers[i].wezel[j].psi<<" , omega: "<<polymers[i].wezel[j].omega<<" fi: "<<polymers[i].wezel[j].fi<<endl;
                                         cout<<"    r: "<<polymers[i].wezel[j].r<<" , v_r: "<<polymers[i].wezel[j].v_r<<endl;
                                         cout<<"  PRZYSPIESZENIA_o:"<<endl;
                                         cout<<"    ax_o: "<<polymers[i].wezel[j].ax_o<<" , ay_o: "<<polymers[i].wezel[j].ay_o<<endl;
                                         cout<<"    axS_o: "<<polymers[i].wezel[j].axS_o<<" , ayS_o: "<<polymers[i].wezel[j].ayS_o<<endl;
                                         */


                               }
     }
}

void set:: f_max_length()
{
     //cout<<"f_max_length"<<endl;
     int i;
     max_length=0;
     for(i=0;i<polymer_num;++i)
     if(polymers[i].beads_num>=max_length) max_length=polymers[i].beads_num;
}

void set:: file_operation()
{
           cout<<"file_operation"<<endl;
           ofstream file;
           char str[30];
           int n=iteration_counter;
           sprintf(str,"k%.0lfpsi%.0lfa%.0lfep%1.0lfg%.0lf",100*polymers[0].wezel[1].k2,100*polymers[0].wezel[1].PSI0,100*skala,polymers[0].wezel[1].epsilon,100*polymers[0].wezel[1].gamma);
           string S1(str),S2(".dat");
           output_dir+=S1+=S2;
           cout<<"output_dir: "<<output_dir<<endl;
           file.open(output_dir.c_str());
           int i,j;
           /*
           for(i=0;i<polymer_num;++i)
                                    {
                                                              file<<polymers[i].Xcm<<" "<<polymers[i].Ycm<<" ";
                                    }
           file<<endl;
           */
           file<<"x          y          d         psi\n"<<endl;
           for(j=0;j<max_length;++j)
           {
                                    for(i=0;i<polymer_num;++i)
                                    {
                                                              if (j<polymers[i].beads_num)
                                                                                          file<<polymers[i].wezel[j].x<<" "<<polymers[i].wezel[j].y<<" "<<polymers[i].wezel[j].d<<" "<<polymers[i].wezel[j].psi;
                                                              else
                                                              {
                                                                  file<<"------- ";
                                                              }
                                    }
                                    if(j<polymer_num) file<<polymers[j].Xcm<<" "<<polymers[j].Ycm<<" ";
                                    file<<endl;
           }
           file.close();
}

void set:: file_operation2()
{
           //cout<<"file_operation"<<endl;
           ofstream file;
           char str[30];
           int n=iteration_counter;
           sprintf(str,"wyniki/final/d%i.dat",NUMBER);
           //cout<<"output_dir: "<<output_dir<<endl;
           file.open(str);
           int i,j;
           /*
           for(i=0;i<polymer_num;++i)
                                    {
                                                              file<<polymers[i].Xcm<<" "<<polymers[i].Ycm<<" ";
                                    }
           file<<endl;
           */
           file<<"x          y          d         psi\n"<<endl;
           for(j=0;j<max_length;++j)
           {
                                    for(i=0;i<polymer_num;++i)
                                    {
                                                              if (j<polymers[i].beads_num)
                                                                                          file<<polymers[i].wezel[j].x<<" "<<polymers[i].wezel[j].y<<" "<<polymers[i].wezel[j].d<<" "<<polymers[i].wezel[j].psi;
                                                              else
                                                              {
                                                                  file<<"------- ";
                                                              }
                                    }
                                    if(j<polymer_num) file<<polymers[j].Xcm<<" "<<polymers[j].Ycm<<" ";
                                    file<<endl;
           }
           file.close();
}

/*
void set:: file_tmp(string path)
{
     if(grafika.ifusing()==true)
     {
                                ofstream tmp;
                                tmp.open(path.c_str());
                                int i,j;

                                //for(i=0;i<polymer_num;++i)
                                //{
                                //                          tmp<<polymers[i].Xcm<<" "<<polymers[i].Ycm<<" ";
                                //}
                                tmp<<endl;

                                for(j=0;j<max_length;++j)
                                {
                                                         for(i=0;i<polymer_num;++i)
                                                         {
                                                                                   if (j<polymers[i].beads_num)
                                                                                   {
                                                                                                               tmp<<polymers[i].wezel[j].x<<" "<<polymers[i].wezel[j].y<<" ";
                                                                                   }
                                                                                   else
                                                                                   {
                                                                                       tmp<<"- - ";
                                                                                   }
                                                         }
                                                         if(j<polymer_num) tmp<<polymers[j].Xcm<<" "<<polymers[j].Ycm<<" ";
                                tmp<<endl;
                                }
                                tmp.close();
     }
}
*/

void set::velocity_dependence(ofstream &plik)
{
     //cout<<"obliczanie wsp. synchronizacji"<<endl;
     int i=0,j,k,dep_d;
     double dep[10],delta_vx,delta_vy,delta(0);
     plik<<time<<"  ";
     for(dep_d=1;dep_d<11;++dep_d)
     {
                  dep[dep_d-1]=0;
                  for(j=0;j<polymers[i].beads_num-dep_d;++j)
                  {
                      //nowa wersja
                      delta=0;
                      delta_vx=polymers[i].wezel[j].v_x*polymers[i].wezel[j+dep_d].v_x;
                      delta_vy=polymers[i].wezel[j].v_y*polymers[i].wezel[j+dep_d].v_y;
                      delta=delta_vx+delta_vy;
                      delta/=sqrt((polymers[i].wezel[j].v_x*polymers[i].wezel[j].v_x+polymers[i].wezel[j].v_y*polymers[i].wezel[j].v_y));
                      delta/=sqrt((polymers[i].wezel[j+dep_d].v_x*polymers[i].wezel[j+dep_d].v_x+polymers[i].wezel[j+dep_d].v_y*polymers[i].wezel[j+dep_d].v_y));
                      delta/=polymers[i].beads_num-dep_d;
                      dep[dep_d-1]+=delta;
                                                            /* stara wersja
                                                            for(k=j;k<=j+dep_d;k+=dep_d)
                                                            norm+=polymers[i].wezel[k].v_x*polymers[i].wezel[k].v_x+polymers[i].wezel[k].v_y*polymers[i].wezel[k].v_y;
                                                            delta_vx=(polymers[i].wezel[j+dep_d].v_x-polymers[i].wezel[j].v_x)*(polymers[i].wezel[j+dep_d].v_x-polymers[i].wezel[j].v_x);
                                                            delta_vy=(polymers[i].wezel[j+dep_d].v_y-polymers[i].wezel[j].v_y)*(polymers[i].wezel[j+dep_d].v_y-polymers[i].wezel[j].v_y);
                                                            dep[dep_d-1]+=(delta_vx+delta_vy);
                                                            */
                  }
                  //dep[dep_d-1]/=norm;
                  //dep[dep_d-1]=1-dep[dep_d-1];
                  plik<<dep[dep_d-1]<<" ";
     }
     plik<<endl;
     //cout<<"\nwspolczynnik: "<<norm<<" "<<delta_vx<<" "<<delta_vy<<endl;
     /*
     double av_kink(0),d_kink(0);
     for(j=1;j<polymers[i].beads_num-1;++j)
     {
                                 av_kink+=polymers[i].wezel[j].kink/(double)(polymers[i].beads_num-2);
                                 if(j<polymers[i].beads_num-2)
                                 d_kink+=fabs(polymers[i].wezel[j+1].kink-polymers[i].wezel[j].kink)/(2.0*(polymers[i].beads_num-3));
     }
     plik<<fabs(av_kink)<<" "<<d_kink;
     plik<<endl;*/
}

void set::reduced_dynamics(ofstream &plik)
{
     //cout<<"obliczanie wsp. dyn. red."<<endl;
     int i=0,j;
     double av_kink(0),d_kink(0),av_rad(0),d_rad(0);
     plik<<time<<"  ";
     for(j=1;j<polymers[i].beads_num-1;++j)
     {
                                 av_kink+=polymers[i].wezel[j].kink/(double)(polymers[i].beads_num-2);
                                 if(j<polymers[i].beads_num-2)
                                 d_kink+=fabs(polymers[i].wezel[j+1].kink-polymers[i].wezel[j].kink)/(2.0*(polymers[i].beads_num-3));
     }
     plik<<fabs(av_kink)<<" "<<d_kink<<" ";
     for(j=1;j<polymers[i].beads_num-1;++j)
     {
                                 av_rad+=polymers[i].wezel[j].rad/(double)(polymers[i].beads_num-2);
                                 if(j<polymers[i].beads_num-2)
                                 d_rad+=fabs(polymers[i].wezel[j+1].rad-polymers[i].wezel[j].rad)/(2.0*(polymers[i].beads_num-3));
     }
     plik<<fabs(av_rad)<<" "<<d_rad;
     plik<<endl;
}

polymer& polymer:: operator=(const polymer &A)
{
          if(&A==this) return *this;
          int i;
          if(beads_num!=A.beads_num)
          {
                                    delete [] wezel;
                                    beads_num=A.beads_num;
                                    wezel=new bead[beads_num];
          }
          if(wezel==NULL)
          {
                         wezel=new bead[beads_num];
          }

          for(i=0;i<beads_num;++i)
          {
                                  if(wezel!=NULL)
                                  {
                                                 wezel[i]=A.wezel[i];
                                  }
          }

        return *this;
}

void set:: delete_matrix(matrix *A)
{
     delete A;
     A=NULL;
}

void set:: set_send(FILE *p)
{
     int i,j;
     for(i=0;i<polymer_num;++i)
     for(j=0;j<polymers[i].beads_num;++j)
     {
                                         fprintf(p,"%lf\n",polymers[i].wezel[j].x);
                                         //fflush(p);
                                         fprintf(p,"%lf\n",polymers[i].wezel[j].y);
                                         //fflush(pin);
     }
     fflush(p);
}

void set:: anim_start(FILE *p)
{
     fprintf(p,"%i\n",polymer_num);
     fflush(p);

     int i,j;
     for(i=0;i<polymer_num;++i)
     {
                               fprintf(p,"%i\n",polymers[i].beads_num);
                               fflush(p);
     }

     for(i=0;i<polymer_num;++i)
     for(j=0;j<polymers[i].beads_num;++j)
     {
                                         fprintf(p,"%lf\n",polymers[i].wezel[j].x);
                                         fflush(p);
                                         fprintf(p,"%lf\n",polymers[i].wezel[j].y);
                                         fflush(p);
                                         fprintf(p,"%lf\n",polymers[i].wezel[j].sigma);
                                         fflush(p);
                                         fprintf(p,"%c\n",'p');
                                         fflush(p);
     }
}

