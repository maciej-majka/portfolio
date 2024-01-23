#include<iostream>
#include<cstdio>
#include<ctime>
#include<string>
#include "plot_anim.h"
using namespace std;

animation::animation(bool ok)
{
                          //cout<<"ANIMATION c-tor"<<endl;
                          string path1("e:\\gnuplot\\gnuplot\\binary\\pgnuplot.exe");
                          string path2("e:\\\\programy\\\\sim\\\\tmp.dat");
                          string path3("e:\\programy\\sim\\tmp.dat");
                          isusing=ok;
                          tmp_path=path2;
                                                    
                           if(isusing==true)
                           {
                                            //cout<<"isusing"<<endl;
                                       pin=popen("\n","r");
                                       pout=popen(path1.c_str(),"w");
                                       waiting(1);
                                       if(pout==NULL)
                                       {
                                                     cout<<"gnuplot: bledna sciezka\n"<<path1<<endl;
                                                     isusing=false;
                                       }
                                       else gplot_path=path1;
                           }
                           tmp_for_gplot=path2;
                           tmp_path=path3;
                           //cout<<"end animation c_tor"<<endl;
}

animation::~animation()
{
                       if(isusing==true)
                       {
                                        pclose(pin);
                                        pclose(pout);
                       }
}

bool animation:: ifusing()
{
     if(isusing==true) return true;
     else return false;
}

void animation::refresh()
{
     if(isusing==true)
     {
                      fprintf(pout,"replot\n");
                      fflush(pout);
     }
}

void animation:: start(int num)
{
     int i;
     if(isusing==true)
     {
                      //fprintf(pout,"set size ratio 1\n");
                      //fflush(pout);
                      //fprintf(pout,"set lmargin 10\n set rmargin 1\n");
                      //fflush(pout);
                      fprintf(pout,"set pointsize 2\n");
                      fflush(pout);
                      fprintf(pout,"set xrange [-70:70]\n set yrange [-70:70]\n");
                      fflush(pout);
                      fprintf(pout,"set pointsize 2\n");
                      fflush(pout);
                      fprintf(pout,"plot ");
                      char c=',';
                      for(i=0;i<2*num;i+=2)
                      {
                                           if(i==2*(num-1)) c='\n';
                                         fprintf(pout,"\"%s\" using %i:%i with linespoints pointtype 6%c ",tmp_for_gplot.c_str(),i+1,i+2,c);
                      }
                      //fprintf(pout,"\"%s\" using %i:%i with points pointtype 14 lt \"black\"\n",tmp_for_gplot.c_str(),num+1,num+2);
                      fflush(pout);
                      waiting(2);
     }
}

void animation:: end()
{
     if(isusing==true)
     {
                      cout<<"KONIEC: nacisnij klawisz"<<endl;
                      //getchar();
                      fprintf(pout,"exit\n");
                      fflush(pout);
     }
}

void animation:: waiting(double t)
{
     time_t begin(clock());
     time_t now(clock());
     while(difftime(now,begin)<(t*1000))
     now=clock();
}

void animation:: get_path()
{
     cout<<"podaj sciezke dostepu do pgnuplot.exe: ";
     cin>>gplot_path;
}

string animation:: file_path()
{
       return tmp_path;
}
