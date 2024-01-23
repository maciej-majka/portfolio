#ifndef PLOT_ANIM_H
#define PLOT_ANIM_H
#include<cstdio>
#include<iostream>
#include<string>
using namespace std;

class animation
{
      public:
             animation(bool);
             ~animation();
             void get_path(); //pobiera lokalizacje gnuplota
             void start(int); //rozpoczyna animacjê
             void end(); //konczy
             void refresh(); //odswieza wykres
             void waiting(double); //przerwa na wyœwietlenie
             bool ifusing();
             string file_path();
      private:
              bool isusing;  //okreœla czy wykonywac animacjê
              string gplot_path; //lokalizacja pgnuplot.exe
              string tmp_path; //lokalizacj¹ tmp.dat
              string tmp_for_gplot;
              FILE *pin,*pout; //pipe do gnuplota
              
};

#endif
