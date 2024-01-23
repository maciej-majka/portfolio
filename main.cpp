//#include<GL/glut.h>
#include<iostream>
#include"matrix.h"
#include"functions.h"
#include"set.h"
#include"animGL.h"
#include<cmath>
#include<ctime>

#include<cstdio>
using namespace std;

int main(int argc, char **argv)
{
    cout<<"*****************"<<endl;
    cout<<"*  NOISE * SIM  *"<<endl;
    cout<<"*****************"<<endl;

    if(true)
    {
        int N,i,j;
        double T_END;

        cout<<"podaj czas trwania symulacji"<<endl;
        cin>>T_END;
        cout<<"podaj dzielnik kroku"<<endl;
        double delta;
        cin>>delta;

        char str[30];
        double a,k2=0.25,psi0=3;
        string katalog;

        sprintf(str,"md wyniki\\k%.0lfpsi%.0lf\\",100*k2,100*psi0);
        katalog=str;
        system(katalog.c_str());
        double sigma=pow(2,1/6.0);
        cout<<"podaj typ korelacji:";
        char c;
        cin>>c;
        double kT2,lambda(0);
        if(c=='c')
        {
              cout<<"podaj lamdba:";
              cin>>lambda;
        }
        cout<<"podaj natezenie szumu:";
        cin>>kT2;
        //uklad(l.lancuchow,l.wezlow,R0,k1,psi0,k2,epsilon,sigma,lambda,gamma,noise,a,margin,delta,typ korelacji,katalog,seed)
        set uklad(1,64,7,7,1.85,2,-1,3,lambda,-20,kT2,20,4,c,katalog,8612345,0,true);

        //FILE *p;
        //uklad.anim_start(p);

        //uklad.print_data();
        //uklad.file_operation();
        double dt=1/delta;

        //uklad.F_stochastic_dynamic();
        time_t t0,t;
        double dif;
        time(&t0);
        uklad.solution(dt,T_END,'R');
        time(&t);
        dif=difftime(t,t0);
        cout<<"czas: "<<dif<<endl;

        //uklad.print_noise();
        //uklad.force_dynamic_file();
        //uklad.print_data();
        //uklad.print_kink_correlations();
        //uklad.print_rad_correlations();
        //uklad.print_cross_correlations();
        uklad.corr_trans_file();
        //uklad.check_force_correlation(100);
    }
    else
    simulation();
    return 0;
}
