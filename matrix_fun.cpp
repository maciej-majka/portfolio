#include<iostream>
#include<cmath>
#include"matrix.h"
#include "functions.h"
using namespace std;

matrix& matrix:: operator=(matrix const &A)
{
        //cout<<"operator="<<endl;
        if(&A==this) return *this;
        int k;
        if(blok!=NULL) usuwanie(blok);
        l_elementow=A.l_elementow;
        list *wskA,*wsk;
        wskA=A.blok;
        blok=copy_struct(A.blok);
        return *this;
}

list* const matrix:: copy_struct() const
{
      return copy_struct(blok);
}

list* const matrix:: copy_struct(list const *wsk) const
{
      int k;
      list *tmp;
      if(wsk==NULL) return NULL;
      tmp=new list();
      tmp->l_zer=wsk->l_zer;
      tmp->dl_bloku=wsk->dl_bloku;
      tmp->tab=new double[wsk->dl_bloku];
      for(k=0;k<wsk->dl_bloku;++k)
      tmp->tab[k]=wsk->tab[k];
      tmp->next=copy_struct(wsk->next);
      return tmp;
}

list * matrix:: oproznij_bufor(double *bufor, int k,list *wsk)
{
     int i;
     wsk->tab=new double[k];
     for(i=0;i<k;++i)
     wsk->tab[i]=bufor[i];
     wsk->dl_bloku=k;
     return wsk;
}

void matrix:: usuwanie(list *ptr)
{
     if(ptr->next!=NULL) usuwanie(ptr->next);
     else 
     {
          delete [] ptr->tab;
          delete ptr;
     }
}

void matrix:: wypisz()
{
     int i,licznik(0);
     if(blok==NULL) cout<<"brak danych"<<endl;
     list *wsk=blok;
     cout<<endl<<endl;
     while(wsk!=NULL)
     {
                           for(i=0;i<wsk->l_zer;++i)
                           {
                                                    cout<<" 0 ";
                                                    ++licznik;
                                                    if(licznik%l_elementow==0)
                                                    {cout<<endl;
                                                    licznik=0;}
                           }
                           if(wsk->tab!=NULL)
                           for(i=0;i<wsk->dl_bloku;++i)
                           {
                                                       cout<<" "<<wsk->tab[i]<<" ";
                                                       ++licznik;
                                                       if(licznik%l_elementow==0)
                                                       {cout<<endl;
                                                       licznik=0;}
                           }
                           wsk=wsk->next;
     }
     cout<<endl;                            
}

void matrix:: wypisz_element(int i, int j)
{
     list *wsk;
     wsk=blok;
     int suma(0);
     int numer=l_elementow*(i-1)+j;
     if(blok==NULL) {cout<<"brak danych"<<endl; return;}
     while(wsk!=NULL)
     {
                           suma+=wsk->l_zer;
                           if(suma>=numer) 
                           {cout<<" 0 "<<endl;
                           break;}
                           suma+=wsk->dl_bloku;
                           if(suma>=numer) 
                           {cout<<" "<<wsk->tab[wsk->dl_bloku-(suma-numer)-1]<<" "<<endl;
                           break;}
                           wsk=wsk->next;
     }
}

double matrix:: element(int i, int j)  //numeracja elemntow jak w matematyce tzn od 11, 12, 13...NN
{
     list *wsk;
     wsk=blok;
     int suma(0);
     int numer=l_elementow*(i-1)+j;
     if(blok==NULL) {cout<<"brak danych"<<endl; return 0;}
     while(wsk!=NULL)
     {
                           suma+=wsk->l_zer;
                           if(suma>=numer) return 0;
                           suma+=wsk->dl_bloku;
                           if(suma>=numer) return wsk->tab[wsk->dl_bloku-(suma-numer)-1];
                           wsk=wsk->next;
     }
}

void matrix:: wypisz_struct()
{
     list *wsk=blok;
     cout<<"\nstruktura:"<<endl;
     while(wsk!=NULL)
     {
                           cout<<"l_zer: "<<wsk->l_zer<<endl;
                           cout<<"dl_bloku: "<<wsk->dl_bloku<<endl;
                           for(int i(0);i<wsk->dl_bloku;++i)
                           cout<<wsk->tab[i]<<", ";
                           cout<<"\nnext"<<endl;
                           wsk=wsk->next;
     }
}


void matrix::transpose()
{
     matrix tmp(transpose_el,*this,*this);
     *this=tmp;
}


void matrix:: solve_vec(double *wektor_x,double *wektor_y)
{
     int i,j(0),licznik(0);
     double suma(0);
     if(blok==NULL) cout<<"brak danych"<<endl;
     list *wsk=blok;
     while(wsk!=NULL)
     {
                           for(i=0;i<wsk->l_zer;++i)
                           {
                                                    if(((licznik+1)%l_elementow==0)&&(licznik!=0))
                                                    {
                                                                              //cout<<"zmiana"<<endl;
                                                                              wektor_y[j]=suma;
                                                                              suma=0;
                                                                              ++j;
                                                    }
                                                    ++licznik;
                           }
                           if(wsk->tab!=NULL)
                           for(i=0;i<wsk->dl_bloku;++i)
                           {
                                                       suma+=wsk->tab[i]*wektor_x[licznik%l_elementow];
                                                       //cout<<"el="<<wsk->tab[i]<<" wek="<<wektor_x[licznik%l_elementow]<<endl;
                                                       if(((licznik+1)%l_elementow==0)&&(licznik!=0))
                                                       {
                                                                                 wektor_y[j]=suma;
                                                                                 suma=0;
                                                                                 //cout<<"zmiana"<<endl;
                                                                                 ++j;
                                                       }
                                                       ++licznik;
                           }
                           wsk=wsk->next;
     }
     //cout<<"koniec"<<endl;
}
