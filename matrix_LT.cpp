#include<iostream>
#include<cmath>
#include "matrix.h"
using namespace std;

//funkcja dokonujaca rozkladu cholesky'ego

double matrix:: LT(int i, int j, matrix &A, list *blok, double *bufor)
{
       int k, pos(0), numer(0), z(1),mk(0);
       double suma(0),l,l1,l2;
       bool buf(false);
       list *wsk;
       wsk=blok;
       //cout<<"bufor "<<bufor[0]<<","<<bufor[1]<<","<<bufor[2]<<","<<bufor[3]<<","<<bufor[4]<<endl;
       if(i==j)
       {
               //cout<<"\ni=j="<<i<<endl;
               if(i==1&&j==1) return sqrt(A.element(1,1));
               for(k=1;k<=i-1;++k)
               {                              
                                 numer=i+(k-1)*A.l_elementow;
                                 //cout<<"nr "<<numer<<endl;
                                 //cout<<"nr "<<numer<<" pos "<<pos<<endl;
                                 while(1)
                                 {
                                 if(buf==false)
                                 {
                                 if(numer>pos&&mk==0) pos+=wsk->l_zer;
                                 if(pos>=numer&&(mk==0||mk==1)&&numer!=0) {mk=1; /*cout<<"zero"<<endl;*/ break;}
                                 if(mk==1) mk=0;
                                 if(wsk->next==NULL&&wsk->tab==NULL) {buf=true; /*cout<<"w buforze! buf="<<buf<<endl;*/ }
                                 if(numer>pos&&mk==0&&buf==false) pos+=wsk->dl_bloku;
                                 if(pos>=numer&&(mk==0||mk==2)&&buf==false)
                                 {
                                                           l=wsk->tab[wsk->dl_bloku-(pos-numer)-1];
                                                           //cout<<"tab l="<<l<<" nr="<<numer<<" pos="<<pos<<" dl_bloku: "<<wsk->dl_bloku<<endl;
                                                           suma+=l*l;
                                                           mk=2;
                                                           break;
                                 }
                                 if(mk==2) mk=0;
                                 if(numer>pos&&buf==false) wsk=wsk->next;
                                 }
                                 
                                 if(buf==true)
                                 {
                                                                    l=bufor[numer-pos-1];
                                                                    //cout<<"buf l="<<l<<" buf.pos="<<numer-1-pos<<" nr="<<numer<<" pos="<<pos<<endl;
                                                                    suma+=l*l; break;
                                 }
                                 }
               }
               return(sqrt(A.element(i,i)-suma));
       }
       //cout<<"\ni="<<i<<" j="<<j<<endl;
       if(i!=j)
       {
               numer=i-(A.l_elementow-(j-i));
               //cout<<"numer: "<<numer<<","<<j-i<<endl;
               for(k=1;k<=2*i-1;++k)
               {
                                   if(k%2==1) numer+=A.l_elementow-(j-i);
                                   else numer+=j-i;
                                   //cout<<"nr:"<<numer<<endl;
                                   while(1)
                                   {
                                   if(buf==false)
                                   {
                                                 
                                                  if(numer>pos&&mk==0) pos+=wsk->l_zer;
                                                  if(pos>=numer&&(mk==0||mk==1)&&numer!=0) {mk=1; l1=0; /*cout<<"zero"<<endl;*/ break;}
                                                  if(mk==1) mk=0;
                                                  if(wsk->next==NULL&&wsk->tab==NULL) {buf=true; /*cout<<"w buforze!"<<endl;*/ }
                                                  if(numer>pos&&mk==0&&buf==false) pos+=wsk->dl_bloku;
                                                  if(pos>=numer&&(mk==0||mk==2)&&buf==false)
                                                  {
                                                                                l=wsk->tab[wsk->dl_bloku-(pos-numer)-1];
                                                                                //cout<<"tab l="<<l<<" nr="<<numer<<" pos="<<pos<<" dl_bloku: "<<wsk->dl_bloku<<endl;
                                                                                if(k%2==1) l1=l;
                                                                                else suma+=l1*l;
                                                                                mk=2;
                                                                                break;
                                                  }
                                                  if(mk==2) mk=0;
                                                  if(numer>pos&&buf==false) {wsk=wsk->next; /*cout<<"next"<<endl;*/ }
                                   }
                                   if(buf==true)
                                   {
                                                l=bufor[numer-pos-1];
                                                //cout<<"buf l="<<l<<" buf.pos="<<numer-1-pos<<" nr="<<numer<<" pos="<<pos<<endl;
                                                if(k%2==1) l1=l;
                                                else suma+=l1*l;
                                                break;
                                   }
                                   }
               }
               return((A.element(i,j)-suma)/l);
               
       }
       
}
