#include <math.h>
#include <Funkcje.h>

//funkcja "init1" inicjalizuje tablice wartosci dla funkcji "funkcja1", jest ona zapisywana w dostarczonym buforze "sol".
void init1(int dimension, double* sol){
  int i;
  for(i=0;i<dimension;i++){
    sol[i] = 3.0;
  }
}

//implementacja "1. Quadratic function" z funkcji testowych do zadań optymalizacji.
double funkcja1(int dimension, double* sol){
  double sum = 0.0;
  int i;
  for(i=2;i<dimension;i++){
    sum += 	100*(pow(sol[i],2.0)+pow(sol[i-1],2.0)) 
		+pow(sol[i-2],2.0); 
  }
  return sum;
}

//funkcja "init2" inicjalizuje tablice wartosci dla funkcji "funkcja2", jest ona zapisywana w dostarczonym buforze "sol".
void init2(int dimension, double* sol){
  int i;
  for(i=0;i<dimension;i++){
    if(i%2 == 1) sol[i] = -1.0;
    else sol[i] = -3.0;
  }
}

//implementacja "2. Woods function" z funkcji testowych do zadań optymalizacji.
double funkcja2(int dimension, double* sol){
  double sum = 0.0;
  int i;
  int suma_do = (int)(dimension/4);
  for(i=1;i<=suma_do;i++){
    sum+= 	100*pow((sol[4*i-2]-pow(sol[4*i-3],2.0)),2.0)
		+pow((1-sol[4*i-3]),2.0)
		+90*pow((sol[4*i]-pow(sol[4*i-1],2.0)),2.0)
		+pow((1-sol[4*i-1]),2.0)
		+10*pow((sol[4*i-2]+sol[4*i]-2),2.0)
		+0.1*pow((sol[4*i-2]-sol[4*i]),2.0); 
  }
  return sum;
}
