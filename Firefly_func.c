#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <Firefly_func.h>
#include <Funkcje.h>

#define PODGLAD 1
#define MAX_D 1000
#define MAX_FFA 1000

int n_local;
int d_local;
int limit_generacji_local;
double alpha_local;
double beta_local;
double gamma_local;

int Index[MAX_FFA]; //swietliki posortowane po jakosci wyniku
double ffa[MAX_FFA][MAX_D]; //swietliki
double ffa_tmp[MAX_FFA][MAX_D]; //nowa generacja
double f[MAX_FFA]; //wartosc funkcji
double I[MAX_FFA]; //intensywnosc swiecenia
double nbest[MAX_FFA]; //najlepsze dotychczasowe rozwiazanie
double lb[MAX_D]; //dolna granica
double ub[MAX_D]; //gorna granica

double fbest; //najlepszy wynik

typedef double (*FunctionalCallback)(int dimension, double sol[MAX_D]);
FunctionalCallback funkcja = &funkcja1;

void inicjalizuj_ffa();
void inicjalizuj_funkcje(int numer_funkcji);
void inicjalizacja_zmiennych(int n, int d, int maxGeneracji, double alpha, double beta, double gamma);
void pokaz_ffa(int numer_generacji);
void sort_ffa();
void move_ffa();
void replace_ffa();
double reduce_alpha(double alpha, int numer_generacji);

void ffa_symulation(int n, int d, int maxGeneracji, double alpha, double beta, double gamma){
  printf("Funkcja otrzymala parametry:\n");
  printf(" n = %d\n", n);
  printf(" d = %d\n", d);
  printf(" ng = %d\n", maxGeneracji);
  printf(" a = %.2f\n", alpha);
  printf(" b = %.2f\n", beta);
  printf(" g = %.2f\n", gamma);

  inicjalizacja_zmiennych(n, d, maxGeneracji, alpha, beta, gamma);
  inicjalizuj_funkcje(1);

  int numer_generacji = 1; //licznik generacji

  srand(1);

  //inicjalizacja roju swietlikow
  inicjalizuj_ffa();
#ifdef PODGLAD
  pokaz_ffa(numer_generacji);
#endif

  int i,j;  
  while(numer_generacji <= limit_generacji_local){
    alpha_local = reduce_alpha(alpha_local, numer_generacji);

    for(i=0;i<n_local;i++){
      f[i] = funkcja(d_local, ffa[i]);
      I[i] = f[i];
    }

    sort_ffa();
    
    replace_ffa();

    for(i=0;i<d_local;i++){
      nbest[i] = ffa[0][i];
    }
    fbest = I[0];

    move_ffa();
  #ifdef PODGLAD
    pokaz_ffa(numer_generacji);
  #endif

    numer_generacji++;
  }

  printf("Koniec optymalizacji. Najlepszy wynik: %.4f\n",fbest);
  
  return;

}

void inicjalizacja_zmiennych(int n, int d, int maxGeneracji, double alpha, double beta, double gamma){
  n_local = n;
  d_local = d;
  limit_generacji_local = maxGeneracji;
  alpha_local = alpha;
  beta_local = beta;
  gamma_local = gamma;
}

void inicjalizuj_funkcje(int numer_funkcji){
  switch(numer_funkcji){
	case 1: funkcja = &funkcja1;
	default: funkcja = &funkcja1;
  }
} 

void inicjalizuj_ffa(){
  int i,j;
  double r;

  //inicjalizacja gornego i dolnedo ograniczenia
  for(i=0;i<d_local;i++){
    lb[i] = 0.0;
    ub[i] = 100.0;
  }

  for(i=0;i<n_local;i++){
    for(j=0;j<d_local;j++){
      r = ((double)rand()/((double)(RAND_MAX)+(double)(1)));
      ffa[i][j] = r * (ub[j] - lb[j]) + lb[j];
    }
    f[i] = 1.0;
    I[i] = f[i];
  }
}

void sort_ffa(){
  int i,j;
  //inicjalizacja indeksow
  for(i=0;i<n_local;i++){
    Index[i] = i;
  }

  //sortowanie babelkowe
  for(i=0;i<n_local;i++){
    for(j=i+1;j<n_local;j++){
      if(I[i] > I[j]){
        double z = I[i]; //zamiana atrakcyjnosci
        I[i] = I[j];
        I[j] = z;
        z = f[i]; //zamiana wartosci funkcji
        f[i] = f[j];
        f[j] = z;
        int k = Index[i]; //zamiana indeksow
        Index[i] = Index[j];
        Index[j] = k;
      }
    }
  }
}

void replace_ffa(){
  int i,j;
  //kopiowanie obecnej populacji do buforu tymczasowego
  for(i=0;i<n_local;i++){
    for(j=0;j<d_local;j++){
      ffa_tmp[i][j] = ffa[i][j];
    }
  }
  //???
  for(i=0;i<n_local;i++){
    for(j=0;j<d_local;j++){
      ffa[i][j] = ffa_tmp[Index[i]][j];
    }
  }
}

void findLimits(int k){
  int i;
  for(i=0;i<d_local;i++){
    if(ffa[k][i] < lb[i]) ffa[k][i] = lb[i];
    if(ffa[k][i] > ub[i]) ffa[k][i] = ub[i]; 
  }
}

void move_ffa(){
  int i,j,k;
  double scale;
  double r,beta;
  unsigned int* seed = (int*)1;

  for(i=0;i<n_local;i++){
    scale = abs(ub[i] - lb[i]);
    for(j=0;j<n_local;j++){
      //oblicz dlugosc r pomiedzy i-tym i j-tym swietlikiek
      r = 0.0;
      for(k=0;k<d_local;k++){
        r += (ffa[i][k] - ffa[j][k]) * (ffa[i][k] - ffa[j][k]);
      }
      r = sqrt(r);
      if(I[i] > I[j]){ //przesun swietlika i w kierunku swietlika j
        //zmodyfikuj atrakcyjność
        double beta0 = 1.0;
        beta = (beta0 - beta_local)*exp(-gamma_local*pow(r, 2.0)) + beta_local;
        //wygeneruj losowy wektor 
        for(k=0;k<d_local;k++){
          r = ((double)rand_r(seed) / ((double)(RAND_MAX) + (double)(1)));
          double tmpf = alpha_local * (r - 0.5) * scale;
          //utworz nowe rozwiazanie
          ffa[i][k] = ffa[i][k] * (1.0 - beta) + ffa_tmp[j][k] * beta + tmpf;
        }
      }
    }
    findLimits(i);
  }
}

double reduce_alpha(double alpha, int numer_generacji){
  double delta;
  delta = 1.0 - pow((pow(10.0, -4.0)/0.9), 1.0/(double) numer_generacji);
  return (1 - delta) * alpha;
}

void pokaz_ffa(int numer_generacji){
  printf("Podglad generacji numer: %d, najlepszy wynik: %.4f\n", numer_generacji, fbest);
}

