#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <Firefly_func.h>
#include <Funkcje.h>

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
double f[MAX_FFA]; //wartosc funkcji
double I[MAX_FFA]; //intensywnosc swiecenia
double nbest[MAX_FFA]; //najlepsze dotychczasowe rozwiazanie

double fbest; //najlepszy wynik obecnej populacji
double global_best = DBL_MAX; //najlepszy wynik globalnie

unsigned int seed = 1;
clock_t start, end; //zmienne do pomiaru czasu
double czas_wykonania;


typedef double (*FunctionalCallback)(int dimension, double sol[MAX_D]);
FunctionalCallback funkcja = &funkcja2;

typedef void (*InitializationCallback)(int dimension, double sol[MAX_D]);
InitializationCallback inicjalizacja_danych = &init2;

//funkcje pomocnicze
void inicjalizuj_ffa();
void inicjalizuj_funkcje(int numer_funkcji);
void inicjalizacja_zmiennych(int n, int d, int maxGeneracji, double alpha, double beta, double gamma);
void pokaz_ffa(int numer_generacji);
void sort_ffa();
void move_ffa();
void replace_ffa();

void ffa_symulation(int n, int d, int maxGeneracji, double alpha, double beta, double gamma, int numer_funkcji){

  inicjalizacja_zmiennych(n, d, maxGeneracji, alpha, beta, gamma);
  inicjalizuj_funkcje(numer_funkcji);

  int numer_generacji = 1; //licznik generacji

  srand(1);

  //inicjalizacja roju swietlikow
  inicjalizuj_ffa();
  
  pokaz_ffa(numer_generacji);

  //glowna petla, wykonywana dla kazdej generacji
  int i,j;  

  start = clock();
  while(numer_generacji <= limit_generacji_local){

    for(i=0;i<n_local;i++){
      f[i] = funkcja(d_local, ffa[i]);
      I[i] = f[i];
    }
    sort_ffa();

    for(i=0;i<d_local;i++){
      nbest[i] = ffa[0][i];
    }
    fbest = I[0];
    if(fbest < global_best)
      global_best = fbest;

    move_ffa();

    pokaz_ffa(numer_generacji);

    numer_generacji++;
  }
  end = clock();
  czas_wykonania = ((double) (end - start))/CLOCKS_PER_SEC;
  printf("Koniec optymalizacji. Najlepszy wynik: %.4f, w czasie: %5.1fms\n",global_best, czas_wykonania * 1000);
  
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
		inicjalizacja_danych = &init1;
		break;
        case 2: funkcja = &funkcja2;
		inicjalizacja_danych = &init2;
		break;
	default: funkcja = &funkcja1;
                 inicjalizacja_danych = &init1;
		break;
  }
} 

void inicjalizuj_ffa(){
  int i,j;
  double r;
  double dane[MAX_D];
  inicjalizacja_danych(d_local, dane);

  for(i=0;i<n_local;i++){
    for(j=0;j<d_local;j++){
      ffa[i][j] = dane[j];
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

void move_ffa(){
  int i,j,k;
  double r,beta;

  for(i=0;i<n_local;i++){
    for(j=0;j<n_local;j++){
      //oblicz dlugosc r pomiedzy i-tym i j-tym swietlikiek
      r = 0.0;
      for(k=0;k<d_local;k++){
        r += (ffa[i][k] - ffa[j][k]) * (ffa[i][k] - ffa[j][k]);
      }
      r = sqrt(r);
      if(I[i] < I[j]){ //przesun swietlika i w kierunku swietlika j

        //zmodyfikuj atrakcyjność
        beta = beta_local*exp(-gamma_local*pow(r, 2.0));

        //wygeneruj losowy wektor 
        for(k=0;k<d_local;k++){
          r = ((double)rand_r(&seed) / ((double)(RAND_MAX) + (double)(1)));
          double u = alpha_local * (r - 0.5);

          //utworz nowe rozwiazanie
          ffa[i][k] = ffa[i][k] + beta * (ffa[j][k] - ffa[i][k]) + u;
        }
      }
      //przesun swietlika z najlepszym rozwiazaniem
      if(f[i] == fbest){
        for(k=0;k<d_local;k++){
          //wygeneruj losowy wektor
          r = ((double)rand_r(&seed) / ((double)(RAND_MAX) + (double)(1)));
          double u = alpha_local * (r - 0.5);
          //utworz nowe rozwiazanie
          ffa[i][k] = ffa[i][k] + u;
        }
      }
    }
  }
}

void pokaz_ffa(int numer_generacji){
  printf("Podglad generacji numer: %d, najlepszy wynik: %.4f\n", numer_generacji, fbest);
}

