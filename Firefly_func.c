#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <Firefly_func.h>
#include <Funkcje.h>

#define PODGLAD 0
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
double ffa_tmp[MAX_FFA][MAX_D]; //bufor pomocniczy
double f[MAX_FFA]; //wartosc funkcji
double I[MAX_FFA]; //intensywnosc swiecenia
double nbest[MAX_FFA]; //najlepsze dotychczasowe rozwiazanie
double lb[MAX_D]; //dolna granica
double ub[MAX_D]; //gorna granica

double fbest; //najlepszy wynik
unsigned int seed = 1;


typedef double (*FunctionalCallback)(int dimension, double sol[MAX_D]);
FunctionalCallback funkcja = &funkcja2;

typedef void (*InitializationCallback)(int dimension, double sol[MAX_D]);
InitializationCallback inicjalizacja_danych = &init2;

void inicjalizuj_ffa();
void inicjalizuj_funkcje(int numer_funkcji);
void inicjalizacja_zmiennych(int n, int d, int maxGeneracji, double alpha, double beta, double gamma);
void pokaz_ffa(int numer_generacji);
void sort_ffa();
void move_ffa();
void replace_ffa();
double reduce_alpha(double alpha, int numer_generacji);

void ffa_symulation(int n, int d, int maxGeneracji, double alpha, double beta, double gamma, int numer_funkcji){

  inicjalizacja_zmiennych(n, d, maxGeneracji, alpha, beta, gamma);
  inicjalizuj_funkcje(numer_funkcji);

  int numer_generacji = 1; //licznik generacji

  srand(1);

  //inicjalizacja roju swietlikow
  inicjalizuj_ffa();
#ifdef PODGLAD
  pokaz_ffa(numer_generacji);
#endif

  //glowna petla, wykonywana dla kazdej generacji
  int i,j;  
  while(numer_generacji <= limit_generacji_local){
    //alpha_local = reduce_alpha(alpha_local, numer_generacji);

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
		inicjalizacja_danych = &init1;
        case 2: funkcja = &funkcja2;
		inicjalizacja_danych = &init2;
	default: funkcja = &funkcja1;
                 inicjalizacja_danych = &init1;
  }
} 

void inicjalizuj_ffa(){
  int i,j;
  double r;
  double dane[MAX_D];
  inicjalizacja_danych(d_local, dane);

  for(i=0;i<d_local;i++){
    lb[i] = -10.0;

    ub[i] = 100.0;
  }

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

  for(i=0;i<n_local;i++){
    scale = abs(ub[i] - lb[i]);
    for(j=0;j<n_local;j++){
      //oblicz dlugosc r pomiedzy i-tym i j-tym swietlikiek
      r = 0.0;
      for(k=0;k<d_local;k++){
        r += (ffa[i][k] - ffa[j][k]) * (ffa[i][k] - ffa[j][k]);
      }
      r = sqrt(r);
      if(I[i] >= I[j]){ //przesun swietlika i w kierunku swietlika j
        //zmodyfikuj atrakcyjność
        double beta0 = 1.0;
        beta = (beta0 - beta_local)*exp(-gamma_local*pow(r, 2.0)) + beta_local;
        //wygeneruj losowy wektor 
        for(k=0;k<d_local;k++){
          r = ((double)rand_r(&seed) / ((double)(RAND_MAX) + (double)(1)));
          double tmpf = alpha_local * (r - 0.5) * scale;
		//printf("ub[i]: %.4f, lb[i]: %.4f, Scale: %.4f, beta: %.4f, r: %.4f, tmp_f: %.4f\n", ub[i], lb[i], scale, beta, r, tmpf);
          //utworz nowe rozwiazanie
          ffa[i][k] = ffa[i][k] + beta * (ffa_tmp[j][k] - ffa[i][k]) + tmpf;
        }
      }
    }
    findLimits(i);
  }
}

double reduce_alpha(double alpha, int numer_generacji){
  double delta;
  delta = 1.0 - pow((pow(10.0, -4.0)/0.9), 1.0/(double) numer_generacji);
	printf("Nowa alpha: %.4f\n", (1-delta)*alpha);
  return (1 - delta) * alpha;
}

void pokaz_ffa(int numer_generacji){
  printf("Podglad generacji numer: %d, najlepszy wynik: %.4f\n", numer_generacji, fbest);
}

