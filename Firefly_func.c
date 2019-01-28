#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <limits.h>
#include <omp.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include <Firefly_func.h>
#include <Funkcje.h>

#define MAX_D 1000
#define MAX_FFA 1000

//lokalne wartosci parametrow
int ilosc_swietlikow;
int wymiar_problemu;
int limit_generacji;
double alpha_local;
double beta_local;
double gamma_local;
bool multithreading_local;

unsigned int seed;

double ffa[MAX_FFA][MAX_D]; //swietliki
double ffa_prev_gen[MAX_FFA][MAX_D]; //zmienna do przechowywania poprzedniej generacji.
double f[MAX_FFA]; //wartosc funkcji
int Index[MAX_FFA]; //index wiazacy wartosc funkcji ze swietlikiem, pozwoli uniknac kopiowania parametrow swietlikow

double fbest; //najlepszy wynik obecnej populacji
double global_best = DBL_MAX; //najlepszy wynik globalnie
double global_best_param[MAX_D];

double start, end; //zmienne do pomiaru czasu
double czas_wykonania;
FILE* plik_wynikowy = NULL;


typedef double (*FunctionalCallback)(int dimension, double sol[MAX_D]);
FunctionalCallback funkcja = &funkcja1;

typedef void (*InitializationCallback)(int dimension, double sol[MAX_D]);
InitializationCallback inicjalizacja_danych = &init1;

//funkcje pomocnicze
void inicjalizuj_ffa();
void inicjalizuj_funkcje(int numer_funkcji);
void inicjalizacja_zmiennych(int n, int d, int maxGeneracji, double alpha, double beta, double gamma, char* nazwa_pliku_wynikowego, bool multithreading);
void pokaz_ffa(int numer_generacji);
void pokaz_parametry_maszyny();
void pokaz_rozwiazanie();
void replace_ffa(double old[MAX_FFA][MAX_D], double new[MAX_FFA][MAX_D]);
void sort_ffa();
void move_ffa();

void ffa_symulation(int n, int d, int g, double alpha, double beta, double gamma, int numer_funkcji, char* nazwa_pliku_wynikowego, bool multithreading){

  inicjalizacja_zmiennych(n, d, g, alpha, beta, gamma, nazwa_pliku_wynikowego, multithreading);
  inicjalizuj_funkcje(numer_funkcji);

  int numer_generacji = 1; //licznik generacji

  //inicjalizacja roju swietlikow
  inicjalizuj_ffa();
  
  //pokaz_ffa(numer_generacji);
  pokaz_parametry_maszyny();

  //glowna petla, wykonywana dla kazdej generacji
  int i,j, index_best;  

  start = omp_get_wtime();
  while(numer_generacji <= limit_generacji){

    #pragma omp pararell for if(multithreading_local)
    for(i=0;i<ilosc_swietlikow;i++){
      f[i] = funkcja(wymiar_problemu, ffa[Index[i]]);
    }
 
    index_best = find_index_of_local_best(f);

    fbest = f[index_best];
    if(fbest < global_best){
      global_best = fbest;
      memcpy(global_best_param, ffa[Index[index_best]], sizeof(double) * MAX_D);
    }

    move_ffa();

    replace_ffa(ffa, ffa_prev_gen); //

    pokaz_ffa(numer_generacji);

    numer_generacji++;
  }
  end = omp_get_wtime();
  czas_wykonania = end - start;
  
  pokaz_rozwiazanie();

  printf("Koniec optymalizacji. Najlepszy wynik: %.4f, w czasie: %5.1fms\n",global_best, czas_wykonania*1000);

  if(plik_wynikowy != NULL){
    fprintf(plik_wynikowy, "Koniec optymalizacji. Najlepszy wynik: %.4f, w czasie: %5.1fms\n",global_best, czas_wykonania*1000);
  }
  fclose(plik_wynikowy);  
  return;

}

void inicjalizacja_zmiennych(int n, int d, int g, double alpha, double beta, double gamma, char* nazwa_pliku_wynikowego, bool multithreading){
  ilosc_swietlikow = n;
  wymiar_problemu = d;
  limit_generacji = g;
  alpha_local = alpha;
  beta_local = beta;
  gamma_local = gamma;
  fbest = DBL_MAX;

  seed = time(NULL) ^ getpid();
  multithreading_local = multithreading;

  if( nazwa_pliku_wynikowego != NULL){
    plik_wynikowy = fopen(nazwa_pliku_wynikowego, "w");
  }
}

void inicjalizuj_funkcje(int numer_funkcji){
  switch(numer_funkcji){
	case 1: funkcja = &funkcja1;
		inicjalizacja_danych = &init1;
		break;
        case 2: funkcja = &funkcja2;
		inicjalizacja_danych = &init2;
		break;
        case 3: funkcja = &funkcja3;
		inicjalizacja_danych = &init3;
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
  inicjalizacja_danych(wymiar_problemu, dane);

  for(i=0;i<ilosc_swietlikow;i++){
    for(j=0;j<wymiar_problemu;j++){
      ffa[i][j] = dane[j];
      ffa_prev_gen[i][j] = dane[j];
    }
    f[i] = 1.0;
    Index[i] = i;
  }
}

int find_index_of_local_best(double f[MAX_FFA]){
  int i, result;
  double local_best = DBL_MAX;
   for(i=0; i<ilosc_swietlikow; i++){
    if(local_best > f[i]){
      local_best = f[i];
      result = i;
    } 
  } 
  return result;
}

void replace_ffa(double old[MAX_FFA][MAX_D], double new[MAX_FFA][MAX_D]){
  int i;  
  for(i=0;i<ilosc_swietlikow;i++){
    memcpy(new[Index[i]], old[Index[i]], sizeof(double) * MAX_D);
  }
}

void move_ffa(){
  int i,j,k;
  double r,beta;

  #pragma omp parallel for private(r, beta, j, k) if(multithreading_local)
  for(i=0;i<ilosc_swietlikow;i++){
    for(j=0;j<ilosc_swietlikow;j++){
      r = 0.0;
      //oblicz wypadkowa dlugosc r pomiedzy i-tym i j-tym swietlikiek
      for(k=0;k<wymiar_problemu;k++){
        r += (ffa_prev_gen[Index[i]][k] - ffa_prev_gen[Index[j]][k]) * (ffa_prev_gen[Index[i]][k] - ffa_prev_gen[Index[j]][k]);
      }
      r = sqrt(r);
      //przesun swietlika z najlepszym rozwiazaniem
      if(f[i] == fbest){
        for(k=0;k<wymiar_problemu;k++){
          //wygeneruj losowy wektor;
          r = ((double)rand_r(&seed) / ((double)(RAND_MAX) + (double)(1)));
          double u = alpha_local * (r - 0.5);
          //utworz nowe rozwiazanie
          ffa[Index[i]][k] = ffa_prev_gen[Index[i]][k] + u;
        }
      } else if(f[i] > f[j]){ //przesun swietlika i w kierunku swietlika j
        //zmodyfikuj atrakcyjność
        beta = beta_local*exp(-gamma_local*pow(r, 2.0));

        //wygeneruj losowy wektor 
        for(k=0;k<wymiar_problemu;k++){
          r = ((double)rand_r(&seed) / ((double)(RAND_MAX) + (double)(1)));
          double u = alpha_local * (r - 0.5);
          //utworz nowe rozwiazanie
          ffa[Index[i]][k] = ffa[Index[i]][k] + beta * (ffa_prev_gen[Index[j]][k] - ffa_prev_gen[Index[i]][k]) + u;
        }
      }
    }
  }
}

void pokaz_ffa(int numer_generacji){
  printf("Podglad generacji numer: %d, najlepszy wynik: %.4f\n", numer_generacji, fbest);
  if(plik_wynikowy != NULL){
    fprintf(plik_wynikowy, "Podglad generacji numer: %d, najlepszy wynik: %.4f\n", numer_generacji, fbest);
  }
}

void pokaz_parametry_maszyny(){
  printf("OMP_DYNAMIC = %d\n",omp_get_dynamic());
  printf("OMP_NESTED = %d\n",omp_get_nested());
  printf("OMP_NUM_THREADS = %d\n",omp_get_max_threads());
  //printf("OMP_SCHEDULE = %s\n",omp_get_schedule());
  printf("OMP_THREAD_LIMIT = %d\n",omp_get_thread_limit());
  printf("OMP_DYNAMIC = %d\n",omp_get_dynamic());
  printf("Multithreading = %d\n", multithreading_local);
  
  if(plik_wynikowy != NULL){
    fprintf(plik_wynikowy, "OMP_DYNAMIC = %d\n",omp_get_dynamic());
    fprintf(plik_wynikowy, "OMP_NESTED = %d\n",omp_get_nested());
    fprintf(plik_wynikowy, "OMP_NUM_THREADS = %d\n",omp_get_max_threads());
    //fprintf(plik_wynikowy, "OMP_SCHEDULE = %s\n",omp_get_schedule());
    fprintf(plik_wynikowy, "OMP_THREAD_LIMIT = %d\n",omp_get_thread_limit());
    fprintf(plik_wynikowy, "OMP_DYNAMIC = %d\n",omp_get_dynamic());
    fprintf(plik_wynikowy, "Multithreading = %d\n", multithreading_local);
  }
}

void pokaz_rozwiazanie(){
  int i, count;
  count = (int)(wymiar_problemu/10);
  printf("Parametry po 10 w wierszu.\n");
  for(i=0;i<wymiar_problemu/10;i++){
    printf("%f, %f, %f, %f, %f, %f, %f, %f, %f, %f.\n",
      global_best_param[i],
      global_best_param[i+1],
      global_best_param[i+2],
      global_best_param[i+3],
      global_best_param[i+4],
      global_best_param[i+5],
      global_best_param[i+6],
      global_best_param[i+7],
      global_best_param[i+8],
      global_best_param[i+9]);

    if(plik_wynikowy != NULL){
      fprintf(plik_wynikowy, "%f, %f, %f, %f, %f, %f, %f, %f, %f, %f.\n",
        global_best_param[i],
        global_best_param[i+1],
        global_best_param[i+2],
        global_best_param[i+3],
        global_best_param[i+4],
        global_best_param[i+5],
        global_best_param[i+6],
        global_best_param[i+7],
        global_best_param[i+8],
        global_best_param[i+9]);
    }  
  }
  for(i=count*10;i<wymiar_problemu;i++){
    printf("%f, ",global_best_param[i]);
    if(plik_wynikowy != NULL){
      fprintf(plik_wynikowy,"%f, ",global_best_param[i]);
    }
  }
  printf("\n");
  if(plik_wynikowy != NULL){
      fprintf(plik_wynikowy, "\n");
  }
 
}

