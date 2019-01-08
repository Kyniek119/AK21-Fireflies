
#include <stdio.h>
#include <Firefly_func.h>

int main(int argc, char* argv[]){
  int n = 100;
  int d = 100;
  int g = 20;
  int funkcja = 1;
  double alpha = 0.5;
  double beta = 1;
  double gamma = 0.01; //duza wartosc powoduje berdzo losowe przeszukiwanie przestrzeni
  char* nazwa_pliku_wynikowego = "Firefly_output.txt";

  int i = 0;
  void pomoc();
 
  //UI
  for(i=1;i<argc;i++){
    if((strncmp(argv[i], "-h", 2) == 0) || (strncmp(argv[i], "-?", 2)) == 0){
      pomoc();
      return 0;
    } else if(strncmp(argv[i], "-n", 2) == 0){ //ilosc swietlikow
      n = atoi(&argv[i][2]);
      if(n < 1 || n > 1000) { printf("Niepoprawna wartość parametru n: %d. Dostepne wartosci to <1, 1000>.\n",n); return -2;}
    } else if(strncmp(argv[i], "-d", 2) == 0){ //wymiar problemu
      d = atoi(&argv[i][2]);
      if(d < 1 || d > 1000) { printf("Niepoprawna wartość parametru d: %d. Dostepne wartosci to <1, 1000>.\n",d); return -2;}
    } else if(strncmp(argv[i], "-g", 2) == 0){ //maksymalna ilosc generacji
      g = atoi(&argv[i][2]);
      if(g < 1) { printf("Niepoprawna wartość parametru g: %d. Dostepne wartosci >1.\n",g); return -2;}
    } else if(strncmp(argv[i], "-a", 2) == 0){ //wspolczynnik alpha
      alpha = atof(&argv[i][2]);
      if(alpha < 0 || alpha > 1)  { printf("Niepoprawna wartość parametru alpha: %f. Dostepne wartosci <0, 1>.\n",alpha); return -2;}
    } else if(strncmp(argv[i], "-b", 2) == 0){ //wspolczynnik beta
      beta = atof(&argv[i][2]);
      if(beta < 0 || beta > 1)  { printf("Niepoprawna wartość parametru beta: %f. Dostepne wartosci <0, 1>.\n",beta); return -2;}
    } else if(strncmp(argv[i], "-c", 2) == 0){ //wspolczynnik gamma
      gamma = atof(&argv[i][2]);
      if(gamma < 0 || gamma > 10)  { printf("Niepoprawna wartość parametru gamma: %f. Dostepne wartosci <0, 10>.\n",gamma); return -2;}
    } else if(strncmp(argv[i], "-f", 2) == 0){ //numer funkcji
      funkcja = atoi(&argv[i][2]);
      if(funkcja < 0 || funkcja > 3)  { printf("Niepoprawna wartość parametru funkcja: %d. Dostepne wartosci (1, 2, 3).\n",funkcja); return -2;}
    } else if(strncmp(argv[i], "-o", 2) == 0){ //nazwa pliku wynikowego
      nazwa_pliku_wynikowego = &argv[i][2];
    } else {
      printf("Error: Niepoprawny parametr: %s .\nProgram zakonczyl dzilanie.\n", argv[i]);
      return -1;
    }
  }

  //inicjalizacja generatora liczb losowych
  ffa_symulation(n, d, g, alpha, beta, gamma, funkcja, nazwa_pliku_wynikowego);
  return(0);
}

void pomoc(){
  printf("Syntax:\n");
  printf(" Firefly [-h|-?] [-n] [-d] [-g] [-a] [-b] [-c]\n");
  printf("  Gdzie: -n = ilosc swietlikow (1-1000)\n");
  printf("         -d = wymiar problemu (1-1000)\n");
  printf("         -g = maksymalna ilosc generacji\n");
  printf("         -f = numer funkcji z zestawu (1-2)\n");
  printf("         -a = wspolczynnik alpha\n");
  printf("         -b = wspolczynnik beta\n");
  printf("         -c = wspolczynnik gamma\n");
  printf("         -o = nazwa pliku wynikowego\n");
}
