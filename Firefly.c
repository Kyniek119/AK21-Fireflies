
#include <stdio.h>
#include <Firefly_func.h>

int main(int argc, char* argv[]){
  int generation_number = 1; //licznik generacji
  int n = 100;
  int d = 1000;
  int MaxGeneracji = 20;
  int funkcja = 1;
  double alpha = 0.5;
  double beta = 0.2;
  double gamma = 1.0;

  int i = 0;
  void pomoc();
 
  //UI
  for(i=1;i<argc;i++){
    if((strncmp(argv[i], "-h", 2) == 0) || (strncmp(argv[i], "-?", 2)) == 0){
      pomoc();
      return 0;
    } else if(strncmp(argv[i], "-n", 2) == 0){ //ilosc swietlikow
      n = atoi(&argv[i][2]);
    } else if(strncmp(argv[i], "-d", 2) == 0){ //wymiar problemu
      d = atoi(&argv[i][2]);
    } else if(strncmp(argv[i], "-g", 2) == 0){ //maksymalna ilosc generacji
      MaxGeneracji = atoi(&argv[i][2]);
    } else if(strncmp(argv[i], "-a", 2) == 0){ //wspolczynnik alpha
      alpha = atof(&argv[i][2]);
    } else if(strncmp(argv[i], "-b", 2) == 0){ //wspolczynnik beta
      beta = atof(&argv[i][2]);
    } else if(strncmp(argv[i], "-c", 2) == 0){ //wspolczynnik gamma
      gamma = atof(&argv[i][2]);
    } else if(strncmp(argv[i], "-f", 2) == 0){ //numer funkcji
      funkcja = atoi(&argv[i][2]);
    } else {
      printf("Error: Niepoprawny parametr: %s .\nProgram zakonczyl dzilanie.\n", argv[i]);
      return -1;
    }
  }

  //inicjalizacja generatora liczb losowych
  ffa_symulation(n, d, MaxGeneracji, alpha, beta, gamma, funkcja);
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
}
