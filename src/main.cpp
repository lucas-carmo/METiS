/* To do: Depois que criar todas as funcoes de leitura separadamente, ver oq da pra transformar em uma funcao comum a elas.
   Ex: Uma funcao checkInput pra verificar se tem o numero certo de elementos pra leitura, se as conversoes pros tipos numericos
   do input deram certo, coisa do tipo. Assim facilita a manutencao, pq se quiser mudar mensagem de erro ou coisa do tipo, mudo
   em um so lugar

   TODO: Talvez fazer dos nodes uma classe a parte, pra deixar os IDs como int e não como double

   TODO: No IO.cpp, adicionar um warning pra caso a keyword seja desconhecida. Provavelmente
   fique melhor trocar aqueles ifs por cases
*/

#include <iostream>
#include <chrono> // For std::chrono functions. It is useful to time the code and verify whether one method or another will be more performant
#include <armadillo> // Linear algebra library with usage similar to MATLAB
#include <string>

#include "IO.h" // IO is a pure static class that manages input/output
#include "FOWT.h"
#include "ENVIR.h"


#ifdef __unix__
# include <unistd.h>
#elif defined _WIN32
# include <windows.h>
#define usleep(x) Sleep((x)/1000)
#endif

int main(int argc, char *argv[])
{
	FOWT fowt;
	ENVIR envir;
    
	if (argc < 2)
	{	
		// Acho melhor throw exception, escrever pro console, escrever pro log file,
		// esperar input do usuario e daí terminar
		std::cout << "Please provide at least one input file.\n";
		return 0;
	}
    
	for (int ii = 1; ii < argc; ++ii)
	{		
		std::cout << "Running METiS with file " << argv[ii] << "\n\n\n";
		IO::setInputFile(argv[ii]);
		IO::readInputFile(fowt, envir);
		IO::print2CheckVariables(fowt, envir);


		std::cout << '\n';
		for (int ii = 1; ii < 100; ++ii)
		{		
			std::cout << ii << "%" << '\r';		
			fflush(stdout);
			usleep(50000);
		}

	}

	std::cin.sync();
	std::cout << "\n\n\nMETiS run completed. Press enter to exit.\n";
	std::cin.get();

	return 0;
}
