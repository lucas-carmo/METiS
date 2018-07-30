/* To do: 
Trocar aquele header do METiS por uma função, dai uso pra printar pro console e pro summary file

Depois que criar todas as funcoes de leitura separadamente, ver oq da pra transformar em uma funcao comum a elas.
   Ex: Uma funcao checkInput pra verificar se tem o numero certo de elementos pra leitura, se as conversoes pros tipos numericos
   do input deram certo, coisa do tipo. Assim facilita a manutencao, pq se quiser mudar mensagem de erro ou coisa do tipo, mudo
   em um so lugar

   TODO: No IO.cpp, adicionar um warning pra caso a keyword seja desconhecida. Provavelmente fica melhor trocar aqueles ifs por cases
*/

#include <iostream>
#include <chrono> // For std::chrono functions. It is useful to time the code and verify whether one method or another will be more performant
#include <armadillo> // Linear algebra library with usage similar to MATLAB
#include <string>
#include <stdexcept> // For std::exception

#include "IO.h" // IO is a pure static class that manages input/output
#include "FOWT.h"
#include "ENVIR.h"



#ifdef __unix__
# include <unistd.h>
#elif defined _WIN32
# include <windows.h>
#define usleep(x) Sleep((x)/1000)
// #define mkdir(x) CreateDirectory((x), NULL)
#endif

int main(int argc, char *argv[])
{
    
	if (argc != 2)
	{	
		// Acho melhor throw exception, escrever pro console, escrever pro log file,
		// esperar input do usuario e daí terminar
		std::cout << "Please provide only one input file.\n";
		return 0;
	}
    
	try{
		FOWT fowt;
		ENVIR envir;

		std::cout << IO::METiS_Header() << std::endl << std::endl;
		std::cout << "Running METiS with file " << argv[1] << "\n\n\n";
		IO::setFiles(argv[1]); // Set paths to input file and output files
		IO::readInputFile(fowt, envir); // Read data from input file to fowt and envir
		IO::printSumFile(fowt, envir);

		std::cout << '\n';
		for (int ii = 1; ii <= 100; ++ii)
		{		
			std::cout << ii << "%" << '\r';		
			fflush(stdout);
			usleep(10000);
		}		
	}
	catch(std::exception &exception)
	{
		IO::writeErrorMessage( "Standard exception: " + std::string(exception.what()) );
	}


	std::cin.sync();
	std::cout << "\n\n\nPress enter to exit.\n";
	std::cin.get();

	return 0;
}
