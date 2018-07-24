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



int main(int argc, char *argv[])
{
	FOWT fowt;
	ENVIR envir;
    
	if (argc == 1)
	{
		// std::string defaultFlNm = "Z:/METiS/Test/mainInputExample.txt";
		std::string defaultFlNm = "Test/mainInputExample.txt";
		
		std::cout << "No input file provided. Running default example located at ./Test/mainInputExample.txt.\n\n\n";
		IO::setInputFile(defaultFlNm);
	}

	else if (argc == 2)
	{
		std::cout << "Running METiS with file " << argv[1] << "\n\n\n";
		IO::setInputFile(argv[1]);
	}

	else {
		std::cout << "Please provide only one input file.\n";
		return 0;
	}
    

    //IO::setInputFile(inFlNm);
    IO::readInputFile(fowt, envir);

    IO::print2CheckVariables(fowt, envir);

	std::cout << "\n\n\nMETiS run completed. Press any key to exit.\n";
	char ii;
	std::cin >> ii;

	return 0;
}
