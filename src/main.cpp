/* To do: Depois que criar todas as funcoes de leitura separadamente, ver oq da pra transformar em uma funcao comum a elas.
   Ex: Uma funcao checkInput pra verificar se tem o numero certo de elementos pra leitura, se as conversoes pros tipos numericos
   do input deram certo, coisa do tipo. Assim facilita a manutencao, pq se quiser mudar mensagem de erro ou coisa do tipo, mudo
   em um so lugar
*/


#include <iostream>
#include <chrono> // For std::chrono functions. It is useful to time the code and verify whether one method or another will be more performant
#include <armadillo> // Linear algebra library with usage similar to MATLAB
#include <string>

#include "FOWT.h"
#include "IO.h" // IO is a pure static class that manages input/output
#include "ENVIR.h"



int main()
{
	FOWT fowt;
	ENVIR envir;
    
    //std::string inFlNm = "Test/mainInputExample.txt";
	std::string inFlNm = "Z:/METiS/Test/mainInputExample.txt";

    IO::setInputFilePath(inFlNm);
    IO::readInputFile(fowt, envir);

	int ii;
	std::cin >> ii;

	return 0;
}
