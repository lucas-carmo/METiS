/* To do: Depois que criar todas as funcoes de leitura separadamente, ver oq da pra transformar em uma funcao comum a elas.
   Ex: Uma funcao checkInput pra verificar se tem o numero certo de elementos pra leitura, se as conversoes pros tipos numericos
   do input deram certo, coisa do tipo. Assim facilita a manutencao, pq se quiser mudar mensagem de erro ou coisa do tipo, mudo
   em um so lugar

   TODO: Criar uma função dentro da classe IO pra leitura de linha, pra colocar o getLine e o update line number juntos e não ter que ficar lembrando disso

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



int main()
{
	FOWT fowt;
	ENVIR envir;
    
    std::string inFlNm = "Test/mainInputExample.txt";
	// std::string inFlNm = "Z:/METiS/Test/mainInputExample.txt";

    IO::setInputFilePath(inFlNm);
    IO::readInputFile(fowt, envir);

	int ii;
	std::cin >> ii;

	return 0;
}
