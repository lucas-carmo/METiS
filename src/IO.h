#pragma once

#include <string>
#include "FOWT.h"
#include "ENVIR.h"


class IO
{
private:
	static std::string m_header;
	static bool m_shouldWriteHeader;

	static std::string m_outputTimeStep; // String with the data that is output at each time step (FOWT position, hydro force components, anything that is a function of time)
	static bool m_shouldWriteOutputTimeStep;

public:
	static void readInputFile(const std::string &inFlNm, FOWT &fowt, ENVIR &envir);
	//static bool checkData(); // Depende do tipo de analise a ser feita

	static void print2outFile(const std::string &outFlNm);
};


