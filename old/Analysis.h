//class Analysis{
//
//public:
//	void timeDomain( FOWT &fowt, ENVIR &envir, NUMP &nump );
//
//	void integrateRK4( FOWT &fowt, ENVIR &envir, NUMP &nump  );
//};

// Essa aqui é a principal, faz a análise completa no dominio do tempo
void timeDomainAnalysis( Fowt &fowt, Envir &envir );

// Printa um arquivo txt com licença, versão e resumo do que foi rodado
void printSummarizedOutput( Fowt &fowt, Envir &envir );

// Pega a FOWT em um time step, faz a integração por runge-kutta e atualiza ela pro próximo time-step
void integrateRK4( Fowt &fowt, Envir &envir );

// Printar um header no txt pra saber o que é cada coluna
void printHeader2txt( Fowt &fowt, Envir &envir );
void print2txt( Fowt &fowt, Envir &envir );



// Depois tem que passar pra um cpp file
// Faz uma classe que faz a análise no domínio do tempo. Dessa forma, posso adicionar outros tipos de análise de forma facil.
// Fazer com que ela seja 'friend' das outras classes pra ter acesso a seus members
// Deixar o método de integração como uma função a parte, talvez como parte da classe, para ser facilmente mudável.
// Fazer uma função pra printar pra um txt que é chamada ao final de cada time step. Talvez fazer uma classe pra poder printar diferentes partes do programa (cabeçalho, finalização, etc)
void timeDomainAnalysis( Fowt &fowt, Envir &envir ){

	printHeader2txt( fowt, envir );
	print2txt( fowt, envir ); // Print instant 0

	FOWT fowt_RK = fowt; // Acho que tenho que dar overload no assignment operator. Ler sobre copy constructores e overloading assignment operator

	while ( envir.getCurrentTime() < envir.getFinalTime() ){
		integrateRK4( fowt_RK, envir ) // Tenho que pensar ainda em como fazer essa aqui.
		print2txt( fowt_RK, envir );
	}

}






