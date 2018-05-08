/* PASSOS

1) Definir como serão os inputs e outputs na forma final (com aerodinamica e etc)
2) Pensar em como vou lidar com arrays, se crio um classe de matrizes ou uso std::array e std::vector. A classe de matrizes é interessante para dar overload nos operadores e usar como se fosse MATLAB
3) Escrever a função/classe ou sei lá oq pra input. 
4) Implementar a parte hidrodinâmica aos poucos. Começar pelo cálculo da velocidade em um vetor de pontos. Depois nas funções de pegar a posição dos nós do cilindro. E assim por diante

*/


// include as parada

int main(){
	Fowt fowt;   // class with the properties of the structural system (floater, tower, rotor, blades, etc)
	Envir envir; // class with the environmental conditions (wave, wind, current) AND numerical parameters (time step, total simulation time, etc)

	bool readOK = dataRead( Fowt &fowt, Envir &envir );
	if (~readOK) {
		// manda um erro aqui
	};


	// Faz uma classe que faz a análise no domínio do tempo. Dessa forma, posso adicionar outros tipos de análise de forma facil.
	// Fazer com que as outras classes sejam 'friends' dela pra ter acesso a seus members
	// Deixar o método de integração como uma função a parte, talvez como parte da classe, para ser facilmente mudável.
	// Fazer uma função pra printar pra um txt que é chamada ao final de cada time step. Talvez fazer uma classe pra poder printar diferentes partes do programa (cabeçalho, finalização, etc)
	timeDomainAnalysis(fowt, envir);


	// Pra printar, posso fazer um string ao qual vou adicionando as coisas que vou imprimir


	return 0;
}
