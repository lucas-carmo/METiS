#pragma once

class Fowt{
private:
	Floater m_floater;
	Tower m_tower;
	Nacelle m_nacelle;
	Rotor m_rotor;

	int nBlade;
    std::vector<Blade, nBlade> m_blade; // fazer um array da classe Blade com tamanho nBlades. Ver como fazer isso direitinho. Talvez valesse a pena deixar um default com 3 pra ganhar tempo de processamento
			
	std::array<double, 3> m_cog; // CoG do conjunto. Por ortogonalidade, isso tem que ser uma função. Mesmo porque não uso esse aqui nas equacoes do movimento.
								 //  Matriz de massa do sistema completo -- Ainda não sei como fazer isso. Por ortogonalidade, isso também tem que ser uma função.
	std::array<double, 6> m_acc, m_vel, m_pos; // fazer essas aqui serem um array de 6 componentes (x, y, z, Rx, Ry, Rz)

public:
	std::array<double, 6> calcHydroForce( Envir &envir ); // fazer com que as funcoes de forca retornem um array de 6 componentes.
	std::array<double, 6> calcAeroForce( Envir &envir );
//	std::array<double, 6> calcMoorForce();
	std::array<double, 6> calcTotalForce( Envir &envir );
};



// vai pro seu próprio header
class Floater {
private:
	int m_numCylp;
	std::vector<Cylp, m_numCylp> m_cylp;

	int m_nRetp;
	std::vector<Retp, m_numRetp> m_retp;

public:
	int getNumCylp() { return m_numCylp; }
	int getNumRetp() { return m_numRetp; }

};

// vai pro seu próprio header
class Cylp {
private:
	Matrix m_Node1pos;
	Matrix m_Node1vel;
	Matrix m_Node1acc;

	Matrix m_Node2pos;
	Matrix m_Node2vel;
	Matrix m_Node2acc;

	int m_numIntPoints;
	Matrix m_IntPoints;

	float cd;
	float cm;
	//etc


public:
	std::array<double, 6> calcMorisonForce( Envir &envir );

	void updateNodes(); // Update node position, velocity and acceleration according to floater position, velocity and acceleration

	Matrix getIntPointsPos(); // calcula a posição dos pontos de integração baseado nas posicoes dos nós
	Matrix getIntPointsVel(); // calcula a velocidade dos pontos de integração baseado nas posicoes dos nós
	Matrix getIntPointsAcc(); // calcula a aceleração dos pontos de integração baseado nas posicoes dos nós

};




/* Vai pro Fowt.cpp
Basicamente, pra cada cylp e cada retp, preciso da velocidade/aceleracao do fluido em seus pontos de integração (que vai vir de uma função da classe Envir),
da velocidade/aceleracao do cilindro em seus nós (o suficiente pra calcular em cada ponto de integração). Depois, calculo a força de Morison,
a força de heave plate e a força hidrostática.

Acho que a velocidade e aceleração estruturais dos cylp/retp podem ser calculados logo quando atualizo a posição, velocidade e aceleração do
flutuador a cada time step.
*/
std::array<double, 6> Fowt::calcHydroForce( Envir &envir ) {
	// Calculate hydro forces on each cylinder
	std::array<double, 6> cylpMorisonForce; // não lembro se preciso inicializar std::array pra evitar lixo
	std::array<double, 6> cylpPlateForce;
	std::array<double, 6> cylpHydrostForce;

	for (int iiCylp = 0; iiCylp < m_floater.getNumCylp(); ++iiCylp) {
		// 1) Calcula a velocidade/aceleração do fluido no referencial global em cada um dos pontos de integração do cilindro. Isso é feito passando a matriz dos nós 
		// do cilindro pra uma função da classe Envir. É algo tipo o abaixo
		Matrix intPointsFlowVel = envir.calcFlowVel( m_floater.m_cylp[ii].getIntPointsPos() );
		Matrix intPointsFlowAcc = envir.calcFlowAcc( m_floater.m_cylp[ii].getIntPointsPos() );

		// 2) Calcula a força de Morison em cada cilindro. Passo apenas a velocidade/aceleração do fluido em cada ponto de integração.
		// Dentro da calcMorisonForce, projeto a velocidade/aceleração nas direções perpendiculares ao cilindro; calculo a velocidade/aceleração dos pontos de integração e projeto tmbm nas
		// direções perpendiculares (funções getIntPointVel e getIntPointAcc); o resto é tranquilo
		cylpMorisonForce += m_floater.cylp[iiCylp].calcMorisonForce( intPointsFlowVel, intPointsFlowAcc );


		// Calcula forca de heave plate
		// Calcula forca hidrostatica
	}


	// Calculate hydro forces on each rectangle
/*	std::array<double, 6> retpForce;
	for (int iiRetp = 0; iiRetp < m_floater.getNumCylp(); ++iiRetp) {
		m_floater.retp[iiRetp].calcPerpVelAcc(); //
		retpForce += m_floater.retp[iiRetp].calcMorisonForce(m_Wave, m_Current);
	}
*/

}
