#pragma once

class Envir{
private:
	// Environmental data
	Wave m_wave;
	Current m_current;
	Wind m_wind;

	// Numerical parameters
	float m_dt;
	float m_time;
	float m_finalTime;

public:
	Matrix calcFlowVel('Referencia constante pra um vetor de pontos nos quais quero calcular velocidade');
	Matrix calcFlowAcc('Referencia constante pra um vetor de pontos nos quais quero calcular velocidade');

};
