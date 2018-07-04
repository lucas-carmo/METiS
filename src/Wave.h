class Wave{
private:
    double m_amplitude;
    double m_period;
    double m_direction;

public:
	/*****************************************************
		Setters
	*****************************************************/
    void setAmplitude();
    void setPeriod();
    void setDirection();

	/*****************************************************
		Getters
	*****************************************************/
	double getAmplitude();
	double getPeriod();
	double getDirection();
};
