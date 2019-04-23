#ifndef GAS_CONTAINER_H
#define GAS_CONTAINER_H

#include "Constants.h"
#include "Molecule.h"

class Container {
private:
	Molecule **gas;
	double length;
	double width;
	double height;
	double pressure;
	double pressureAccumulator;
	double surface;
	double volume;
	int countInVolume;
	
//	bool isInVolume(const sf::Vector3<double> *coordinates);
	bool collisionTest(Molecule *m1, Molecule *m2);
	void collide(Molecule *m1, Molecule *m2);

public:
	Container(double length, double width, double height, double temperature);
	~Container();
	
	void update(long tick);
	
	void getVelocityDistribution(sf::Vector3<double> *data);
	double getPressure();
	int getCountInVolume();
	
	void setTemperature(double temperature);
};

#endif //GAS_CONTAINER_H
