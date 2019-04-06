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
	double tempeture;
	double pressure;
	double pressureAccumulator;
	double surface;
	double volume;

public:
	Container(double length, double width, double height, double temperature);
	~Container();
	
	void update(long tick);
	
	void getVelocityDistribution(sf::Vector3<double> *data);
	double getPressure();
};

#endif //GAS_CONTAINER_H
