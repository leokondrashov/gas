#ifndef GAS_MOLECULE_H
#define GAS_MOLECULE_H

#include "Constants.h"
#include <SFML/System/Vector3.hpp>

class Molecule {
private:
	sf::Vector3<double> *r;
	sf::Vector3<double> *v;
	double mass;
//	double radius;

public:
	Molecule(double rx, double ry, double rz, double vx, double vy, double vz, double mass);
	~Molecule();
	
	void update(double delta);
	
	const sf::Vector3<double> *getCoordinates();
	const sf::Vector3<double> *getVelocity();
	double getMass();
	
	void setVelocity(double vx, double vy, double vz);
};

#endif
