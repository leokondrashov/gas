#ifndef GAS_MOLECULE_H
#define GAS_MOLECULE_H

#include <SFML/System/Vector3.hpp>

class Molecule {
private:
	sf::Vector3<double> *r;
	sf::Vector3<double> *v;
	double mass;
	double radius;
	bool collided;

public:
//	Molecule(double rx, double ry, double rz, double vx, double vy, double vz, double mass);
	Molecule(double mass, double radius);
	~Molecule();
	
	void update(double delta);
	
	const sf::Vector3<double> *getCoordinates();
	const sf::Vector3<double> *getVelocity();
	double getMass();
	double getRadius();
	double getEnergy();
	
	void setVelocity(double vx, double vy, double vz);
	void setCoordinates(double rx, double ry, double rz);
	void setMass(double mass);
	
	
	
	void setCollided();
	bool isCollided();
};

#endif
