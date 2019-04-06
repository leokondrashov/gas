#include "Molecule.h"
#include "Constants.h"

Molecule::Molecule(double rx, double ry, double rz, double vx, double vy, double vz, double mass) {
	r = new sf::Vector3<double>();
	r->x = rx;
	r->y = ry;
	r->z = rz;
	
	v = new sf::Vector3<double>();
	v->x = vx;
	v->y = vy;
	v->z = vz;
	
	this->mass = mass;
}

Molecule::~Molecule() {
	delete r;
	delete v;
}

void Molecule::update(double delta) {
	*r += delta * *v;
}

const sf::Vector3<double> *Molecule::getCoordinates() {
	return r;
}

const sf::Vector3<double> *Molecule::getVelocity() {
	return v;
}

void Molecule::setVelocity(double vx, double vy, double vz) {
	v->x = vx;
	v->y = vy;
	v->z = vz;
}
double Molecule::getMass() {
	return mass;
}

