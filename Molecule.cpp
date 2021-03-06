#include <assert.h>
#include "Molecule.h"

//Molecule::Molecule(double rx, double ry, double rz, double vx, double vy, double vz, double mass) {
//	r = new sf::Vector3<double>();
//	r->x = rx;
//	r->y = ry;
//	r->z = rz;
//
//	v = new sf::Vector3<double>();
//	v->x = vx;
//	v->y = vy;
//	v->z = vz;
//
//	this->mass = mass;
//}

Molecule::Molecule(double mass, double radius) {
	r = new sf::Vector3<double>();
	r->x = 0;
	r->y = 0;
	r->z = 0;
	
	v = new sf::Vector3<double>();
	v->x = 0;
	v->y = 0;
	v->z = 0;
	
	this->mass = mass;
	this->radius = radius;
	
	collided = false;
}

Molecule::~Molecule() {
	delete r;
	delete v;
}

void Molecule::update(double delta) {
	*r += delta * *v;
	collided = false;
}

const sf::Vector3<double> *Molecule::getCoordinates() {
	assert(this);
	
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

void Molecule::setCoordinates(double rx, double ry, double rz) {
	r->x = rx;
	r->y = ry;
	r->z = rz;
}

double Molecule::getMass() {
	return mass;
}

void Molecule::setCollided() {
	collided = true;
}

bool Molecule::isCollided() {
	return collided;
}

double Molecule::getRadius() {
	return radius;
}

void Molecule::setMass(double mass) {
	this->mass = mass;
}

double Molecule::getEnergy() {
	double velocitySqr = v->x * v->x + v->y * v->y + v->z * v->z;
	return mass * velocitySqr / 2;
}

