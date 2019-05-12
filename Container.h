#ifndef GAS_CONTAINER_H
#define GAS_CONTAINER_H

#include "Molecule.h"
#include "OctTree.h"

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
	
	OctTree *tree;
	
//	bool isInVolume(const sf::Vector3<double> *coordinates);
	bool collisionTest(Molecule *m1, Molecule *m2);
	void collide(Molecule *m1, Molecule *m2);

public:
	Container(double length, double width, double height, double temperature);
	Container(double length, double width, double height);
	~Container();
	
	void update(long tick);
	
	void getVelocityDistribution(sf::Vector3<double> *data);
	void getEnergyDistribution(double *data);
	double getPressure();
	int getCountInVolume();
	double getEnergy();
	
	void setTemperature(double temperature);
	void setRandomDistribution();
	
	void testCollisions(const OctTreeNode *node, Molecule *molecule);
	
	void addMassiveMolecules(int num);
	void addFastMolecules(int num);
	void addSlowMolecules(int num);
};

#endif //GAS_CONTAINER_H
