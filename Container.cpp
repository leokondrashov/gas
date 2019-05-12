#include "Container.h"
#include "Constants.h"
#include <random>
#include <assert.h>

#define isInVolume(coordinates) ((coordinates->x >= VOLUME_X) && (coordinates->x <= VOLUME_X + VOLUME_LENGTH) && \
(coordinates->y >= VOLUME_Y) && (coordinates->y <= VOLUME_Y + VOLUME_WIDTH) && \
(coordinates->z >= VOLUME_Z) && (coordinates->z <= VOLUME_Z + VOLUME_HEIGHT))

Container::Container(double length, double width, double height, double temperature) {
	this->length = length;
	this->width = width;
	this->height = height;
	pressure = 0;
	surface = 2 * length * width + 2 * length * height + 2 * width * height;
	volume = length * width * height;
	pressureAccumulator = 0;
	countInVolume = 0;
	
	tree = new OctTree(length, width, height);
	
	gas = new Molecule*[MOLECULES_NUM];
	for (int i = 0; i < MOLECULES_NUM; i++) {
		this->gas[i] = new Molecule(MOLECULE_MASS, MOLECULE_RADIUS);
	}
	setTemperature(temperature);
}

Container::Container(double length, double width, double height) {
	this->length = length;
	this->width = width;
	this->height = height;
	pressure = 0;
	surface = 2 * length * width + 2 * length * height + 2 * width * height;
	volume = length * width * height;
	pressureAccumulator = 0;
	countInVolume = 0;
	
	tree = new OctTree(length, width, height);
	
	gas = new Molecule*[MOLECULES_NUM];
	for (int i = 0; i < MOLECULES_NUM; i++) {
		this->gas[i] = new Molecule(MOLECULE_MASS, MOLECULE_RADIUS);
	}
}

Container::~Container() {
	for (int i = 0; i < MOLECULES_NUM; i++) {
		delete gas[i];
	}
	delete[] gas;
	
	delete tree;
}

void Container::getVelocityDistribution(sf::Vector3<double> *data) {
	assert(data);
	
	const sf::Vector3<double> *curVelocity = nullptr;
	for (int i = 0; i < MOLECULES_NUM; i++) {
		curVelocity = gas[i]->getVelocity();
		data[i].x = curVelocity->x;
		data[i].y = curVelocity->y;
		data[i].z = curVelocity->z;
	}
}

void Container::testCollisions(const OctTreeNode *node, Molecule *molecule) {
	assert(molecule);
	if (node == nullptr) {
		return;
	}
	
	if (node->isNear(molecule)) {
		Molecule *m = node->getMolecule();
		if (m != nullptr) {
			if (m->isCollided() || m == molecule)
				return;
			collide(m, molecule);
			return;
		}
		
		for (int i = 0; i < 8; i++) {
			testCollisions(node->getChild(i), molecule);
		}
	}
}

void Container::update(long tick) {
	int countTmp = 0;
	double accumulator = 0;
	for (int i = 0; i < MOLECULES_NUM; i++) {
		
		const sf::Vector3<double> *curPosition = nullptr;
		gas[i]->update(DELTA);
		curPosition = gas[i]->getCoordinates();
		
		if (curPosition->x > length || curPosition->x < 0) {
			const sf::Vector3<double> *curVelocity = gas[i]->getVelocity();
			gas[i]->setVelocity(-curVelocity->x, curVelocity->y, curVelocity->z);
			accumulator += 2 * gas[i]->getMass() * fabs(curVelocity->x);
		}
		
		if (curPosition->y > width || curPosition->y < 0) {
			const sf::Vector3<double> *curVelocity = gas[i]->getVelocity();
			gas[i]->setVelocity(curVelocity->x, -curVelocity->y, curVelocity->z);
			accumulator += 2 * gas[i]->getMass() * fabs(curVelocity->y);
		}
		
		if (curPosition->z > height || curPosition->z < 0) {
			const sf::Vector3<double> *curVelocity = gas[i]->getVelocity();
			gas[i]->setVelocity(curVelocity->x, curVelocity->y, -curVelocity->z);
			accumulator += 2 * gas[i]->getMass() * fabs(curVelocity->z);
		}
		
		if (isInVolume(curPosition)) {
			countTmp++;
		}
	}
	
	tree->clear();
	
	for (int i = 0; i < MOLECULES_NUM; i++) {
		tree->AddMolecule(gas[i]);
	}
	
//	tree->dump();
	
	for (int i = 0; i < MOLECULES_NUM; i++) {
		testCollisions(tree->gerRoot(), gas[i]);
		gas[i]->setCollided();
	}
	
	pressureAccumulator += accumulator;
	countInVolume = countTmp;
	
	if (tick % PRESSURE_TICKS_AVERAGE == 0) {
		pressure = pressureAccumulator / DELTA / PRESSURE_TICKS_AVERAGE / surface;
		pressureAccumulator = 0;
	}
}

double Container::getPressure() {
	return pressure;
}

int Container::getCountInVolume() {
	return countInVolume;
}

void Container::setTemperature(double temperature) {
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_real_distribution<> xDistribution(0.0, length);
	std::uniform_real_distribution<> yDistribution(0.0, width);
	std::uniform_real_distribution<> zDistribution(0.0, length);
	std::normal_distribution<> vDistribution(0.0, sqrt(kBoltzmann * temperature / MOLECULE_MASS));
	
	for (int i = 0; i < MOLECULES_NUM; i++) {
		this->gas[i]->setCoordinates(xDistribution(gen), yDistribution(gen), xDistribution(gen));
		this->gas[i]->setVelocity(vDistribution(gen), vDistribution(gen), vDistribution(gen));
	}
}

double dotProduct(const sf::Vector3<double> *v1, const sf::Vector3<double> *v2) {
	assert(v1);
	assert(v2);
	return v1->x * v2->x + v1->y * v2->y + v1->z * v2->z;
}

bool Container::collisionTest(Molecule *m1, Molecule *m2) {
	assert(m1);
	assert(m2);
	sf::Vector3<double> dr = *(m1->getCoordinates()) - *(m2->getCoordinates());
	return dotProduct(&dr, &dr) < (m1->getRadius() + m2->getRadius()) * (m1->getRadius() + m2->getRadius());
}

void Container::collide(Molecule *m1, Molecule *m2) {
//	printf("collision\n");
	assert(m1);
	assert(m2);
	assert(m1 != m2);

	if (!collisionTest(m1, m2))
		return;
	
	const sf::Vector3<double> *r1 = (m1->getCoordinates());
	const sf::Vector3<double> *r2 = (m2->getCoordinates());
	sf::Vector3<double> dr = *r1 - *r2;
	const sf::Vector3<double> *v1 = (m1->getVelocity());
	const sf::Vector3<double> *v2 = (m2->getVelocity());
	
	sf::Vector3<double> u = *v2 - *v1;
	sf::Vector3<double> uPerp = dr * dotProduct(&u, &dr) / dotProduct(&dr, &dr);
	
	sf::Vector3<double> newV1 = uPerp / m1->getMass() * m2->getMass() + *v1;
	sf::Vector3<double> newV2 = *v2 - uPerp;
	
	m1->setVelocity(newV1.x, newV1.y, newV1.z);
	m2->setVelocity(newV2.x, newV2.y, newV2.z);
}


void Container::setRandomDistribution() {
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_real_distribution<> xDistribution(0.0, length);
	std::uniform_real_distribution<> yDistribution(0.0, width);
	std::uniform_real_distribution<> zDistribution(0.0, length);
	std::uniform_real_distribution<> vDistribution(-MAX_VELOCITY / sqrt(3) / 2, MAX_VELOCITY / sqrt(3) / 2);
	
	for (int i = 0; i < MOLECULES_NUM; i++) {
		this->gas[i]->setCoordinates(xDistribution(gen), yDistribution(gen), xDistribution(gen));
		this->gas[i]->setVelocity(vDistribution(gen), vDistribution(gen), vDistribution(gen));
	}
}


double Container::getEnergy() {
	double energy = 0;
	for (int i = 0; i < MOLECULES_NUM; i++) {
		const sf::Vector3<double> *velocity = gas[i]->getVelocity();
		double velocitySqr = dotProduct(velocity, velocity);
		energy += gas[i]->getMass() * velocitySqr / 2;
	}
	
	return energy;
}

void Container::addMassiveMolecules(int num) {
	for (; num > 0; num--) {
		gas[MOLECULES_NUM - num]->setMass(MASSIVE_MOLECULE_MASS);
		const sf::Vector3<double> *velocity = gas[MOLECULES_NUM - num]->getVelocity();
		double coefficient = sqrt(MASSIVE_MOLECULE_MASS / MOLECULE_MASS);
		gas[MOLECULES_NUM - num]->setVelocity(velocity->x / coefficient, velocity->y / coefficient, velocity->z / coefficient);
	}
}

void Container::addFastMolecules(int num) {
	for (; num > 0; num--) {
		const sf::Vector3<double> *velocity = gas[MOLECULES_NUM - num]->getVelocity();
		gas[MOLECULES_NUM - num]->setVelocity(velocity->x * 10, velocity->y * 10, velocity->z * 10);
	}
}

void Container::getEnergyDistribution(double *data) {
	assert(data);
	
	for (int i = 0; i < MOLECULES_NUM; i++) {
		data[i] = gas[i]->getEnergy();
	}
}

void Container::addSlowMolecules(int num) {
	for (; num > 0; num--) {
		gas[MOLECULES_NUM - num]->setVelocity(0, 0, 0);
	}
}
