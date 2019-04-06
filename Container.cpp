#include "Container.h"
#include <random>

Container::Container(double length, double width, double height, double temperature) {
	this->length = length;
	this->width = width;
	this->height = height;
	this->tempeture = temperature;
	pressure = 0;
	surface = 2 * length * width + 2 * length * height + 2 * width * height;
	volume = length * width * height;
	
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_real_distribution<> xDistribution(0.0, length);
	std::uniform_real_distribution<> yDistribution(0.0, width);
	std::uniform_real_distribution<> zDistribution(0.0, length);
	std::normal_distribution<> vDistribution(0.0, sqrt(kBoltzmann * temperature / MOLECULE_MASS));
	
	gas = new Molecule*[MOLECULES_NUM];
	for (int i = 0; i < MOLECULES_NUM; i++) {
		this->gas[i] = new Molecule(xDistribution(gen), yDistribution(gen), xDistribution(gen),
				vDistribution(gen), vDistribution(gen), vDistribution(gen), MOLECULE_MASS);
	}
}

Container::~Container() {
	for (int i = 0; i < MOLECULES_NUM; i++) {
		delete gas[i];
	}
	delete gas;
}

void Container::getVelocityDistribution(sf::Vector3<double> *data) {
	const sf::Vector3<double> *curVelocity = nullptr;
	for (int i = 0; i < MOLECULES_NUM; i++) {
		curVelocity = gas[i]->getVelocity();
		data[i].x = curVelocity->x;
		data[i].y = curVelocity->y;
		data[i].z = curVelocity->z;
	}
}

void Container::update(long tick) {
	for (int i = 0; i < MOLECULES_NUM; i++) {
		const sf::Vector3<double> *curPosition = nullptr;
		gas[i]->update(DELTA);
		curPosition = gas[i]->getCoordinates();
		
		if (curPosition->x > length || curPosition->x < 0) {
			const sf::Vector3<double> *curVelocity = gas[i]->getVelocity();
			gas[i]->setVelocity(-curVelocity->x, curVelocity->y, curVelocity->z);
			pressureAccumulator += 2 * gas[i]->getMass() * fabs(curVelocity->x);
		}
		
		if (curPosition->y > width || curPosition->y < 0) {
			const sf::Vector3<double> *curVelocity = gas[i]->getVelocity();
			gas[i]->setVelocity(curVelocity->x, -curVelocity->y, curVelocity->z);
			pressureAccumulator += 2 * gas[i]->getMass() * fabs(curVelocity->y);
		}
		
		if (curPosition->z > height || curPosition->z < 0) {
			const sf::Vector3<double> *curVelocity = gas[i]->getVelocity();
			gas[i]->setVelocity(curVelocity->x, curVelocity->y, -curVelocity->z);
			pressureAccumulator += 2 * gas[i]->getMass() * fabs(curVelocity->z);
		}
	}
	
	if (tick % PRESSURE_TICKS_AVERAGE == 0) {
		pressure = pressureAccumulator / DELTA / PRESSURE_TICKS_AVERAGE / surface;
		pressureAccumulator = 0;
	}
}

double Container::getPressure() {
	return pressure;
}