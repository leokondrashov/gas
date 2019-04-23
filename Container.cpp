#include "Container.h"
#include <random>

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
	
	gas = new Molecule*[MOLECULES_NUM];
	for (int i = 0; i < MOLECULES_NUM; i++) {
		this->gas[i] = new Molecule(MOLECULE_MASS, MOLECULE_RADIUS);
	}
	setTemperature(temperature);
}

Container::~Container() {
	for (int i = 0; i < MOLECULES_NUM; i++) {
		delete gas[i];
	}
	delete[] gas;
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
	
	for (int i = 0; i < MOLECULES_NUM; i++) {
		for (int j = i + 1; j < MOLECULES_NUM; j++) {
			if (collisionTest(gas[i], gas[j])) {
				collide(gas[i], gas[j]);
			}
		}
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
	return v1->x * v2->x + v1->y * v2->y + v1->z * v2->z;
}

bool Container::collisionTest(Molecule *m1, Molecule *m2) {
	sf::Vector3<double> dr = *(m1->getCoordinates()) - *(m2->getCoordinates());
	return dotProduct(&dr, &dr) < (MOLECULE_RADIUS * 2) * (MOLECULE_RADIUS * 2);
}

void Container::collide(Molecule *m1, Molecule *m2) {
//	printf("collision\n");
	const sf::Vector3<double> *r1 = (m1->getCoordinates());
	const sf::Vector3<double> *r2 = (m2->getCoordinates());
	sf::Vector3<double> dr = *r1 - *r2;
	const sf::Vector3<double> *v1 = (m1->getVelocity());
	const sf::Vector3<double> *v2 = (m2->getVelocity());
	sf::Vector3<double> u = *v2 - *v1;
	sf::Vector3<double> uPerp = dr * dotProduct(&u, &dr) / dotProduct(&dr, &dr);
	
	sf::Vector3<double> newV1 = uPerp + *v1;
	sf::Vector3<double> newV2 = *v2 - uPerp;
	
	m1->setVelocity(newV1.x, newV1.y, newV1.z);
	m2->setVelocity(newV2.x, newV2.y, newV2.z);
}
