#include <iostream>
#include <cmath>
#include "Container.h"
#include "Constants.h"

void printVelocityDistribution(long time, sf::Vector3<double> *data, FILE *velFile);
void initFiles(FILE *velFile, FILE *presFile);
void printPressure(long time, double pressure, FILE *presFile);

int main() {
	Container container(CONTAINER_LENGTH, CONTAINER_WIDTH, CONTAINER_HEIGHT, GAS_TEMPERATURE);
	
	FILE *velocities = fopen(VELOCITY_FILE, "wb");
	FILE *pressures = fopen(PRESSURE_FILE, "wb");
	initFiles(velocities, pressures);
	
	sf::Vector3<double> *data = new sf::Vector3<double>[MOLECULES_NUM];
	
	for (long t = 0; t <= (long) 1e5; t++) {
		if (t % PRESSURE_TICKS_AVERAGE == 0) {
			container.getVelocityDistribution(data);
			printVelocityDistribution(t, data, velocities);
			printPressure(t, container.getPressure(), pressures);
		}
		container.update(t);
	}
	
	fclose(velocities);
	fclose(pressures);
	
	delete[] data;
	
	return 0;
}

void printVelocityDistribution(long time, sf::Vector3<double> *data, FILE *velFile) {
	int *count = new int[(int) (MAX_VELOCITY / VELOCITY_PRECISION) + 1] { 0 };
	
	for (int i = 0; i < MOLECULES_NUM; i++) {
		double v = sqrt(data[i].x * data[i].x + data[i].y * data[i].y + data[i].z * data[i].z);
		if (v > MAX_VELOCITY) {
			count[(int) (MAX_VELOCITY / VELOCITY_PRECISION)]++;
		} else {
			count[(int) (v / VELOCITY_PRECISION)]++;
		}
	}
	
	fprintf(velFile, "%12ld ", time);
	
	for (int i = 0; i < MAX_VELOCITY / VELOCITY_PRECISION; i++) {
		fprintf(velFile, "%6d ", count[i]);
	}
	
	fprintf(velFile, "\n");
	
	delete[] count;
}

void initFiles(FILE *velFile, FILE *presFiles) {
	fprintf(velFile, "t, %7g s", DELTA);
	for (int i = 0; i < MAX_VELOCITY / VELOCITY_PRECISION; i++) {
		fprintf(velFile, "%6.0f ", i * VELOCITY_PRECISION);
	}
	fprintf(velFile, "\n");
	
	fprintf(presFiles, "t, %7g s Pressure\n", DELTA);
}

void printPressure(long time, double pressure, FILE *presFile) {
	fprintf(presFile, "%12ld %12g\n", time, pressure);
}