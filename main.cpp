#include <iostream>
#include <cmath>
#include "Container.h"
#include "Constants.h"

void initFiles(FILE *velFile, FILE *presFile, FILE *pressVsTempFile, FILE *fluctsVsTempFile);
void printVelocityDistribution(long time, sf::Vector3<double> *data, FILE *velFile);
void printPressure(long time, double pressure, FILE *presFile);
void printPressureVsTemperature(double temp, double press, FILE *presVsTempFile);
void printFluctuations(double temp, double fluctuation, double mean, FILE *fluctsVsTempFile);
double fluctuations(int *data);
double mean(int *data);

int main() {
	Container container(CONTAINER_LENGTH, CONTAINER_WIDTH, CONTAINER_HEIGHT, GAS_TEMPERATURE_START);
	
	FILE *velocities = fopen(VELOCITY_FILE, "wb");
	FILE *pressures = fopen(PRESSURE_FILE, "wb");
	FILE *pressVsTemp = fopen(PRESS_VS_TEMP_FILE, "wb");
	FILE *fluctsVsTemp = fopen(FLUCTUATIONS_VS_TEMP, "wb");
	
	initFiles(velocities, pressures, pressVsTemp, fluctsVsTemp);
	
//	sf::Vector3<double> *data = new sf::Vector3<double>[MOLECULES_NUM];
	
/*	for (long t = 0; t <= (long) 1e5; t++) {
		container.update(t);
		if (t % TICKS_AVERAGE == 0) {
			container.getVelocityDistribution(data);
			printVelocityDistribution(t, data, velocities);
			printPressure(t, container.getPressure(), pressures);
		}
	}*/

	int *data = new int[TICKS_AVERAGE + 1];

	for (double temperature = GAS_TEMPERATURE_START; temperature <= GAS_TEMPERATURE_END; temperature += GAS_TEMPERATURE_STEP, container.setTemperature(temperature)) {
		for (int t = 0; t <= TICKS_AVERAGE; t++) {
			container.update(t);
			data[t] = container.getCountInVolume();
		}
		printFluctuations(temperature, fluctuations(data), mean(data), fluctsVsTemp);
		printPressureVsTemperature(temperature, container.getPressure(), pressVsTemp);
	}
	
	fclose(velocities);
	fclose(pressures);
	fclose(pressVsTemp);
	
	delete[] data;
	
	return 0;
}

void initFiles(FILE *velFile, FILE *presFiles, FILE *pressVsTempFile, FILE *fluctsVsTempFile) {
	fprintf(velFile, "t (%6g s), ", DELTA);
	for (int i = 0; i < MAX_VELOCITY / VELOCITY_PRECISION; i++) {
		fprintf(velFile, "%6.0f, ", i * VELOCITY_PRECISION);
	}
	fprintf(velFile, "\n");
	
	fprintf(presFiles, "t (%6g s), Pressure\n", DELTA);
	
	fprintf(pressVsTempFile, "T(1K)\t\t, Pressure\n");
	
	fprintf(fluctsVsTempFile, "T(1K)\t\t, Fluctuations, Mean\t\t\n");
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
	
	fprintf(velFile, "%12ld, ", time);
	
	for (int i = 0; i < MAX_VELOCITY / VELOCITY_PRECISION; i++) {
		fprintf(velFile, "%6d, ", count[i]);
	}
	
	fprintf(velFile, "\n");
	
	delete[] count;
}

void printPressure(long time, double pressure, FILE *presFile) {
	fprintf(presFile, "%12ld, %12g\n", time, pressure);
}

void printPressureVsTemperature(double temp, double press, FILE *presVsTempFile) {
	fprintf(presVsTempFile, "%12g, %12g\n", temp, press);
}


void printFluctuations(double temp, double fluctuation, double mean, FILE *fluctsVsTempFile) {
	fprintf(fluctsVsTempFile, "%12g, %12g, %12g\n", temp, fluctuation, mean);
}

double mean(const int *data) {
	double meanN = 0;
	for (int i = 0; i <= TICKS_AVERAGE; i++) {
		meanN += data[i];
	}
	meanN /= TICKS_AVERAGE + 1;
	
	return meanN;
}

double fluctuations(const int *data) {
	double meanN = mean(data);
	
	double variation = 0;
	
	for (int i = 0; i <= TICKS_AVERAGE; i++) {
		variation += (data[i] - meanN) * (data[i] - meanN);
	}
	
	return sqrt(variation / (TICKS_AVERAGE + 1));
}