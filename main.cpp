#include <iostream>
#include <cmath>
#include <assert.h>
#include <cstring>
#include "Container.h"
#include "Constants.h"

#define TIME_STAMP_LENGTH 12

char timeStamp[TIME_STAMP_LENGTH + 1] = "";

void initTime();

FILE *creatVelocities();
FILE *creatPressures();
FILE *creatTimes();
FILE *creatEnergies();
FILE *creatEnergyDistr();
FILE *creatPressVsTemp();
FILE *creatFluctsVsTemp();
FILE *creatPressFlucts();

void printConsts();

void printVelocityDistribution(long time, sf::Vector3<double> *data, FILE *velFile);
void printPressure(long time, double pressure, FILE *presFile);
void printTime(double temperature, long time, FILE *timeFile);
void printEnergy(long time, double energy, FILE *enrgFile);
void printEnergyDistribution(long time, double *data, FILE *energyDistrFile);
void printPressureVsTemperature(double temp, double press, FILE *presVsTempFile);
void printFluctuations(double temp, double fluctuation, double mean, FILE *fluctsVsTempFile);
void printPressureFluctuations(double temp, double fluctuation, double mean, FILE *pressFluctsFile);

double mean(const int *data, int n);
double mean(const double *data, int n);

double fluctuations(const int *data, int n);
double fluctuations(const double *data, int n);

double maxDifference(const double *data, const double *prevData, int n);
int maxDifference(const int *data, const int *prevData, int n);

void getEnergyDistribution(const double *data, int n, int *dest);

int main() {
	Container container(CONTAINER_LENGTH, CONTAINER_WIDTH, CONTAINER_HEIGHT, TEMPERATURE);
//	container.setRandomDistribution();
//	container.addMassiveMolecules(MOLECULES_NUM / 2);
	container.addFastMolecules(MOLECULES_NUM / 2);
	
	initTime();
	
	printConsts();

//	FILE *velocities = creatVelocities();
//	FILE *pressures = creatPressures();
//	FILE *energies = creatEnergies();
	FILE *energyDistr = creatEnergyDistr();
	FILE *timeFile = creatTimes();
//	FILE *pressVsTemp = creatPressVsTemp();
//	FILE *fluctsVsTemp = creatFluctsVsTemp();
//	FILE *pressFlucts = creatPressFlucts();

//	double *data = new double[MOLECULES_NUM];
//
//	int *count = new int[(int) (MAX_ENERGY / ENERGY_PRECISION) + 1] { 0 };
//	int *prevCount = new int[(int) (MAX_ENERGY / ENERGY_PRECISION) + 1] { 0 };
//	
//	for (long t = 0; t <= EXPERIMENT_LENGTH; t++) {
//		if (t % TICKS_AVERAGE == 0) {
//			container.getEnergyDistribution(data);
//			printEnergyDistribution(t, data, energyDistr);
//			
//			int *tmp = count;
//			count = prevCount;
//			prevCount = tmp;
//			
//			memset(count, 0, sizeof(int) * ((int) (MAX_ENERGY / ENERGY_PRECISION) + 1));
//			getEnergyDistribution(data, MOLECULES_NUM, count);
//			
//			if ((double) maxDifference(count, prevCount, (int) (MAX_ENERGY / ENERGY_PRECISION) + 1)
//					/ MOLECULES_NUM < EPSILON) {
//				printf("%ld", t);
//				break;
//			}
//			
//			//			printPressure(t, container.getPressure(), pressures);
//			//			printEnergy(t, container.getEnergy(), energies);
//		}
//		
//		container.update(t);
//	}
//	
	sf::Vector3<double> *data = new sf::Vector3<double>[MOLECULES_NUM];
	double *energyData = new double[MOLECULES_NUM];
	int *count = new int[(int) (MAX_ENERGY / ENERGY_PRECISION) + 1] { 0 };
	int *prevCount = new int[(int) (MAX_ENERGY / ENERGY_PRECISION) + 1] { 0 };
	for (double temp = GAS_TEMPERATURE_START; temp <= GAS_TEMPERATURE_END; temp += GAS_TEMPERATURE_STEP) {
		for (int i = 0; i < 30; i ++) {
			container.setTemperature(temp);
			container.addSlowMolecules(MOLECULES_NUM / 2);
			for (long t = 0; t <= EXPERIMENT_LENGTH; t++) {
				if (t % TICKS_AVERAGE == 0) {
					container.getEnergyDistribution(energyData);
					printEnergyDistribution(t, energyData, energyDistr);

					int *tmp = count;
					count = prevCount;
					prevCount = tmp;

					memset(count, 0, sizeof(int) * ((int) (MAX_ENERGY / ENERGY_PRECISION) + 1));
					getEnergyDistribution(energyData, MOLECULES_NUM, count);

					if ((double) maxDifference(count, prevCount, (int) (MAX_ENERGY / ENERGY_PRECISION) + 1)
							/ MOLECULES_NUM < EPSILON) {
						printTime(temp, t, timeFile);
						printf("%g\n", temp);
						break;
					}

					//			printPressure(t, container.getPressure(), pressures);
					//			printEnergy(t, container.getEnergy(), energies);
				}

				container.update(t);
			}
		}
	}



//	int *count = new int[TICKS_AVERAGE];
//	double *pressData = new double[TICKS_AVERAGE / PRESSURE_TICKS_AVERAGE];
//	double pressure = 0;
//
//	for (double temperature = GAS_TEMPERATURE_START; temperature <= GAS_TEMPERATURE_END; temperature += GAS_TEMPERATURE_STEP, container.setTemperature(temperature)) {
//		for (int t = 1; t <= TICKS_AVERAGE; t++) {
//			container.update(t);
//			count[t] = container.getCountInVolume();
//			if (t % PRESSURE_TICKS_AVERAGE == 0) {
//				pressData[t / PRESSURE_TICKS_AVERAGE - 1] = container.getPressure();
//				printPressure(t, pressData[t / PRESSURE_TICKS_AVERAGE - 1], pressures);
//			}
//		}
//		pressure = mean(pressData, TICKS_AVERAGE / PRESSURE_TICKS_AVERAGE);
//		printFluctuations(temperature, fluctuations(count, TICKS_AVERAGE), mean(count, TICKS_AVERAGE), fluctsVsTemp);
//		printPressureFluctuations(temperature, fluctuations(pressData, TICKS_AVERAGE / PRESSURE_TICKS_AVERAGE),
//				pressure, pressFlucts);
//		printPressureVsTemperature(temperature, pressure, pressVsTemp);
//	}
	
//	fclose(velocities);
//	fclose(pressures);
//	fclose(energies);
	fclose(energyDistr);
	fclose(timeFile);
//	fclose(pressVsTemp);
//	fclose(fluctsVsTemp);
//	fclose(pressFlucts);
	
//	delete[] count;
//	delete[] pressData;
	delete[] data;
	
	return 0;
}

void initTime() {
	time_t tp = time(nullptr);
	strftime(timeStamp, TIME_STAMP_LENGTH + 1, "%y%m%d%H%M%S", localtime(&tp));
}

FILE *creatVelocities() {
	char name[sizeof(VELOCITY_FILE) - 2 + TIME_STAMP_LENGTH] = "";
	sprintf(name, VELOCITY_FILE, timeStamp);
	
	FILE *velFile = fopen(name, "wb");
	
	fprintf(velFile, "t (%6g s)", DELTA);
	for (int i = 0; i <= MAX_VELOCITY / VELOCITY_PRECISION; i++) {
		fprintf(velFile, ", %6.0f", i * VELOCITY_PRECISION);
	}
	fprintf(velFile, "\n");
	
	return velFile;
}

FILE *creatPressures() {
	char name[sizeof(PRESSURE_FILE) - 2 + TIME_STAMP_LENGTH] = "";
	sprintf(name, PRESSURE_FILE, timeStamp);
	
	FILE *pressFile = fopen(name, "wb");
	
	fprintf(pressFile, "t (%6g s), Pressure\n", DELTA);
	
	return pressFile;
}

FILE *creatTimes() {
	char name[sizeof(TIME_FILE) - 2 + TIME_STAMP_LENGTH] = "";
	sprintf(name, TIME_FILE, timeStamp);
	
	FILE *timeFile = fopen(name, "wb");
	
	fprintf(timeFile, " temperature, time\n");
	
	return timeFile;
}

FILE *creatEnergies() {
	char name[sizeof(ENERGY_FILE) - 2 + TIME_STAMP_LENGTH] = "";
	sprintf(name, ENERGY_FILE, timeStamp);
	
	FILE *enrgFile = fopen(name, "wb");
	
	fprintf(enrgFile, "t (%6g s), Energy\n", DELTA);
	
	return enrgFile;
}

FILE *creatEnergyDistr() {
	char name[sizeof(ENERGY_DISTR_FILE) - 2 + TIME_STAMP_LENGTH] = "";
	sprintf(name, ENERGY_DISTR_FILE, timeStamp);
	
	FILE *energyDistrFile = fopen(name, "wb");
	
	fprintf(energyDistrFile, "t (%6g s)", DELTA);
	for (int i = 0; i <= MAX_ENERGY / ENERGY_PRECISION; i++) {
		fprintf(energyDistrFile, ", %8.2g", i * ENERGY_PRECISION);
	}
	fprintf(energyDistrFile, "\n");
	
	return energyDistrFile;
}

FILE *creatPressVsTemp() {
	char name[sizeof(PRESS_VS_TEMP_FILE) - 2 + TIME_STAMP_LENGTH] = "";
	sprintf(name, PRESS_VS_TEMP_FILE, timeStamp);
	
	FILE *pressVsTempFile = fopen(name, "wb");
	
	fprintf(pressVsTempFile, "T(1K)\t\t, Pressure\n");
	
	return pressVsTempFile;
}

FILE *creatFluctsVsTemp() {
	char name[sizeof(FLUCTUATIONS_VS_TEMP) - 2 + TIME_STAMP_LENGTH] = "";
	sprintf(name, FLUCTUATIONS_VS_TEMP, timeStamp);
	
	FILE *fluctsVsTempFile = fopen(name, "wb");
	
	fprintf(fluctsVsTempFile, "T(1K)\t\t, Fluctuations, Mean\t\t\n");
	
	return fluctsVsTempFile;
}

FILE *creatPressFlucts() {
	char name[sizeof(PRESS_FLUCTS_FILE) - 2 + TIME_STAMP_LENGTH] = "";
	sprintf(name, PRESS_FLUCTS_FILE, timeStamp);
	
	FILE *pressFluctsFile= fopen(name, "wb");
	
	fprintf(pressFluctsFile, "T(1K)\t\t, Fluctuations, Mean\t\t\n");
	
	return pressFluctsFile;
}

void printConsts() {
	char name[sizeof(CONSTS_FILE) - 2 + TIME_STAMP_LENGTH] = "";
	sprintf(name, CONSTS_FILE, timeStamp);
	FILE *consts = fopen(name, "wb");
	
	fprintf(consts, "MOLECULES_NUM, %d\n", MOLECULES_NUM);
	fprintf(consts, "MOLECULE_MASS, %g\n", MOLECULE_MASS);
	fprintf(consts, "MASSIVE_MOLECULE_MASS, %g\n", MASSIVE_MOLECULE_MASS);
	fprintf(consts, "MOLECULE_RADIUS, %g\n", MOLECULE_RADIUS);
	fprintf(consts, "\n");
	
	fprintf(consts, "TEMPERATURE, %g\n", TEMPERATURE);
	
	fprintf(consts, "CONTAINER_LENGTH, %g\n", CONTAINER_LENGTH);
	fprintf(consts, "CONTAINER_WIDTH, %g\n", CONTAINER_WIDTH);
	fprintf(consts, "CONTAINER_HEIGHT, %g\n", CONTAINER_HEIGHT);
	fprintf(consts, "\n");
	
	fprintf(consts, "VOLUME_LENGTH, %g\n", VOLUME_LENGTH);
	fprintf(consts, "VOLUME_WIDTH, %g\n", VOLUME_WIDTH);
	fprintf(consts, "VOLUME_HEIGHT, %g\n", VOLUME_HEIGHT);
	fprintf(consts, "\n");
	
	fprintf(consts, "VOLUME_X, %g\n", VOLUME_X);
	fprintf(consts, "VOLUME_Y, %g\n", VOLUME_Y);
	fprintf(consts, "VOLUME_Z, %g\n", VOLUME_Z);
	fprintf(consts, "\n");
	
	fprintf(consts, "DELTA, %g\n", DELTA);
	fprintf(consts, "\n");
	
	fprintf(consts, "TICKS_AVERAGE, %i\n", TICKS_AVERAGE);
	fprintf(consts, "PRESSURE_TICKS_AVERAGE, %i\n", PRESSURE_TICKS_AVERAGE);
	fprintf(consts, "\n");
	
	fclose(consts);
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
	
	fprintf(velFile, "%12ld", time);
	
	for (int i = 0; i <= MAX_VELOCITY / VELOCITY_PRECISION; i++) {
		fprintf(velFile, ", %6d", count[i]);
	}
	
	fprintf(velFile, "\n");
	
	delete[] count;
}

void printPressure(long time, double pressure, FILE *presFile) {
	fprintf(presFile, "%12ld, %12g\n", time, pressure);
}

void printTime(double temperature, long time, FILE *timeFile) {
	fprintf(timeFile, "%12g, %12ld\n", temperature, time);
}

void printEnergy(long time, double energy, FILE *enrgFile) {
	fprintf(enrgFile, "%12ld, %12g\n", time, energy);
}

void printEnergyDistribution(long time, double *data, FILE *velFile) {
	int *count = new int[(int) (MAX_ENERGY / ENERGY_PRECISION) + 1] { 0 };
	
	for (int i = 0; i < MOLECULES_NUM; i++) {
		if (data[i] > MAX_ENERGY) {
			count[(int) (MAX_ENERGY / ENERGY_PRECISION)]++;
		} else {
			count[(int) (data[i] / ENERGY_PRECISION)]++;
		}
	}
	
	fprintf(velFile, "%12ld", time);
	
	for (int i = 0; i <= MAX_ENERGY / ENERGY_PRECISION; i++) {
		fprintf(velFile, ", %8d", count[i]);
	}
	
	fprintf(velFile, "\n");
	
	delete[] count;
}

void printPressureVsTemperature(double temp, double press, FILE *presVsTempFile) {
	fprintf(presVsTempFile, "%12g, %12g\n", temp, press);
}

void printFluctuations(double temp, double fluctuation, double mean, FILE *fluctsVsTempFile) {
	fprintf(fluctsVsTempFile, "%12g, %12g, %12g\n", temp, fluctuation, mean);
}

void printPressureFluctuations(double temp, double fluctuation, double mean, FILE *pressFluctsFile) {
	fprintf(pressFluctsFile, "%12g, %12g, %12g\n", temp, fluctuation, mean);
}

double mean(const int *data, int n) {
	assert(data);
	
	double meanN = 0;
	for (int i = 0; i < n; i++) {
		meanN += data[i];
	}
	meanN /= n;
	
	return meanN;
}

double mean(const double *data, int n) {
	assert(data);
	
	double meanN = 0;
	for (int i = 0; i < n; i++) {
		meanN += data[i];
	}
	meanN /= n;
	
	return meanN;
}

double fluctuations(const int *data, int n) {
	assert(data);
	
	double sum = 0;
	double sumSqr = 0;
	
	for (int i = 0; i < n; i++) {
		sum += data[i];
		sumSqr += data[i] * data[i];
	}
	
	return sqrt((sumSqr - sum * sum / n) / n);
}

double fluctuations(const double *data, int n) {
	assert(data);
	
	double sum = 0;
	double sumSqr = 0;
	
	for (int i = 0; i < n; i++) {
		sum += data[i];
		sumSqr += data[i] * data[i];
	}
	
	return sqrt((sumSqr - sum * sum / n) / n);
}

double maxDifference(const double *data, const double *prevData, int n) {
	assert(data);
	assert(prevData);
	
	double max = 0;
	
	for (int i = 0; i < n; i++) {
		if (fabs(data[i] - prevData[i]) > max)
			max = fabs(data[i] - prevData[i]);
	}
	
	return max;
}

int maxDifference(const int *data, const int *prevData, int n) {
	assert(data);
	assert(prevData);
	
	int max = 0;
	
	for (int i = 0; i < n; i++) {
		if (abs(data[i] - prevData[i]) > max)
			max = abs(data[i] - prevData[i]);
	}
	
	return max;
}

void getEnergyDistribution(const double *data, int n, int *dest) {
	assert(data);
	assert(dest);
	
	for (int i = 0; i < n; i++) {
		if (data[i] > MAX_ENERGY) {
			dest[(int) (MAX_ENERGY / ENERGY_PRECISION)]++;
		} else {
			dest[(int) (data[i] / ENERGY_PRECISION)]++;
		}
	}
}
