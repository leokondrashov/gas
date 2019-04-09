#ifndef GAS_CONSTANTS_H
#define GAS_CONSTANTS_H

#define kBoltzmann 1.38e-23

#define MOLECULES_NUM 1000
#define MOLECULE_MASS (28e-3/6.02e23)

#define CONTAINER_LENGTH 1.0e-6
#define CONTAINER_WIDTH 1.0e-6
#define CONTAINER_HEIGHT 1.0e-6

#define VOLUME_LENGTH 1.0e-8
#define VOLUME_WIDTH 1.0e-6
#define VOLUME_HEIGHT 1.0e-6

#define VOLUME_X 0.5e-6
#define VOLUME_Y 0.5e-6
#define VOLUME_Z 0.5e-6

#define GAS_TEMPERATURE_START 300.0
#define GAS_TEMPERATURE_STEP 20.0
#define GAS_TEMPERATURE_END 500.0

#define DELTA 1e-12

#define VELOCITY_PRECISION 20.0
#define MAX_VELOCITY 1000.0

#define VELOCITY_FILE "velocities.csv"
#define PRESSURE_FILE "pressure.csv"
#define PRESS_VS_TEMP_FILE "pressVsTemp.csv"
#define FLUCTUATIONS_VS_TEMP "flutcsVsTemp.csv"

#define TICKS_AVERAGE 10000

#endif //GAS_CONSTANTS_H
