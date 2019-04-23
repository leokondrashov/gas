#ifndef GAS_CONSTANTS_H
#define GAS_CONSTANTS_H

#define kBoltzmann 1.38e-23

#define MOLECULES_NUM 1000
#define MOLECULE_MASS (28e-3/6.02e23)
#define MOLECULE_RADIUS 1.6e-10

#define CONTAINER_LENGTH 1.0e-8
#define CONTAINER_WIDTH 1.0e-8
#define CONTAINER_HEIGHT 1.0e-8

#define VOLUME_LENGTH 1.0e-7
#define VOLUME_WIDTH 1.0e-7
#define VOLUME_HEIGHT 1.0e-7

#define VOLUME_X 0.5e-6
#define VOLUME_Y 0.5e-6
#define VOLUME_Z 0.5e-6

#define GAS_TEMPERATURE_START 300.0
#define GAS_TEMPERATURE_STEP 20.0
#define GAS_TEMPERATURE_END 500.0

#define DELTA 1e-13

#define VELOCITY_PRECISION 20.0
#define MAX_VELOCITY 1000.0

#define CONSTS_FILE "log/consts%s.csv"
#define VELOCITY_FILE "log/velocities%s.csv"
#define PRESSURE_FILE "log/pressure%s.csv"
#define PRESS_VS_TEMP_FILE "log/pressVsTemp%s.csv"
#define PRESS_FLUCTS_FILE "log/pressFlucts%s.csv"
#define FLUCTUATIONS_VS_TEMP "log/flutcsVsTemp%s.csv"

#define EXPERIMENT_LENGTH 100000
#define TICKS_AVERAGE 1000
#define PRESSURE_TICKS_AVERAGE 1000

#endif //GAS_CONSTANTS_H
