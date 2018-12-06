
#pragma once

#define M_PI       3.14159265358979323846   // pi

#include <cmath>

#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <iostream>

#ifndef uint
typedef unsigned int uint;
#endif

#ifndef byte
typedef unsigned char byte;
#endif

// ***********
// * Logging

// error log -- silent; will not stop program
std::ostream &err();

#ifndef ENSURE
#define ENSURE(STATEMENT) { \
    if(!(STATEMENT)) { \
        std::cerr << "ENSURE FAILED at " \
                  << __FILE__ << ", line #" << __LINE__ << ":\n" \
                  << "    " << #STATEMENT << std::endl; \
        err()     << "ENSURE FAILED at " \
                  << __FILE__ << ", line #" << __LINE__ << ":\n" \
                  << "    " << #STATEMENT << std::endl; \
        exit(1); \
    } \
}
#endif // ENSURE

// Use ERROR to print an error message tagged with the given file/line #
#ifndef CORK_ERROR
#define CORK_ERROR(message) { \
    std::cerr << "error at " \
              << __FILE__ << ", line #" << __LINE__ << ": " \
              << (message) << std::endl; \
    err()     << "error at " \
              << __FILE__ << ", line #" << __LINE__ << ": " \
              << (message) << std::endl; \
}
#endif // CORK_ERROR

// Use MARKER for debugging to create a trace of control flow...
#ifndef MARKER
#define MARKER(message) { \
    std::cout << "marker at " \
              << __FILE__ << ", line #" << __LINE__ << ": " \
              << (message) << std::endl; \
}
#endif // MARKER

// ***********
// * Assorted

// snap the value a into the specified range
inline double clamp(double a, double mina, double maxa) {
	return std::min(maxa, std::max(mina, a));
}
inline float  clamp(float  a, float  mina, float  maxa) {
	return std::min(maxa, std::max(mina, a));
}

// modulo the value a into the specified range
inline double wrap(double a, double mina, double maxa) {
	double val = std::fmod(a - mina, maxa - mina);
	if (val < 0.0) val += maxa - mina;
	return val + mina;
}
inline float  wrap(float  a, float  mina, float  maxa) {
	float val = std::fmod(a - mina, maxa - mina);
	if (val < 0.0) val += maxa - mina;
	return val + mina;
}

inline double deg2rad(double deg) {
	return (M_PI / 180.0) * deg;
}
inline double rad2deg(double rad) {
	return (180.0 / M_PI) * rad;
}

// **********
// * Timing

#ifdef _WIN32
#include <winsock.h>
#endif

class Timer {
public:
	Timer(); // automatically start timer on creation
	~Timer();
public:
	void start();
	// returns lap time in milliseconds
	double lap();
	// returns total ellapsed time in milliseconds
	double stop();
	// re-retrieve values
	// if this operation is meaningless, you will receive 0.0 instead
	double lastLap() const;
	double ellapsed() const;
private:
	timeval init;
	timeval prev;
	double last_lap_duration;
	double ellapsed_time_duration;
	bool isRunning;
};

// ***********
// * Random

// need functions that allow for the random source to be made more deterministic for replays...

inline void initRand() {
	// currently none!  Should seed using clock
	srand(uint(time(0)));
}

inline double drand(double min, double max) {
	const double invMAX = 1.0 / double(RAND_MAX);
	double rand0to1 = double(std::rand())*invMAX;
	return (max - min)*rand0to1 + min;
}

inline uint randMod(uint range) {
	return std::rand() % range;
}
