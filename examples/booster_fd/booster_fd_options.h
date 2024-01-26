#ifndef BOOSTER_FD_OPTIONS_H_
#define BOOSTER_FD_OPTIONS_H_
#include <string>

struct Booster_fd_options {
    int gridx = 32;
    int gridy = 32;
    int gridz = 32;
    // at end of pip2_realistic_l6aperture
    double real_particles = 82716049382.71603;
    int num_particles = 5837;
    int turns = 100;
    int steps = 24 * 6; // 6 steps per period
    double pipesizex = 0.14986;
    double pipesizey = 0.077216;
    double pipesizez = 474.202752 / 84;
    int commsize;
};
#endif /* BOOSTER_FD_OPTIONS_H_ */
