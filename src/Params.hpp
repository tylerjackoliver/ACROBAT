#ifndef __PARAMS_H__
#define __PARAMS_H__

namespace PARAMS
{
    double ECC   = 0.95;
    double INC   = 30*std::atan(1.0)*4.0 / 180.;
    double LONGTD= 0.0;
    double M     = 0.0;
    double EPOCH = 24023040.;
    double targetGM = 1.327e11; // E.g. Earth
    double hostGM = 1.e13; // E.g. Sun
    double R = 6367.; // Planetary radius, in km
    double RS = 10000.; // Sphere of influence
}

#endif