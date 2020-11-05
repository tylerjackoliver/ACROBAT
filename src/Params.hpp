#ifndef __PARAMS_H__
#define __PARAMS_H__

#include <cmath>

namespace PARAMS
{
    /* ~~~~~~ User-defined parameters first ~~~~~~ */
    /* Eccentricity of the ballistic capture orbit */ 
    const double ECC = 0.95;
    /* Inclination of the ballistic capture plane */
    const double INC   = 45*std::atan(1.0)*4.0 / 180.;
    /* Longitude of the ballistic capture orbit */
    const double LONGTD= 0.0;
    /* Epoch of the transfer */
    const double EPOCH = 24023040.;
    /* Common identifier or SPK ID of the TARGET */
    std::string TARGET="Earth";
    /* Common identifier or SPK ID of the HOST */
    std::string HOST = "Sun";
    
    /* ~~~~~~ Derived parameters ~~~~~~ */
    /* GM of TARGET */
    double targetGM = 0.0; // E.g. Earth
    /* GM of the HOST */
    double hostGM = 0.0; // E.g. Sun
    /* Planetary radius of the TARGET */
    double R = 0.0; // Planetary radius, in km
    /* Sphere of influence of the TARGET */
    double RS = 0.0; // Sphere of influence
    /* Mean anomaly of TARGET around HOST */
    double M     = 0.0;
    /* Maximum integration time for determining a trajectory to be acrobatic */
    double maxT = 0.0;
}

#endif