#ifndef __PARAMS_H__
#define __PARAMS_H__

#include <cmath>
#include <vector>
#include <string>
extern "C"
{
    #include <SpiceUsr.h>
}

namespace PARAMS
{
    /* ~~~~~~ User-defined parameters first ~~~~~~ */
    /* Eccentricity of the ballistic capture orbit */ 
    const double ECC = 0.95;
    /* Inclination of the ballistic capture plane */
    const double INC   = 45.04*std::atan(1.0)*4.0 / 180.;
    /* Longitude of the ballistic capture orbit */
    const double LONGTD= 202.5 * 4.0 * std::atan(1.0) / 180.;
    /* Epoch of the transfer */
    const double EPOCH = 634780800;
    /* Common identifier or SPK ID of the TARGET */
    std::string TARGET="Mercury";
    /* Common identifier or SPK ID of the HOST */
    std::string HOST = "Sun";
    /* Additional planets to consider */
    std::vector<std::string> additionalPlanets = {"Venus", "Jupiter", "Saturn"};
    
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
    /* GMs of other planets */
    std::vector<double> additionalPlanetsGM;  
}

#endif
