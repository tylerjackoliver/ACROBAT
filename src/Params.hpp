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
    // const double INC = 45.00 * std::atan(1.0) * 4.0 / 180.;
    /* Longitude of the ballistic capture orbit */
    const double LONGTD= 202.5 * 4.0 * std::atan(1.0) / 180.;
    // const double LONGTD = 233.82 * std::atan(1.0) * 4.0 / 180.; 
    /* Epoch of the transfer */
    const double EPOCH = 634798080; // Mercury
    // const double EPOCH = 631341216; //Europa
    // const double EPOCH = 634506048;// Earth
    /* Common identifier or SPK ID of the TARGET */
    std::string TARGET="Mercury";
    /* Common identifier or SPK ID of the HOST */
    std::string HOST = "Sun";
    /* Additional planets to consider */
    std::vector<std::string> additionalPlanets = {"Venus", "Jupiter", "Saturn", "Io", "Europa", "Ganymede", "Callisto", "Titan", "Titania", "Charon"};
    // std::vector<std::string> additionalPlanets = {"Sun",
    // "Saturn", "Io", "Ganymede", "Callisto"}; // Europa study
    // std::vector<std::string> additionalPlanets = {"Moon", "Jupiter", "Saturn"};

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
