#ifndef __EQUATIONS_OF_MOTION_H__
#define __EQUATIONS_OF_MOTION_H__

#include <eigen3/Core>
#include "OEs.hpp"
#include "Params.hpp"

extern "C"
{
    #include <SpiceUsr.h>
}

template <typename Type>
using Vector = Eigen::Vector<Type, 6>;

/* Computes the direction vector to the Sun */
template <typename Type>
void getSunVector(Vector r, const Type t)
{
    Vector rotationVector;
    double trueAnomaly;

    // Get orbital elements of host about target in the EME2000 frame (actually the J2000 frame in SPICE!)
    SCROTAL::OEs planetOEs = getOEs(const Type t, PARAMS::host, PARAMS::target); // Host, target

    // Construct rotation vector
    getRotationMatrix(&rotationMatrix);

    // Get true anomaly
    meanToTrue(&planetOEs.M, &planetOEs.ecc, &trueAnomaly);

    // Compute vector
    double scalar = planetOEs.a * (1 - planetOEs.e) / (1 + planetOEs.ecc * std::cos(trueAnomaly));

    r = scalar * rotationVector;
}

/* Gets the orbital elements at time t of an object target about a host host */
template <typename Type>
void getOEs(const Type t, std::string host, std::string target)
{

}

#endif