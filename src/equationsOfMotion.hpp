#ifndef __EQUATIONS_OF_MOTION_H__
#define __EQUATIONS_OF_MOTION_H__

#include <eigen3/Core>
#include "OEs.hpp"
#include "Params.hpp"
#include <cmath>
#include <limits>

extern "C"
{
    #include <SpiceUsr.h>
}

/* Computes the direction vector to the Sun */
template <typename Type>
void getSunVector(Eigen::Vector<Type, 3> &r, const Type t)
{
    Eigen::Vector3<Type, 3> rotationVector;
    double trueAnomaly;

    // Get orbital elements of host about target in the EME2000 frame (actually the J2000 frame in SPICE!)
    SCROTAL::OEs sunOEs = getOEs(t, PARAMS::host, PARAMS::target); // Need orientation of the Sun wrt planet we're seeking ballistic capture orbits for
    SCROTAL::OEs planetOEs = getOEs(t, PARAMS::target, PARAMS::host); // Need a, e of planet we're constructing ballistic orbits for

    // Get the true anomaly of the Sun about the target planet
    meanToTrue(sunOEs.M, sunOEs.ecc, &trueAnomaly);

    // Construct rotation vector
    getRotationMatrix(&rotationVector, &sunOEs, &trueAnomaly);

    // Compute vector
    double a = planetOEs.rp / (1. - planetOEs.ecc);
    double scalar = a * (1 - planetOEs.e) / (1 + planetOEs.ecc * std::cos(trueAnomaly));

    r = scalar * rotationVector;
}

/* C++ wrapper to the C implementation */
SCROTAL::OEs getOEs(double t, std::string body, std::string obs)
{
    SCROTAL::OEs ret;
    extern "C"
    {
        ConstSpiceChar body[] = "Sun";
        ConstSpiceChar ref[] = "J2000";
        ConstSpiceChar abcorr[] = "NONE";
        ConstSpiceChar obs[] = target.c_str();
        SpiceDouble state[6];
        SpiceDouble elts[8];
        SpiceDouble et = t;
        SpiceDouble lt = 0.0;

        spkezr_c(body, et, ref, abcorr, obs, starg, &lt);
        oscelt_c(state, et, PARAMS::GM, elts);
    }
    // Copy elts into vector
    ret.rp = elts[0];
    ret.ecc = elts[1];
    ret.inc = elts[2];
    ret.longtd = elts[3];
    ret.omega = elts[4];
    ret.M = elts[5];
    ret.epoch = elts[6];

    return ret;
}

template <typename Type>
void getRotationMatrix(Eigen::Vector<Type, 3> &rotationVector, SCROTAL::OEs &sunOEs, double &trueAnomaly)
{
    double theta = sunOEs.omega + trueAnomaly;
    
    rotationVector(0) = std::sin(theta) * std::sin(sunOEs.longtd) * std::cos(sunOEs.inc) - std::cos(theta) * std::cos(sunOEs.longtd);
    rotationVector(1) = - std::cos(theta) * std::sin(sunOEs.longtd) - std::sin(theta) * std::cos(sunOEs.longtd) * std::cos(sunOEs.inc);
    rotationVector(2) = -std::sin(theta) * std::sin(sunOEs.inc);
}

template <typename Type>
void forceFunction(Eigen::Vector<Type, 6> &x, Eigen::Vector<Type, 6> &dx, const double t)
{
    // First, the trivial derivatives
    dx(0) = x(3);
    dx(1) = x(4);
    dx(2) = x(5);

    // Now the not-so-trivial
    Eigen::Vector<Type, 3> sunVector, positionVector, velocityVector, solarTerm;

    positionVector(0) = x(0);
    positionVector(1) = x(1);
    positionVector(2) = x(2);

    getSunVector(&sunVector, t);
    posDifference = positionVector - sunVector;

    solarTerm = PARAMS::hostGM * (sunVector / ( std::pow(sunVector.norm(), 3) ) + ( posDifference / ( std::pow(posDifference.norm(), 3) ) ));
    velocityVector = -PARAMS::targetGM / positionVector.norm() - solarTerm;

    dx(3) = velocityVector(0);
    dx(4) = velocityVector(1);
    dx(5) = velocityVector(2);
}

/* Controls the ingration; returns various status integers depending on whether conditions are satisfied.

Integer
~~~~~~~
0: Do nothing!
1: Trajectory has crashed (R < Rmin)
2: Trajectory has escaped (R > Rs)
3: Trajectory is weakly stable (wacky geometrics)
4: Trajectory is acrobatic (nothing happens!)

*/
template <typename Type>
int integrationController(Eigen::Vector<Type, 6> &x, Eigen::Vector<Type, 6> &x0, const double t)
{
    int flag = 0;
    Eigen::Vector<Type, 3> r = x(Eigen::seq(0,2));
    Eigen::Vector<Type, 3> r0 = x0(Eigen::seq(0,2));
    Eigen::Vector<Type, 3> v = x(Eigen::seq(3,5));
    Eigen::Vector<Type, 3> v0 = x0(Eigen::seq(3,5));

    // First check for a crash
    // Note: we have not yet regularised the dynamics, so for now we are checking that R < Rrad (not 1)
    Type rMag = r.norm();
    if (rMag <= PARAMS::R) return 1;    // Premature exit to prevent computing unnecessary branches

    // Next, check for an escape
    Type vMag = v.norm();
    Type keplerEnergy = (vMag * vMag / 2. - 1./rMag);
    if (rMag >= PARAMS::RS && H > 0) return 2;

    // Now, check for weakly stable
    Eigen::Vector<Type, 3> angMomentum0 = r0.cross(v0);
    Type conditionOne = r.dot(angMomentum0.cross(v0));
    Type conditionTwo = r.dot(r0);
    Type conditionThree = v.dot(v0) * v0.dot(v0); // Supposed to be v.dot(v0) * v(k-1).dot(v0), but only doing one at a time here

    if (conditionOne < std::numeric_limits<Type>::epsilon() * 1000 && conditionTwo > 0 && conditionThree > 0) return 3;

    // Lastly, check time
    double pi = 4.0 * std::atan(1.0);
    if (t >= 8.0 * pi * std::pow(PARAMS::RS, 1.5) ) return 4;

    // If nothing else
    return 0;
}
#endif