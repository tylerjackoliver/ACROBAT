#ifndef __EQUATIONS_OF_MOTION_H__
#define __EQUATIONS_OF_MOTION_H__

#include <eigen3/Eigen/Core>
#include "OEs.hpp"
#include "Params.hpp"
#include <cmath>
#include <limits>
#include <boost/numeric/odeint.hpp>
#include <stdexcept>

extern "C"
{
    #include <SpiceUsr.h>
}

/* @brief Computes the direction vector to the Sun from the spacecraft.
 * @param[in] r: Vector to the Sun; Eigen::Matrix<Type, Size>
 * @param[in] t: Time at which this vector is desired; const Type
 */
template <typename Type>
void getSunVector(Eigen::Matrix<Type, 3, 1> &r, const Type t)
{
    Eigen::Matrix3<Type, 3> rotationVector;
    double trueAnomaly;

    // Get orbital elements of host about target in the EME2000 frame (actually the J2000 frame in SPICE!)
    ACROBAT::OEs sunOEs = getOEs(t, PARAMS::host, PARAMS::target, PARAMS::targetGM);      // Need orientation of the Sun wrt planet we're seeking ballistic capture orbits for
    ACROBAT::OEs planetOEs = getOEs(t, PARAMS::target, PARAMS::host, PARAMS::hostGM);   // Need a, e of planet we're constructing ballistic orbits for

    // Get the true anomaly of the Sun about the target planet
    meanToTrue(sunOEs.M, sunOEs.ecc, &trueAnomaly);

    // Construct rotation vector
    getRotationMatrix(&rotationVector, &sunOEs, &trueAnomaly);

    // Compute vector
    double a = planetOEs.rp / (1. - planetOEs.ecc);
    double scalar = a * (1 - planetOEs.e) / (1 + planetOEs.ecc * std::cos(trueAnomaly));

    r = scalar * rotationVector;
}

/* @brief C++ wrapper to the SPICE orbital element routines.
 * @param[in] t: Time (ephemeris seconds) at which the orbital elements are required.
 * @param[in] body: String identifying the body the orbital elements are sought for.
 * @param[in] obs: String identifying the body the orbital elements are defined about/
 * @returns ACROBAT::OEs struct, holding the orbital elements of body about obs at time t.
 */
ACROBAT::OEs getOEs(double t, std::string body, std::string obs, double &gm)
{
    ACROBAT::OEs ret;

    std::string ref = "IAU_" + obs;     // Reference frame for the BME field is the IAU frame about the target
    
    ConstSpiceChar abcorr[] = "NONE";
    
    SpiceDouble state[6];
    SpiceDouble elts[8];
    SpiceDouble et = t;
    SpiceDouble lt = 0.0;

    spkezr_c(body.c_str(), et, ref.c_str(), abcorr, obs.c_str(), state, &lt);
    oscelt_c(state, et, gm, elts);

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

/* @brief Get the rotation matrix to move from the EME2000 to point at the Sun
 * @param[inout] rotationVector: Eigen::Matrix<Type, 3, 1> containing the rotation vector
 * @param[in] sunOEs: Reference to an OEs struct containing the orbital elements of the Sun
 * @param[in] trueAnomaly: Current true anomaly of planet around the Sun.
 */
template <typename Type>
void getRotationMatrix(Eigen::Matrix<Type, 3, 1> &rotationVector, ACROBAT::OEs &sunOEs, double &trueAnomaly)
{
    double theta = sunOEs.omega + trueAnomaly;
    
    rotationVector(0) = std::sin(theta) * std::sin(sunOEs.longtd) * std::cos(sunOEs.inc) - std::cos(theta) * std::cos(sunOEs.longtd);
    rotationVector(1) = - std::cos(theta) * std::sin(sunOEs.longtd) - std::sin(theta) * std::cos(sunOEs.longtd) * std::cos(sunOEs.inc);
    rotationVector(2) = -std::sin(theta) * std::sin(sunOEs.inc);
}

/* @brief Compute the derivative vector of the current state for numerical integration routines.
*  @param[in] x Eigen::Matrix<Type, 6> of the initial [r, v] vector for the particle
*  @param[out] dx Eigen::Matrix<Type, 6> corresponding to the derivative of x
*  @param[in] t Current integration time-step
*/
template <typename Type>
void forceFunction(Eigen::Matrix<Type, 6, 1> &x, Eigen::Matrix<Type, 6, 1> &dx, const double t)
{
    // First, the trivial derivatives
    dx(0) = x(3);
    dx(1) = x(4);
    dx(2) = x(5);

    // Now the not-so-trivial
    Eigen::Matrix<Type, 3, 1> sunVector, positionVector, velocityVector, solarTerm;

    positionVector(0) = x(0);
    positionVector(1) = x(1);
    positionVector(2) = x(2);

    double currentTime = PARAMS::EPOCH + t; // t is measured in seconds

    getSunVector(&sunVector, currentTime);
    posDifference = positionVector - sunVector;

    solarTerm = PARAMS::hostGM * (sunVector / ( std::pow(sunVector.norm(), 3) ) + ( posDifference / ( std::pow(posDifference.norm(), 3) ) ));
    velocityVector = -PARAMS::targetGM * positionVector / ( std::pow(positionVector.norm(), 3) ) - solarTerm;

    dx(3) = velocityVector(0);
    dx(4) = velocityVector(1);
    dx(5) = velocityVector(2);
}

/* @brief Custom integration observer function; returns various status integers depending on whether stopping conditions have been triggered.
   @param[in] x Eigen::Matrix<Type, 6> containing the state vector of the particle at the given point.
   @param[in] x0 Eigen::Matrix<Type, 6> containing the initial position of the particle on the given trajectory
   @param[in] t Const double containing the current time-step of the integration (assumed seconds)
   @returns An integer corresponding to one of four separate events.

    Integer Return Codes
    ~~~~~~~~~~~~~~~~~~~~
    0: Do nothing!
    1: Trajectory has crashed (R < Rmin)
    2: Trajectory has escaped (R > Rs)
    3: Trajectory is weakly stable (wacky geometrics)
    4: Trajectory is acrobatic (nothing happens!)
*/
template <typename Type>
int integrationController(Eigen::Matrix<Type, 6, 1> &x, Eigen::Matrix<Type, 6, 1> &x0, const double t)
{
    int flag = 0;
    Eigen::Matrix<Type, 3, 1> r = x(Eigen::seq(0,2));
    Eigen::Matrix<Type, 3, 1> r0 = x0(Eigen::seq(0,2));
    Eigen::Matrix<Type, 3, 1> v = x(Eigen::seq(3,5));
    Eigen::Matrix<Type, 3, 1> v0 = x0(Eigen::seq(3,5));

    // First check for a crash
    // Note: we have not yet regularised the dynamics, so for now we are checking that R < Rrad (not 1)
    Type rMag = r.norm();
    if (rMag <= PARAMS::R) return 1;    // Premature exit to prevent computing unnecessary branches

    // Next, check for an escape
    Type vMag = v.norm();
    Type keplerEnergy = (vMag * vMag / 2. - 1./rMag);
    if (rMag >= PARAMS::RS && H > 0) return 2;

    // Now, check for weakly stable
    Eigen::Matrix<Type, 3, 1> angMomentum0 = r0.cross(v0);
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
