#ifndef __EQUATIONS_OF_MOTION_H__
#define __EQUATIONS_OF_MOTION_H__

#include <eigen3/Eigen/Core>
#include "OEs.hpp"
#include "Params.hpp"
#include <cmath>
#include <limits>
#include <boost/numeric/odeint.hpp>
#include <stdexcept>
#include "laguerreConway.hpp"
#include <cmath>
#include <string>

extern "C"
{
    #include <SpiceUsr.h>
}

/* @brief Computes and returns the dimensional time from the non-dimensional time.
 * @param[in] nonDim Non-dimensional time to convert
 * @returns The dimensional time corresponding to nonDim
 */
template <typename Type> 
Type getDimensionalTime(const Type &nonDim)
{
	if ( std::abs(PARAMS::targetGM) < 1e-6 ) throw std::domain_error("Error: target GM is zero in getDimensionalTime.");
    return nonDim * std::sqrt( std::pow(PARAMS::R, 3) / PARAMS::targetGM );
}

/* @brief Computes and returns the the non-dimensional time from the dimensional time
 * @param[in] dim Dimensional time to convert
 * @returns The non-dimensional time corresponding to dim
 */
template <typename Type> 
Type getNonDimensionalTime(const Type &dim)
{

	if ( std::abs(PARAMS::R) <  1e-06) throw std::domain_error("Error: R is zero in getNonDimensionalTime.");
    return dim * std::sqrt( PARAMS::targetGM / std::pow(PARAMS::R, 3) );
}

/* @brief Computes the non-dimensional position relating to the dimensional position
   @param[in] dim Dimensional position to convert
   @param[out] nonDim Non-dimensional position corresponding to dim
*/
template <typename Type>
void getNonDimPosition(Type& dim, Type& nonDim)
{
    nonDim = dim;
    if ( std::abs(PARAMS::R) < 1e-06 ) throw std::domain_error("Error: R is ~0 in getNonDimPosition.");
    for (size_t idx = 0; idx < nonDim.size(); ++idx) nonDim[idx] = nonDim[idx] / PARAMS::R;
}

/* @brief Computes the dimensional position relating to the non-dimensional position
   @param[in] nonDim Non-dimensional position to convert
   @param[out] dim Dimensional position corresponding to dim
*/
template <typename Type>
void getDimPosition(Type& nonDim, Type& dim)
{
    dim = nonDim;
    for (size_t idx = 0; idx < dim.size(); ++idx) dim[idx] = dim[idx] * PARAMS::R;
}

/* @brief Computes the dimensional velocity corresponding to the non-dimensional velocity
   @params[in] nonDim The non-dimensional velocity to convert.
   @params[out] dim The dimensional velocity corresponding to nonDim
*/
template <typename Type>
void getDimVelocity(Type& nonDim, Type& dim)
{
    dim = nonDim;
    if ( std::abs(PARAMS::R) < 1e-06 ) throw std::domain_error("Error: R is ~0 in getDimVelocity.");
    double coeff = std::sqrt(PARAMS::targetGM / PARAMS::R);
    for (size_t idx = 0; idx < dim.size(); ++idx) dim[idx] = dim[idx] * coeff;
}

/* @brief Computes the non-dimensional velocity corresponding to the dimensional velocity
   @params[in] dim The dimensional velocity to convert.
   @params[out] nonDim The non-dimensional velocity corresponding to dim
*/
template <typename Type> 
void getNonDimVelocity(Type& dim, Type& nonDim)
{
    nonDim = dim;
    if ( std::abs(PARAMS::targetGM) < 1e-06 ) throw std::domain_error("Error: targetGM is ~0 in getDimVelocity.");
    double coeff = std::sqrt(PARAMS::R / PARAMS::targetGM);
    for (size_t idx = 0; idx < nonDim.size(); ++idx) nonDim[idx] = nonDim[idx] * coeff;
}

/* @brief Makes a given input state non-dimensional
   @param[in] dim Dimensional state to non-dimensionalise
   @param[out] nonDim Non-dimensional state corresponding to dim
*/
template <typename Type>
void getNonDimState(std::vector<Type>& dim, std::vector<Type>& nonDim)
{
    std::vector<Type> pos = std::vector<Type>(dim.begin(), dim.begin() + 3);
    std::vector<Type> vel = std::vector<Type>(dim.end()-3, dim.end());

    std::vector<Type> nonDimPos;
    std::vector<Type> nonDimVel;

    getNonDimPosition(pos, nonDimPos);
    getNonDimVelocity(vel, nonDimVel);

    nonDim.insert(nonDim.begin(), nonDimPos.begin(), nonDimPos.end());
    nonDim.insert(nonDim.end(), nonDimVel.begin(), nonDimVel.end());
}

/* @brief Makes a given input state dimensional
   @param[in] nonDim Non-dimensional state to convert
   @param[out] dim Dimensional state corresponding to nonDim
*/
template <typename Type>
void getDimState(std::vector<Type>& nonDim, std::vector<Type>& dim)
{
    std::vector<Type> pos = std::vector<Type>(nonDim.begin(), nonDim.begin() + 3);
    std::vector<Type> vel = std::vector<Type>(nonDim.end()-3, nonDim.end());

    std::vector<Type> dimPos;
    std::vector<Type> dimVel;

    getDimPosition(pos, dimPos);
    getDimVelocity(vel, dimVel);

    dim.insert(dim.begin(), dimPos.begin(), dimPos.end());
    dim.insert(dim.end(), dimVel.begin(), dimVel.end());
}

/* @brief Computes the direction vector to the Sun from the spacecraft.
 * @param[in] r: Vector to the Sun; Eigen::Matrix<Type, Size>
 * @param[in] t: Time at which this vector is desired; const Type
 */
template <typename Type>
void getSunVector(Eigen::Matrix<Type, 3, 1> &r, const Type t)
{
    Eigen::Matrix<Type, 3, 1> rotationVector;
    double trueAnomaly;

    // Get orbital elements of host about target in the EME2000 frame (actually the J2000 frame in SPICE!)
    ACROBAT::OEs sunOEs = getOEs(t, PARAMS::HOST, PARAMS::TARGET, PARAMS::hostGM);      // Need orientation of the Sun wrt planet we're seeking ballistic capture orbits for
    ACROBAT::OEs planetOEs = getOEs(t, PARAMS::TARGET, PARAMS::HOST, PARAMS::hostGM);   // Need a, e of planet we're constructing ballistic orbits for

    // Get the true anomaly of the Sun about the target planet
    meanToTrue(sunOEs.M, sunOEs.ecc, trueAnomaly);

    // Construct rotation vector
    getRotationMatrix(rotationVector, sunOEs, trueAnomaly);

    // Compute vector
    double a = planetOEs.rp / (PARAMS::R * (1. - planetOEs.ecc)); // PARAMS::R normalizes
    double scalar = a * (1 - planetOEs.ecc) / (1 + planetOEs.ecc * std::cos(trueAnomaly));

    r = scalar * rotationVector;
}

/* @brief Get the vector to a given planet from the orbital ephemerides
   @param[in] t Dimensional time at which the vector is sought
   @param[inout] r Normalised position vector between the TARGET and the HOST
*/
template <typename Type>
void getPlanetVectorEphemeris(std::string& targ, Eigen::Matrix<Type, 3, 1> &r, const Type t)
{
    ConstSpiceChar abcorr[] = "NONE";
    
    SpiceDouble state[3];
    SpiceDouble et = t;
    SpiceDouble lt = 0.0;

    // Get state
    spkpos_c(targ.c_str(), et, "J2000", abcorr, PARAMS::TARGET.c_str(), state, &lt);

    for (size_t idx = 0; idx < 3; ++idx) r(idx) = state[idx] / PARAMS::R; // Copy position into output vector, normalised
}

/* @brief Get a dimensional version of the vector to a given planet from the orbital ephemerides
   @param[in] t Dimensional time at which the vector is sought
   @param[inout] r Normalised position vector between the TARGET and the HOST
*/
template <typename Type>
void getDimensionalPlanetVectorEphemeris(std::string& targ, Eigen::Matrix<Type, 3, 1> &r, const Type t)
{
    ConstSpiceChar abcorr[] = "NONE";
    
    SpiceDouble state[3];
    SpiceDouble et = t;
    SpiceDouble lt = 0.0;

    // Get state
    spkpos_c(targ.c_str(), et, "J2000", abcorr, PARAMS::TARGET.c_str(), state, &lt);

    for (size_t idx = 0; idx < 3; ++idx) r(idx) = state[idx]; // Copy position into output vector, normalised
}

/* @brief Get the dimensional state of a given planet, measured from the target
 * @param[in] targ A std::string containing the planet for which to get the vector to
 * @param[inout] r Normalised vector between the target and the host
 */
template <typename Type>
void getDimensionalPlanetState(std::string& targ, Eigen::Matrix<Type, 6, 1> &r, const Type t)
{
    ConstSpiceChar abcorr[] = "NONE";
    SpiceDouble state[6];
    SpiceDouble lt = 0.0;
    SpiceDouble et = t;

    // Get the state
    spkezr_c(targ.c_str(), et, "J2000", abcorr, PARAMS::TARGET.c_str(), state, &lt);

    for (size_t idx = 0; idx < 6; ++idx) r(idx) = state[idx];
}

/* @brief Get the dimensional velocity of a given planet, measured from the target
 * @param[in] targ A std::string containing the planet for which to get the vector to
 * @param[inout] r Normalised vector between the target and the host
 */
template <typename Type>
void getDimensionalPlanetVelocity(std::string& targ, Eigen::Matrix<Type, 3, 1>& r, const Type t)
{
    ConstSpiceChar abcorr[] = "NONE";
    SpiceDouble state[6];
    SpiceDouble lt = 0.0;
    SpiceDouble et = t;

    // Get the state
    spkezr_c(targ.c_str(), et, "J2000", abcorr, PARAMS::TARGET.c_str(), state, &lt);

    for (size_t idx = 0; idx < 3; ++idx) r(idx) = state[idx + 3];
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

    ConstSpiceChar abcorr[] = "NONE";
    
    SpiceDouble state[6];
    SpiceDouble elts[8];
    SpiceDouble et = t;
    SpiceDouble lt = 0.0;

    spkezr_c(body.c_str(), et, "J2000", abcorr, obs.c_str(), state, &lt);
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
*  @param[in] x std::vector<Type> of the initial [r, v] vector for the particle
*  @param[out] dx std::vector<Type> corresponding to the derivative of x
*  @param[in] t Current integration time-step
*/
void forceFunction(const std::vector<double> &x, std::vector<double> &dx, const double t)
{
    Eigen::Matrix<double, 3, 1> sunVector, positionVector, velocityVector, solarTerm;
    // First, the trivial derivatives
    for (size_t idx = 0; idx < 3; ++idx)
    {
        positionVector(idx) = x[idx];
        dx[idx] = x[idx+3];
    }

    // Now the not-so-trivial

    double dimensionalTime = getDimensionalTime(t);
    double currentTime = PARAMS::EPOCH + dimensionalTime; // t is measured in seconds

    getPlanetVectorEphemeris(PARAMS::HOST, sunVector, currentTime);
    Eigen::Matrix<double, 3, 1> posDifference = positionVector - sunVector;

    double normalHostGM = PARAMS::hostGM / PARAMS::targetGM; // Normalise
    double normalTargetGM = 1.0;

    solarTerm = normalHostGM * (sunVector / ( std::pow(sunVector.norm(), 3) ) + ( posDifference / ( std::pow(posDifference.norm(), 3) ) ));
    velocityVector = -normalTargetGM * positionVector / ( std::pow(positionVector.norm(), 3) );
    /* Add on any additional perturbing bodies */
    for (size_t idx = 0; idx < PARAMS::additionalPlanets.size(); ++idx)
    {
        double tmpNormalGM = PARAMS::additionalPlanetsGM[idx] / PARAMS::targetGM;
        Eigen::Matrix<double, 3, 1> planetVector, planetParticleDifference;
        // Get vector to the planet
        getPlanetVectorEphemeris(PARAMS::additionalPlanets[idx], planetVector, currentTime);
        planetParticleDifference = positionVector - planetVector;
        solarTerm += tmpNormalGM * ( (planetVector / std::pow(planetVector.norm(), 3) ) + ( planetParticleDifference / std::pow(planetParticleDifference.norm(), 3)) );
    }
    
    velocityVector -= solarTerm;
    for (size_t i = 0; i < 3; ++i) dx[i+3] = velocityVector(i);
}


/* @brief Test version of the force function using a fully dimensionalised state
*/
void dimensionalForceFunction(const std::vector<double> &x, std::vector<double> &dx, const double t)
{
    Eigen::Matrix<double, 3, 1> sunVector, positionVector, velocityVector, solarTerm;

    // First, the trivial derivatives
    for (size_t idx = 0; idx < 3; ++idx)        // Unrolls for -O2+
    {
        positionVector(idx) = x[idx];
        dx[idx] = x[idx+3];
    }

    // Get current integration time
    double currentTime = PARAMS::EPOCH + t; // EPOCH of the study (seconds) + the integration time so far (seconds)

    // Get Mercury -> Sun vector
    getDimensionalPlanetVectorEphemeris(PARAMS::HOST, sunVector, currentTime);  // Get ephemeris position from target planet (Mercury) to host planet (Sun)
    // Eigen::Matrix<double, 3, 1> posDifference = positionVector - sunVector;     // Position vector between spacecraft and Sun
    Eigen::Matrix<double, 3, 1> posDifference = sunVector - positionVector;

    // Compute solar contributions -> mu * [ (rSun / rSun^3) + ( rSatellite - rSun ) / (rSatellite - rSun) ^ 3 ]
    solarTerm = PARAMS::hostGM * (sunVector / ( std::pow(sunVector.norm(), 3) ) + ( posDifference / ( std::pow(posDifference.norm(), 3) ) ));
    
    // Add on additional perturbing bodies
    for (size_t idx = 0; idx < PARAMS::additionalPlanets.size(); ++idx)     // PARAMS::additionalPlanets contains std::strings of other planets to consider
    {
        double tmpNormalGM = PARAMS::additionalPlanetsGM[idx];              // GM of the other target planet (computed via ephemeris). Normally normalised but not here
        Eigen::Matrix<double, 3, 1> planetVector, planetParticleDifference; // Holding vectors for Mercury <=> Planet vector and (Mercury - Planet) vector
        // Get vector to the planet
        getDimensionalPlanetVectorEphemeris(PARAMS::additionalPlanets[idx], planetVector, currentTime);     // Position vector from host (Mercury) -> Planet
        planetParticleDifference = positionVector - planetVector;                                           // Position vector rSatellite - rPlanet
        // Add this contribution to the overall force on the body
        solarTerm += tmpNormalGM * ( (planetVector / std::pow(planetVector.norm(), 3) ) + ( planetParticleDifference / std::pow(planetParticleDifference.norm(), 3)) );
    }

    /* Get the kinetic energy of the system */
    // Get masses of target and host
    const double G = 6.6743015e-011;
    double massTarget = PARAMS::targetGM * 1000 / G;
    double massHost = PARAMS::hostGM * 1000 / G;

    double T = 0.0;
    double V = 0.0;
    double energy = 0.0;

    // Add sun kinetic energy
    Eigen::Matrix<double, 3, 1> sunVelocity;
    getDimensionalPlanetVelocity(PARAMS::TARGET, sunVelocity, currentTime);

    energy += 0.5 * sunVelocity.norm() * massHost - G / 2.0 * massTarget * massHost / ( sunVector.norm() );
    std::cout << energy << std::endl;
    
    // Compute overall derivative (ignore naming!)
    velocityVector = -PARAMS::targetGM * positionVector / ( std::pow(positionVector.norm(), 3) ) + PARAMS::hostGM * ( (posDifference) / (std::pow(posDifference.norm(), 3)) - sunVector / std::pow(sunVector.norm(), 3)  );

    for (size_t idx = 0; idx < 3; ++idx) dx[idx+3] = velocityVector(idx);     // Reassign
}

void nbody(const std::vector<double> &x, std::vector<double> &dx, const double t)
{
    std::vector<double> spacecraftPos(3), sunPos(3), mercuryPos(3), spacecraftSunVector(3);
    double spacecraftPosNorm = 0;
    double sunPosNorm = 0;
    double spacecraftSunVectorNorm = 0;
    for (size_t idx = 0; idx < x.size(); ++idx)
    {
        spacecraftPos[idx] = x[idx];
        spacecraftPosNorm += x[idx] * x[idx];
        sunPos[idx] = x[idx+3];
        sunPosNorm += sunPos[idx] * sunPos[idx];
        spacecraftSunVector[idx] = spacecraftPos[idx] - sunPos[idx];
        spacecraftSunVectorNorm += spacecraftSunVector[idx] * spacecraftSunVector[idx];
    }

    spacecraftPosNorm = std::pow( std::sqrt(spacecraftPosNorm), 3);
    sunPosNorm = std::pow( std::sqrt(sunPosNorm), 3);
    spacecraftSunVectorNorm = std::pow( std::sqrt(sunPosNorm), 3);

    // EoM for the s/c about Mercury
    for (size_t idx = 0; idx < 3; ++idx)
    {
        dx[idx] =  -PARAMS::targetGM * (spacecraftPos[idx] / spacecraftPosNorm) - PARAMS::hostGM * (spacecraftSunVector[idx] / spacecraftSunVectorNorm + sunPos[idx] / sunPosNorm);
        dx[idx+3] = PARAMS::targetGM;
    }


}

/* @brief Custom integration observer function; returns various status integers depending on whether stopping conditions have been triggered.
   @param[in] x Eigen::Matrix<Type, 6> containing the state vector of the particle at the given point.
   @param[in] x0 Eigen::Matrix<Type, 6> containing the initial position of the particle on the given trajectory
   @param[in] t Const double containing the current time-step of the integration (assumed seconds)
   @param[inout] prevCondOne Holds the previous value of condition one to check for sign reversals
   @returns An integer corresponding to one of four separate events.

    Integer Return Codes
    ~~~~~~~~~~~~~~~~~~~~
    0: Do nothing!  
    1: Trajectory has crashed (R < Rmin)
    2: Trajectory has escaped (R > Rs)
    3: Trajectory is weakly stable (wacky geometrics)
    4: Trajectory is acrobatic (nothing happens!)
*/
// template <typename Type>
int integrationController(std::vector<double> &x, std::vector<double>& xkm1, std::vector<double>& x0, const double t, double & prevCondOne)
{
    Eigen::Matrix<double, 3, 1> r, r0, v0, v, rkm1, vkm1;
    for (size_t idx = 0; idx < 3; ++idx)
    {
        r(idx) = x[idx];
        r0(idx) = x0[idx];
        v(idx) = x[idx+3];
        v0(idx) = x0[idx+3];
        rkm1(idx) = xkm1[idx];
        vkm1(idx) = xkm1[idx+3];
    }

    // First check for a crash - note regularised dynamics
    double rMag = r.norm();
    if (rMag <= 1)
    {
        return 1;    // Premature exit to prevent computing unnecessary branches
    }

    // Next, check for an escape
    double vMag = v.norm();
    double keplerEnergy = (vMag * vMag / 2. - 1./rMag);

    if (rMag > PARAMS::RS/PARAMS::R && keplerEnergy > 0) 
    {
        return 2;
    }

    // Now, check for weakly stable
    Eigen::Matrix<double, 3, 1> angMomentum0 = r0.cross(v0);
    double conditionOne = r.dot(angMomentum0.cross(r0));
    bool signChange = std::signbit(conditionOne * prevCondOne); // Returns True if negative

    double conditionTwo = r.dot(r0);
    double conditionThree = v.dot(v0) * vkm1.dot(v0); // Supposed to be v.dot(v0) * v(k-1).dot(v0), but only doing one at a time here

    if ( signChange && conditionTwo > 0 && conditionThree > 0) return 3;

    // Lastly, check time
    if (fabs(t) >= PARAMS::maxT)  // maxT is defined in interface.hpp as 4 orbits about the Hill region for the given body
    {
        return 4;
    }

    // If nothing else
    prevCondOne = conditionOne;
    return 0;
}

#endif
