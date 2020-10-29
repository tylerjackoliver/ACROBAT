#ifndef __EQUATIONS_OF_MOTION_H__
#define __EQUATIONS_OF_MOTION_H__

#include <eigen3/Core>
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
 * @param[in] r: Vector to the Sun; Eigen::Vector<Type, Size>
 * @param[in] t: Time at which this vector is desired; const Type
 */
template <typename Type>
void getSunVector(Eigen::Vector<Type, 3> &r, const Type t)
{
    Eigen::Vector3<Type, 3> rotationVector;
    double trueAnomaly;

    // Get orbital elements of host about target in the EME2000 frame (actually the J2000 frame in SPICE!)
    ACROBAT::OEs sunOEs = getOEs(t, PARAMS::host, PARAMS::target);      // Need orientation of the Sun wrt planet we're seeking ballistic capture orbits for
    ACROBAT::OEs planetOEs = getOEs(t, PARAMS::target, PARAMS::host);   // Need a, e of planet we're constructing ballistic orbits for

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
ACROBAT::OEs getOEs(double t, std::string body, std::string obs)
{
    ACROBAT::OEs ret;
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

/* @brief Get the rotation matrix to move from the EME2000 to point at the Sun
 * @param[inout] rotationVector: Eigen::Vector<Type, 3> containing the rotation vector
 * @param[in] sunOEs: Reference to an OEs struct containing the orbital elements of the Sun
 * @param[in] trueAnomaly: Current true anomaly of planet around the Sun.
 */
template <typename Type>
void getRotationMatrix(Eigen::Vector<Type, 3> &rotationVector, ACROBAT::OEs &sunOEs, double &trueAnomaly)
{
    double theta = sunOEs.omega + trueAnomaly;
    
    rotationVector(0) = std::sin(theta) * std::sin(sunOEs.longtd) * std::cos(sunOEs.inc) - std::cos(theta) * std::cos(sunOEs.longtd);
    rotationVector(1) = - std::cos(theta) * std::sin(sunOEs.longtd) - std::sin(theta) * std::cos(sunOEs.longtd) * std::cos(sunOEs.inc);
    rotationVector(2) = -std::sin(theta) * std::sin(sunOEs.inc);
}

/* @brief Compute the derivative vector of the current state for numerical integration routines.
*  @param[in] x Eigen::Vector<Type, 6> of the initial [r, v] vector for the particle
*  @param[out] dx Eigen::Vector<Type, 6> corresponding to the derivative of x
*  @param[in] t Current integration time-step
*/
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
   @param[in] x Eigen::Vector<Type, 6> containing the state vector of the particle at the given point.
   @param[in] x0 Eigen::Vector<Type, 6> containing the initial position of the particle on the given trajectory
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

/* @brief Obtains the points in the set for a given stability index stabNum across the whole EMEJ2000 field.
   @param[in] stabNum The stability index
   @param[in] field EME2000 Field containing the domain to integrate
   @param[out] points std::vector<> of Points containing the indices of points in the set
*/
template <typename integerType, typename fieldType, typename vectorType>
void getSet(integerType stabNum, ACROBAT::emeField<fieldType>& field, std::vector<Point<vectorType>> &points)
{
    // Initialise integrator
    typedef Eigen::Vector<double, 6> stateType;
    boost::numeric::odeint::runge_kutta_fehlberg78<stateType> method;
    auto stepper = boost::numeric::odeint::make_controlled(/*reltol*/ 1e-011, /*absTol*/ 1e-011);
    
    // Statuses
    unsigned long numAcrobatic = 0;
    unsigned long numCrash = 0;
    unsigned long numEscape = 0;
    unsigned long numStable = 0;

    // Get the integration direction (+ve or -ve time?)
    int directionTime = sgn(field.getFinalTime() - field.getInitialTime());

    // Iterate through the conditions on the EME2000 field
    #pragma omp parallel
    {
    #pragma omp for
    for (unsigned i = 0; i < field.getXExtent(); ++i)
    {
        // Private copies of statuses
        std::vector<Point<int>> privatePoints;
        unsigned long localNumAcrobatic = 0;
        unsigned long localNumCrash = 0;
        unsigned long localNumEscape = 0;
        unsigned long localNumStable = 0;

        for (unsigned j = 0; j < field.getYExtent(); ++j)
        {
            for (unsigned k = 0; k < field.getZExtent(); ++k)
            {
                // Get an initial position
                Point<double> tmp = field.getValue(i, j, k);

                // Fill an initial condition vector
                stateType x0, x;
                for (unsigned idx = 0; idx < 6; ++i) x0[idx] = tmp[idx];
                x = x0;

                // Initialise current time
                double currentTime = field.getInitialTime();
                double finalTime = field.getFinalTime();
                double dt = 0.01 * direction; // Multiply by the direction of integration (sign)

                // Step success?
                boost::numeric::odeint::controlled_step_result result;

                while ( fabs(currentTime) < fabs(finalTime) )
                {
                    // Try and make a step
                    result = stepper.try_step(forceFunction, x, currentTime, dt);

                    // If within tolerance requirements
                    if (result == boost::numeric::odeint::success)
                    {
                        /* Call the integrator observer function */
                        int status = integrationController(x, x0, currentTime);

                        /* Check the return code - do most likely first */ // Switch-case?
                        if (status == 0) continue;
                        if (status == 1) // Crash
                        {
                            localNumCrash++;
                            currentTime = 1.e6; // Force while-loop break
                        }
                        if (status == 2) // Escaped
                        {
                            localNumEscape++;
                            currentTime = 1.e6;
                        }
                        if (status == 3) // Weakly stable - what we want!
                        {
                            localNumStable++;
                            // Add indices to privatePoints
                            Point<int> tmp;
                            tmp[0] = i; tmp[1] = j; tmp[2] = k;
                            privatePoints.push_back(tmp);
                            currentTime = 1.e6;
                        }
                        if (status == 4) // Acrobatic
                        {
                            localNumAcrobatic++;
                            currentTime = 1.e6;
                        }
                    }
                }
            }
        }
    }
    // Globalise all the local values and append everything to points for a clean exit
    #pragma omp critical
    {
        numCrash += localNumCrash;
        numEscape += localNumEscape;
        numStable += localNumStable;
        numAcrobatic += localNumAcrobatic;

        points.insert(points.end(), privatePoints.begin(), privatePoints.end());
    }
    } // end of parallel section
    std::cout << "For a stability index of stabNum, the statistics are as follows:" << std::endl;
    std::cout << "\t Number of crashes:   " << numCrash << std::endl;
    std::cout << "\t Number of escapes:   " << numEscapes << std::endl;
    std::cout << "\t Number of stable:    " << numStable << std::endl;
    std::cout << "\t Number of acrobatic: " << numAcrobatic << std::endl;
} // function

/* @brief Obtains the capture set from a given set of indices corresponding to points in an EME2000 field.
   @param[in] stabNum The stability index the set is sought for
   @param[in] domain The vector of points containing indices for which the sets are to be based off of
   @param[in] field The domain for which the points in domain correspond to
   @param[out] points A vector containing the indices of the corresponding capture set
*/ 
template <typename integerType, typename fieldType, typename vectorType>
void getSetFromPoints(const integerType& stabNum, const std::vector<Point<vectorType>> &domain, const ACROBAT::emeField<fieldType>& field, std::vector<Point<vectorType>> &points)
{
    // Initialise integrator
    typedef Eigen::Vector<double, 6> stateType;
    boost::numeric::odeint::runge_kutta_fehlberg78<stateType> method;
    auto stepper = boost::numeric::odeint::make_controlled(/*reltol*/ 1e-011, /*absTol*/ 1e-011);
    
    // Statuses
    unsigned long numAcrobatic = 0;
    unsigned long numCrash = 0;
    unsigned long numEscape = 0;
    unsigned long numStable = 0;

    // Get the integration direction (+ve or -ve time?)
    int directionTime = sgn(field.getFinalTime() - field.getInitialTime());

    // Iterate through the conditions on the EME2000 field
    #pragma omp parallel
    {
    #pragma omp for
    for (unsigned i = 0; i < domain.size(); ++i)
    {
        // Private copies of statuses
        std::vector<Point<int>> privatePoints;
        std::unordered_map<std::string, unsigned long long> setStatistics;

        // Populate map
        setStatistics["a"] = 0;  // acrobatic
        setStatistics["c"] = 0;  // crash
        setStatistics["e"] = 0;  // escape
        setStatistics["s"] = 0;  // stable for one more revolution

        unsigned long localNumAcrobatic = 0;
        unsigned long localNumCrash = 0;
        unsigned long localNumEscape = 0;
        unsigned long localNumStable = 0;

        // Get an initial position
        Point<vectorType> indices = domain[i];
        Point<double> tmp = field.getValue(indices[0], indices[1], indices[2]);

        // Fill an initial condition vector
        stateType x0, x;
        for (unsigned idx = 0; idx < 6; ++i) x0[idx] = tmp[idx];
        x = x0;

        // Initialise current time
        double currentTime = field.getInitialTime();
        double finalTime = field.getFinalTime();
        double dt = 0.01 * direction; // Multiply by the direction of integration (sign)

        // Step success?
        boost::numeric::odeint::controlled_step_result result;

        while ( fabs(currentTime) < fabs(finalTime) )
        {
            /* Make a step using the given solver & force function */
            make_step(stepper, forceFunction, x, currentTime, dt);

            /* Call the integrator observer function */
            int status = integrationController(x, x0, currentTime);

            /* Check the return code - do most likely first */ // Switch-case?
            if (status == 0) continue;
            if (status == 1) // Crash
            {
                localNumCrash++;
                currentTime = 1.e6; // Force while-loop break
            }
            if (status == 2) // Escaped
            {
                localNumEscape++;
                currentTime = 1.e6;
            }
            if (status == 3) // Weakly stable - what we want!
            {
                localNumStable++;
                // Add indices to privatePoints
                Point<int> tmp;
                tmp[0] = i; tmp[1] = j; tmp[2] = k;
                privatePoints.push_back(tmp);
                currentTime = 1.e6;
            }
            if (status == 4) // Acrobatic
            {
                localNumAcrobatic++;
                currentTime = 1.e6;
            }
        }
    }
    // Globalise all the local values and append everything to points for a clean exit
    #pragma omp critical
    {
        numCrash += localNumCrash;
        numEscape += localNumEscape;
        numStable += localNumStable;
        numAcrobatic += localNumAcrobatic;

        points.insert(points.end(), privatePoints.begin(), privatePoints.end());
    }
    } // end of parallel section
    std::cout << "For a stability index of stabNum, the statistics are as follows:" << std::endl;
    std::cout << "\t Number of crashes:   " << numCrash << std::endl;
    std::cout << "\t Number of escapes:   " << numEscapes << std::endl;
    std::cout << "\t Number of stable:    " << numStable << std::endl;
    std::cout << "\t Number of acrobatic: " << numAcrobatic << std::endl;
} // function

/* @brief Performs a step using a given boost stepper; returns the current time and the new state if successful
 * @param[in] stepper Boost::odeint::numeric object corresponding to a controlled stepper
 * @param[inout] currentTime On input, it contains the time of integration at the start of step. On exit, it contains the new time (i.e. after one step)
 * @param[inout] x On input, it contains the state of the particle at the previous timestep. On exit, it contains the new state of the particle (i.e. after one step)
 * @param[inout] dt On input, it contains a guess for the time-step to be taken. On exit, it contains the actual time-step taken.
 * @param[in] f The force function to integrate.
 */
template <typename stepperType, typename timeType, typename stateType, typename function>
void make_step(stepperType& stepper, stateType &x, timeType& currentTime, timeType& dt, function& f
{
    boost::numeric::odeint::controlled_step_result result = boost::numeric::odeint::fail;
    
    /* We only want this to return back to the calling function when the step was successful.
       Therefore, keep retrying this until the stepper returns a value that was within tolerance.
    */
    while(result == boost::numeric::odeint::fail)
    {
        result = stepper.try_step(f, x, currentTime, dt);
    }
}

#endif
