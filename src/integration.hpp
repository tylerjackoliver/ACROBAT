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
    // Statuses
    std::vector<unsigned long> setStatistics(4);

    // Get the integration direction (+ve or -ve time?)
    int directionTime = sgn(field.getFinalTime() - field.getInitialTime());

    // Iterate through the conditions on the EME2000 field
    #pragma omp parallel
    {
    #pragma omp for
    for (unsigned i = 0; i < field.getXExtent(); ++i)
    {
        for (unsigned j = 0; j < field.getYExtent(); ++j)
        {
            for (unsigned k = 0; k < field.getZExtent(); ++k)
            {
                // Get current point
                Point<double> currentPoint = field.getValue(indices[0], indices[1], indices[2]);

                // Get the status of this points behaviour
                int status = getStatus(currentPoint, field.getInitialTime(), directionTime);

                // Increment the correct counter in the setStatistics vector (ordered s.t. indexes is status-1)
                #pragma omp atomic
                setStatistics[status-1]++;

                // If weakly stable, append its indices to the points vector
                if (status == 2)
                {
                    #pragma omp critical
                    {
                        points.inset(points.end(), currentPoint);
                    }
                }
            }
        }
    }
    } // parallel
    printStatistics(stabNum, setStatistics);
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
    // Statuses - ordered s.t. index related to condition is (status code - 1) - crash, escape, stable, acrobatic
    std::vector<unsigned long> setStatistics(4);

    // Get the integration direction (+ve or -ve time?)
    int directionTime = sgn(field.getFinalTime() - field.getInitialTime());

    // Iterate through the conditions on the EME2000 field
    #pragma omp parallel
    {
    #pragma omp for
    for (unsigned i = 0; i < domain.size(); ++i)
    {
        // Get an initial position
        Point<vectorType> indices = domain[i];
        Point<double> currentPoint = field.getValue(indices[0], indices[1], indices[2]);

        // Get the status of this points behaviour
        int status = getStatus(currentPoint, field.getInitialTime(), directionTime);
        
        // Increment the correct counter in the setStatistics vector (ordered s.t. index is status-1)
        #pragma omp atomic
        setStatistics[status-1]++;

        // If weakly stable, append its indices to the points vector
        if (status == 2)
        {
            #pragma omp critical
            {
                points.inset(points.end(), currentPoint);
            }
        }
    }
    } // end of parallel section
    printStatistics(stabNum, setStatistics);
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

/* @brief Prints the statistics for a given ballistic set computation
   @param[in] stabNum The number of the sability index for the current computation
   @param[in] setStatistics std::vector<> containing the statistics defined as in the set determination routines.
*/
template <typename integerType, typename vectorType>
void printStatistics(integerType stabNum, std::vector<vectorType>& setStatistics)
{
    std::cout << "For a stability index of " stabNum << "the statistics are as follows:" << std::endl;
    std::cout << "\t Number of crashes:   " << setStatistics[0] << std::endl;
    std::cout << "\t Number of escapes:   " << setStatistics[1] << std::endl;
    std::cout << "\t Number of stable:    " << setStatistics[2] << std::endl;
    std::cout << "\t Number of acrobatic: " << setStatistics[3] << std::endl;
}

/* @brief Computes whether a given point is acrobatic, weakly stable, crashes, or escapes.
*  @param[in] point The initial conditions (in a Point<> structure) to be tested
*  @param[in] initTime The initial time for the integration of the trajectory
*  @param[in] direction The direction for the timespan of the integration (+1/-1)
*  @returns Status code corresponding to the behaviour of the trajectory
*/
template <typename pointType, typename doubleType, typename integerType>
int getStatus(Point<PointType>& point, doubleType& initTime, integerType& direction)
{
    // Initialise the stepper to be used
    typedef Eigen::Vector<double, 6> stateType;
    boost::numeric::odeint::runge_kutta_fehlberg78<stateType> method;
    auto stepper = boost::numeric::odeint::make_controlled(/*reltol*/ 1e-011, /*absTol*/ 1e-011);

    // Fill an initial condition vector
    stateType x0, x;
    for (unsigned idx = 0; idx < 6; ++i) 
    {
        x0[idx] = point[idx];
        x[idx] = point[idx];
    }

    // Initialise current time
    double currentTime = initTime;
    double dt = 0.01 * currentTime * direction;

    // Integration status
    int status = 0;

    while (status == 0) // While none of the stopping conditions have been verified
    {
        /* Make a step using the given solver & force function */
        make_step(stepper, forceFunction, x, currentTime, dt);

        /* Call the integrator observer function */
        status = integrationController(x, x0, currentTime);
    }
    return status; // Return status when it doesn't correspond to 'keep going'
}

#endif
