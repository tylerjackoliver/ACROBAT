#ifndef __ROOT_FINDING_H__
#define __ROOT_FINDING_H__

#include "integration.hpp"
#include <boost/numeric/odeint/stepper/generation/make_controlled.hpp>
#include <boost/numeric/odeint/stepper/generation/make_dense_output.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta_dopri5.hpp>
#include <stdexcept>
#include <vector>
#include <boost/numeric/odeint.hpp>
#include <boost/math/tools/roots.hpp>

// /* @brief Computes the cross product of two size-3 std::vectors.
//  * @param[in] a The first vector in the cross product
//  * @param[in] b The second vector in the cross product
//  * @param[out] c The last vector in the cross product. Space is reserved if it does not already exist.
//  */
// void cross3(std::vector<double> &a, std::vector<double>& b, std::vector<double>& c)
// {
//     if (a.size() != 3 or b.size() != 3)
//     {
//         throw std::domain_error("One of the input vectors to cross3 are not of size 3.");
//     }
//     // Reserve space in c
//     c.reserve(3);
//     // Compute
//     c[0] = a[1] * b[2] - a[2] * b[1];
//     c[1] = a[2] * b[0] - a[0] * b[2];
//     c[2] = a[0] * b[1] - a[1] * b[0];
// }

// /* @brief Computes the dot product of two std::vectors. The vectors must be the same size.
//  * @param[in] a The first vector in the dot product
//  * @param[in] b The second vector in the dot product
//  * @returns dot The result of the dot product.
//  */
// template <typename Type>
// Type dotProduct(std::vector<Type>& a, std::vector<Type>& b)
// {
//     Type dot = 0;                                                                               // Output
//     if (a.size() != b.size())                                                                   // Check valid inputs
//     {
//         throw std::domain_error("The input vectors to dotProduct are not the same size.");
//     }
//     for (size_t idx = 0; idx < a.size(); ++idx)
//     {
//         dot += a[idx] * b[idx];                                                                 // Intel compiler will auto-unroll this for -O2
//     }
//     return dot;
// }

/* @brief Obtain the value of condition one - the intersection of a trajectory with a plane - for a given state state.
 * @param[in] state State at which to compute the value of condition one for
 * @param[in] x0 The initial condition for the given trajectory
 * @returns The value of condition one at the given point.
 */
double getConditionOne(std::vector<double>& state, std::vector<double>& x0)
{
    Eigen::Matrix<double, 3, 1> r, r0, v0, v, angMomentum0;
    for (size_t idx = 0; idx < 3; ++idx)
    {
        r(idx) = state[idx];
        r0(idx) = x0[idx];
        v(idx) = state[idx+3];
        v0(idx) = x0[idx+3];
    }
    angMomentum0 = r0.cross(v0);
    return r.dot(angMomentum0.cross(r0));
}

// double getConditionOneIntrinsic(std::vector<double>& state, std::vector<double>& x0)
// {
//     std::vector<double> r0(3), v0(3), r(3), v(3), angularMomentum0(3), angMomentumCrossR0;
//     for (int idx = 0; idx < 3; ++idx)
//     {
//         r0[idx] = x0[idx];
//         v0[idx] = x0[idx+3];
//         r[idx] = state[idx];
//         v[idx] = state[idx + 3];
//     }
//     cross3(r0, v0, angularMomentum0);
//     cross3(angularMomentum0, r0, angMomentumCrossR0);
//     return dotProduct(r, angMomentumCrossR0);
// }

/* @brief Does what it says on the tin
*/
bool getConditionTwoAndThree(std::vector<double>& currentPosition, std::vector<double>& x0ThisOrbit, std::vector<double>& x0Overall)
{
    Eigen::Matrix<double, 3, 1> r, r0, v0, v, rkm1, vkm1;
    for (size_t idx = 0; idx < 3; ++idx)
    {
        r(idx) = currentPosition[idx];
        r0(idx) = x0Overall[idx];
        v(idx) = currentPosition[idx+3];
        v0(idx) = x0Overall[idx+3];
        rkm1(idx) = x0ThisOrbit[idx];
        vkm1(idx) = x0ThisOrbit[idx+3];
    }
    bool conditionTwo = r.dot(r0) > 0;
    bool conditionThree = v.dot(v0) * vkm1.dot(v0) > 0; // Supposed to be v.dot(v0) * v(k-1).dot(v0), but only doing one at a time here
    return conditionTwo && conditionThree;
}


/* @brief Objective function for the interval bisection.
*/
template <typename stepper>
double costFunction(const double& initTime, stepper& step, std::vector<double>& x0)
{
    std::vector<double> thisState(6);
    step.calc_state(initTime, thisState);
    return getConditionOne(thisState, x0);
}

template <typename T, typename stepper>
bool getTolerance(T& lowerTime, T& upperTime, stepper& step, std::vector<double>& x0)
{
    double c1 = costFunction(lowerTime, step, x0);
    double c2 = costFunction(upperTime, step, x0);
    if ( fabs(c2 - c1) <= 1e-010 ) return true;
    return false;
}

/* @brief Obtains the zero of the intersection of a trajectory with a given reference plane.
 * @param[in] x0 The initial conditon of the trajectory at t = 0.
 * @param[inout] xGuess On function entry, this is a guess for the state immediately _after_ the root. On exit, this is the state _at_ the root.
 * @param[inout] t0Guess On function entry, this is a guess for the time immediately _after_ the root. On exit, this is the time _at_ the root.
 * @param[in] tol Optional. The tolerance of the interval bisection method. Defaults to 1e-06.
 *
 * This function will receive as input the step immediately after the root. The strategy to regain the root is thus:
 *      - Integrate the trajectory backwards using make_dense_output until the value of the current step and the value of the future step
 *        exist on either side of the root. At this point, the entire root is bracketed and no further integration need take place. In theory,
 *        this should occur within only a few steps (and, therefore, function evaluations)
 *      - An interval bisection method is then used in this root bracket to identify, within a tolerance tol, where the root lies.
 *      - The state at this time is then computed and assigned to the input xGuess. The time at the root is assigned to the input t0Guess.
 */
void obtainZero(std::vector<double>& x0, std::vector<double>& xGuess, double& t0Guess, double tol=1e-12)
{
    
    typedef boost::numeric::odeint::result_of::make_dense_output                                // (Unattractively) instantiate stepper types in the integration
                                               <boost::numeric::odeint::runge_kutta_dopri5
                                               <std::vector<double>>>::type dense_stepper_type;
    dense_stepper_type stepper = make_dense_output(1.e-013, 1.e-013, 
                                                   boost::numeric::odeint::runge_kutta_dopri5
                                                   <std::vector<double>>());                    // Makes a dense stepper from the RK5 solver
    typedef dense_stepper_type::time_type time_type;                                            // Define later units in the same datatype the integrator is using 
    stepper.initialize(xGuess, t0Guess, -.001);                                                 // Initialise solver with the initial values of the ODE

    double conditionOnePrevStep = getConditionOne(xGuess, x0);
    // std::cout << "Condition one on entry " << conditionOnePrevStep << std::endl;
    double conditionOne = conditionOnePrevStep;                                                  // Condition one on this step; check for sign change
    bool signChange = false;                                                                    // Have we reached the root?
    std::pair<time_type, time_type> currentTime;                                                // Contains (t, t+dt) on exit from do_step

    while (!signChange)                                                                         // While we aren't bracketing the root
    {
        std::vector<double> currentState(6);
        conditionOnePrevStep = conditionOne;                                                    // Save the previous value of the objective function 
        currentTime = stepper.do_step(forceFunction);                                           // Make a step
        stepper.calc_state(currentTime.second, currentState);                                   // Calculate the state at the later time (t + dt)
        conditionOne = getConditionOne(currentState, x0);                                       // Get the value of the objective function at this point
        signChange = conditionOne * conditionOnePrevStep < 0;                                   // Check if we're bracketing the root here
    }
    /*
     Now we know we're bracketing the root. Get the bracket and use the TOMS748 algorithm from
     Boost to find the root
    */
    
    time_type lower, upper, mid;
    lower = std::min(currentTime.first, currentTime.second);
    upper = std::max(currentTime.first, currentTime.second);

    auto func = [&stepper, &x0](const time_type& t){return costFunction(t, stepper, x0);};
    auto tolFunct = [&stepper, &x0](const double& lower, const double& upper){return getTolerance(lower, upper, stepper, x0);};
    boost::uintmax_t max_iter = 150, initMaxIter = 150;
    std::pair<double, double> rootBracket = boost::math::tools::toms748_solve(func, lower, upper, tolFunct, max_iter);
    
    if (max_iter < initMaxIter)
    {
        // std::cout << "Bracketing root " << rootBracket.first << " " << rootBracket.second << std::endl;
        t0Guess = rootBracket.first;
        stepper.calc_state(t0Guess, xGuess);
        // std::cout << "Lower condition " << getConditionOne(xGuess, x0) << std::endl;
        t0Guess = rootBracket.first + (rootBracket.second - rootBracket.first) / 2.;
        stepper.calc_state(t0Guess, xGuess);
        std::cout << "Upper condition " << getConditionOne(xGuess, x0) << std::endl;
    } else
    {
        throw std::runtime_error("Maximum number of iterations exceeded in finding the root.");
    }
}
#endif
