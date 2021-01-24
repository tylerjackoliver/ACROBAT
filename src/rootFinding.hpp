#ifndef __ROOT_FINDING_H__
#define __ROOT_FINDING_H__

#include "integration.hpp"
#include <boost/numeric/odeint/stepper/generation/make_controlled.hpp>
#include <boost/numeric/odeint/stepper/generation/make_dense_output.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta_dopri5.hpp>
#include <stdexcept>
#include <vector>
#include <boost/numeric/odeint.hpp>

/* @brief Computes the cross product of two size-3 std::vectors.
 * @param[in] a The first vector in the cross product
 * @param[in] b The second vector in the cross product
 * @param[out] c The last vector in the cross product. Space is reserved if it does not already exist.
 */
void cross3(std::vector<double> &a, std::vector<double>& b, std::vector<double>& c)
{
    if (a.size() != 3 or b.size() != 3)
    {
        throw std::domain_error("One of the input vectors to cross3 are not of size 3.");
    }

    // Reserve space in c
    c.reserve(3);
    // Compute
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2]; 
    c[2] = a[0] * b[1] - a[1] * b[0];
}

/* @brief Computes the dot product of two std::vectors. The vectors must be the same size.
 * @param[in] a The first vector in the dot product
 * @param[in] b The second vector in the dot product
 * @returns dot The result of the dot product.
 */
template <typename Type>
Type dotProduct(std::vector<Type>& a, std::vector<Type>& b)
{
    Type dot = 0;                                                                               // Output
    if (a.size() != b.size())                                                                   // Check valid inputs
    {
        throw std::domain_error("The input vectors to dotProduct are not the same size.");
    }
    for (size_t idx = 0; idx < a.size(); ++idx)
    {
        dot += a[idx] * b[idx];                                                                 // Intel compiler will auto-unroll this for -O2
    }
    return dot;
}

/* @brief Obtain the value of condition one - the intersection of a trajectory with a plane - for a given state state.
 * @param[in] state State at which to compute the value of condition one for
 * @param[in] x0 The initial condition for the given trajectory
 * @returns The value of condition one at the given point.
 */
double getConditionOne(std::vector<double>& state, std::vector<double>& x0)
{
    std::vector<double> angularMomentum(3), r0(3), v0(3), r(3), v(3), angMomentumCrossR0(3);
    for (size_t idx = 0; idx < state.size() / 2; idx++)                                         // Assign values from state vector to position, position derivative
    {
        r[idx] = state[idx];
        r0[idx] = x0[idx];
        v[idx] = state[idx+3]; 
        v0[idx] = x0[idx+3];
    }
    cross3(r0, v0, angularMomentum);                                                            // Get the angular momentum vector
    cross3(angularMomentum, r0, angMomentumCrossR0);                                            // Cross angular momentum vector with the initial position 
    double dot = dotProduct(r, angMomentumCrossR0);                                             // Dot product of the above yields condition one
    return dot;
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
void obtainZero(std::vector<double>& x0, std::vector<double>& xGuess, double& t0Guess, double tol=1e-06)
{
    
    typedef boost::numeric::odeint::result_of::make_dense_output                                // (Unattractively) instantiate stepper types in the integration
                                               <boost::numeric::odeint::runge_kutta_dopri5
                                               <std::vector<double>>>::type dense_stepper_type;
    dense_stepper_type stepper = make_dense_output(1.e-012, 1.e-012, 
                                                   boost::numeric::odeint::runge_kutta_dopri5
                                                   <std::vector<double>>());                    // Makes a dense stepper from the RK5 solver
    typedef dense_stepper_type::time_type time_type;                                            // Define later units in the same datatype the integrator is using 
    stepper.initialize(xGuess, t0Guess, -.001);                                                 // Initialise solver with the initial values of the ODE

    double conditionOnePrevStep = 1.0;                                                          // Condition one on the previous step 
    double conditionOne = 1.0;                                                                  // Condition one on this step; check for sign change
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

    // Now we know we're bracketing the root. Start an interval bisection
    time_type lower, upper, mid;
    lower = std::min(currentTime.first, currentTime.second);
    upper = std::max(currentTime.first, currentTime.second);
    const int maxIterations = 100;                                                              // Guards against stale while loop
    int numberOfIterations = 0;

    while (lower < upper && numberOfIterations < maxIterations)
    {
        mid = lower + (upper - lower) / 2.;                                                     // Interval bisection - compute midpoint
        std::vector<double> state(6);
        stepper.calc_state(mid, state);                                                         // Extract state at the midpoint
        conditionOne = getConditionOne(state, x0);
        
        if ( fabs(conditionOne) < tol)                                                          // We've found the root
        {
            xGuess = state;                                                                     // Set state to the output
            t0Guess = mid;                                                                      // Set time to the output
            break;
        }
        else if (conditionOne < 0)                                                              // Lower is too low
        {
            lower = mid;
        } else                                                                                  // Upper is too high
        {
            upper = mid;
        }
        numberOfIterations++;                                                                // Increment counter
    }
    if (numberOfIterations == maxIterations)
    {
        throw std::runtime_error("Maximum number of iterations exceeded in finding the root.");
    }
}
#endif
