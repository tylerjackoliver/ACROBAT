#ifndef __TOPS_748_H__
#define __TOPS_748_H__

#include "integration.hpp"
#include <boost/numeric/odeint/stepper/generation/make_controlled.hpp>
#include <boost/numeric/odeint/stepper/generation/make_dense_output.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta_dopri5.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta_fehlberg78.hpp>
#include <stdexcept>
#include <vector>
#include <boost/math/tools/roots.hpp> // Bracket_and_solve_root
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
    c[1] = a[2] * b[0] - a[0] * b[2]; // Sign changed d/t determinant
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
    Type dot;
    if (a.size() != b.size())
    {
        throw std::domain_error("The input vectors to dotProduct are not the same size.");
    }
    dot = 0;
    for (size_t idx = 0; idx < a.size(); ++idx)
    {
        dot += a[idx] * b[idx]; // Intel compiler will auto-unroll this for -O2
    }
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
template <typename stepperType, typename vectorType>
void obtainZero(std::vector<vectorType>& x0, std::vector<vectorType>& xGuess, double& t0Guess, double tol=1e-06)
{
    // (Unattractively) instantiate a stepper for use in the integration
    typedef boost::numeric::odeint::result_of::make_dense_output
                                               <boost::numeric::odeint::runge_kutta_dopri5
                                               <std::vector<double>>>::type dense_stepper_type;
    dense_stepper_type stepper = make_dense_output(1.e-012, 1.e-012, 
                                                   boost::numeric::odeint::runge_kutta_dopri5
                                                   <std::vector<double>>());

    // Initialise the integration - going backwards...
    typedef dense_stepper_type::time_type time_type;
    double conditionOnePrevStep = 1.0; // Condition one on the previous step 
    double conditionOne = 1.0; // Condition one on this step; check for sign change
    bool signChange = false;
    std::pair<time_type, time_type> currentTime = {1, 1};
    stepper.initialize(xGuess, t0Guess, -.001);

    while (!signChange)
    {
        conditionOnePrevStep = conditionOne;
        currentTime = stepper.do_step(forceFunction);
        std::vector<double> currentState(6);
        stepper.calc_state(currentTime.second, currentState);
        conditionOne = getConditionOne(currentState, x0);
        signChange = conditionOne * conditionOnePrevStep < 0;
    }

    // Now we know we're bracketing the root. Start an interval bisection
    time_type lower, upper, mid;
    lower = currentTime.first;
    upper = currentTime.second;

    while (lower <= upper)
    {
        mid = lower + (upper - lower) / 2.;
        std::vector<double> state(6);
        stepper.calc_state(mid, state);
        
        if ( fabs(conditionOne) < tol)  // We've found the root
        {
            xGuess = state;
            t0Guess = mid;
            break;
        }
        else if (conditionOne < 0) // Lower is too low
        {
            lower = mid;
        } else  // Upper is too high
        {
            upper = mid;
        }
    }
}

/* @brief Obtain the value of condition one - the intersection of a trajectory with a plane - for a given state state.
 * @param[in] state State at which to compute the value of condition one for
 * @param[in] x0 The initial condition for the given trajectory
 * @returns The value of condition one at the given point.
 */
double getConditionOne(std::vector<double>& state, std::vector<double>& x0)
{
    std::vector<double> angularMomentum(3), r0(3), v0(3), r(3), v(3), angMomentumCrossR0(3);
    for (size_t idx = 0; idx < state.size() / 2; idx++)
    {
        r[idx] = state[idx];
        r0[idx] = x0[idx];
        v[idx] = state[idx+3]; // Assuming 6-dimensional vector
        v0[idx] = x0[idx+3];
    }
    // Get the initial angular momentum vector
    cross3(r0, v0, angularMomentum);
    cross3(angularMomentum, r0, angMomentumCrossR0);
    double dot = dotProduct(r, angMomentumCrossR0);
    return dot;
}

#endif
