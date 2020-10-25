#ifndef __LAGUERRE_CONWAY_H__
#define __LAGUERRE_CONWAY_H__

/* @brief Implements the Laguerre-Conway method for solving Kepler's equation. Convergence is generally good within 5 iterations.
   @param[in] E0: Initial guess for the eccentric anomaly
   @param[out] Ef: Final iterate; solution of Kepler's equation.
   @param[in] ecc: Eccentricity for Kepler's equation.
   @param[in] M: Mean anomaly for which the solution is desired.
   @param[in]: eps: Tolerance for the solver
   @param[in]: f: Function that returns the error in Kepler's equation
   @param[in]: fp: Function that returns the error in first derivative of Kepler's equation
   @param[in]: fpp: Function that returns the error in second derivative of Kepler's equation
*/
template <typename Type, typename Function>
void laguerreConway(Type &E0, Type &Ef, Type &ecc, Type &M, Type &eps, Function f, Function fp, Function fpp)
{
    const int n = 5;        // Tuning parameter - 5 for now, as per paper
    Type tolerance = 100.   // Stopping tolerance
    Type xi = E0;
    unsigned num_iters = 0;

    while (tolerance >= eps && num_iters < 10) // Iteration should be within 4 iterations for majority of E, ecc
    {   
        // Pre-compute function evaluations and derivatives
        Type fval = f(xi, ecc, M);
        Type deriv = fp(xi, ecc);
        Type dDeriv = fpp(xi, ecc);

        // Compute numerator, square root
        Type numerator = - n * f(E0);
        Type root = std::sqrt( fabs((n-1) * (n-1) * (fp * fp) - n * (n-1) * fval * dDeriv) );

        // Denominator is such that absolute value is maximised
        Type denominator = std::max( fabs(fp + root), fabs(fp - root) );

        // Compute update to iterate
        Type delta_n1 = numerator / denominator;
        xi += delta_n1;

        tolerance = delta_n1; // Change in latest iterate
        num_iters++;
    };
    //
    Ef = xi;
}

/* @brief Computes the error in Kepler's equation.
   @param[in] E The value of eccentric anomaly at the current step
   @param[in] ecc The value of eccentricity at the current step
   @param[in] M The value of Mean anomaly at the current step
   @returns The error in Kepler's equation.
*/
template <typename Type>
Type keplersEquation(const Type E, const Type ecc, const Type M)
{
    return E - ecc * std::sin(E) - M;
}

/* @brief Computes the error in the first derivative of Kepler's equation.
   @param[in] E The value of eccentric anomaly at the current step
   @param[in] ecc The value of eccentricity at the current step
   @returns The error in the first derivative of Kepler's equation.
*/
template <typename Type>
Type dMdE(const Type E, const Type ecc)
{
    return 1 - ecc * std::cos(E);
}

/* @brief Computes the error in the second derivative of Kepler's equation.
   @param[in] E The value of eccentric anomaly at the current time-step
   @param[in] ecc The value of eccentricity at the current step.
   @returns The error in the second derivative of Kepler's equation.
 */
template <typename Type>
Type dMMdEE(const Type E, const Type ecc)
{
    return ecc * std::sin(E);
}

/* @brief Wrapper function to convert from Mean anomaly to Eccentric Anomaly
   @param[in] M: Mean anomaly to solver for
   @param[in] ecc: Eccentricity to solve for
   @param[in] eps: Tolerance to be used in the solver.
   @returns The eccentric anomaly at that point.
*/
template <typename Type>
Type meanToEccentric(Type &M, Type &ecc, Type &eps)
{
    Type E;
    laguerreConway(M, E, ecc, M, eps, keplersEquation, dMdE, dMMdEE);
    return E;
}

/* @brief Wrapper function to convert from eccentric to true anomaly
   @param[in] E: Eccentric anomaly to convert
   @param[in] ecc: Eccentricity of orbit to convert
   @returns The true anomaly corresponding to the eccentric anomaly.
*/
template <typename Type>
Type eccentricToTrue(Type &E, Type &ecc)
{
    double pi = 4.0 * atan(1.0);
    Type root = sqrt::( (1-ecc) / (1+ecc) );
    Type inv = root * std::tan(E / 2.0);
    return ( 2.0 * std::atan(inv) ) % 2. * pi;
}

/* Wrapper: Converts from Mean to True anomaly */
/* @brief Wrapper function to convert from mean anomaly to true anomaly.
   @param[in] M Mean anomaly to convert
   @param[in] ecc Eccentricity of orbit for which to convert
   @returns The true anomaly corresponding to the mean anomaly and eccentricity.
*/
template <typename Type>
Type meanToTrue(Type &M, Type &ecc, Type &f)
{
    double eps = 1.e-012;
    Type E =  meanToEccentric(M, ecc, eps);
    f = eccentricToTrue(E, ecc);
}

#endif