#ifndef __LAGUERRE_CONWAY_H__
#define __LAGUERRE_CONWAY_H__

/* Implements the Laguerre-Conway method for solving Kepler's equation */
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
    Ef = xi;
}

/* Computes Kepler's equation for a given input value X */
template <typename Type>
Type keplersEquation(const Type E, const Type ecc, const Type M)
{
    return E - ecc * std::sin(E) - M;
}

/* Computes the derivative of Kepler's equation for a given input value X */
template <typename Type>
Type dMdE(const Type E, const Type ecc)
{
    return 1 - ecc * std::cos(E);
}

/* Computes the second derivative of Kepler's equation */
template <typename Type>
Type dMMdEE(const Type E, const Type ecc)
{
    return ecc * std::sin(E);
}

/* Wrapper: Converts from Mean anomaly to Eccentric anomaly */
template <typename Type>
Type meanToEccentric(Type &M, Type &ecc, Type &eps)
{
    laguerreConway(M, E, ecc, M, eps, keplersEquation, dMdE, dMMdEE);
}

/* Converts from Eccentric anomaly to true anomaly */
template <typename Type>
Type eccentricToTrue(Type &E, Type &ecc)
{
    double pi = 4.0 * atan(1.0);
    Type root = sqrt::( (1-ecc) / (1+ecc) );
    Type inv = root * std::tan(E / 2.0);
    return ( 2.0 * std::atan(inv) ) % 2. * pi;
}

/* Wrapper: Converts from Mean to True anomaly */
template <typename Type>
Type meanToTrue(Type &M, Type &ecc, Type &f)
{
    double eps = 1.e-012;
    Type E =  meanToEccentric(M, ecc, eps);
    f = eccentricToTrue(E, ecc);
}

#endif