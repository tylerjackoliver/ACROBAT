#ifndef __LAGUERRE_CONWAY_H__
#define __LAGUERRE_CONWAY_H__

template <typename Type, typename Function>
void laguerreConway(Type &x0, Type &xf, Type &eps, Function f, Function fp, Function fpp)
{
    const int n = 5;        // Tuning parameter - 5 for now, as per paper
    Type tolerance = 100.   // Stopping tolerance
    Type xi = x0;

    while (tolerance >= eps)
    {   
        // Pre-compute function evaluations and derivatives
        Type fval = f(xi);
        Type deriv = fp(xi);
        Type dDeriv = fpp(xi);
        // Compute numerator, square root
        Type numerator = - n * f(x0);
        Type root = std::sqrt( fabs((n-1) * (n-1) * (fp * fp) - n * (n-1) * fval * dDeriv) );
        // Denominator is such that absolute value is maximised
        Type denominator = std::max( fabs(fp + root), fabs(fp - root) );
        // Compute update to iterate
        Type delta_n1 = numerator / denominator;
        xi += delta_n1;
        tolerance = delta_n1; // Change in latest iterate
    };
    xf = xi;
}

#endif