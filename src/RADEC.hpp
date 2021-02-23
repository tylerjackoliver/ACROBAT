#ifndef __RADEC_H__
#define __RADEC_H__

#include <cmath>
#include <string>
#include <functional>
#include <map>

/* @brief Namespace that contains definitions of functions that determine the RA and DEC of a planet for the BME->EME conversion.
   Abstracted behind a namespace to prevent scope cluttering.
*/
namespace RADEC
{
    const double secondsToCenturies = 1./(86400. * 36525.);
    const double secondsToDays = 1./(86400.); 
    const double deg2rad = 4.0 * std::atan(1.0) / 180.;

    void earthRADEC(const double& epoch, double &alpha, double &delta)
    {
        double alpha0 = 0.;
        double delta0 = 90.0 * deg2rad;

        double timeDelta = epoch * secondsToCenturies;
        
        double alphaPrec = -0.641 * deg2rad;
        alpha = alpha0 + alphaPrec * timeDelta;

        double deltaNut = -0.557 * deg2rad;
        delta = delta0 - deltaNut * timeDelta;
    }

    void marsRADEC(const double &epoch, double &alpha, double &delta)
    {
        double alpha0 = 317.68143 * deg2rad;
        double delta0 = 52.88650 * deg2rad;

        double timeDelta = epoch * secondsToCenturies;

        double alphaPrec = -0.1061 * deg2rad;
        double deltaNut = -0.0609 * deg2rad;
        
        alpha = alpha0 + alphaPrec * timeDelta;
        delta = delta0 + deltaNut * timeDelta;
    }

    void mercuryRADEC(const double& epoch, double& alpha, double& delta)
    {
        double alpha0 = 281.0097 * deg2rad;
        double delta0 = 61.4143 * deg2rad;

        double timeDelta = epoch * secondsToCenturies;

        double alphaPrec = -0.0328 * deg2rad;
        double deltaNut = -.0049 * deg2rad;

        alpha = alpha0 + alphaPrec * timeDelta;
        delta = delta0 + deltaNut * timeDelta;
    }

    void europaRADEC(const double& epoch, double& alpha, double& delta)
    {
        const double timeDelta = epoch * secondsToCenturies;

        const double J4Term = (355.80 + 1191.3 * timeDelta) * deg2rad;
        const double J5Term = (119.90 + 262.1 * timeDelta) * deg2rad;
        const double J6Term = (229.80 + 64.3 * timeDelta) * deg2rad;
        const double J7Term = (352.25 + 2382.6 * timeDelta) * deg2rad;

        alpha = (268.20 - 0.009 * timeDelta + 1.086 * std::sin(J4Term) + 0.06 * std::sin(J5Term)
                 + 0.015 * std::sin(J6Term) + 0.009 * std::sin(J7Term)) * deg2rad;
        delta = (64.51 + 0.003 * timeDelta + 0.468 * std::cos(J4Term) + 0.026 * std::cos(J5Term)
                 + 0.007 * std::cos(J6Term) + 0.002 * std::cos(J7Term)) * deg2rad;
    }

    void europaDerivatives(const double& epoch, double& dAlpha, double& dDelta)
    {
        const double timeDelta = epoch * secondsToCenturies;

        const double J4Term = (355.80 + 1191.3 * timeDelta) * deg2rad;
        const double J5Term = (119.90 + 262.1 * timeDelta) * deg2rad;
        const double J6Term = (229.80 + 64.3 * timeDelta) * deg2rad;
        const double J7Term = (352.25 + 2382.6 * timeDelta) * deg2rad;

        dAlpha = (-0.009/secondsToCenturies + 1.086 * 1191.3 * std::cos(J4Term) + 0.060 * 262.1 * std::cos(J5Term)
                  + 0.015 * 64.3 * std::cos(J6Term) + 0.009 * 2382.6 * std::cos(J7Term)) * deg2rad * secondsToCenturies;
        dDelta = (0.003/secondsToCenturies - 0.468 * std::sin(J4Term) * 1191.3 - 0.026 * 262.1 * std::sin(J5Term)
                  - 0.007 * 64.3 * std::sin(J6Term) - 0.002 * 2382.6 * std::sin(J7Term)) * deg2rad * secondsToCenturies;
    }

    void getAlphaDelta(const double& epoch, double& alpha, double& delta)
    {
        mercuryRADEC(epoch, alpha, delta);
    }

    void mercuryDerivatives(const double& epoch, double& dAlpha, double& dDelta)
    {
    	dAlpha = -.0328 * deg2rad * secondsToCenturies;
    	dDelta = -.0049 * deg2rad * secondsToCenturies;
    }

    void alphaDeltaDerivatives(const double& epoch, double& dAlpha, double &dDelta)
    {
    	mercuryDerivatives(epoch, dAlpha, dDelta);
    }
}
#endif
