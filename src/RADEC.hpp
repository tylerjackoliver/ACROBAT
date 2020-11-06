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

    void getAlphaDelta(const double& epoch, double& alpha, double& delta)
    {
        marsRADEC(epoch, alpha, delta);
    }
}
#endif
