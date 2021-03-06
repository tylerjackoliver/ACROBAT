#ifndef __OES_H__
#define __OES_H__

#include "field.hpp"
#include <ostream>
#include "Params.hpp"
extern "C"
{
    #include "SpiceUsr.h"
}

namespace ACROBAT
{

    struct OEs
    {
        public:
            SpiceDouble rp = 0.0;
            SpiceDouble ecc = 0.0;
            SpiceDouble inc = 0.0;
            SpiceDouble longtd = 0.0;
            SpiceDouble omega = 0.0;
            SpiceDouble M = 0;
            SpiceDouble epoch = 0.0;

        SpiceDouble operator()(unsigned i)
        {
            SpiceDouble ret;
            switch(i)
            {
                case 0:
                    ret = rp;
                    break;
                case 1:
                    ret = ecc;
                    break;
                case 2:
                    ret = inc;
                    break;
                case 3:
                    ret = longtd;
                    break;
                case 4:
                    ret = omega;
                    break;
                case 5:
                    ret = M;
                    break;
                case 6:
                    ret = epoch;
                    break;
            }
            return ret;
        }
    };

    std::ostream& operator<<(std::ostream& os, const OEs& in)
    {
        return (os << "OEs: (" << in.rp << ", " << in.ecc << ", " << in.inc << ", " << in.longtd << ", " << in.omega << ", " << in.M << ", " << in.epoch << ")");
    }

    class oeField : public field2D<ACROBAT::OEs>
    {
        public:
            /* Constructor for initialising the position field
            * @param[in] nx The number of points in the \f$x\f$-direction
            * @param[in] ny The number of points in the \f$y\f$-direction
            * @param[in] nz The number of points in the \f$z\f$-direction
            */
            oeField(unsigned nx, unsigned ny): ACROBAT::field2D<ACROBAT::OEs>(nx, ny)
            {};
 
            void initialiseField(double rMin, double rMax, double omegaMin, double omegaMax)
            {
                #pragma omp parallel for
                for (unsigned i = 0; i < this->getXExtent(); ++i)
                {
                    for (unsigned j = 0; j < this->getYExtent(); ++j)
                    {
                        ACROBAT::OEs temp;
                        SpiceDouble r0, omega;

                        r0 = rMin + (rMax - rMin) / (this->nx_-1) * i;
                        omega = omegaMin + (omegaMax - omegaMin) / (this->ny_-1) * j;

                        temp.rp = r0;
                        temp.omega = omega;

                        temp.ecc = PARAMS::ECC;
                        temp.inc = PARAMS::INC;
                        temp.longtd = PARAMS::LONGTD;
                        temp.M = 0.0;                   // All trajectories start at the periapsis of their orbit
                        temp.epoch = PARAMS::EPOCH;

                        this->setValue(temp, i, j);
                    }
                }
            };
    };
}
#endif