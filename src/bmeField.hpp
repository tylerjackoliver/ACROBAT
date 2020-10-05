#ifndef __BME_FIELD_H__
#define __BME_FIELD_H__

#include "field.hpp"
#include "SpiceUsr.h"
#include <ostream>
#include "Params.hpp"
#include "OEs.hpp"

namespace SCROTAL
{

    template <class Type>
    class bmeField : public field2D<Point<Type>>
    {
        public:

            bmeField(int nx, int ny, int nz) : field2D<Point<Type>>(nx, ny)
            {};

            void initialiseField(SCROTAL::oeField &input)
            {
                #pragma omp parallel for
                for (unsigned int i = 0; i < this->getXExtent(); ++i)
                {
                    for (unsigned int j = 0; j < this->getYExtent(); ++j)
                    {
                        SCROTAL::OEs tempOE = input.getValue(i, j);
                        Point<Type> tempPoint;
                        this->OEsToState(tempOE, tempPoint);
                        this->setValue(tempPoint, i, j);
                    }
                }
            }

        private:

            template <typename Type>
            void OEsToState(SCROTAL::OEs &OE, Point<Type> &stateOut)
            {
                // Convert to SpiceDouble for SPICE library
                ConstSpiceDouble elts[8] = {OE.rp, OE.ecc, OE.inc, OE.longtd, OE.omega, OE.M, OE.epoch, PARAMS::GM};
                
                SpiceDouble et = OE.epoch;

                SpiceDouble state[6];

                // Call converter
                conics_c(elts, et, state);

                // Reassign
                for (unsigned i = 0; i < 5; ++i) stateOut[i] = state[i];

            };

    };

};

#endif