#ifndef __BME_FIELD_H__
#define __BME_FIELD_H__

#include "field.hpp"
#include <ostream>
#include "Params.hpp"
#include "OEs.hpp"
#include <stdexcept>
#include <stdio.h>
extern "C"
{
#include "SpiceUsr.h"
}

namespace ACROBAT
{
    template <class Type>
    class bmeField : public field2D<Point<Type>>
    {
        public:
            bmeField(int nx, int ny) : field2D<Point<Type>>(nx, ny)
            {};

            void initialiseField(ACROBAT::oeField &input)
            {
                // Check input sizes
                if ( (input.getXExtent() != this->getXExtent() ) || (input.getYExtent() != this->getYExtent()) ){
                    throw std::out_of_range("OE Field and BME Field not the same size.");
                }

                for (unsigned int i = 0; i < this->getXExtent(); ++i)
                {
                    for (unsigned int j = 0; j < this->getYExtent(); ++j)
                    {
                        ACROBAT::OEs tempOE = input.getValue(i, j);
                        Point<Type> tempPoint;
                        OEsToState(tempOE, tempPoint);
                        this->setValue(tempPoint, i, j);
                    }
                }
            }
    };

};

extern "C" void OEsToState(ACROBAT::OEs &OE, Point<double> &stateOut)
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
#endif