#ifndef __COORDINATE_TRANSFORMS_H__
#define __COORDINATE_TRANSFORMS_H__

#include "field.hpp"
#include <eigen3/Eigen/Core>
#include "Params.hpp"
#include "bmeField.hpp"
#include "Point.hpp"
#include "Params.hpp"
#include <stdexcept>
extern "C"
{
    #include <SpiceUsr.h>
}

/* @brief Converts an orbital element field to a BME field.
   @param[in] oeField: ACROBAT::oeField of ACROBAT::OEs to convert
   @param[in] bmeField: ACROBAT::bmeField of Point<Type> to store the conversion.
*/
template <typename Type>
void OEstoBME(ACROBAT::oeField &oeField, ACROBAT::bmeField<Point<Type>>)
{
    for (unsigned int i = 0; i < bmeField.getXExtent(); ++i)
    {
        for (unsigned int j = 0; j < bmeField.getYExtent(); ++j)
        {
            for (unsigned int k = 0; k < bmeField.getZExtent(); ++k)
            {
                ACROBAT::OEs tempOE = oeField.getValue(i, j, k);
                Point<Type> tempPoint;
                OEsToState(&tempOE, &tempPoint);
                bmeField.setValue(i, j, k, &tempPoint);
            }
        }
    }
}
#endif
