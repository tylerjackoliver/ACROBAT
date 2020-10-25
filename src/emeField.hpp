#ifndef __EME_FIELD_H__
#define __EME_FIELD_H__

#include "field.hpp"
#include <ostream>
#include "Params.hpp"
#include "OEs.hpp"
#include "coordinateTransforms.hpp"

extern "C"
{
    #include "SpiceUsr.h"
}

namespace ACROBAT
{

    template <class Type>
    class emeField : public field3D<Type>
    {
        public:
            emeField(int nx, int ny, int nz) : field3D<Type>(nx, ny, nz)
            {};

            void initialiseField(ACROBAT::oeField &input)
            {
                OEstoEME(input, this);
            }
    };

};

#endif