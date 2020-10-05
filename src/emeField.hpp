#ifndef __EME_FIELD_H__
#define __EME_FIELD_H__

#include "field.hpp"
#include "SpiceUsr.h"
#include <ostream>
#include "Params.hpp"
#include "OEs.hpp"
#include "coordinateTransforms.hpp"

namespace SCROTAL
{

    template <class Type>
    class emeField : public field3D<Type>
    {
        public:

            emeField(int nx, int ny, int nz) : field3D<Type>(nx, ny, nz)
            {};

            void initialiseField(SCROTAL::oeField &input)
            {
                OEstoEME(input, this);
            }

    };

};

#endif