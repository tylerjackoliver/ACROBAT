#ifndef __EME_FIELD_H__
#define __EME_FIELD_H__

#include "field.hpp"
#include <ostream>
#include "Params.hpp"
#include "OEs.hpp"
#include "coordinateTransforms.hpp"
#include "integration.hpp"

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

            /* @brief Initialises the elements in the field using an OEField.
               @params[in] ACROBAT::oeField containing the orbital elements defining the domain
            */
            void initialiseField(ACROBAT::oeField &input)
            {
                OEstoEME(input, this);
            }

            /* @brief Returns the final time for the trajectory integration
               @returns The final time for the integration.
            */
            double getFinalTime() const
            {
                return this->_finalTime;
            }

            /* @brief Returns the initial time for the trajectory integration
               @returns The initial time for the integration
            */
            double getInitialTime() const
            {
                return this->_initialTime;
            }

            /* @brief Sets the final time for the trajectory integration
               @param[in] finalTime The final time for the trajectory integration
            */
            void setFinalTime(const double finalTime)
            {
                this->_finalTime = finalTime;
            }

            /* @brief Sets the initial time for the trajectory integration
               @param[in] initialTime The initial time for the trajectory integration
            */
            void setinitialTime(const double initialTime)
            {
                this->_initialTime = initialTime;
            }

            void getStableSet(unsigned int, std::unordered_map<int, std::vector<Point<Type>>>);
        private:
            double _finalTime;      // Final time for the trajectory integration
            double _initialTime;    // Initial time for the trajectory integration
    };
};


/* @brief Obtains the stabNumber-revolution stable set for a given EME2000 field. Note than stabNumber entries in the map are returned.
 * @param[in] stabNumber The number of revolutions a particle should complete in order to be classed as stable.
 * @param[out] std::unordered_map The keys are the number of rotations, and the values are a std::vector<> containing the points comprising the set.
 */
template<typename Type>
void emeField<Type>::getStableSet(unsigned int stabNumber, std::unordered_map<int, std::vector<Point<Type>>>& out);
{

    for (unsigned revs = 1; revs <= stabNumber; ++revs)  // Must complete at least one orbit about the host planet
    {
        std::vector<Point<Type>> points;
        getSet(revs, this, points);
        out[revs] = points;
    }
}

#endif