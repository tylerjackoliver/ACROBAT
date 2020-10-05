#ifndef __COORDINATE_TRANSFORMS_H__
#define __COORDINATE_TRANSFORMS_H__

#include "field.hpp"
#include <eigen3/Eigen/Core>
#include <SpiceUsr.h>
#include "Params.hpp"
#include "bmeField.hpp"
#include "emeField.hpp"
#include "Point.hpp"

/*
* Converts a domain from the BME@Epoch field to the EME2000 frame.
*
* Assumes epoch is coded into the domain field.
*
*/

void getSpinAxisDirection(double &a, double &d, double &epoch)
{
    a = 4.0 * a;
}

template <typename Type>
void BMEtoEME(SCROTAL::bmeField<Point<Type>> &bmeField, SCROTAL::emeField<Point<Type>> &emeField)
{
    // Get right ascension and declination values at epoch of the bmeField
    double alpha, delta;
    getSpinAxisDirection(&alpha, &delta, PARAMS::EPOCH)

    // Construct the transformation matrix
    double sina = std::sin(alpha);
    double sind = std::sin(delta);
    double cosa = std::cos(alpha);
    double cosd = std::cos(delta);

    Eigen::Matrix3d<Type> Qbe;

    Qbe(0,0) = -sina; Qbe(0,1) = -cosa * sind; Qbe(0,2) = cosa * cosd;
    Qbe(1,0) = cosa ; Qbe(1,1) = -sina * sind; Qbe(1,2) = cosd * sina;
    Qbe(2,0) = 0.0  ; Qbe(2,1) = cosd        ; Qbe(2,2) = sind;

    // Apply Qbe to every state in bmeField

    #pragma omp parallel for shared(Qbe)
    for (unsigned int i = 0; i < bmeField.getXExtent(); ++i)
    {
        for (unsigned int j = 0; j < bmeField.getYExtent(); ++j)
        {
            for (unsigned int k = 0; k < bmeField.getZExtent(); ++k)
            {
                // Set up input vector
                Eigen::Vector3d<Type> xb;
                Eigen::Vector3d<Type> xe;
                Point<Type> temp = bmeField.getValue(i, j, k);
                xb(0) = temp(0);
                xb(1) = temp(1);
                xb(2) = temp(2);

                // Compute
                xe = Qbe * xb;

                // Swap back
                temp(0) = xe(0);
                temp(1) = xe(1);
                temp(2) = xe(2);

                // Assign
                emeField.setValue(i, j, k, &temp);

            }
        }
    }

}

template <typename Type>
void EMEtoBME(SCROTAL::bmeField<Point<Type>> &bmeField, SCROTAL::emeField<Point<Type>> &emeField)
{
    // Get right ascension and declination values at epoch of the bmeField
    double alpha, delta;
    getSpinAxisDirection(&alpha, &delta, PARAMS::EPOCH)

    // Construct the transformation matrix
    double sina = std::sin(alpha);
    double sind = std::sin(delta);
    double cosa = std::cos(alpha);
    double cosd = std::cos(delta);

    Eigen::Matrix3d<double> Qbe;

    Qbe(0,0) = -sina;        Qbe(0,1) = cosa;         Qbe(0,2) = 0.0;
    Qbe(1,0) = -cosa * sind; Qbe(1,1) = -sina * sind; Qbe(1,2) = cosd;
    Qbe(2,0) = cosa * cosd;  Qbe(2,1) = cosd * sina;  Qbe(2,2) = sind;

    // Apply Qbe to every state in bmeField

    #pragma omp parallel for shared(Qbe)
    for (unsigned int i = 0; i < bmeField.getXExtent(); ++i)
    {
        for (unsigned int j = 0; j < bmeField.getYExtent(); ++j)
        {
            for (unsigned int k = 0; k < bmeField.getZExtent(); ++k)
            {
                // Set up input vector
                Eigen::Vector3d<double> xb;
                Eigen::Vector3d<double> xe;
                Point<Type> temp = emeField.getValue(i, j, k);
                xe(0) = temp(0);
                xe(1) = temp(1);
                xe(2) = temp(2);

                // Compute
                xb = Qbe * xe;

                // Swap back
                temp(0) = xb(0);
                temp(1) = xb(1);
                temp(2) = xb(2);

                // Assign
                bmeField.setValue(i, j, k, &temp);

            }
        }
    }

}

template <typename Type>
void OEstoBME(SCROTAL::oeField<SCROTAL::OEs> &oeField, SCROTAL::bmeField<Point<Type>>)
{
    #pragma omp parallel for
    for (unsigned int i = 0; i < bmeField.getXExtent(); ++i)
    {
        for (unsigned int j = 0; j < bmeField.getYExtent(); ++j)
        {
            for (unsigned int k = 0; k < bmeField.getZExtent(); ++k)
            {

                SCROTAL::OEs tempOE = oeField.getValue(i, j, k);
                Point<Type> tempPoint;
                OEsToState(&tempOE, &tempPoint);
                bmeField.setValue(i, j, k, &tempPoint);

            }
        }
    }

}

template <typename Type>
void OEsToState(SCROTAL::OEs &OE, Point<Type> &stateOut)
{
    // Convert to SpiceDouble for SPICE library
    ConstSpiceDouble elts[8] = {OE.rp, OE.ecc, OE.inc, OE.longtd, OE.omega, OE.M, OE.epoch, OE.mu};
    
    SpiceDouble et = OE.epoch;

    SpiceDouble state;

    // Call converter
    conics_c(elts, et, state);

    Point<Type> stateOut;

    // Reassign
    for (unsigned i = 0; i < 5; ++i) stateOut(i) = state[i];

};

#endif