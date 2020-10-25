#ifndef __COORDINATE_TRANSFORMS_H__
#define __COORDINATE_TRANSFORMS_H__

#include "field.hpp"
#include <eigen3/Eigen/Core>
#include "Params.hpp"
#include "bmeField.hpp"
#include "emeField.hpp"
#include "Point.hpp"
#include <stdexcept>
extern "C"
{
    #include <SpiceUsr.h>
}

/* @brief Obtains the right ascension $\alpha$ and declination $\delta$ for the PARAMS::TARGET body at a given epoch.
*  @param[out] a The right ascension
*  @param[out] d The declination
*  @param[in] epoch The epoch of the RA and DEC pair.
*/
void getSpinAxisDirection(double &a, double &d, double &epoch)
{
    // Convert target string to integer ID
    extern "C"
    {
        SpiceInt code;
        SpiceBoolean found;
        SpiceDouble RA, DEC, et, w, lambda;
        et = epoch;
        ConstSpiceChar name[] = PARAMS::TARGET;

        // Call the conversion routine
        bods2c_c(name, &code, &found);

        // Check if it was found or not
        if (!found) throw std::invalid_argument("Bad TARGET name string in getSpinAxisDirection.");

        // Now compute RA and DEC for the given identifier
        bodeul_(&code, &et, &RA, &DEC, &w, &lambda);
        a = RA;
        d = DEC;
    }
}

/* @brief Converts a BME@Epoch field to an EME2000 field.
   @param[in] bmeField: bmeField to convert
   @param[out] emeField: emeField to store the conversion in.
*/
template <typename Type>
void BMEtoEME(ACROBAT::bmeField<Point<Type>> &bmeField, ACROBAT::emeField<Point<Type>> &emeField)
{
    // Get right ascension and declination values at epoch of the bmeField
    double alpha, delta;

    // Compute the direction of the spin axis
    getSpinAxisDirection(&alpha, &delta, PARAMS::EPOCH)

    // Construct the transformation matrix
    double sina = std::sin(alpha);
    double sind = std::sin(delta);
    double cosa = std::cos(alpha);
    double cosd = std::cos(delta);

    Eigen::Matrix3d<Type> Qbe;

    // Create the transformation matrix
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

/* @brief Converts an EME2000 fiel to a BME@Epoch field
   @param[in] emeField: emeField to convert
   @param[out] bmeField: bmeField to store the conversion in.
*/
template <typename Type>
void EMEtoBME(ACROBAT::bmeField<Point<Type>> &bmeField, ACROBAT::emeField<Point<Type>> &emeField)
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

/* @brief Converts an orbital element field to a BME field.
   @param[in] oeField: ACROBAT::oeField of ACROBAT::OEs to convert
   @param[in] bmeField: ACROBAT::bmeField of Point<Type> to store the conversion.
*/
template <typename Type>
void OEstoBME(ACROBAT::oeField<ACROBAT::OEs> &oeField, ACROBAT::bmeField<Point<Type>>)
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
