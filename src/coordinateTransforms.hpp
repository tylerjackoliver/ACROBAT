#ifndef __COORDINATE_TRANSFORMS_H__
#define __COORDINATE_TRANSFORMS_H__

#include "field.hpp"
#include <eigen3/Eigen/Core>
#include "Params.hpp"
#include "bmeField.hpp"
#include "emeField.hpp"
#include "Point.hpp"
#include "Params.hpp"
#include <stdexcept>
extern "C"
{
    #include <SpiceUsr.h>
}

/* @brief Obtains the rotation matrix for transforming from the BME frame to the EME frame at a given epoch
*  @param[in] epoch The epoch of the transformation in ephemeris seconds past J2000
*  @param[out] rot The rotation matrix to transform from the BME to the EME frame (6x6)
*/
template <typename matrixType>
void getBMEtoEMERotationMatrix(double &epoch, Eigen::Matrix<MatrixType,6,6> &rot)
{
    // Convert target string to integer ID
    SpiceDouble rotationMatrix[6][6];

    // Construct the name of the reference frame
    std::string refFrame = "IAU_"+PARAMS::TARGET;
    std::string desFrame = "J2000";

    // Call the rotation matrix generator
    sxform_c(refFrame.c_str(), desFrame.c_str(), epoch, rotationMatrix);

    // Copy the rotation matrix into the output vector
    for (unsigned i = 0; i < 6; ++i)
    {
        for (unsigned j = 0; j < 6; ++j)
        {
            rot(i,j) = rotationMatrix[i][j];
        }
    }
}

/* @brief Obtains the rotation matrix for transforming from the EME frame to the BME frame at a given epoch
*  @param[in] epoch The epoch of the transformation in ephemeris seconds past J2000
*  @param[out] rot The rotation matrix to transform from the EME to the BME frame (6x6)
*/
template <typename matrixType>
void getEMEtoBMERotationMatrix(double &epoch, Eigen::Matrix<MatrixType,6,6> &rot)
{
    // Convert target string to integer ID
    SpiceDouble rotationMatrix[6][6];

    // Construct the name of the reference frame
    std::string desFrame = "IAU_"+PARAMS::TARGET;
    std::string refFrame = "J2000";

    // Call the rotation matrix generator
    sxform_c(refFrame.c_str(), desFrame.c_str(), epoch, rotationMatrix);

    // Copy the rotation matrix into the output vector
    for (unsigned i = 0; i < 6; ++i)
    {
        for (unsigned j = 0; j < 6; ++j)
        {
            rot(i,j) = rotationMatrix[i][j];
        }
    }
}

/* @brief Converts a BME@Epoch field to an EME2000 field.
   @param[in] bmeField: bmeField to convert
   @param[out] emeField: emeField to store the conversion in.
*/
template <typename Type>
void BMEtoEME(ACROBAT::bmeField<Point<Type>> &bmeField, ACROBAT::emeField<Point<Type>> &emeField)
{
    Eigen::Matrix<Type, 6, 6> Qbe;

    // Create the transformation matrix
    getBMEtoEMERotationMatrix(PARAMS::EPOCH, Qbe);

    // Apply Qbe to every state in bmeField
    #pragma omp parallel for shared(Qbe)
    for (unsigned int i = 0; i < bmeField.getXExtent(); ++i)
    {
        for (unsigned int j = 0; j < bmeField.getYExtent(); ++j)
        {
            for (unsigned int k = 0; k < bmeField.getZExtent(); ++k)
            {
                // Set up input vector
                Eigen::Matrix<Type, 6, 1> xb, xe;
                Point<Type> temp = bmeField.getValue(i, j, k);
                
                // Assign
                for(unsigned idx = 0; idx < 6; ++idx) xb(i) = temp[i];

                // Compute
                xe = Qbe * xb;

                // Swap back
                for (unsigned idx = 0; idx < 6; ++idx) temp[i] = xe(i);

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
void EMEtoBME(ACROBAT::emeField<Point<Type>> &emeField, ACROBAT::bmeField<Point<Type>> &bmeField)
{
     Eigen::Matrix<Type, 6, 6> Qeb;

    // Create the transformation matrix
    getEMEtoBMERotationMatrix(PARAMS::EPOCH, Qeb);

    // Apply Qbe to every state in bmeField
    #pragma omp parallel for shared(Qbe)
    for (unsigned int i = 0; i < bmeField.getXExtent(); ++i)
    {
        for (unsigned int j = 0; j < bmeField.getYExtent(); ++j)
        {
            for (unsigned int k = 0; k < bmeField.getZExtent(); ++k)
            {
                // Set up input vector
                Eigen::Matrix<Type, 6, 1> xb, xe;
                Point<Type> temp = emeField.getValue(i, j, k);
                
                // Assign
                for(unsigned idx = 0; idx < 6; ++idx) xe(i) = temp[i];

                // Compute
                xb = Qbe * xe;

                // Swap back
                for (unsigned idx = 0; idx < 6; ++idx) temp[i] = xb(i);

                // Assign
                emeField.setValue(i, j, k, &temp);
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
