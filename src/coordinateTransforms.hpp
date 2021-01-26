#ifndef __COORDINATE_TRANSFORMS_H__
#define __COORDINATE_TRANSFORMS_H__

#include "field.hpp"
#include <eigen3/Eigen/Core>
#include "Params.hpp"
#include "bmeField.hpp"
#include "Point.hpp"
#include "Params.hpp"
#include <stdexcept>
#include "RADEC.hpp"
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

/* @brief Obtains the rotation matrix for transforming from the BME frame to the EME frame at a given epoch
*  @param[in] epoch The epoch of the transformation in ephemeris seconds past J2000
*  @param[out] rot The rotation matrix to transform from the BME to the EME frame (6x6)
*/
template <typename matrixType>
void getBMEtoEMERotationMatrix(const double &epoch, Eigen::Matrix<matrixType,3,3> &rot)
{
    // Get the current values of alpha and delta at this time
    double alpha, delta;
    RADEC::getAlphaDelta(epoch, alpha, delta);

    /* Construct the matrix. In theory, we would need to also construct the derivative
       of the rotation matrix in order to get the velocity rotated into the frame, but
       since the derivative of the rotation matrix is of magnitude 1e-10...we ignore it.
    */
    double sina = std::sin(alpha);
    double sind = std::sin(delta);

    double cosa = std::cos(alpha);
    double cosd = std::cos(delta);

    rot(0, 0) = -sina; rot(0,1) = -cosa * sind; rot(0,2) = cosa * cosd;
    rot(1, 0) = cosa;  rot(1,1) = -sina * sind; rot(1,2) = cosd * sina;
    rot(2, 0) = 0.0;   rot(2,1) = cosd;         rot(2,2) = sind;
}


/* @brief Obtains the rotation matrix and its derivative for transforming from the BME frame to the EME frame at a given epoch.
 * @param[in] epoch The epoch of the transformation in ephemeris seconds past J2000
 * @param[out] rot The rotation matrix to transform from the BME to EME frame
 * @param[out] dRot The derivative of the rotation matrix rot
 */
//template <typename epochType, typename matrixType>
void getBMEtoEMERotationMatrix(const double& epoch, Eigen::Matrix<double, 3, 3>& rot, Eigen::Matrix<double, 3, 3>& dRot)
{
	double alpha, delta, alphaDeriv, deltaDeriv;
	RADEC::getAlphaDelta(epoch, alpha, delta);
	RADEC::alphaDeltaDerivatives(alphaDeriv, deltaDeriv); // Derivatives with respect to time

	/* Construct the rotation matrix */
	double sina = std::sin(alpha);
	double sind = std::sin(delta);
	double cosa = std::cos(alpha);
	double cosd = std::cos(delta);

    rot(0, 0) = -sina; rot(0,1) = -cosa * sind; rot(0,2) = cosa * cosd;
    rot(1, 0) = cosa;  rot(1,1) = -sina * sind; rot(1,2) = cosd * sina;
    rot(2, 0) = 0.0;   rot(2,1) = cosd;         rot(2,2) = sind;
    /* And now the derivative of the dRotation matrix */
    dRot(0, 0) = -cosa; dRot(0, 1) = sina * sind * alphaDeriv - cosa * cosd * deltaDeriv; dRot(0, 2) = -(sina * cosd * alphaDeriv + deltaDeriv * sind * cosd);
    dRot(1, 0) = -sina; dRot(1, 1) = -(cosa * sind * alphaDeriv + deltaDeriv * sina * cosd); dRot(1, 2) = -(sind * sina * deltaDeriv - cosd * cosa * alphaDeriv);
    dRot(2, 0) = 0.0; dRot(2, 1) = -sind; dRot(2,2) = cosd;
}


/* @brief Obtains the rotation matrix for transforming from the EME frame to the BME frame at a given epoch
*  @param[in] epoch The epoch of the transformation in ephemeris seconds past J2000
*  @param[out] rot The rotation matrix to transform from the EME to the BME frame (6x6)
*/
template <typename matrixType>
void getEMEtoBMERotationMatrix(double &epoch, Eigen::Matrix<matrixType,6,6> &rot)
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

#endif
