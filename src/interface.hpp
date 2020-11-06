#ifndef __INTERFACE_H__
#define __INTERFACE_H__

#include "Params.hpp"
#include <cmath>
#include <mpi.h>

extern "C"
{
    #include "SpiceUsr.h"
}

struct opts
{
    // Domain settings
    double xDomainMin = 0.0;
    double xDomainMax = 0.0;

    // Domain settings
    double yDomainMin = 0.0;
    double yDomainMax = 0.0;

    // Discretisation in x- and y-
    unsigned long nX = 0;
    unsigned long nY = 0;

    // Mass parameter
    double mu = 0.;
};

opts OPTIONS;


/* @brief Computes derived problem parameters from the independent problem parameters given in Params.hpp
   
   The independent problem variables are the target and host planets, defined by std::strings of their IAU
   or common identifier, the epoch of the problem, given in ephemeris seconds past J2000, and the eccentricity,
   longitude and inclination of the target ballistic capture orbit.

   The remaining variables - TARGET planetary radii, sphere of influence and graviational parameter - are derived
   from the IAU common string via SPICE kernels. The HOST gravivational parameter is also derived using SPICE kernels.
*/
void initialiseParams()
{
    // Initialise temporary outputs
    double hostRadii[3];
    double avgHostRadii = 0.0;
    double targetRadii[3];
    double avgTargetRadii = 0.0;
    int itemsReturned;
    double positionVector[6];
    double oes[8];
    double lt; // Temp lighttime entry
    double euclideanNorm = 0.0;

    SpiceDouble tempTargetGM;
    SpiceDouble tempHostGM;

    /* Get GM for the HOST and TARGET planets */
    bodvrd_c(PARAMS::TARGET.c_str(), "GM", 1, &itemsReturned, &PARAMS::targetGM); // (body, value, max items returned, actual items returned, where to store)
    bodvrd_c(PARAMS::HOST.c_str(), "GM", 1, &itemsReturned, &PARAMS::hostGM);

    /* Get the radii for the TARGET planet */
    bodvrd_c(PARAMS::TARGET.c_str(), "RADII", 3, &itemsReturned, targetRadii);

    /* Compute average Radii from the given directions */
    for (size_t idx = 0; idx < itemsReturned; ++idx)
    {
        avgTargetRadii += targetRadii[idx];
    }
    avgTargetRadii /= itemsReturned;
    PARAMS::R = avgTargetRadii;

    /* Compute the sphere of influence of the given planet - assumed to be a perfect sphere.
       For this, we need the distance between the HOST and the TARGET. We'll therefore use SPKPOS_C
       to get the position of the TARGET around the HOST (although the inverse would be perfecrtly
       legitimate.) We can then compute the L2-norm of this vector to give the position, and then
       compute the sphere of influence from that. */
    spkezr_c(PARAMS::TARGET.c_str(), PARAMS::EPOCH, "J2000", "NONE", PARAMS::HOST.c_str(), positionVector, &lt);
    oscelt_c(positionVector, PARAMS::EPOCH, PARAMS::hostGM, oes);

    for (size_t idx = 0; idx < 3; ++idx)
    {
        euclideanNorm += positionVector[idx] * positionVector[idx];
    }
    euclideanNorm = std::sqrt(euclideanNorm);
    
    // Finally, can compute the sphere of influence (R_SOI = (mTarg/mHost)^(2/5) * R_distance)
    PARAMS::RS = std::pow(PARAMS::targetGM / PARAMS::hostGM, 2./5.) * (oes[0] / (1-oes[1]));

    /* Now, we can get the Mean anomaly of the TARGET around the HOST using the orbital elements of TARGET about HOST
      and the given computational epoch */
    PARAMS::M = oes[5];

    /* And set the maximum time for the trajectory to be acrobatic in non-dimensional units */
    PARAMS::maxT = 8.0 * (4.0 * std::atan(1.0)) * std::pow(PARAMS::RS/PARAMS::R, 1.5);
}

void welcomeMessage()
{
    int rank, status;
    status = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0)
    {
        std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
        std::cout << "             WELCOME TO ACROBAT              " << std::endl;
        std::cout << "             ~~~~~~~~~~~~~~~~~~              " << std::endl;
        std::cout << "  The balliStic CaptuRe OrbiT Analysis tooL  " << std::endl;
        std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
        std::cout << std::endl;
        std::cout << "Author: Jack Tyler. jack.tyler@soton.ac.uk   " << std::endl;
        std::cout << "Version: 0.0.1, Oct 1 2020                   " << std::endl;
        std::cout << std::endl;
        std::cout << "USER PARAMETERS                              " << std::endl;
        std::cout << "~~~~~~~~~~~~~~~                              " << std::endl;
        std::cout << "X-Domain: [" << OPTIONS.xDomainMin << ", " << OPTIONS.xDomainMax << "]" << std::endl;
        std::cout << "Y-Domain: [" << OPTIONS.yDomainMin << ", " << OPTIONS.yDomainMax << "]" << std::endl;
        std::cout << "Discretisations in X: " << OPTIONS.nX          << std::endl;
        std::cout << "Discretisations in Y: " << OPTIONS.nY          << std::endl;
        std::cout << "System mass parameter: " << OPTIONS.mu         << std::endl;
    }
    // Load the ephemeris files here
    furnsh_c("../spice/metakernel.tm");
    // Initialise parameters
    initialiseParams();
};

#endif
