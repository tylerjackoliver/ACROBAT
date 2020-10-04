#ifndef __INTERFACE_H__
#define __INTERFACE_H__

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

void welcomeMessage()
{

    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    std::cout << "             WELCOME TO SCROTAL              " << std::endl;
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

};

#endif
