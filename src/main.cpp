#include <iostream>
#include "interface.hpp"
#include "field.hpp"
#include "OEs.hpp"
#include "bmeField.hpp"
// #include "coordinateTransforms.hpp"

int main(void)
{
    welcomeMessage();

    // Initialise domain
    SCROTAL::field2D<double> domainBME(500, 500);
    SCROTAL::oeField oeDomain(500, 500);
    SCROTAL::bmeField<double> bmeDomain(500, 500);

    oeDomain.initialiseField(1., 2., 1., 2.);
    bmeDomain.initialiseField(oeDomain);

    std::cout << bmeDomain.getValue(50, 50);
}
