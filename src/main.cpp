#include <iostream>
#include "interface.hpp"
#include "field.hpp"
#include "OEs.hpp"
#include "bmeField.hpp"
#include "emeField.hpp"
#include "coordinateTransforms.hpp"

int main(void)
{
    welcomeMessage();

    // Initialise domain
    ACROBAT::field2D<double> domainBME(500, 500);
    ACROBAT::oeField oeDomain(500, 500);
    ACROBAT::bmeField<Point<double>> bmeDomain(500, 500);

    oeDomain.initialiseField(1., 2., 1., 2.);
    bmeDomain.initialiseField(oeDomain);

    ACROBAT::emeField<Point<double>> emeDomain(500,500);
    emeDomain.initialiseField(bmeDomain);
    std::cout << emeDomain.getValue(2, 2) << std::endl;
}
