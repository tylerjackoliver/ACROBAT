#include <iostream>
#include "interface.hpp"
#include "field.hpp"
#include "OEs.hpp"
#include "bmeField.hpp"
#include "emeField.hpp"
#include "coordinateTransforms.hpp"
#include <unordered_map>
#include <chrono>

int main(void)
{
    welcomeMessage();
    int n = 1;
    int rank;
    int mpiStatus;

    mpiStatus = MPI_Init(NULL, NULL);

    // Initialise domain
    ACROBAT::field2D<double> domainBME(500, 500);
    ACROBAT::oeField oeDomain(500, 500);
    ACROBAT::bmeField<Point<double>> bmeDomain(500, 500);

    oeDomain.initialiseField(PARAMS::R + 100, PARAMS::RS-100, 0., 8.*std::atan(1.0)); // R, omega
    bmeDomain.initialiseField(oeDomain);

    ACROBAT::emeField<Point<double>> emeDomain(500,500);
    emeDomain.initialiseField(bmeDomain);
    emeDomain.setInitialTime(0.0);
    std::unordered_map<int, std::vector<Point<int>>> results;

    auto start = std::chrono::system_clock::now();
    emeDomain.getStableSet(1, results);
    auto end = std::chrono::system_clock::now();
    // emeDomain.outputElementsWrite(n, oeDomain, results);
    if (rank == 0){
        emeDomain.outputWrite(n, results);
        std::cout << "Completed obtaining the stable sets. The time required was " << std::chrono::duration_cast<std::chrono::seconds>(end-start).count() << " seconds.";
    }
    mpiStatus = MPI_Barrier(MPI_COMM_WORLD);
    mpiStatus = MPI_Finalize()
}
