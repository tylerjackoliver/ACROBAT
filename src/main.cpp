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
    int n = 6;
    int rank;
    int mpiStatus;

    mpiStatus = MPI_Init(NULL, NULL);
    welcomeMessage();
    mpiStatus = MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::vector<double> x = {1, 2, 3, 4, 5, 6};
    std::vector<double> dx = x;

    // Initialise domain
    ACROBAT::field2D<double> domainBME(100,100);
    ACROBAT::oeField oeDomain(100,100);
    ACROBAT::bmeField<Point<double>> bmeDomain(100,100);

    oeDomain.initialiseField(PARAMS::R, PARAMS::RS, 0., 8.*std::atan(1.0)); // R, omega
    bmeDomain.initialiseField(oeDomain);

    ACROBAT::emeField<Point<double>> emeDomain(100,100);
    emeDomain.initialiseField(bmeDomain);
    emeDomain.setInitialTime(0.0);
    std::unordered_map<int, std::vector<Point<int>>> results;

    auto start = std::chrono::system_clock::now();
    emeDomain.getStableSet(n, results);
    auto end = std::chrono::system_clock::now();
    // emeDomain.outputElementsWrite(n, oeDomain, results);
    if (rank == 0){
        emeDomain.outputWrite(n, results);
        std::cout << "Completed obtaining the stable sets. The time required was " << std::chrono::duration_cast<std::chrono::seconds>(end-start).count() << " seconds.";
    }
    mpiStatus = MPI_Barrier(MPI_COMM_WORLD);
    mpiStatus = MPI_Finalize();
}
