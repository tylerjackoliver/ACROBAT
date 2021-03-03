#include <iostream>
#include "interface.hpp"
#define EIGEN_USE_MKL_ALL
#include "field.hpp"
#include "OEs.hpp"
#include "bmeField.hpp"
#include "emeField.hpp"
#include "coordinateTransforms.hpp"
#include <unordered_map>

int main(void)
{
    int n = 7;
    int rank;
    int mpiStatus;

    mpiStatus = MPI_Init(NULL, NULL);
    welcomeMessage();
    mpiStatus = MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::vector<double> x = {1, 2, 3, 4, 5, 6};
    std::vector<double> dx = x;

    // Initialise domain
    ACROBAT::field2D<double> domainBME(548, 360);
    ACROBAT::oeField oeDomain(548, 360);
    ACROBAT::bmeField<Point<double>> bmeDomain(548, 360);

    std::vector<double> startingPoint = {1270.563, -8042.341, 1822.657, 0.7, 0.215, 0.460};
    std::pair<int, int> status = ACROBAT::getStatus(n, startingPoint, false);
    std::cout << "First: " << status.first << " Second: " << status.second << std::endl;
    // exit(-1);

    oeDomain.initialiseField(PARAMS::R, PARAMS::RS, 0., 8.*std::atan(1.0)); // R, omega
    bmeDomain.initialiseField(oeDomain);

    ACROBAT::emeField<Point<double>> emeDomain(548, 360);
    emeDomain.initialiseField(bmeDomain);
    emeDomain.setInitialTime(0.0);
    emeDomain.setFinalTime(10.0);

    if (rank == 0){
        std::ofstream testOutput, testOutput1;
        testOutput.open("../emeDomain");
        testOutput1.open("../bmeDomain");
        for (size_t i = 0; i < 548; ++i)
        {
            for (size_t j = 0; j < 360; ++j)
            {
                Point<double> thisPoint = emeDomain.getValue(i, j);
                for (size_t dimension = 0; dimension < 6; ++dimension)
                {
                    testOutput << thisPoint.state[dimension] << " ";
                }
                testOutput << '\n';
                thisPoint = bmeDomain.getValue(i, j);
                for (size_t dimension = 0; dimension < 6; ++dimension)
                {
                    testOutput1 << thisPoint.state[dimension] << " ";
                }
                testOutput1 << '\n';
            }
        }
        testOutput.close();
        testOutput1.close();
    }

    MPI_Barrier(MPI_COMM_WORLD);
    std::unordered_map<int, std::vector<std::vector<std::pair<int, int>>>> forwardSet, backwardSet;
    ACROBAT::initialiseSetMap(n, -1, forwardSet);
    ACROBAT::initialiseSetMap(n, -1, backwardSet);
    std::vector<std::vector<double>> points;
    ACROBAT::getSets(n, emeDomain, forwardSet, backwardSet);

    std::cout << "The number of points in W_1 is " << forwardSet[0][2].size() << std::endl;
    std::cout << "The number of points in W_-1 is " << backwardSet[0][2].size() << std::endl;

    mpiStatus = MPI_Barrier(MPI_COMM_WORLD);
    for (uint16_t i = 0; i < 4; ++i)
    {
        std::cout << rank << " forward " << i << " " << forwardSet[0][i].size() << std::endl;
        std::cout << rank << " backward " << i << " " << backwardSet[0][i].size() << std::endl;
    }

    mpiStatus = MPI_Barrier(MPI_COMM_WORLD);
    for (size_t i = 0; i <= n; ++i) std::cout << "rank " << rank << " i " << i << " " << forwardSet[i][2].size() << std::endl;
    
    for (unsigned i = 0; i <= n; ++i)
    {
        points.clear();
        ACROBAT::extractStableSet(i, forwardSet, backwardSet, bmeDomain, points);
        mpiStatus = MPI_Barrier(MPI_COMM_WORLD);
        broadcastMapping(points);
        MPI_Barrier(MPI_COMM_WORLD);
        std::cout << "rank " << rank << " final points for n=" << i<< " " << points.size() << std::endl;
    }



    if (rank == 0)
    {
        ACROBAT::writeSetMapping(n, points);
        ACROBAT::writeAllSets(n, forwardSet, backwardSet, bmeDomain);
    }
    // // mpiStatus = MPI_Barrier(MPI_COMM_WORLD);

    MPI_Finalize();
    return 0;
}
