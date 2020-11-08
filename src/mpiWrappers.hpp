#ifndef __MPI_WRAPPERS_H__
#define __MPI_WRAPPERS_H__

#include <mpi.h>
#include "Point.hpp"

void broadcastVector(std::vector<Point<int>>& vectorToBroadcast, int& rank)
{
    MPI_Status mpiStatus;
    int returnCode;
    int numToBroadcast = vectorToBroadcast.size();
    /* First, send the number of broadcasts we're going to have to do */
    if (rank != 0) vectorToBroadcast.clear();
    returnCode = MPI_Bcast(&numToBroadcast, 1, MPI_INT, 0, MPI_COMM_WORLD);
    for (size_t idx = 0; idx < numToBroadcast; ++idx)
    {
        Point<int> tmpPoint;
        for (size_t dimension = 0; dimension < 6; ++dimension)
        {
            int tmp = vectorToBroadcast[idx].state[dimension];
            returnCode = MPI_Bcast(&tmp, 1, MPI_INT, 0, MPI_COMM_WORLD);
            tmpPoint.state[dimension] = tmp;
        }
        if (rank != 0) vectorToBroadcast.push_back(tmpPoint);
    }
}

void reduceVector(std::vector<Point<int>>& vectorToReduce, std::vector<Point<int>>& vectorToStore, int& rank, int& poolSize)
{
    MPI_Status mpiStatus;
    int returnCode;

    if (rank == 0)
    {
        int numToReceive;
        for (int worker = 1; worker < poolSize; ++worker)
        {
            returnCode = MPI_Recv(&numToReceive, 1, MPI_INT, worker, worker, MPI_COMM_WORLD, &mpiStatus);
            Point<int> tmpPoint;
            for (int toReceive = 0; toReceive < numToReceive; ++toReceive)
            {
                std::vector<int> tmpState(6);
                for (int dimension = 0; dimension < 6; ++dimension)
                {
                    int tmpValue;
                    returnCode = MPI_Recv(&tmpValue, 1, MPI_INT, worker, worker, MPI_COMM_WORLD, &mpiStatus);
                    tmpState[dimension] = tmpValue;
                }
                tmpPoint.state = tmpState;
                vectorToReduce.push_back(tmpPoint);
            }
        }
    } else
    {
        int numToSend = vectorToReduce.size();
        returnCode = MPI_Send(&numToSend, 1, MPI_INT, 0, rank, MPI_COMM_WORLD);
        for (int toSend = 0; toSend < numToSend; ++toSend)
        {
            for (int dimension = 0; dimension < 6; ++dimension)
            {
                int tmpValue = vectorToReduce[toSend].state[dimension];
                returnCode = MPI_Send(&tmpValue, 1, MPI_INT, 0, rank, MPI_COMM_WORLD);
            }
        }
    }
}

void reduceCount(int& myValue, int& reducedValue)
{
    MPI_Status mpiStatus;
    int returnCode;
    int rank;
    int poolSize;
    returnCode = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    returnCode = MPI_Comm_size(MPI_COMM_WORLD, &poolSize);
    if (rank == 0)
    {
        for (int worker = 1; worker < poolSize; ++worker)
        {
            int tmpValue;
            returnCode = MPI_Recv(&tmpValue, 1, MPI_INT, worker, worker, MPI_COMM_WORLD, &mpiStatus);
            reducedValue += tmpValue;
        }
    } else
    {
        returnCode = MPI_Send(&myValue, 1, MPI_INT, 0, rank, MPI_COMM_WORLD);
    }
}

void broadcastCount(int& countToBroadcast)
{
    int poolSize;
    int rank;
    int returnCode;

    returnCode = MPI_Bcast(&countToBroadcast, 1, MPI_INT, 0, MPI_COMM_WORLD);
}




#endif
