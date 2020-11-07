#ifndef __MPI_WRAPPERS_H__
#define __MPI_WRAPPERS_H__

#include <mpi.h>
#include "Point.hpp"

void broadcastVector(std::vector<Point<int>>& vectorToBroadcast, int& rank)
{
    MPI_Status mpiStatus;
    int returnCode;
    int numToBroadcast;
    /* First, send the number of broadcasts we're going to have to do */
    if (rank != 0) vectorToBroadcast.clear();
    returnCode = MPI_Bcast(numToBroadcast, 1, MPI_INT, 0, MPI_COMM_WORLD);
    for (size_t idx = 0; idx < numToBroadcast; ++idx)
    {
        std::vector<int> destVec(tmpVec.size());
        Point<int> tmpPoint;
        for (size_t dimension = 0; dimension < 6; ++dimension)
        {
            int tmp = vectorToBroadcast[idx].state[dimension];
            returnCode = MPI_Bcast(tmp, 1, MPI_INT, 0, MPI_COMM_WORLD);
            tmpPoint.state[dimension] = tmp;
        }
        if (rank != 0) vectorToBroadcast.push_back(tmpPoint);
    }
}

void reduceVector(std::vector<Point<int>>& vectorToReduce, std::vector<Point<int>>& vectorToStore, int& rank, int& poolSize)
{
    MPI_Status mpiStatus;
    int returnCode;
    int numToReduce = 0;

    if (rank == 0)
    {
        int numToReceive;
        for (int worker = 0; worker < poolSize; ++worker)
        {
            returnCode = MPI_Recv(numToReceive, 1, MPI_INT, worker, worker, MPI_COMM_WORLD, &mpiStatus);
            Point<int> tmpPoint;
            std::vector<int> tmpState;
            for (int toReceive = 0; toReceive < numToReceive, ++toReceive)
            {
                for (int dimension = 0; dimension < 6; ++dimension)
                {
                    int tmpValue;
                    returnCode = MPI_Recv(tmpValue, 1, MPI_INT, worker, worker, MPI_COMM_WORLD, &mpiStatus);
                    tmpState[dimenson] = tmpValue;
                }
                tmpPoint.state = tmpState;
                vectorToReduce.push_back(tmpPoint);
            }
        }
    } else
    {
        int numToSend = vectorToReduce.size();
        returnCode = MPI_Send(numToSend, 1, MPI_INT, 0, rank, MPI_COMM_WORLD);
        for (int toSend = 0; toSend < numToSend; ++toSend)
        {
            for (int dimension = 0; dimension < 6; ++dimension)
            {
                int tmpValue = vectorToReduce[toSend].state[dimension];
                returnCode = MPI_Send(tmpValue, 1, MPI_INT, 0, worker, MPI_COMM_WORLD);
            }
        }
    }
}

void reduceCount(int& myValue, int& reducedValue)
{
    MPI_Status mpiStatus;
    int returnCode;
    int numToReduce = 0;
    reducedValue = 0;

    if (rank == 0)
    {
        int numToReceive;
        for (int worker = 0; worker < poolSize; ++worker)
        {
            int tmpValue;
            returnCode = MPI_Recv(tmpValue, 1, MPI_INT, worker, worker, MPI_COMM_WORLD, &mpiStatus);
            reducedValue += tmpValue;
        }
    } else
    {
        returnCode = MPI_Send(myValue, 1, MPI_INT, 0, worker, MPI_COMM_WORLD);
    }
}

#endif