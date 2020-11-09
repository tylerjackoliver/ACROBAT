#ifndef __MPI_WRAPPER_H__
#define __MPI_WRAPPER_H__

#include <mpi.h>
#include <vector>
#include <boost/multi_array.hpp>
#include "Point.hpp"

/* @brief Scatters the data included in a std::vector across all workers in the pool via intermediary arrays
 * @param[in] std::vector of Points to scatter
 * @param[out] Scattered vector-of-Points
 */ 
template <typename Type>
void scatterData(std::vector<Point<Type>>& toCopy, std::vector<Point<Type>>& toReceive)
{
    // First, get the MPI ID and worker sizes for this pool
    int status, rank, poolSize;
    status = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    poolSize = MPI_Comm_size(MPI_COMM_WORLD, &poolSize);
    MPI_Status mpiStatus;

    // Now, if we're rank zero, we need to extract the data from the only populated array we have
    if (rank == 0)
    {
        // Get the number of entries in toCopy
        size_t vectorSize = toCopy.size();
        // Determine how many values each worker is going to have to take, plus some remainder
        unsigned long valuesPerWorker = vectorSize / poolSize; // Integer division
        unsigned long remainder = vectorSize - valuesPerWorker * poolSize; // Remainder
        unsigned long valuesForMe; // My values

        // The first (remainder-th)-cores are going to have valuesPerWorker+1 of entries
        for (size_t rank = 1; rank < remainder; rank++)
        {
            status = MPI_Send(valuesPerWorker+1, 1, MPI_INTEGER, rank, 0, MPI_COMM_WORLD);
        }
        for (size_t rank = remainder; rank < poolSize; rank++)
        {
            status = MPI_Send(valuesPerWorker, 1, MPI_INTEGER, rank, 0, MPI_COMM_WORLD);
        }

        // Now determine the size required for this array
        if (remainder != 0) valuesForMe = valuesPerWorker + 1;

        // Statically allocate size arrays of size vectorSize to store the intermediary points
        double *x = new double[vectorSize];
        double *y = new double[vectorSize];
        double *z = new double[vectorSize];
        double *xdot = new double [vectorSize];
        double *ydot = new double [vectorSize];
        double *zdot = new double [vectorSize];

        // Assign values into these vectors
        for (size_t idx = 0; idx < vectorSize; ++idx)
        {
            Point<Type> tmp = toCopy[idx];
            x[idx] = tmp.state[0];
            y[idx] = tmp.state[1];
            z[idx] = tmp.state[2];
            xdot[idx] = tmp.state[3];
            ydot[idx] = tmp.state[4];
            zdot[idx] = tmp.state[5];
        }

        // Start sending values. Send the first chunk to all processes...
        for (size_t chunkNum = 0; chunkNum < (valuesPerWorker * poolSize); chunk+=poolSize)
        {
            for (size_t rankToSend = 0; rankToSend < poolSize; ++rankToSend)
            {
                if (rankToSend == 0) // Sending to myself - we handle this later
                {
                    Point tmp;
                    tmp[0] = x[rankToSend+chunkNum];
                    tmp[1] = y[rankToSend+chunkNum];
                    tmp[2] = z[rankToSend+chunkNum];
                    tmp[3] = xdot[rankToSend+chunkNum];
                    tmp[4] = ydot[rankToSend+chunkNum];
                    tmp[5] = zdot[rankToSend+chunkNum];
                    toReceive.push_back(tmp);
                } else
                {
                    status = MPI_Send(x[rankToSend+chunkNum], 1, MPI_DOUBLE, rankToSend, 0, MPI_COMM_WORLD);
                    status = MPI_Send(y[rankToSend+chunkNum], 1, MPI_DOUBLE, rankToSend, 0, MPI_COMM_WORLD);
                    status = MPI_Send(z[rankToSend+chunkNum], 1, MPI_DOUBLE, rankToSend, 0, MPI_COMM_WORLD);
                    status = MPI_Send(xdot[rankToSend+chunkNum], 1, MPI_DOUBLE, rankToSend, 0, MPI_COMM_WORLD);
                    status = MPI_Send(ydot[rankToSend+chunkNum], 1, MPI_DOUBLE, rankToSend, 0, MPI_COMM_WORLD);
                    status = MPI_Send(zdot[rankToSend+chunkNum], 1, MPI_DOUBLE, rankToSend, 0, MPI_COMM_WORLD);
                }
            }
        }
        // Now handle the remainder
        for (size_t rankToSend = 0; rankToSend < remainder; ++rankToSend)
        {
            if (rankToSend == 0) // Sending to myself - we handle this later
            {
                Point tmp;
                tmp[0] = x[rankToSend+chunkNum];
                tmp[1] = y[rankToSend+chunkNum];
                tmp[2] = z[rankToSend+chunkNum];
                tmp[3] = xdot[rankToSend+chunkNum];
                tmp[4] = ydot[rankToSend+chunkNum];
                tmp[5] = zdot[rankToSend+chunkNum];
                toReceive.push_back(tmp);
            } else
            {
                status = MPI_Send(x[rankToSend+chunkNum], 1, MPI_DOUBLE, rankToSend, 0, MPI_COMM_WORLD);
                status = MPI_Send(y[rankToSend+chunkNum], 1, MPI_DOUBLE, rankToSend, 0, MPI_COMM_WORLD);
                status = MPI_Send(z[rankToSend+chunkNum], 1, MPI_DOUBLE, rankToSend, 0, MPI_COMM_WORLD);
                status = MPI_Send(xdot[rankToSend+chunkNum], 1, MPI_DOUBLE, rankToSend, 0, MPI_COMM_WORLD);
                status = MPI_Send(ydot[rankToSend+chunkNum], 1, MPI_DOUBLE, rankToSend, 0, MPI_COMM_WORLD);
                status = MPI_Send(zdot[rankToSend+chunkNum], 1, MPI_DOUBLE, rankToSend, 0, MPI_COMM_WORLD);
            }
        }

    } else
    {
        // Receive the number of values this worker is going to have to do
        int valuesForMe;
        status = MPI_Recv(valuesForMe, 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, &mpiStatus);
        toReceive.reserve(valuesForMe);

        // Receive the sent values
        for (size_t idx = 0; idx < valuesForMe; ++idx)
        {
            double x, y, z, xdot, ydot, zdot;
            status = MPI_Recv(x, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &mpiStatus);
            status = MPI_Recv(y, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &mpiStatus);
            status = MPI_Recv(z, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &mpiStatus);
            status = MPI_Recv(xdot, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &mpiStatus);
            status = MPI_Recv(ydot, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &mpiStatus);
            status = MPI_Recv(zdot, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &mpiStatus);

            Point<Type> tmp;
            tmp.state[0] = x; tmp.state[1] = y; tmp.state[2] = z; 
            tmp.state[3] = xdot; tmp.state[4] = ydot; tmp.state[5] = zdot;
            toReceive[idx] = tmp;
        }
    }    
}

#endif
