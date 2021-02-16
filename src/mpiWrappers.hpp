#ifndef __MPI_WRAPPERS_H__
#define __MPI_WRAPPERS_H__

#include <mpi.h>
#include "Point.hpp"


/* @brief Obtains the rank of a given worker, and the total number of workers in the pool.
 * @param[inout] myRank The rank of the worker calling the function
 * @param[inout] poolSize The total number of workers in the pool
 */
void getWorkerInfo(int& myRank, int& poolSize)
{
	MPI_Status status;
	int returnCode;
	returnCode = MPI_Comm_size(MPI_COMM_WORLD, &poolSize);
	returnCode = MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
}

/* @brief Broadcasts a vector of integer points from the master worker onto all other workers.
   @param[in] vectorToBroadcast Vector, fully defined on rank zero, to be broadcast to all other workers.
   @param[in] rank Rank of this worker in the pool.
*/
void broadcastVector(std::vector<Point<int>>& vectorToBroadcast, int& rank)
{
    int returnCode;
    int numToBroadcast;
    std::vector<Point<int>> tmpVector;
    if (rank == 0) numToBroadcast = vectorToBroadcast.size();
    /* First, send the number of broadcasts we're going to have to do */
    returnCode = MPI_Bcast(&numToBroadcast, 1, MPI_INT, 0, MPI_COMM_WORLD);
    for (size_t idx = 0; idx < numToBroadcast; ++idx)
    {
        Point<int> tmpPoint;
        for (size_t dimension = 0; dimension < 6; ++dimension)
        {
            int tmp;
            if (rank == 0) tmp = vectorToBroadcast[idx].state[dimension];
            returnCode = MPI_Bcast(&tmp, 1, MPI_INT, 0, MPI_COMM_WORLD);
            tmpPoint.state[dimension] = tmp;
        }
        if (rank != 0) tmpVector.push_back(tmpPoint);
    }
    if (rank != 0) vectorToBroadcast = tmpVector;
}

/* @brief Broadcasts a vector of double-precision points from the master worker onto all other workers.
   @param[in] vectorToBroadcast Vector, fully defined on rank zero, to be broadcast to all other workers.
   @param[in] rank Rank of this worker in the pool.
   */
void broadcastVector(std::vector<Point<double>>& vectorToBroadcast, int& rank)
{
    int returnCode;
    int numToBroadcast;
    std::vector<Point<double>> tmpVector;
    if (rank == 0) numToBroadcast = vectorToBroadcast.size();
    /* First, send the number of broadcasts we're going to have to do */
    returnCode = MPI_Bcast(&numToBroadcast, 1, MPI_INT, 0, MPI_COMM_WORLD);
    for (size_t idx = 0; idx < numToBroadcast; ++idx)
    {
        Point<double> tmpPoint;
        for (size_t dimension = 0; dimension < 6; ++dimension)
        {
            double tmp;
            if (rank == 0) tmp = vectorToBroadcast[idx].state[dimension];
            returnCode = MPI_Bcast(&tmp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            tmpPoint.state[dimension] = tmp;
        }
        if (rank != 0) tmpVector.push_back(tmpPoint);
    }
    if (rank != 0) vectorToBroadcast = tmpVector;
}

/* @brief Reduces a std::vector<Point<double>> from all workers onto the master worker (rank zero.)
   @param[in] vectorToReduce Local vector, defined on all workers, to reduce
   @param[out] vectorToStore Target for the vector reduction. Valid only on rank zero.
   @param[in] rank The rank of this worker in the overall pool.
   @param[in] poolSize The number of workers in the pool.
   */
void reduceVector(std::vector<Point<double>>& vectorToReduce, std::vector<Point<double>>& vectorToStore, int& rank, int& poolSize)
{
    MPI_Status mpiStatus;
    int returnCode;

    if (rank == 0)
    {
        int numToReceive;
        // Copy vectorToReduce into vectorToStore on rank 0
        vectorToStore.insert(vectorToStore.end(), vectorToReduce.begin(), vectorToReduce.end());
        for (int worker = 1; worker < poolSize; ++worker)
        {
            returnCode = MPI_Recv(&numToReceive, 1, MPI_INT, worker, worker, MPI_COMM_WORLD, &mpiStatus);
            Point<double> tmpPoint;
            for (int toReceive = 0; toReceive < numToReceive; ++toReceive)
            {
                std::vector<double> tmpState(6);
                for (int dimension = 0; dimension < 6; ++dimension)
                {
                    double tmpValue;
                    returnCode = MPI_Recv(&tmpValue, 1, MPI_DOUBLE, worker, worker, MPI_COMM_WORLD, &mpiStatus);
                    tmpState[dimension] = tmpValue;
                }
                tmpPoint.state = tmpState;
                vectorToStore.push_back(tmpPoint);
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
                double tmpValue = vectorToReduce[toSend].state[dimension];
                returnCode = MPI_Send(&tmpValue, 1, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD);
            }
        }
    }
}

/* @brief Reduces a std::vector<Point<int>> from all workers onto the master worker (rank zero.)
   @param[in] vectorToReduce Local vector, defined on all workers, to reduce
   @param[out] vectorToStore Target for the vector reduction. Valid only on rank zero.
   @param[in] rank The rank of this worker in the overall pool.
   @param[in] poolSize The number of workers in the pool.
*/
void reduceVector(std::vector<Point<int>>& vectorToReduce, std::vector<Point<int>>& vectorToStore, int& rank, int& poolSize)
{
    MPI_Status mpiStatus;
    int returnCode;

    if (rank == 0)
    {
        int numToReceive;
        // Copy vectorToReduce into vectorToStore on rank 0
        vectorToStore.insert(vectorToStore.end(), vectorToReduce.begin(), vectorToReduce.end());
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
                vectorToStore.push_back(tmpPoint);
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

/* @brief Reduces (sums) all values of a given integer onto the master worker (rank zero.)
   @param[in] myValue Local value of the integer to reduce
   @param[out] reducedValue The reduced integer; valid only on rank zero.
*/
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

/* @brief Broadcasts the count of a given variable from the master worker to all other workers in the pool.
   @param[in] countToBroadcast Integer to broadcast
*/
void broadcastCount(int& countToBroadcast)
{
    int returnCode;
    returnCode = MPI_Bcast(&countToBroadcast, 1, MPI_INT, 0, MPI_COMM_WORLD);
}

/* @brief Reduces and broadcasts a vector of initial conditions to each worker.
 * @param[inout] ics. On entry: initial conditions valid on each worker. On exit: initial conditions
 */
void broadcastMapping(std::vector<std::vector<double>>& ics)
{
	MPI_Status mpiStatus;
	int rank, poolSize;
	getWorkerInfo(rank, poolSize);
	int returnCode;
	/* If master, need to receive everything first. */
	if (rank == 0)
	{
        std::cout << "here in the rank zero portion " << std::endl;
		/* Iterate through workers */
		for (int worker = 1; worker < poolSize; ++worker)
		{
            std::cout << "here2 in the rank zero portion " << std::endl;
			int numberOfICsToReceive;
			returnCode = MPI_Recv(&numberOfICsToReceive, 1, MPI_INT, worker, worker, MPI_COMM_WORLD, &mpiStatus);

            std::cout << "here3 in the rank zero portion " << std::endl;
			for (int ICNum = 0; ICNum < numberOfICsToReceive; ++ICNum)
			{
				std::vector<double> thisCondition(6);
				for (int dimension = 0; dimension < thisCondition.size(); ++dimension)
				{
					double tmp;

                    std::cout << "here4 in the rank zero portion " << std::endl;
					returnCode = MPI_Recv(&tmp, 1, MPI_DOUBLE, worker, worker, MPI_COMM_WORLD, &mpiStatus);

                    std::cout << "here5 in the rank zero portion " << std::endl;
					thisCondition[dimension] = tmp;
				}
				ics.push_back(thisCondition);
			}
		}
		/* Now all the data has been amassed on the master worker, send it to all other workers. */
		int numberOfICsToSend = ics.size();

		for (int worker = 1; worker < poolSize; ++worker)
		{
			returnCode = MPI_Send(&numberOfICsToSend, 1, MPI_INT, worker, worker, MPI_COMM_WORLD);
            std::cout << "here6 in the rank zero portion " << std::endl;
			for (int ICNum = 0; ICNum < numberOfICsToSend; ++ICNum)
			{
				std::vector<double> thisIC = ics[ICNum];
				for (int dimension = 0; dimension < thisIC.size(); ++dimension)
				{
					double tmp = thisIC[dimension];
                    std::cout << "here7 in the rank zero portion " << rank << std::endl;
					returnCode = MPI_Send(&tmp, 1, MPI_DOUBLE, worker, worker, MPI_COMM_WORLD);
                    std::cout << "here8 in the rank zero portion " << rank << std::endl;
				}
			}
		}
	} else
	{
        std::cout << "here in the rank n portion " << std::endl;
		/* Send the number of things this worker has to send to the master core */
		int numberOfICsToSend = ics.size();
		returnCode = MPI_Send(&numberOfICsToSend, 1, MPI_INT, 0, rank, MPI_COMM_WORLD);
        std::cout << "here2 in the rank n portion " << std::endl;
		/* Now send the actual data */
		for (int ICNum = 0; ICNum < numberOfICsToSend; ++ICNum)
		{
			std::vector<double> thisCondition = ics[ICNum];
			for (int dimension = 0; dimension < thisCondition.size(); ++dimension)
			{
				double tmp = thisCondition[dimension];
                std::cout << "here3 in the rank n portion " << rank << std::endl;
				returnCode = MPI_Send(&tmp, 1, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD);
                std::cout << "here4 in the rank n portion " << rank << std::endl;
			}
		}
		/* Now the data has been sent, clear our set ready to receive from the master */
		ics.clear();
		/* Receive from master - first get the number to receive */
		int numberOfICsToReceive;
        std::cout << "here5 in the rank n portion " << std::endl;
		returnCode = MPI_Recv(&numberOfICsToReceive, 1, MPI_INT, 0, rank, MPI_COMM_WORLD, &mpiStatus);
        std::cout << "here6 in the rank n portion " << std::endl;
		/* Now receive the conditions */
		for (int ICNum = 0; ICNum < numberOfICsToReceive; ++ICNum)
		{
			std::vector<double> icToReceive(6);
			for (int dimension = 0; dimension < icToReceive.size(); ++dimension)
			{
				double tmp;
                std::cout << "here7 in the rank n portion " << std::endl;
				returnCode = MPI_Recv(&tmp, 1, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD, &mpiStatus);
                std::cout << "here8 in the rank n portion " << std::endl;
				icToReceive[dimension] = tmp;
			}
			ics.push_back(icToReceive);
		}
	}
    std::cout << "here9 in the rank n portion " << std::endl;
	/* Wait for all to complete - not needed in theory, safe in practice. */
	MPI_Barrier(MPI_COMM_WORLD);
    std::cout << "here10 in the rank n portion " << std::endl;
}


#endif
