#ifndef __EME_FIELD_H__
#define __EME_FIELD_H__

#include "field.hpp"
#include <ostream>
#include "Params.hpp"
#include "OEs.hpp"
#include "integration.hpp"
#include "bmeField.hpp"
#include <unordered_map>
#include "RADEC.hpp"
#include <ostream>
#include <algorithm>
#include "mpiWrappers.hpp"
#include <string>
#include "rootFinding.hpp"
#include <mpi.h>
#include "coordinateTransforms.hpp"

extern "C"
{
    #include "SpiceUsr.h"
}


/* @brief Returns the sign of a given numeric input
@params[in] Numeric type for which the sign is desired, passed by reference
@returns An integer giving the sign of the expression ( {1, 0, -1} )
*/
template <typename numericType>
int sgn(numericType val)
{
    return (val > 0) - (val < 0);
}

template <typename pointType>
void broadcastMap(unsigned int stabNum, std::unordered_map<int, std::vector<Point<pointType>>>&map)
{
    int mpiReturn;
    int rank, size;
    mpiReturn = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    mpiReturn = MPI_Comm_size(MPI_COMM_WORLD, &size);
    for (int revNumber = 1; revNumber <= stabNum; ++revNumber)
    {
        std::vector<Point<pointType>> thisValue = map[revNumber];
        std::vector<Point<pointType>> toBroadcast;
        reduceVector(thisValue, toBroadcast, rank, size);
        broadcastVector(toBroadcast, rank);
        map[revNumber] = toBroadcast;
    }
    int revNumber = -1;
    std::vector<Point<pointType>> thisValue = map[revNumber];
    std::vector<Point<pointType>> toBroadcast;
    reduceVector(thisValue, toBroadcast, rank, size);
    broadcastVector(toBroadcast, rank);
    map[revNumber] = toBroadcast;
}


namespace ACROBAT
{
    template <class Type>
    class emeField : public field2D<Type>
    {
        public:
            emeField(int nx, int ny) : field2D<Type>(nx, ny)
            {};

            /* @brief Initialises the elements in the field using a BME field.
               @params[in] ACROBAT::bmeField containing the points defining the domain
            */
            void initialiseField(ACROBAT::bmeField<Type> &input)
            {
                BMEtoEME(input, *this);
            }

            /* @brief Returns the final time for the trajectory integration
               @returns The final time for the integration.
            */
            double getFinalTime() const
            {
                return this->_finalTime;
            }

            /* @brief Returns the initial time for the trajectory integration
               @returns The initial time for the integration
            */
            double getInitialTime() const
            {
                return this->_initialTime;
            }

            /* @brief Sets the final time for the trajectory integration
               @param[in] finalTime The final time for the trajectory integration
            */
            void setFinalTime(const double finalTime)
            {
                this->_finalTime = finalTime;
            }

            /* @brief Sets the initial time for the trajectory integration
               @param[in] initialTime The initial time for the trajectory integration
            */
            void setInitialTime(const double initialTime)
            {
                this->_initialTime = initialTime;
            }

            /* @brief Writes the results of the set generation to an output file
               @param[in] stabNumber The number of revolutions a particle must complete.
               @param[in] map Map of Point<int> indexed by the number of revolutions
            */
            void outputWrite(int& stabNumber, std::unordered_map<int, std::vector<Point<int>>>& map)
            {
                // Define n=-1 set
                std::vector<Point<int>> backwardsSet = map[-1];
                // Define the n-th set
                std::vector<Point<int>> nthSet = map[stabNumber];

                // Open output file
                std::string fname = "stableSet_n=" + std::to_string(stabNumber);
                std::ofstream output;
                output.open(fname);

                // The final set is those in nthSet that are also in backwardsSet
                for(auto val : nthSet)
                {
                    if (std::find(backwardsSet.begin(), backwardsSet.end(), val) != backwardsSet.end()) // val is in backwardsSet
                    {
                        int idx1, idx2;
                        idx1 = val.state[0];
                        idx2 = val.state[1];
                        Point<double> initialCondition = this->getValue(idx1, idx2);
                        double x = initialCondition.state[0];
                        double y = initialCondition.state[1];
                        double z = initialCondition.state[2];
                        output << x << "," << y << "," << z << std::endl;
                    }
                }
                output.close();
            }

            /* @brief Writes the results of the set generation to an output file
               @param[in] stabNumber The number of revolutions a particle must complete.
               @param[in] map Map of Point<int> indexed by the number of revolutions
            */
            void outputWriteBMEFrame(int& stabNumber, std::unordered_map<int, std::vector<Point<int>>>& map, ACROBAT::bmeField<Point<double>>& field)
            {
                // Define n=-1 set
                std::vector<Point<int>> backwardsSet = map[-1];
                // Define the n-th set
                std::vector<Point<int>> nthSet = map[stabNumber];

                // Open output file
                std::string fname = "BMEField_stableSet_n=" + std::to_string(stabNumber);
                std::ofstream output;
                output.open(fname);

                // The final set is those in nthSet that are also in backwardsSet
                for(auto val : nthSet)
                {
                    if (std::find(backwardsSet.begin(), backwardsSet.end(), val) != backwardsSet.end()) // val is in backwardsSet
                    {
                        int idx1, idx2;
                        idx1 = val.state[0];
                        idx2 = val.state[1];
                        Point<double> initialCondition = field.getValue(idx1, idx2);
                        double x = initialCondition.state[0];
                        double y = initialCondition.state[1];
                        double z = initialCondition.state[2];
                        output << x << "," << y << "," << z << std::endl;
                    }
                }
                output.close();
            }

            /* @brief Writes the results of the set generation to an output file - writes r0 and w0
               @param[in] stabNumber The number of revolutions a particle must complete.
               @param[in] map Map of Point<int> indexed by the number of revolutions
            */
            void outputElementsWrite(int& stabNumber, ACROBAT::oeField& field, std::unordered_map<int, std::vector<Point<int>>>& map)
            {
                // Define n=-1 set
                std::vector<Point<int>> backwardsSet = map[-1];
                // Define the n-th set
                std::vector<Point<int>> nthSet = map[stabNumber];

                // Open output file
                std::string fname = "stableSet_n=" + std::to_string(stabNumber);
                std::ofstream output;
                output.open(fname);

                // The final set is those in nthSet that are also in backwardsSet
                for(auto val : nthSet)
                {
                    if (std::find(backwardsSet.begin(), backwardsSet.end(), val) != backwardsSet.end()) // val is in backwardsSet
                    {
                        int idx1, idx2;
                        idx1 = val.state[0];
                        idx2 = val.state[1];
                        ACROBAT::OEs initialCondition = field.getValue(idx1, idx2);
                        double r0 = initialCondition.rp;
                        double omega = initialCondition.omega;
                        output << r0 << "," << omega << std::endl;
                    }
                }
                output.close();
            }
        private:
            double _finalTime;      // Final time for the trajectory integration
            double _initialTime;    // Initial time for the trajectory integration
    };

    /* @brief Prints the statistics for a given ballistic set computation
        @param[in] stabNum The number of the sability index for the current computation
        @param[in] setStatistics std::vector<> containing the statistics defined as in the set determination routines.
    */
    template <typename integerType, typename vectorType>
    void printStatistics(const integerType stabNum, const std::vector<vectorType>& setStatistics)
    {
        std::cout << "For a stability index of " << stabNum << "the statistics are as follows:" << std::endl;
        std::cout << "\t Number of crashes:   " << setStatistics[0] << std::endl;
        std::cout << "\t Number of escapes:   " << setStatistics[1] << std::endl;
        std::cout << "\t Number of stable:    " << setStatistics[2] << std::endl;
        std::cout << "\t Number of acrobatic: " << setStatistics[3] << std::endl;
    }

    /*  @brief Converts a BME@Epoch field to an EME2000 field.
        @param[in] bmeField: bmeField to convert
        @param[out] emeField: emeField to store the conversion in.
    */
    template <typename Type>
    void BMEtoEME(ACROBAT::bmeField<Point<Type>> &bmeField, ACROBAT::emeField<Point<Type>> &emeField)
    {
        Eigen::Matrix<Type, 3, 3> Qbe;
        Eigen::Matrix<Type, 3, 3> QbeDeriv;

        // Create the transformation matrix
        getBMEtoEMERotationMatrix(PARAMS::EPOCH, Qbe, QbeDeriv);

        // Apply Qbe to every state in bmeField
        #pragma omp parallel for shared(Qbe, QbeDeriv)
        for (unsigned int i = 0; i < bmeField.getXExtent(); ++i)
        {
            for (unsigned int j = 0; j < bmeField.getYExtent(); ++j)
            {
                // Set up input vector
                Eigen::Matrix<Type, 3, 1> xb, xe, vb, ve;
                Point<Type> temp = bmeField.getValue(i, j);
                
                // Assign
                for(unsigned idx = 0; idx < 3; ++idx) 
                {
                    xb(idx) = temp[idx];
                    vb(idx) = temp[idx+3];
                }

                // Compute
                xe = Qbe * xb;
                ve = QbeDeriv * xb + Qbe * vb;

                // Swap back
                for (unsigned idx = 0; idx < 3; ++idx) 
                {
                    temp[idx] = xe(idx);
                    temp[idx+3] = ve(idx);
                }

                // Assign
                emeField.setValue(temp, i, j);
            }
        }
    }

    /*  @brief Converts an EME2000 field to a BME@Epoch field
        @param[in] emeField: emeField to convert
        @param[out] bmeField: bmeField to store the conversion in.
    */
    template <typename Type>
    void EMEtoBME(const ACROBAT::emeField<Point<Type>> &emeField, ACROBAT::bmeField<Point<Type>> &bmeField)
    {
        Eigen::Matrix<Type, 6, 6> Qeb;

        // Create the transformation matrix
        getEMEtoBMERotationMatrix(PARAMS::EPOCH, Qeb);

        // Apply Qbe to every state in bmeField
        #pragma omp parallel for shared(Qbe)
        for (unsigned int i = 0; i < bmeField.getXExtent(); ++i)
        {
            for (unsigned int j = 0; j < bmeField.getYExtent(); ++j)
            {
                // Set up input vector
                Eigen::Matrix<Type, 6, 1> xb, xe;
                Point<Type> temp = emeField.getValue(i, j);
                
                // Assign
                for(unsigned idx = 0; idx < 6; ++idx) xe(i) = temp[i];

                // Compute
                xb = Qbe * xe;

                // Swap back//    if (rank == 0){
                //        emeDomain.outputWrite(n, results);
                //        std::cout << "Completed obtaining the stable sets. The time required was " << std::chrono::duration_cast<std::chrono::seconds>(end-start).count() << " seconds.";
                //    }
                for (unsigned idx = 0; idx < 6; ++idx) temp[i] = xb(i);

                // Assign
                emeField.setValue(temp, i, j);
            }
        }
    }

    /* @brief Performs a single step of the ODE integration.
     * @param[in] Stepper Boost::numeric::odeint::stepper to use.
     * @param[inout] x On entry: point to advect forward. On exit: advected point.
     * @param[inout] t On entry: previous time step. On exit: new time step.
     * @param[inout] dt On entry: guess for the time step needed. On exit: time step taken.
     */
     template <typename stepperType, typename vectorType>
     void makeStep(stepperType& stepper, std::vector<vectorType>& x, vectorType& t, vectorType& dt)
     {
    	 boost::numeric::odeint::controlled_step_result status = boost::numeric::odeint::fail;

    	 /* Return to the calling function only when the step was successful. */
    	 while (status == boost::numeric::odeint::fail)
    	 {
    		 status = stepper.try_step(forceFunction, x, currentTime, dt); // Force function defined in integration.hpp
    	 }
     }


     /* @brief Obtain the sets of solutions from a given field.
      * @param[in] stabNum Stability number for which the sets should be obtained
      * @param[inout] field The initial conditions field over which the sets should be obtained
      * @param[inout] sets unordered_map, mapping orbit numbers to a vector of vectors containing the Points in the given set
      */
     template <typename integerType, typename fieldType, typename vectorType>
     void getSets(const integerType& stabNum, ACROBAT::emeField<fieldType>& field,
     		std::unordered_map<int, std::vector<std::vector<std::pair<int, int>>>>& sets)
     {
     	/* Going to iterate through the domain using all the workers in the MPI pool, so get the
     	 * size of the MPI pool and the position of each worker in the pool here
     	 */
     	int myRank, poolSize;
     	getWorkerInfo(myRank, poolSize); // Defined in mpiWrappers.hpp

     	/* How big is the domain we're studying? */
     	int domainExtentX = field.getXExtent(), domainExtentY = field.getYExtent();

     	/* Now, iterate through the domain.
     	 * At each point, we call the getStatus function to obtain the its status for the n-th revolution. The result is returned
     	 * as a pair: (number of revolutions, status on n-th revolutions). This allows us to populate the relevant vectors in the mapping
     	 * for its starting point.
     	 */
     	for (long i = myRank; i < domainExtentX; i += poolSize)
     	{
     		for (long j = 0; j < domainExtentY; ++j)
     		{
     			/* Extract this starting point from the domain */
     			std::vector<fieldType> startingPoint = emeField.getValue(i, j).state; // emeField.getValue returns Point<fieldType>
     			/* Integrate it forward and get the status */
     			std::pair<int, int> forwardStatus = getStatus(stabNumber, startingPoint);
     			/* Integrate it backward and get the status */
     			std::pair<int, int> backwardStatus = getStatus(-1, startingPoint);
     			/* Add the relevant positions to the final set mapping
     			 * The first value in the pair goes up to the final revolution number;
     			 * the second is the status at the final revolution.
     			 *
     			 * E.g. <0, 1> => on the zero-th revolution, the particle crashed.
     			 * 		<2, 3> => on the second revolution, the particle crashed. It successfully orbited for n = 2.
     			 */
     			std::pair pairToPush = std::make_pair(i, j);
     			for (int n  = 0; n < forwardStatus.first; ++n)
     			{
     				sets[n][2].push_back(pairToPush); // Set -> vector -> vector -> pair
     			}
     			sets[forwardStatus.first][forwardStatus.second].push_back(pairToPush); // Set -> vector -> vector -> pair
     			sets[backwardStatus.first][backwardStatus.second].push_back(pairToPush); // Set -> vector -> vector -> pair
     		}
     	}
     	/* Send mapping to all the workers */
     	broadcastMapping(sets);
     	/* Need to write a function to broadcast each value of the mapping from all cores to the worker and back again */
     }

     /* @brief Extract the stable set for a given stability number. The stable set is formed of n forward, -1 backward.
      * @param[in] setMapping Contains the set statistics for all points tested.
      * @param[in] stabNumber The desired forward stability number
      * @param[in] field Desired field for which to extract the initial conditions.
      * @param[out] points vector<point<type>> containing the initial conditions for the set in the desired frame.
      */
      template <typename pointType, typename fieldType>
      void extractStableSet(const int stabNumber,
     		               const std::unordered_map<int, std::vector<std::vector<std::pair<int, int>>>>& setMapping,
     		 	 	 	   const fieldType& field, std::vector<std::vector<pointType>>& points)
      {
     	 /* First, extract the n-forward points */
     	 std::vector< std::pair<int, int> > nForward = setMapping[n][2]; // n-Forward, full revs
     	 /* Now the -1-backwards */
     	 std::vector< std::pair<int, int> > nBackward = setMapping[-1][1]; // n-Backward, escape
     	 /* Add initial conditions to points */
     	 points.clear();
     	 for (std::pair<int, int>& forwardPoint : nForward)
     	 {
     		 /* Check if this point in nBackward */
     		 for (std::pair<int, int>& backwardPoint : nBackward)
     		 {
     			 if (forwardPoint == backwardPoint)
     			 {
     				 Point<pointType> tmpPoint = field.getValue(forwardPoint.first, forwardPoint.second);
     				 std::vector<pointType> tmpVector = tmpPoint.state; // Underlying vector storage
     				 points.push_back(tmpVector);
     			 }
     		 }
     	 }
      }

      /* @brief Computes the status of a single trajectory. For stabNumber > 0, integration is in forward time. For stabNumber < 0, backward.
       * @param[in] stabNumber The number of orbits up to (and including) which to test
       * @param[in] startingPoint C++ vector containing the initial conditions of the point
       * @param[in] targetCondition Integer status code for which the final behaviour is desired.
       * @returns std::pair<int, int> of the final number of orbits completed, and the status on the n-th revolution
       */
      template <typename integerType, typename vectorType>
      std::pair<int, int> getStatus(const integerType stabNumber, const std::vector<vectorType>& startingPoint)
      {
      	/* Conditions in EMEField need to be non-dimensionalised for the integration to take place correctly.
      	 * The normalisation prevents errors occurring when the motion is close to the target body.
      	 *
      	 * startingPoint contains the original position before *any* integration. nonDimInitialCondition
      	 * stores the starting point of the particle on a given orbit (so if n = 1 revolutions occur, it is the initial condition
      	 * after one revolution. currentPoint is the...current position.
      	 */
      	std:vector<vectorType> nonDimStartingPoint(6), nonDimCurrentPoint(6), nonDimInitialCondition(6);

      	getNonDimState(startingPoint, nonDimStartingPoint);
      	nonDimCurrentPoint = nonDimStartingPoint; // On the first orbit, all three are equal
      	nonDimInitialCondition = nonDimStartingPoint;

      	/* Initialise integration stepper */
      	boost::numeric::odeint::runge_kutta_fehlberg78<std::vector<vectorType>> method;
      	auto stepper = boost::numeric::odeint::make_controlled(1.e-012, 1.e-012, method); // relTol, absTol, method
      	int integrationDirection = sgn(stabNumber); // Returns (1, 0, -1) depending on sign
      	double dt = 0.01 * integrationDirection;

      	/* Initialise output result */
      	std::pair<int, int> result = std::make_pair(-1, -1); // Initialise to "impossible" number to check assignment later

      	int n;
      	for(n = 0; std::abs(n) < std::abs(stabnumber); n += integrationDirection) // Abs and integrationDirection allow for -ve or +ve time
      	{
      		int integrationStatus = 0;
      		double previousConditionOne = 1.0 * integrationDirection;
      		/* Begin running through revolutions */
      		double currentIntegrationTime = 0.0;

      		while (integrationStatus == 0)
      		{
      			/* Make a step on the trajectory */
      			make_step(stepper, nonDimCurrentPoint, currentTime, dt);
      			/* Obtain the status of the trajectory.
      			 * 0 => nothing to report, continue
      			 * 1 => crashed
      			 * 2 => escaped
      			 * 3 => full revolution
      			 * 4 => acrobatic
      			 */
      			integrationStatus = integrationController(nonDimCurrentPoint, nonDimStartingPoint, nonDimInitialCondition,
      													  currentIntegrationTime, previousConditionOne);
      		}
      		/* Now, check why the integration failed. If it wasn't because it completed a revolution - stop. */
      		if ( integrationStatus == 3) )
      		{
      			result.first = n; result.second = integrationStatus; // How far it got, what it failed with
      		} else
      		{
      			/* If it went round on its orbit, we're going to continue the integration from here. We need to find the
      			 * exact (ish) point it completed its revolution, and reset the definitions of currentPoint etc. accordingly
      			 */
      			if ( integrationDirection > 0 ) obtainZero(nonDimInitialCondition, nonDimCurrentPoint, currentIntegrationTime);
      		}
      	}
      	/* If we dropped out of the loop only ever orbiting, result won't have been populated. Check that here and override
      	 * if necessary.
      	 */
      	if (result.second < 0)
      	{
      		result.first = n; result.second = 3;
      	}
      	return result;
      }
}
#endif
