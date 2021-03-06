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
        // #pragma omp parallel for shared(Qbe, QbeDeriv)
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
                ve = Qbe * vb; 
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
            status = stepper.try_step(forceFunction, x, t, dt); // Force function defined in integration.hpp
        }
    }

    /* @brief Computes the status of a single trajectory. For stabNumber > 0, integration is in forward time. For stabNumber < 0, backward.
    * @param[in] stabNumber The number of orbits up to (and including) which to test
    * @param[in] startingPoint C++ vector containing the initial conditions of the point
    * @param[in] targetCondition Integer status code for which the final behaviour is desired.
    * @returns std::pair<int, int> of the final number of orbits completed, and the status on the n-th revolution
    */
    template <typename integerType, typename vectorType>
    std::pair<int, int> getStatus(const integerType stabNumber, std::vector<vectorType>& startingPoint, bool flag)
    {
        /* Conditions in EMEField need to be non-dimensionalised for the integration to take place correctly.
            * The normalisation prevents errors occurring when the motion is close to the target body.
            *
            * startingPoint contains the original position before *any* integration. nonDimInitialCondition
            * stores the starting point of the particle on a given orbit (so if n = 1 revolutions occur, it is the initial condition
            * after one revolution. currentPoint is the...current position.
            */
        std::vector<vectorType> nonDimStartingPoint, nonDimCurrentPoint, nonDimInitialCondition;

        getNonDimState(startingPoint, nonDimStartingPoint);
        nonDimCurrentPoint = nonDimStartingPoint; // On the first orbit, all three are equal
        nonDimInitialCondition = nonDimStartingPoint;
        std::ofstream output, planeintersections;
        std::cout << nonDimStartingPoint[0] << " " << nonDimStartingPoint[1] << " " << nonDimStartingPoint[2] << " " << nonDimStartingPoint[3] << " " << nonDimStartingPoint[4] << " " << nonDimStartingPoint[5] << '\n';
        output.open("../testIntegrationStates");
        planeintersections.open("../intersections");
        // std::cout << "starting point " << nonDimCurrentPoint[0] << " " << nonDimCurrentPoint[1] << " " << nonDimCurrentPoint[2] << '\n';    

        /* Initialise output result */
        std::pair<int, int> result = std::make_pair(-1, -1); // Initialise to "impossible" number to check assignment later
        int n = 0;
        double currentIntegrationTime = 0.0;
        while ( std::abs(n) < std::abs(stabNumber) )
        {
            // std::cout << n << std::endl;
            /* Initialise integration stepper */
            boost::numeric::odeint::runge_kutta_fehlberg78<std::vector<vectorType>> method;
            auto stepper = boost::numeric::odeint::make_controlled(1.e-014, 1.e-014, method); // relTol, absTol, method
            int integrationDirection = sgn(stabNumber); // Returns (1, 0, -1) depending on sign
            double dt = 0.01 * integrationDirection;
            int integrationStatus = 0;

            /* Make first step to force update */
            makeStep(stepper, nonDimCurrentPoint, currentIntegrationTime, dt);
            double previousConditionOne = getConditionOneIntrinsic(nonDimCurrentPoint, nonDimStartingPoint);    

            /* Begin running through revolutions */
            while (integrationStatus == 0)
            {
                /* Make a step on the trajectory */
                makeStep(stepper, nonDimCurrentPoint, currentIntegrationTime, dt);
                /* Obtain the status of the trajectory.
                    * 0 => nothing to report, continue
                    * 1 => crashed
                    * 2 => escaped                                     
                    * 3 => _potentially_ a full revolution
                    * 4 => acrobatic
                    */
                integrationStatus = integrationController(n, nonDimCurrentPoint, nonDimInitialCondition, nonDimStartingPoint,
                                                          currentIntegrationTime, previousConditionOne);
                output << nonDimCurrentPoint[0] << " " << nonDimCurrentPoint[1] << " " << nonDimCurrentPoint[2] << std::endl;
            }
            if ( integrationStatus != 3)
            {
                result.first = n; result.second = integrationStatus; // How far it got, what it failed with
                break;
            } else 
            {
                if ( integrationDirection > 0 )
                {
                    /* If it went round on its orbit, we're going to continue the integration from here. We need to find the
                        * exact (ish) point it completed its revolution, and reset the definitions of currentPoint etc. accordingly
                        */
                    obtainZero(nonDimStartingPoint, nonDimCurrentPoint, currentIntegrationTime);
                    Eigen::Matrix<double, 3, 1> r, r0, v0, v, rkm1, vkm1;
                    for (size_t idx = 0; idx < 3; ++idx)
                    {
                        r(idx) = nonDimCurrentPoint[idx];
                        r0(idx) = nonDimStartingPoint[idx];
                        v(idx) = nonDimCurrentPoint[idx+3];
                        v0(idx) = nonDimStartingPoint[idx+3];
                        rkm1(idx) = nonDimInitialCondition[idx];
                        vkm1(idx) = nonDimInitialCondition[idx+3];
                    }
                    planeintersections << r << '\n' << v << std::endl;
                    std::cout << "Condition One is " << getConditionOneIntrinsic(nonDimCurrentPoint, nonDimStartingPoint) << " Condition Two and Three are " << r.dot(r0) << " " << v.dot(v0) * vkm1.dot(v0) << std::endl;
                    if ( getConditionTwoAndThree(nonDimCurrentPoint, nonDimInitialCondition, nonDimStartingPoint) ) n += integrationDirection;
                    if ( r.dot(r0) > 0) nonDimInitialCondition = nonDimCurrentPoint;
                } else
                {
                    n += integrationDirection;
                }
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

    /* @brief Initialise the mapping used to store set indices.
     * @param[in] nForward Number of forward revolutions.
     * @param[in] nBackward Number of backward revolutions.
     * @param[inout] unordered_map<int, vector<vector<pair<int>>>>
     */
    void initialiseSetMap(const int nForward, const int nBackward, std::unordered_map<int, std::vector<std::vector<std::pair<int, int>>>>& sets)
    {
        for (int n = nBackward; n <= nForward; n++)
        {
            std::vector<std::vector<std::pair<int, int>>> toAssign(4, std::vector<std::pair<int, int>>());
            sets[n] = toAssign;
        }
    }

    /* @brief Obtain the sets of solutions from a given field.
    * @param[in] stabNum Stability number for which the sets should be obtained
    * @param[inout] field The initial conditions field over which the sets should be obtained
    * @param[inout] sets unordered_map, mapping orbit numbers to a vector of vectors containing the Points in the given set
    */
    template <typename fieldType>
    void getSets(const int stabNum, ACROBAT::emeField<fieldType>& field,
     		std::unordered_map<int, std::vector<std::vector<std::pair<int, int>>>>& forwardSet,
            std::unordered_map<int, std::vector<std::vector<std::pair<int, int>>>>& backwardSet)
    {
     	/* Going to iterate through the domain using all the workers in the MPI pool, so get the
     	 * size of the MPI pool and the position of each worker in the pool here
     	 */
     	int myRank, poolSize;
     	getWorkerInfo(myRank, poolSize); // Defined in mpiWrappers.hpp

     	/* How big is the domain we're studying? */
     	int domainExtentX = field.getXExtent(), domainExtentY = field.getYExtent();
     	int progressCounter = 0;
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
                Point<double> tmpPoint = field.getValue(i, j);
     			std::vector<double> startingPoint = tmpPoint.state; // emeField.getValue returns Point<fieldType>
     			/* Integrate it forward and get the status */
     			std::pair<int, int> forwardStatus = getStatus(stabNum, startingPoint, j == 19);
                std::cout << "Forward status: " << forwardStatus.first << " " << forwardStatus.second << '\n';
     			/* Integrate it backward and get the status */
     			std::pair<int, int> backwardStatus = getStatus(-1, startingPoint, j == 19);
     			/* Add the relevant positions to the final set mapping
     			 * The first value in the pair goes up to the final revolution number;
     			 * the second is the status at the final revolution.
     			 *
     			 * E.g. <0, 1> => on the zero-th revolution, the particle crashed.
     			 * 		<2, 3> => on the second revolution, the particle crashed. It successfully orbited for n = 2.
     			 */
     			std::pair<int, int> pairToPush = std::make_pair(i, j);
     			for (int n  = 0; n < forwardStatus.first; ++n)
     			{
     				forwardSet[n][2].push_back(pairToPush); // Set -> vector -> vector -> pair
     			}
                for (int n = 0; n < backwardStatus.first; ++n)
                {
                    backwardSet[n][2].push_back(pairToPush);
                }
     			forwardSet[forwardStatus.first][forwardStatus.second-1].push_back(pairToPush); // Set -> vector -> vector -> pair
     			backwardSet[backwardStatus.first][backwardStatus.second-1].push_back(pairToPush); // Set -> vector -> vector -> pair
            }
     		if (myRank == 0) std::cout << "Finished integrating " << static_cast<double>(i) / domainExtentX * 100 << "% of trajectories." << "\n";
     	}
    }

    /* @brief Outputs the set mapping generated previously
     * @param[in] stabNumber Desired stability number
     * @param[in] points vector<vector<type>> containing ICs to write out
     */
    template <typename pointType>
    void writeSetMapping(const int stabNumber, std::vector<std::vector<pointType>>& points)
    {
        std::ofstream output; 
        output.open("../stableSet_n=" + std::to_string(stabNumber));

        for (int i = 0; i < points.size(); ++i)
        {
            for (int j = 0; j < points[0].size(); ++j)
            {
                output << points[i][j] << ", "; 
            }
            output << '\n';
        }
        output.close();
    }

    /* @brief Extract the stable set for a given stability number. The stable set is formed of n forward, -1 backward.
    * @param[in] setMapping Contains the set statistics for all points tested.
    * @param[in] stabNumber The desired forward stability number
    * @param[in] field Desired field for which to extract the initial conditions.
    * @param[out] points vector<point<type>> containing the initial conditions for the set in the desired frame.
    */
    template <typename pointType, typename fieldType>
    void extractStableSet(const int stabNumber,
                        std::unordered_map<int, std::vector<std::vector<std::pair<int, int>>>>& forwardSet,
                        std::unordered_map<int, std::vector<std::vector<std::pair<int, int>>>>& backwardSet,
                        fieldType& field, std::vector<std::vector<pointType>>& points)
    {
        /* First, extract the n-forward points */
        std::vector<std::pair<int, int> > nForward = forwardSet[stabNumber][2];
        /* Now the -1-backwards */
        std::vector<std::pair<int, int> > nBackward = backwardSet[0][1];
        int progressCounter = 0, myRank, poolSize;
        getWorkerInfo(myRank, poolSize);
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
            progressCounter++;
            if (myRank == 0)
            {
                // std::cout << "Completed extracting " << static_cast<double>(progressCounter) / nForward.size() * 100 << "%." << "\n";
            }
        }
    }

    /* Write all the sets associated with a field */
    template <typename fieldType>
    void writeAllSets(const int stabNum, std::unordered_map<int, std::vector<std::vector<std::pair<int, int>>>>& forwardSet,
                      std::unordered_map<int, std::vector<std::vector<std::pair<int, int>>>>& backwardSet,
                      fieldType& field)
    {
        std::unordered_map<int, std::string> setNames = { {0, "Crashed"},
                                                          {1, "Escaped"},
                                                          {2, "Stable"},
                                                          {3, "Acrobatic"} };
        for (unsigned i = 0; i <= stabNum; ++i)
        {
            for (unsigned statusNum = 0; statusNum <= 3; ++statusNum)
            {
                std::ofstream output;
                output.open("../results/"+setNames[statusNum] + "_" + std::to_string(i));
                std::vector<std::pair<int, int>> theseConditions = forwardSet[i][statusNum];
                for (std::pair<int, int>& pair : theseConditions)
                {
                    Point<double> tmpPoint = field.getValue(pair.first, pair.second);
                    for (size_t dimension = 0; dimension < tmpPoint.state.size(); ++dimension)
                    {
                        output << tmpPoint.state[dimension] << ",";
                    }
                    output << '\n';
                }
                output.close();
            }
        }
        for (unsigned i = 0; i <= 1; ++i)
        {
            for (unsigned statusNum = 0; statusNum <= 3; ++statusNum)
            {
                std::ofstream output;
                output.open("../results/backward"+setNames[statusNum] + "_" + std::to_string(i));
                std::vector<std::pair<int, int>> theseConditions = backwardSet[i][statusNum];
                for (std::pair<int, int>& pair : theseConditions)
                {
                    Point<double> tmpPoint = field.getValue(pair.first, pair.second);
                    for (size_t dimension = 0; dimension < tmpPoint.state.size(); ++dimension)
                    {
                        output << tmpPoint.state[dimension] << ",";
                    }
                    output << '\n';
                }
                output.close();
            }
        }
    }

}
#endif
