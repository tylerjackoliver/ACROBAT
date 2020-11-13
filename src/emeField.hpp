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

            /*  @brief Obtains the capture set from a given set of indices corresponding to points in an EME2000 field.
                @param[in] stabNum The stability index the set is sought for
                @param[in] domain The vector of points containing indices for which the sets are to be based off of
                @param[in] field The domain for which the points in domain correspond to
                @param[out] points A vector containing the indices of the corresponding capture set
            */ 
            template <typename integerType, typename fieldType, typename vectorType>
            void getSetFromPoints(const integerType& stabNum, std::vector<Point<vectorType>> &domain, ACROBAT::emeField<fieldType>& field, std::vector<Point<vectorType>> &points, std::vector<Point<double>>& initialConditions)
            {    
                // Statuses - ordered s.t. index related to condition is (status code - 1) - crash, escape, stable, acrobatic
                std::vector<int> setStatistics(4);
                std::vector<Point<vectorType>> privatePoints;
                std::vector<Point<double>> privateCoordinates;

                // Get the integration direction (+ve or -ve time?)
                int directionTime = sgn(stabNum);
                unsigned prog = 0;

                // Iterate through the conditions on the EME2000 field - do jumps of our rank and pool size
                int rank, poolSize, returnStatus;
                MPI_Status mpiStatus;
                returnStatus = MPI_Comm_size(MPI_COMM_WORLD, &poolSize);
                returnStatus = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
                for (unsigned i = rank; i < domain.size(); i += poolSize)
                {
                    // Get an initial position
                    Point<vectorType> indices = domain[i];
                    int idx1 = indices.state[0];
                    int idx2 = indices.state[1];
                    Point<double> currentPoint = initialConditions[i];
                    std::vector<double> stateAtCrossing(6);

                    // Get the status of this points behaviour
                    int status = getStatus(currentPoint, field.getInitialTime(), directionTime, stateAtCrossing);
                    // Increment the correct counter in the setStatistics vector (ordered s.t. index is status-1)
                    setStatistics[status-1]++;
                    // If weakly stable, append its indices to the points vector
                    if ( (directionTime > 0 && status == 3) || (directionTime < 0 && status == 2) )  // Want n-revolution points or backwards escape
                    {
                        Point<int> thisPoint; thisPoint[0] = idx1; thisPoint[1] = idx2;
                        std::vector<double> stateToInsert;
                        getDimState(stateAtCrossing, stateToInsert);
                        Point<double> pointToInsert(stateToInsert);
                        privatePoints.insert(privatePoints.end(), thisPoint);
                        privateCoordinates.insert(privateCoordinates.end(), pointToInsert);
                    }
                    prog += poolSize;
                    if (rank == 0) std::cout << "Completed integration " << prog << " of " << domain.size() << "." << std::endl;
                }
                for (int dimension = 0; dimension < setStatistics.size(); ++dimension)
                {
                    reduceCount(setStatistics[dimension], setStatistics[dimension]);
                    broadcastCount(setStatistics[dimension]);
                }
                initialConditions.clear();
                reduceVector(privatePoints, points, rank, poolSize);
                reduceVector(privateCoordinates, initialConditions, rank, poolSize);
                broadcastVector(initialConditions, rank);
                broadcastVector(points, rank);
                if (rank == 0) printStatistics(stabNum, setStatistics);
            } // function


            /* @brief Obtains the stabNumber-revolution stable set for a given EME2000 field. Note than stabNumber entries in the map are returned.
            * @param[in] stabNumber The number of revolutions a particle should complete in order to be classed as stable.
            * @param[out] std::unordered_map The keys are the number of rotations, and the values are a std::vector<> containing the points comprising the set.
            */
            template <typename integerType, typename pointType>
            void getStableSet(integerType stabNumber, std::unordered_map<int, std::vector<Point<pointType>>>& outIndices)
            {
                std::unordered_map<int, std::vector<Point<double>>> initConds;
                std::vector<Point<pointType>> points;
                std::vector<Point<double>> initialConditions;
                // Get the set defined solely from the initial domain (i.e. 1-stable)
                getSet(1, *this, points, initialConditions);
                outIndices[1] = points;
                initConds[1] = initialConditions;
                // Get the backwards set
                points.clear();
                initialConditions.clear();
                getSet(-1, *this, points, initialConditions);
                outIndices[-1] = points;
                initConds[-1] = initialConditions;
                // Now get the rest 
                for (unsigned revs = 2; revs <= stabNumber; ++revs)  // Must complete at least one orbit about the host planet
                {
                    // Clear points from the previous run
                    points.clear();
                    initialConditions.clear();
                    initialConditions = initConds[revs-1];
                    getSetFromPoints(revs, outIndices[revs-1], *this, points, initialConditions);
                    outIndices[revs] = points;
                    initConds[revs] = initialConditions;
                }
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

    /* @brief Obtains the rotation matrix for transforming from the BME frame to the EME frame at a given epoch
    *  @param[in] epoch The epoch of the transformation in ephemeris seconds past J2000
    *  @param[out] rot The rotation matrix to transform from the BME to the EME frame (6x6)
    */
    template <typename matrixType>
    void getBMEtoEMERotationMatrix(const double &epoch, Eigen::Matrix<matrixType,3,3> &rot)
    {
        // Get the current values of alpha and delta at this time
        double alpha, delta;
        RADEC::getAlphaDelta(epoch, alpha, delta);

        /* Construct the matrix. In theory, we would need to also construct the derivative 
           of the rotation matrix in order to get the velocity rotated into the frame, but
           since the derivative of the rotation matrix is of magnitude 1e-10...we ignore it.
        */
        double sina = std::sin(alpha);
        double sind = std::sin(delta);

        double cosa = std::cos(alpha);
        double cosd = std::cos(delta);

        rot(0, 0) = -sina; rot(0,1) = -cosa * sind; rot(0,2) = cosa * cosd;
        rot(1, 0) = cosa;  rot(1,1) = -sina * sind; rot(1,2) = cosd * sina;
        rot(2, 0) = 0.0;   rot(2,1) = cosd;         rot(2,2) = sind;
    }

    /* @brief Obtains the rotation matrix for transforming from the EME frame to the BME frame at a given epoch
    *  @param[in] epoch The epoch of the transformation in ephemeris seconds past J2000
    *  @param[out] rot The rotation matrix to transform from the EME to the BME frame (6x6)
    */
    template <typename matrixType>
    void getEMEtoBMERotationMatrix(double &epoch, Eigen::Matrix<matrixType,6,6> &rot)
    {
        // Convert target string to integer ID
        SpiceDouble rotationMatrix[6][6];

        // Construct the name of the reference frame
        std::string desFrame = "IAU_"+PARAMS::TARGET;
        std::string refFrame = "J2000";

        // Call the rotation matrix generator
        sxform_c(refFrame.c_str(), desFrame.c_str(), epoch, rotationMatrix);

        // Copy the rotation matrix into the output vector
        for (unsigned i = 0; i < 6; ++i)
        {
            for (unsigned j = 0; j < 6; ++j)
            {
                rot(i,j) = rotationMatrix[i][j];
            }
        }
    }

    /*  @brief Converts a BME@Epoch field to an EME2000 field.
        @param[in] bmeField: bmeField to convert
        @param[out] emeField: emeField to store the conversion in.
    */
    template <typename Type>
    void BMEtoEME(ACROBAT::bmeField<Point<Type>> &bmeField, ACROBAT::emeField<Point<Type>> &emeField)
    {
        Eigen::Matrix<Type, 3, 3> Qbe;

        // Create the transformation matrix
        getBMEtoEMERotationMatrix(PARAMS::EPOCH, Qbe);

        // Apply Qbe to every state in bmeField
        #pragma omp parallel for shared(Qbe)
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

                // Swap back
                for (unsigned idx = 0; idx < 6; ++idx) temp[i] = xb(i);

                // Assign
                emeField.setValue(temp, i, j);
            }
        }
    }
        
    /*  @brief Obtains the points in the set for a given stability index stabNum across the whole EMEJ2000 field.
        @param[in] stabNum The stability index
        @param[in] field EME2000 Field containing the domain to integrate
        @param[out] points std::vector<> of Points containing the indices of points in the set
    */
    template <typename integerType, typename fieldType>
    void getSet(const integerType stabNum, ACROBAT::emeField<fieldType>& field, std::vector<Point<int>> &pointIndices, std::vector<Point<double>>& pointCoordinates)
    {
        // Statuses
        std::vector<int> setStatistics(4);
        std::vector<Point<int>> privatePoints;
        std::vector<Point<double>> privateCoordinates;
        // Get the integration direction (+ve or -ve time?)
        int directionTime = sgn(stabNum);
        unsigned prog = 0;
        // Iterate through the conditions on the EME2000 field - do jumps of our rank and pool size
        int rank, poolSize, returnStatus;
        MPI_Status mpiStatus;
        returnStatus = MPI_Comm_size(MPI_COMM_WORLD, &poolSize);
        returnStatus = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        for (unsigned i = rank; i < field.getXExtent(); i += poolSize)
        {
            for (unsigned j = 0; j < field.getYExtent(); ++j)
            {
                // Get current point
                Point<double> currentPoint = field.getValue(i, j);
                std::vector<double> stateAtCrossing(6);
 
                
                // Get the status of this points behaviour
                int status = getStatus(currentPoint, field.getInitialTime(), directionTime, stateAtCrossing);
                setStatistics[status-1]++;
                // If weakly stable, append its indices to the points vector
                if ( (directionTime > 0 && status == 3) || (directionTime < 0 && status == 2) )  // Want n-revolution points or backwards escape
                {
                    Point<int> thisPoint; thisPoint[0] = i; thisPoint[1] = j;
                    std::vector<double> dimensionalisedStateAtCrossing;
                    getDimState(stateAtCrossing, dimensionalisedStateAtCrossing);
                    Point<double> coordinateToAdd(dimensionalisedStateAtCrossing);
                    privatePoints.insert(privatePoints.end(), thisPoint);
                    privateCoordinates.insert(privateCoordinates.end(), coordinateToAdd);
                }
                prog += poolSize;
            }
            if (rank == 0) std::cout << "Completed integration " << prog << " of " << (field.getXExtent() * field.getYExtent()) << "." << std::endl;
        }
        // Now every worker has completed their chunk, their results must all be sent to the rank 0 processor
        reduceVector(privatePoints, pointIndices, rank, poolSize);
        reduceVector(privateCoordinates, pointCoordinates, rank, poolSize);
        broadcastVector(pointIndices, rank);
        broadcastVector(pointCoordinates, rank);
        for (unsigned i = 0; i < setStatistics.size(); ++i)
        {
            reduceCount(setStatistics[i], setStatistics[i]); // Send from all cores to receive at master
            broadcastCount(setStatistics[i]);  // Send from master to all cores
        }
        if (rank == 0) printStatistics(stabNum, setStatistics);
    } // function


    /*  @brief Performs a step using a given boost stepper; returns the current time and the new state if successful
    *   @param[in] stepper Boost::odeint::numeric object corresponding to a controlled stepper
    *   @param[inout] currentTime On input, it contains the time of integration at the start of step. On exit, it contains the new time (i.e. after one step)
    *   @param[inout] x On input, it contains the state of the particle at the previous timestep. On exit, it contains the new state of the particle (i.e. after one step)
    *   @param[inout] dt On input, it contains a guess for the time-step to be taken. On exit, it contains the actual time-step taken.
    *   @param[in] f The force function to integrate.
    */
    template <typename stepperType>
    void make_step(stepperType& stepper, std::vector<double>& x, double& currentTime, double& dt)
    {
        boost::numeric::odeint::controlled_step_result result = boost::numeric::odeint::fail;
        
        /* We only want this to return back to the calling function when the step was successful.
        Therefore, keep retrying this until the stepper returns a value that was within tolerance.
        */
        while(result == boost::numeric::odeint::fail)
        {
            result = stepper.try_step(forceFunction, x, currentTime, dt);
        }
    }

    /* @brief Computes whether a given point is acrobatic, weakly stable, crashes, or escapes.
    *  @param[in] point The initial conditions (in a Point<> structure) to be tested
    *  @param[in] initTime The initial time for the integration of the trajectory
    *  @param[in] direction The direction for the timespan of the integration (+1/-1)
    *  @returns Status code corresponding to the behaviour of the trajectory
    */
    template <typename pointType, typename doubleType, typename integerType>
    int getStatus(Point<pointType>& point, const doubleType& initTime, const integerType& direction, std::vector<double>& pointAtCrossing)
    {
        // Initialise the stepper to be used
        typedef std::vector<double> stateType;
        boost::numeric::odeint::runge_kutta_fehlberg78<stateType> method;
        auto stepper = boost::numeric::odeint::make_controlled(/*reltol*/ 1e-012, /*absTol*/ 1e-012, method);

        // Fill an initial condition vector and normalize on entry
        stateType x0Dim(6), xDim(6), x0, x;
        for (unsigned idx = 0; idx < 6; ++idx) 
        {
            x0Dim[idx] = point[idx];
            xDim[idx] = point[idx];
        }

        // Normalise
        getNonDimState(x0Dim, x0);
        getNonDimState(xDim, x);

        // Initialise current time
        double currentTime = initTime;
        double dt = 0.01 * direction;

        // Integration status
        int status = 0;
        // Initialiser for the previous value of condition one
        double prevCondOne = 1.0 * direction;

        while (status == 0) // While none of the stopping conditions have been verified
        {
            /* Make a step using the given solver & force function */
            make_step(stepper, x, currentTime, dt);

            /* Call the integrator observer function */
            status = integrationController(x, x0, currentTime, prevCondOne);
        }
        // Obtain the exact point of the zero if status == 3
        if (status == 3)
        {
            obtainZero(x0, x, currentTime);
            pointAtCrossing = x;
        } else 
        {
            pointAtCrossing = x;
        }
        return status; // Return status when it doesn't correspond to 'keep going'
    }
}
#endif
