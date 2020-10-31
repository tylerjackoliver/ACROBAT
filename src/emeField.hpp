#ifndef __EME_FIELD_H__
#define __EME_FIELD_H__

#include "field.hpp"
#include <ostream>
#include "Params.hpp"
#include "OEs.hpp"
#include "integration.hpp"
#include "bmeField.hpp"
#include <unordered_map>

extern "C"
{
    #include "SpiceUsr.h"
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
            void initialiseField(ACROBAT::oeField &input)
            {
                BMEtoEME(input, this);
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
            void setinitialTime(const double initialTime)
            {
                this->_initialTime = initialTime;
            }

            /* @brief Obtains the stabNumber-revolution stable set for a given EME2000 field. Note than stabNumber entries in the map are returned.
            * @param[in] stabNumber The number of revolutions a particle should complete in order to be classed as stable.
            * @param[out] std::unordered_map The keys are the number of rotations, and the values are a std::vector<> containing the points comprising the set.
            */
            template <typename integerType, typename vectorType>
            void getStableSet(const integerType, std::unordered_map<int, std::vector<Point<vectorType>>>);
                    
            /* @brief Obtains the rotation matrix for transforming from the BME frame to the EME frame at a given epoch
            *  @param[in] epoch The epoch of the transformation in ephemeris seconds past J2000
            *  @param[out] rot The rotation matrix to transform from the BME to the EME frame (6x6)
            */
            template <typename matrixType>
            void getBMEtoEMERotationMatrix(const double&, Eigen::Matrix<matrixType, 6, 6>&);

            /* @brief Obtains the rotation matrix for transforming from the EME frame to the BME frame at a given epoch
            *  @param[in] epoch The epoch of the transformation in ephemeris seconds past J2000
            *  @param[out] rot The rotation matrix to transform from the EME to the BME frame (6x6)
            */
            template <typename matrixType>
            void getEMEtoBMERotationMatrix(double&, Eigen::Matrix<matrixType, 6, 6>&);

            /*  @brief Converts a BME@Epoch field to an EME2000 field.
                @param[in] bmeField: bmeField to convert
                @param[out] emeField: emeField to store the conversion in.
            */
            template <typename pointType>
            void BMEtoEME(const ACROBAT::bmeField<Point<Type>>, ACROBAT::emeField<Point<Type>>&);

            /* @brief Converts an EME2000 field to a BME@Epoch field
                @param[in] emeField: emeField to convert
                @param[out] bmeField: bmeField to store the conversion in.
            */
            template <typename pointType>
            void EMEtoBME(const ACROBAT::emeField<Point<pointType>>&, ACROBAT::bmeField<Point<pointType>>&);
        
            /*  @brief Obtains the points in the set for a given stability index stabNum across the whole EMEJ2000 field.
                @param[in] stabNum The stability index
                @param[in] field EME2000 Field containing the domain to integrate
                @param[out] points std::vector<> of Points containing the indices of points in the set
            */
            template <typename integerType, typename fieldType, typename vectorType>
            void getSet(const integerType&, const ACROBAT::emeField<fieldType>&, std::vector<Point<vectorType>>&);
        
            /*  @brief Obtains the capture set from a given set of indices corresponding to points in an EME2000 field.
                @param[in] stabNum The stability index the set is sought for
                @param[in] domain The vector of points containing indices for which the sets are to be based off of
                @param[in] field The domain for which the points in domain correspond to
                @param[out] points A vector containing the indices of the corresponding capture set
            */ 
           template <typename integerType, typename fieldType, typename vectorType>
           void getSetFromPoints(const integerType&, const std::vector<Point<vectorType>>&, 
                                 const ACROBAT::emeField<fieldType>&, std::vector<Point<vectorType>>&);
                
            /*  @brief Performs a step using a given boost stepper; returns the current time and the new state if successful
            *   @param[in] stepper Boost::odeint::numeric object corresponding to a controlled stepper
            *   @param[inout] currentTime On input, it contains the time of integration at the start of step. On exit, it contains the new time (i.e. after one step)
            *   @param[inout] x On input, it contains the state of the particle at the previous timestep. On exit, it contains the new state of the particle (i.e. after one step)
            *   @param[inout] dt On input, it contains a guess for the time-step to be taken. On exit, it contains the actual time-step taken.
            *   @param[in] f The force function to integrate.
            */
           template <typename stepperType, typename timeType, typename stateType, typename function>
           void make_step(stepperType&, stateType&, timeType&, timeType&, function&);

            /* @brief Prints the statistics for a given ballistic set computation
               @param[in] stabNum The number of the sability index for the current computation
               @param[in] setStatistics std::vector<> containing the statistics defined as in the set determination routines.
            */
           template <typename integerType, typename vectorType>
           void printStatistics(const integerType stabNum, std::vector<vectorType>&) const;

            /* @brief Computes whether a given point is acrobatic, weakly stable, crashes, or escapes.
            *  @param[in] point The initial conditions (in a Point<> structure) to be tested
            *  @param[in] initTime The initial time for the integration of the trajectory
            *  @param[in] direction The direction for the timespan of the integration (+1/-1)
            *  @returns Status code corresponding to the behaviour of the trajectory
            */
            template <typename pointType, typename doubleType, typename integerType>
            int getStatus(Point<PointType>&, const doubleType&, const integerType&)

        private:
            double _finalTime;      // Final time for the trajectory integration
            double _initialTime;    // Initial time for the trajectory integration
    };

    template<class classType, typename integerType, typename vectorType>
    void emeField<classType>::getStableSet(const integerType stabNumber, std::unordered_map<int, std::vector<Point<vectorType>>>& out);
    {
        std::vector<Point<vectorType>> points;
        // Get the set defined solely from the initial domain (i.e. 1-stable)
        getSet(1, this, points);
        out[1] = points;

        for (unsigned revs = 2; revs <= stabNumber; ++revs)  // Must complete at least one orbit about the host planet
        {
            // Clear points from the previous run
            points.clear();
            getSetFromPoints(revs, this, points);
            out[revs] = points;
        }
    }

    template <typename matrixType>
    void getBMEtoEMERotationMatrix(const double &epoch, Eigen::Matrix<matrixType,6,6> &rot)
    {
        // Convert target string to integer ID
        SpiceDouble rotationMatrix[6][6];

        // Construct the name of the reference frame
        std::string refFrame = "IAU_"+PARAMS::TARGET;
        std::string desFrame = "J2000";

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

    template <typename Type>
    void BMEtoEME(const ACROBAT::bmeField<Point<Type>> &bmeField, ACROBAT::emeField<Point<Type>> &emeField)
    {
        Eigen::Matrix<Type, 6, 6> Qbe;

        // Create the transformation matrix
        getBMEtoEMERotationMatrix(PARAMS::EPOCH, Qbe);

        // Apply Qbe to every state in bmeField
        #pragma omp parallel for shared(Qbe)
        for (unsigned int i = 0; i < bmeField.getXExtent(); ++i)
        {
            for (unsigned int j = 0; j < bmeField.getYExtent(); ++j)
            {
                for (unsigned int k = 0; k < bmeField.getZExtent(); ++k)
                {
                    // Set up input vector
                    Eigen::Matrix<Type, 6, 1> xb, xe;
                    Point<Type> temp = bmeField.getValue(i, j, k);
                    
                    // Assign
                    for(unsigned idx = 0; idx < 6; ++idx) xb(i) = temp[i];

                    // Compute
                    xe = Qbe * xb;

                    // Swap back
                    for (unsigned idx = 0; idx < 6; ++idx) temp[i] = xe(i);

                    // Assign
                    emeField.setValue(i, j, k, &temp);
                }
            }
        }
    }

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
                for (unsigned int k = 0; k < bmeField.getZExtent(); ++k)
                {
                    // Set up input vector
                    Eigen::Matrix<Type, 6, 1> xb, xe;
                    Point<Type> temp = emeField.getValue(i, j, k);
                    
                    // Assign
                    for(unsigned idx = 0; idx < 6; ++idx) xe(i) = temp[i];

                    // Compute
                    xb = Qbe * xe;

                    // Swap back
                    for (unsigned idx = 0; idx < 6; ++idx) temp[i] = xb(i);

                    // Assign
                    emeField.setValue(i, j, k, &temp);
                }
            }
        }
    }

    template <typename integerType, typename fieldType, typename vectorType>
    void getSet(const integerType stabNum, const ACROBAT::emeField<fieldType>& field, std::vector<Point<vectorType>> &points)
    {
        // Statuses
        std::vector<unsigned long> setStatistics(4);

        // Get the integration direction (+ve or -ve time?)
        int directionTime = sgn(field.getFinalTime() - field.getInitialTime());

        // Iterate through the conditions on the EME2000 field
        #pragma omp parallel for shared(directionTime, setStatistics, points)
        for (unsigned i = 0; i < field.getXExtent(); ++i)
        {
            for (unsigned j = 0; j < field.getYExtent(); ++j)
            {
                for (unsigned k = 0; k < field.getZExtent(); ++k)
                {
                    // Get current point
                    Point<double> currentPoint = field.getValue(indices[0], indices[1], indices[2]);

                    // Get the status of this points behaviour
                    int status = getStatus(currentPoint, field.getInitialTime(), directionTime);

                    // Increment the correct counter in the setStatistics vector (ordered s.t. indexes is status-1)
                    #pragma omp atomic
                    setStatistics[status-1]++;

                    // If weakly stable, append its indices to the points vector
                    if (status == 2)
                    {
                        #pragma omp critical
                        {
                            points.inset(points.end(), currentPoint);
                        }
                    }
                }
            }
        }
        printStatistics(stabNum, setStatistics);
    } // function

    template <typename integerType, typename fieldType, typename vectorType>
    void getSetFromPoints(const integerType& stabNum, const std::vector<Point<vectorType>> &domain, const ACROBAT::emeField<fieldType>& field, std::vector<Point<vectorType>> &points)
    {    
        // Statuses - ordered s.t. index related to condition is (status code - 1) - crash, escape, stable, acrobatic
        std::vector<unsigned long> setStatistics(4);

        // Get the integration direction (+ve or -ve time?)
        int directionTime = sgn(field.getFinalTime() - field.getInitialTime());

        // Iterate through the conditions on the EME2000 field
        #pragma omp parallel for shared(directionTime)
        for (unsigned i = 0; i < domain.size(); ++i)
        {
            // Get an initial position
            Point<vectorType> indices = domain[i];
            Point<double> currentPoint = field.getValue(indices[0], indices[1], indices[2]);

            // Get the status of this points behaviour
            int status = getStatus(currentPoint, field.getInitialTime(), directionTime);
            
            // Increment the correct counter in the setStatistics vector (ordered s.t. index is status-1)
            #pragma omp atomic
            setStatistics[status-1]++;

            // If weakly stable, append its indices to the points vector
            if (status == 2)
            {
                #pragma omp critical
                {
                    points.inset(points.end(), currentPoint);
                }
            }
        }
        } // end of parallel section
        printStatistics(stabNum, setStatistics);
    } // function

    template <typename stepperType, typename timeType, typename stateType, typename function>
    void make_step(stepperType& stepper, stateType &x, timeType& currentTime, timeType& dt, function& f)
    {
        boost::numeric::odeint::controlled_step_result result = boost::numeric::odeint::fail;
        
        /* We only want this to return back to the calling function when the step was successful.
        Therefore, keep retrying this until the stepper returns a value that was within tolerance.
        */
        while(result == boost::numeric::odeint::fail)
        {
            result = stepper.try_step(f, x, currentTime, dt);
        }
    }

    template <typename integerType, typename vectorType>
    void printStatistics(const integerType stabNum, const std::vector<vectorType>& setStatistics) const
    {
        std::cout << "For a stability index of " stabNum << "the statistics are as follows:" << std::endl;
        std::cout << "\t Number of crashes:   " << setStatistics[0] << std::endl;
        std::cout << "\t Number of escapes:   " << setStatistics[1] << std::endl;
        std::cout << "\t Number of stable:    " << setStatistics[2] << std::endl;
        std::cout << "\t Number of acrobatic: " << setStatistics[3] << std::endl;
    }

    template <typename pointType, typename doubleType, typename integerType>
    int getStatus(Point<PointType>& point, const doubleType& initTime, const integerType& direction)
    {
        // Initialise the stepper to be used
        typedef Eigen::Matrix<double, 6, 1> stateType;
        boost::numeric::odeint::runge_kutta_fehlberg78<stateType> method;
        auto stepper = boost::numeric::odeint::make_controlled(/*reltol*/ 1e-011, /*absTol*/ 1e-011);

        // Fill an initial condition vector
        stateType x0, x;
        for (unsigned idx = 0; idx < 6; ++i) 
        {
            x0[idx] = point[idx];
            x[idx] = point[idx];
        }

        // Initialise current time
        double currentTime = initTime;
        double dt = 0.01 * currentTime * direction;

        // Integration status
        int status = 0;

        while (status == 0) // While none of the stopping conditions have been verified
        {
            /* Make a step using the given solver & force function */
            make_step(stepper, forceFunction, x, currentTime, dt);

            /* Call the integrator observer function */
            status = integrationController(x, x0, currentTime);
        }
        return status; // Return status when it doesn't correspond to 'keep going'
    }
}
#endif
