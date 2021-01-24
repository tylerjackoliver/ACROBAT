#include "../src/rootFinding.hpp"
#include <cmath>
#include <vector>
#include <Eigen/eigen3/core>
#include <random>

/* @brief Contains functions for testing what can be tested programmatically and simply
 * in rootFinding.hpp.
 *
 * cross3 and dotProduct are tested here.
 * The remaining functions - findConditionOne, obtainZero - were tested at write-time using
 * 'manual' methods. Unfortunately, since these solutions are prone to numerical error, testing them
 * using the GTest framework is difficult.
 */

/* Function for generating reference solutions */
Eigen::Vector3d generateReferenceCrossProduct(std::vector<double>& a, std::vector<double>& b)
{
	if (a.size() != b.size()) throw std::runtime_error("Inconsistent sizes: generateReferenceCrossProduct");
	Eigen::Matrix<double, a.size(), 1> testA, testB;
	for (int i = 0; i < a.size(); ++i)
	{
		testA(i) = a[i];
		testB(i) = b[i];
	}
	return testA.cross(testB);
}

double generateReferenceDotProduct(std::vector<double> &a, std::vector<double>& b)
{
	if (a.size() != b.size()) throw std::runtime_error("Inconsistent sizes: generateReferenceDotProduct");
	Eigen::Matrix<double, a.size(), 1> testA, testB;
	for (int i = 0; i < a.size(); ++i)
	{
		testA(i) = a[i];
		testB(i) = b[i];
	}
	return testA.dot(testB);
}

TEST( cross3, EigenCompare )
{
	std::minstd_rand rng(42); // Seed
	int numReps = 1e4;
	// Perform a bunch of randomly-generated cross products
	for (int rep = 0; rep < numReps; ++rep)
	{
		std::vector<double> a(3), b(3), crossProduct(3);
		for (int dimension = 0; dimension < 3; ++dimension)
		{
			a[i] = rng();
			b[i] = rng();
		}
		Eigen::Vector3d reference = generateReferenceCrossProduct(a, b);
		cross3(a, b, crossProduct); // One we're testing
		for (int dimension = 0; dimension < 3; ++dimension)
		{
			ASSERT_DOUBLE( reference(i), crossProduct[i] );
		}
	}
}

TEST( dotProduct, EigenCompare )
{
	std::minstd_rand(42); // Seed
	int numReps = 1e4;
	// Perform randomly-generated dot products
	for (int rep = 0; rep < numReps; ++rep)
	{
		std::vector<double> a(3), b(3);
		double dotProduct;
		for (int dimension = 0; dimension < 3; ++dimension)
		{
			a[i] = rng();
			b[i] = rng();
		}
		dotProduct = dotProduct(a, b);
		double referenceSolution = generateReferenceDotProduct(a, b);
		ASSERT_DOUBLE( referenceSolution, dotProduct );
	}
}

TEST( cross3, checkDimensionOK )
{
	std::vector<double> a(3), b(3);
	EXPECT_NO_THROW( cross3(a, b) );
}

TEST( cross3, checkDimensionNotOK)
{
	std::vector<double> a(3), b(3), c(2), d;
	try {
		cross3(a, c, d);
		FAIL() << "Expected std::domain_error";
	}
	catch(std::out_of_range const & err) {
		EXPECT_EQ(err.what(), std::string("One of the input vectors to cross3 are not of size 3."));
	}
	catch(...) {
		FAIL() << "Expected std::domain_error";
	}

}
