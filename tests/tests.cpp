#include <gtest/gtest.h>
#include "../src/laguerreConway.hpp"
#include <iostream>
#include <cmath>

double referenceSolutions(double &E, double &ecc)
{
	return E - ecc * std::sin(E);
}

TEST (laguerreConwayTest, meanToEccentric)
{
	/* Compute a bunch of mean anomalies manually, and check that the iterative scheme can recover
	 * the answer at each step */
	double solverTolerance = 1e-013;
	for (double ecc = 0.0; ecc < 1.0; ecc += .05)
	{
		for (double E = 0.0; E < 360.; E+= 10.) // 0 -> 360 degrees
		{
			// Get the reference solution
			double meanAnomaly = referenceSolutions(E, ecc);
			// Attempt to recover the value of E from the Laguerre-conway method
			ASSERT_NEAR( meanToEccentric(meanAnomaly, ecc, solverTolerance), E, solverTolerance*10); // meanAnomaly, eccentricity, solver epsilon
		}
	}
}

int main(int argc, char **argv)
{
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}

