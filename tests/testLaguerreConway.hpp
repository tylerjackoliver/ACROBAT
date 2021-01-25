#include "../src/laguerreConway.hpp"
#include <cmath>

double referenceMeanAnomaly(double &E, double &ecc)
{
	return E - ecc * std::sin(E);
}

TEST (laguerreConwayTest, meanToEccentric)
{
	/* Compute a bunch of mean anomalies manually, and check that the iterative scheme can recover
	 * the answer at each step */
	double solverTolerance = 1e-013;
	for (double ecc = 0.0; ecc < 1.; ecc += .05)
	{
		for (double E = 0.0; E <= 360.; E+= 5.) // 0 -> 360 degrees
		{
			// Get the reference solution
			double meanAnomaly = referenceMeanAnomaly(E, ecc);
			// Attempt to recover the value of E from the Laguerre-conway method
			ASSERT_NEAR( meanToEccentric(meanAnomaly, ecc, solverTolerance), E, solverTolerance*10); // meanAnomaly, eccentricity, solver epsilon
		}
	}
}
