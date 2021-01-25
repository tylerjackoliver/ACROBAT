#include "../src/integration.hpp"
#include <cmath>
#include "../src/Params.hpp"
#include <random>

double _ASTRONOMICAL_UNIT_ = 1.495978707e8;

TEST( getDimensionalTime, badTargetGM )
{
	double hold = PARAMS::targetGM;
	PARAMS::targetGM = 0.0;
	double testTime = 1.0;

	try {
		getDimensionalTime(testTime);
		FAIL() << "Expected std::domain_error; no error thrown";
	} catch (std::domain_error const & err)
	{
		EXPECT_EQ(err.what(), std::string( "Error: target GM is zero in getDimensionalTime.") );
	} catch( ... )
	{
		FAIL() << "Expected std::domain_error; other error thrown";
	}
	PARAMS::targetGM = hold;
}

TEST( getDimensionalTime, checkValue )
{
	PARAMS::targetGM = 1.0;
	PARAMS::R = 1.0;
	std::minstd_rand rng(42);
	int numReps = 1e3;
	for (int i = 0; i < numReps; ++i)
	{
		double testTime = static_cast<double>( rng() );
		double result = getDimensionalTime( testTime );
		double reference = testTime * std::sqrt( std::pow( PARAMS::R, 3 ) / PARAMS::targetGM );
		EXPECT_DOUBLE_EQ(reference, result);
	}
	PARAMS::targetGM = 0.0;
	PARAMS::R = 0.0;
}

TEST( getNonDimensionalTime, badTargetGM )
{
	double hold = PARAMS::R;
	PARAMS::R = 0.0;
	double testTime = 1.0;

	try {
		getNonDimensionalTime(testTime);
		FAIL() << "Expected std::domain_error; no error thrown";
	} catch (std::domain_error const & err)
	{
		EXPECT_EQ(err.what(), std::string("Error: R is zero in getNonDimensionalTime.") );
	} catch( ... )
	{
		FAIL() << "Expected std::domain_error; other error thrown";
	}
	PARAMS::R = hold;
}

TEST( getNonDimensionalTime, checkValue )
{
	std::minstd_rand rng(42);
	int numReps = 1e3;
	for (int i = 0; i < numReps; ++i)
	{
		double testTime = static_cast<double>( rng() );
		PARAMS::R = static_cast<double>( rng() );
		PARAMS::targetGM = static_cast<double>( rng() );
		double result = getNonDimensionalTime( testTime );
		double reference = testTime * std::sqrt( PARAMS::targetGM / std::pow( PARAMS::R, 3 ) );
		EXPECT_DOUBLE_EQ(reference, result);
	}
	PARAMS::R = 0.0;
	PARAMS::targetGM = 0.0;
}
