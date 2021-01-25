#include <gtest/gtest.h>
#include "testLaguerreConway.hpp"
#include "testRootFinding.hpp"
#include "testIntegration.hpp"


int main(int argc, char **argv)
{
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}

