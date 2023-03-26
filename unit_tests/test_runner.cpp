// TODO: Do this with a test module (Gtest, boost test)

//g++ -std=c++14 utility_func.cpp  unit_tests/test_runner.cpp -o test -I/opt/homebrew/Cellar/boost/1.81.0_1/include -lm

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>


using namespace std;

void runAllTests();

int main(int argc, char** argv){
    runAllTests();
    // run_test(TEST_readEDSFile);
    // run_test(TEST_EDSToMatrix);
    // run_test(TEST_findMaxScoringPaths_01);
    // run_test(TEST_findMaxScoringPaths_02);
    // run_test(TEST_findMaxScoringPaths_03);
    // run_test(TEST_findMaxScoringPaths_04);
    // run_test(TEST_findMaxScoringPaths_05);
    //printEDSwithWeights();

    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}