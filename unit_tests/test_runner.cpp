// TODO: Do this with a test module (Gtest, boost test)

//g++ -std=c++14 utility_func.cpp  unit_tests/test_runner.cpp -o test -I/opt/homebrew/Cellar/boost/1.81.0_1/include -lm

#include "tests.cpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

void run_test(bool (*function)()) {
    cout << (function() ? "passed\n" : "FAILED!\n");
}

int main(){
    cout << "Running test_runner" << endl;

    run_test(TEST_readEDSFile);
    run_test(TEST_EDSToMatrix);
    printEDSwithWeights();
}