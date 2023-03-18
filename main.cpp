//g++ -std=c++14 utility_func.cpp  main.cpp -o main -I/opt/homebrew/Cellar/boost/1.81.0_1/include -lm

#include "utility_func.hpp"

#include <boost/algorithm/string/trim.hpp>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

using namespace std;

int main() {
    string EDS = readEDSFile("input_04.txt");
    // Remove leading and trailing white spaces. - unnecessary
    boost::trim(EDS);
    cout << EDS << endl;

    eds_matrix eds_segments = EDSToMatrix(EDS);
    weight_matrix weights = getGCContentWeights(eds_segments);

    score_matrix scores = initScoreMatrix(weights);
    score_matrix choices = initScoreMatrix(weights);

    int result = findMaxScoringPaths(eds_segments, weights, scores, choices, 2);
    cout << "RESULT: " << result << endl;
    return 0;
}
