// g++ -std=c++17 utility_func.cpp  main.cpp -o main
// clang++ -std=c++17 utility_func.cpp  main.cpp -o main

#include <iostream>

#include "utility_func.hpp"

using namespace std;

int main(int argc, char* argv[]) {
    string EDS = readEDSFile(argc > 1 ? argv[1] : "input_01.txt");

    eds_matrix eds_segments = EDSToMatrix(EDS);
    weight_matrix weights = getGCContentWeights(eds_segments);

    score_matrix scores = initScoreMatrix(weights);
    score_matrix choices = initScoreMatrix(weights);

    int result = findMaxScoringPaths(eds_segments, weights, scores, choices, 2);
    auto paths = getPaths(eds_segments, scores, choices);
    cout << "Score: " << result << endl;
    printPaths(paths);
    return 0;
}
