// cd build
// make
// ./main

#include <iostream>
#include <iomanip>

#include "utility_func.hpp"

using namespace std;

int main(int argc, char* argv[]) {
    string EDS = readEDSFile(argc > 1 ? argv[1] : "../unit_tests/test_inputs/input_01.txt");

    eds_matrix eds_segments = EDSToMatrix(EDS);
    //cout << "Loaded the graph" << endl;
    weight_matrix weights = getGCContentWeights(eds_segments, 1, -2);
    //cout << "Assigned weights" << endl;

    score_matrix scores = initScoreMatrix(weights);
    score_matrix choices = initScoreMatrix(weights);

    int result = findMaxScoringPaths(eds_segments, weights, scores, choices, 10);
    //cout << "Found paths and calculated max score" << endl;
    cout << "Score: " << result << endl;
    auto paths = getPaths(eds_segments, scores, choices);
    //cout << "Finished getting the paths" << endl;
    cout << "Number of found paths: " << paths.size() << endl;
    cout << setprecision(2) << fixed;
    cout << "Paths cover the " << pathCoverPercentage(eds_segments, paths) << "\% of the graph\n" ;
    cout << "Average length of paths is: " << pathsAverageLength(paths) << endl;
    cout << result << "\t\t" << paths.size() << "\t\t" << pathCoverPercentage(eds_segments, paths) << "%\t\t" << pathsAverageLength(paths) << endl;
    //printPaths(paths);
    return 0;
}
