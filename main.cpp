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

void printEDSwithWeights() {
    string EDS = readEDSFile("input_04.txt");
    cout << EDS << endl;

    eds_matrix eds_segments = EDSToMatrix(EDS);
    weight_matrix weights = getGCContentWeights(eds_segments);
    for (int i = 0; i < weights.size(); i++) {
        for (int j = 0; j < weights[i].size(); j++) {
            cout << eds_segments[i][j] << endl;
            for (int k = 0; k < weights[i][j].size(); k++) {
                cout << weights[i][j][k];
            }
            cout << endl;
        }
    }
}

void printChoices(const score_matrix &choices) {
    cout << "Choices" << endl;
    for (int i = 0; i < choices.size(); i++) {
        for (int j = 0; j < choices[i].size(); j++) {
            for (int k = 0; k < choices[i][j].size(); k++) {
                for (int l = 0; l < choices[i][j][k].size(); l++) {
                    for (int m = 0; m < choices[i][j][k][l].size(); m++) {
                        cout << choices[i][j][k][l][m];
                    }
                }
                cout <<" ";
            }
            cout << endl;
        }
    }
}

int main() {
    string EDS = readEDSFile("input_01.txt");
    // Remove leading and trailing white spaces. - unnecessary
    boost::trim(EDS);
    cout << EDS << endl;

    eds_matrix eds_segments = EDSToMatrix(EDS);
    weight_matrix weights = getGCContentWeights(eds_segments);

    score_matrix scores = initScoreMatrix(weights);
    score_matrix choices = initScoreMatrix(weights);

    int result = findMaxScoringPaths(eds_segments, weights, scores, choices, 2);
    auto paths = getPaths(eds_segments, scores, choices);
    cout << "RESULT: " << result << endl;
    printPaths(paths);
    //printChoices(choices);
    return 0;
}
