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
    string EDS = readEDSFile("test_input.txt");
    // Remove leading and trailing white spaces.
    boost::trim(EDS);
    cout << EDS << endl;

    vector<vector<string>> eds_segments = EDSToMatrix(EDS);
    vector<vector<vector<int>>> weights = getGCContentWeights(eds_segments);

    for (int i = 0; i < weights.size(); i++) {
        for (int j = 0; j < weights[i].size(); j++) {
            cout << eds_segments[i][j] << endl;
            for (int k = 0; k < weights[i][j].size(); k++) {
                cout << weights[i][j][k];
            }
            cout << endl;
        }
    }

    score_matrix scores = initScoreMatrix(weights);
    cout << scores.size() << endl;
    // cout << scores[1].size() << endl;
    // cout << scores[1][1].size() << endl;
    // cout << scores[1][1][0].size() << endl;
    // cout << scores[1][1][0][0].size() << endl;

    int result = findMaxScoringPaths(eds_segments, weights, scores, 10);
    cout << "RESULT: " << result << endl;
    cout << "END" << endl;
    return 0;
}
