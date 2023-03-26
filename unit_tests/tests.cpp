// TODO: Do this with a test module (Gtest, boost test)

#include "../utility_func.hpp"
#include <iostream>

using namespace std;

bool TEST_readEDSFile() {
    cout << __func__ << ": ";
    return readEDSFile(
               "/Users/evi/Documents/magister/diplomovka/main_program/"
               "unit_tests/input_test_01.txt") ==
           "_GG{AGAA,GGGA,,ACCCCC}{AG,G}AGG{A,G}{C,}{A,AG}G{A,GA,CCC}{,A}_";
}

bool TEST_EDSToMatrix() {
    cout << __func__ << ": ";
    string EDS = "_GG{AGAA,GGGA,,ACCCCC}{AG,G}AGG{A,G}{C,}{A,AG}G{A,GA,CCC}{,A}_";
    eds_matrix eds_segments = EDSToMatrix(EDS);

    eds_matrix expected;
    expected.emplace_back(vector<string>{"_GG"});
    expected.emplace_back(vector<string>{"AGAA","GGGA","_","ACCCCC"});
    expected.emplace_back(vector<string>{"_"});
    expected.emplace_back(vector<string>{"AG","G"});
    expected.emplace_back(vector<string>{"AGG"});
    expected.emplace_back(vector<string>{"A","G"});
    expected.emplace_back(vector<string>{"_"});
    expected.emplace_back(vector<string>{"C","_"});
    expected.emplace_back(vector<string>{"_"});
    expected.emplace_back(vector<string>{"A","AG"});
    expected.emplace_back(vector<string>{"G"});
    expected.emplace_back(vector<string>{"A","GA","CCC"});
    expected.emplace_back(vector<string>{"_"});
    expected.emplace_back(vector<string>{"_","A"});
    expected.emplace_back(vector<string>{"_"});

    return eds_segments == expected;
}

int run_findMaxScoringPaths_01(string EDS) {
    eds_matrix eds_segments = EDSToMatrix(EDS);
    weight_matrix weights = getGCContentWeights(eds_segments, 1, -1);

    score_matrix scores = initScoreMatrix(weights);
    score_matrix choices = initScoreMatrix(weights);

    return findMaxScoringPaths(eds_segments, weights, scores, choices, 2);
}

bool TEST_findMaxScoringPaths_01() {
    cout << __func__ << ": ";

    string EDS = "_GGG{CCG,AGGGA}A_";
    return run_findMaxScoringPaths_01(EDS) == 5;
}

bool TEST_findMaxScoringPaths_02() {
    cout << __func__ << ": ";

    string EDS = "_GGG{ATT,ATA}AGCGC_";
    return run_findMaxScoringPaths_01(EDS) == 3;
}

bool TEST_findMaxScoringPaths_03() {
    cout << __func__ << ": ";

    string EDS = "_GGG{ATT,ATA,}AGCGC_";
    return run_findMaxScoringPaths_01(EDS) == 4;
}

bool TEST_findMaxScoringPaths_04() {
    cout << __func__ << ": ";

    string EDS = "_AG{GGG,,CCC}{AG,GCGG,AA}A{A,G}{G,CC}{AAAA,}_";
    return run_findMaxScoringPaths_01(EDS) == 9;
}

bool TEST_findMaxScoringPaths_05() {
    cout << __func__ << ": ";

    string EDS = "_GG{AGAA,GGGA,,ACCCCC}{AG,G}AGG{A,G}{C,}{A,AG}G{A,GA,CCC}{,A}_";
    return run_findMaxScoringPaths_01(EDS) == 14;
}

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
