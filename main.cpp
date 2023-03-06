#include "utility_func.hpp"

#include <boost/algorithm/string/trim.hpp>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

using namespace std;

// // Returns the one line elastic degenerate string file as a string.
// string readEDSFile(const string &file_path) {
//     ifstream input_stream(file_path);

//     if (!input_stream.good()) {
//         cout << "Reading of file " << file_path << " failed." << endl;
//         return "";
//     }
//     assert(input_stream.good());
//     return static_cast<stringstream const &>(stringstream()
//                                              << input_stream.rdbuf())
//         .str();
// }

// // Store the elastic degenerate string `text` in a 2D matrix.
// vector<vector<string>> EDSToMatrix(const string &EDS) {
//     vector<vector<string>> eds_segments;
//     vector<string> current_segment;
//     string current_string;
//     bool in_nondet_segment = false;

//     for (char const &c : EDS) {
//         // Inside a string or segment.
//         if (c != '{' && c != '}') {
//             // Commas can be only inside segments.
//             assert(in_nondet_segment || c != ',');

//             if (c == ',') {
//                 current_segment.emplace_back(current_string);
//                 current_string.clear();
//             } else {
//                 current_string += c;
//             }
//         }
//         // Start of a non-deterministic segment.
//         else if (c == '{') {
//             assert(!in_nondet_segment && current_segment.empty());
//             in_nondet_segment = true;
//             // Start of a non-deterministic segment after a deterministic
//             // segment.
//             if (!current_string.empty()) {
//                 eds_segments.emplace_back(vector<string>{move(current_string)});
//                 current_string.clear();
//             }
//         }
//         // End of a segment.
//         else {
//             assert(c == '}' && in_nondet_segment && !current_segment.empty());
//             current_segment.emplace_back(current_string);
//             current_string.clear();

//             eds_segments.emplace_back(current_segment);
//             current_segment.clear();

//             in_nondet_segment = false;
//         }
//     }

//     // EDS ended with a deterministic string.
//     if (!current_string.empty()) {
//         assert(!in_nondet_segment && current_segment.empty());
//         eds_segments.emplace_back(vector<string>{move(current_string)});
//     }
//     return eds_segments;
// }

// // Assign score based on GC content:
// // - bases G and C get score 1
// // - empty non-deterministic segment parts get score 0
// // - every other character gets score 0.
// vector<vector<vector<int>>> getGCContentWeights(
//     const vector<vector<string>> &eds_segments) {
//     vector<vector<vector<int>>> weights;
//     weights.reserve(eds_segments.size());
//     for (const auto &segment : eds_segments) {
//         vector<vector<int>> w_segment;
//         w_segment.reserve(segment.size());
//         for (const auto &str : segment) {
//             vector<int> w_str(str.size(), 0);

//             // To an empty string assign value 0.
//             if (str.empty()) w_str.push_back(0);

//             for (int i = 0; i < str.size(); i++) {
//                 w_str[i] = (str[i] == 'G' || str[i] == 'C') ? 1 : 0;
//             }
//             w_segment.emplace_back(w_str);
//         }
//         weights.emplace_back(w_segment);
//     }
//     return weights;
// }

// // Find maximum-scoring paths.

// struct vertex {
//     int segment;
//     int layer;
//     int index;
// };

// bool operator==(const vertex &a, const vertex &b) {
//     tie(a.segment, a.layer, a.index) == tie(b.segment, b.layer, b.index);
// }

// #define SELECTED true

// // Helper functions for vertex type.
// bool isLayerVertex(vertex v, const vector<vector<string>> &eds_segments) {
//     return eds_segments[v.segment].size() > 1;
// }

// bool isJVertex(vertex v, const vector<vector<string>> &eds_segments) {
//     return !isLayerVertex(v, eds_segments) &&
//            (v.index == 0 && v.segment > 0 &&
//             eds_segments[v.segment - 1].size() > 1);
// }

// bool isNVertex(vertex v, const vector<vector<string>> &eds_segments) {
//     return !isLayerVertex(v, eds_segments) && !isJVertex(v, eds_segments);
// }

// bool isBaseLayerVertex(vertex v, const vector<vector<string>> &eds_segments) {
//     assert(isBaseLayerVertex(v, eds_segments));
//     return v.layer == 0;
// }

// bool hasPredecessorVertex(vertex v) { return v.segment != 0 && v.index != 0; }

// vertex getPredecessorVertex(const vector<vector<string>> &eds_segments,
//                             vertex v, int layer = 0) {
//     assert(hasPredecessorVertex(v));
//     if (v.index != 0) {
//         return {.segment = v.segment, .layer = v.layer, .index = v.index - 1};
//     }
//     return {.segment = v.segment - 1,
//             .layer = layer,
//             .index = static_cast<int>(
//                 eds_segments[v.segment - 1][layer].size() - 1)};
// }

// int getPredecessorScore(const vector<vector<string>> &eds_segments,
//                         vector<vector<vector<vector<vector<int>>>>> &scores,
//                         vertex v, bool selected, int layer = 0) {
//     if (!hasPredecessorVertex(v)) {
//         return 0;
//     }
//     vertex predecessor = getPredecessorVertex(eds_segments, v, layer);
//     return scores[v.segment][v.layer][v.index][selected][layer];
// }

// void storePath(vertex store_for, vertex closest_predecessor_on_path,
//                unordered_map<vertex, vector<vertex>> &paths,
//                unordered_map<vertex, vertex> &belonging_path_start) {
//     vertex first_on_path;
//     // The vertex we store the path for is the beginning of the path.
//     // first_on_path = (store_for == closest_predecessor_on_path)
//     //                     ? store_for
//     //                     : belonging_path_start[closest_predecessor_on_path];
//     // belonging_path_start[store_for] = first_on_path;
//     // paths[first_on_path].emplace_back(store_for);
// }

// void findMaxScoringPaths(const vector<vector<string>> &eds_segments,
//                          const vector<vector<vector<int>>> &weights,
//                          vector<vector<vector<vector<vector<int>>>>> &scores,
//                          unordered_map<vertex, vector<vertex>> &paths,
//                          unordered_map<vertex, vertex> &belonging_path_start,
//                          int penalty = -1) {
//     for (int segment = 0; segment < eds_segments.size(); segment++) {
//         for (int layer = 0; layer < eds_segments[segment].size(); layer++) {
//             for (int index = 0; index < eds_segments[segment][layer].size();
//                  index++) {
//                 vertex a{.segment = segment, .layer = layer, .index = index};
//                 // a has type N vertex.
//                 if (isNVertex(a, eds_segments)) {
//                     assert(layer == 0);
//                     int not_selected_pred_score =
//                         getPredecessorScore(eds_segments, scores, a, !SELECTED,
//                                             layer) -
//                         penalty;
//                     int selected_pred_score = getPredecessorScore(
//                         eds_segments, scores, a, SELECTED, layer);

//                     scores[segment][layer][index][SELECTED][layer] =
//                         weights[segment][layer][index];
//                     scores[segment][layer][index][!SELECTED][layer] = max(1, 2);
//                 }
//             }
//         }
//     }
// }

int main() {
    string EDS = readEDSFile("input.txt");
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

    return 0;
}
