#include "utility_func.hpp"

#include <boost/algorithm/string/trim.hpp>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

// Returns the one line elastic-degenerate string file as a string.
string readEDSFile(const string &file_path) {
    ifstream input_stream(file_path);

    if (!input_stream.good()) {
        cout << "Reading of file " << file_path << " failed." << endl;
        return "";
    }
    assert(input_stream.good());
    return static_cast<stringstream const &>(stringstream()
                                             << input_stream.rdbuf())
        .str();
}

// Store the elastic degenerate string `text` in a 2D matrix.
// eds_matrix[s][l]
eds_matrix EDSToMatrix(const string &EDS) {
    eds_matrix eds_segments;
    vector<string> current_segment;
    string current_string;
    bool in_nondet_segment = false;

    for (char const &c : EDS) {
        // Inside a string or segment.
        if (c != '{' && c != '}') {
            // Commas can be only inside segments.
            assert(in_nondet_segment || c != ',');

            if (c == ',') {
                current_segment.emplace_back(current_string);
                current_string.clear();
            } else {
                current_string += c;
            }
        }
        // Start of a non-deterministic segment.
        else if (c == '{') {
            assert(!in_nondet_segment && current_segment.empty());
            in_nondet_segment = true;
            // Start of a non-deterministic segment after a deterministic
            // segment.
            if (!current_string.empty()) {
                eds_segments.emplace_back(vector<string>{move(current_string)});
                current_string.clear();
            }
        }
        // End of a segment.
        else {
            assert(c == '}' && in_nondet_segment && !current_segment.empty());
            current_segment.emplace_back(current_string);
            current_string.clear();

            eds_segments.emplace_back(current_segment);
            current_segment.clear();

            in_nondet_segment = false;
        }
    }

    // EDS ended with a deterministic string.
    if (!current_string.empty()) {
        assert(!in_nondet_segment && current_segment.empty());
        eds_segments.emplace_back(vector<string>{move(current_string)});
    }
    return eds_segments;
}

// Assign score based on GC content:
// - bases G and C get score 1
// - empty non-deterministic segment parts get score 0
// - every other character gets score 0.
vector<vector<vector<int>>> getGCContentWeights(
    const eds_matrix &eds_segments) {
    vector<vector<vector<int>>> weights;
    weights.reserve(eds_segments.size());
    for (const auto &segment : eds_segments) {
        vector<vector<int>> w_segment;
        w_segment.reserve(segment.size());
        for (const auto &str : segment) {
            vector<int> w_str(str.size(), 0);

            // To an empty string assign value 0.
            if (str.empty()) w_str.push_back(0);

            for (int i = 0; i < str.size(); i++) {
                w_str[i] = (str[i] == 'G' || str[i] == 'C') ? 1 : 0;
            }
            w_segment.emplace_back(w_str);
        }
        weights.emplace_back(w_segment);
    }
    return weights;
}

// Find maximum-scoring paths.

bool operator==(const Vertex &a, const Vertex &b) {
    return tie(a.segment, a.layer, a.index) == tie(b.segment, b.layer, b.index);
}

// Helper functions for Vertex type.
bool isLayerVertex(Vertex v, const eds_matrix &eds_segments) {
    return eds_segments[v.segment].size() > 1;
}

bool isJVertex(Vertex v, const eds_matrix &eds_segments) {
    return !isLayerVertex(v, eds_segments) &&
           (v.index == 0 && v.segment > 0 &&
            eds_segments[v.segment - 1].size() > 1);
}

bool isNVertex(Vertex v, const eds_matrix &eds_segments) {
    return !isLayerVertex(v, eds_segments) && !isJVertex(v, eds_segments);
}

bool isBaseLayerVertex(Vertex v, const eds_matrix &eds_segments) {
    assert(isBaseLayerVertex(v, eds_segments));
    return v.layer == 0;
}

bool hasPredecessorVertex(Vertex v) { return v.segment != 0 && v.index != 0; }

Vertex getPredecessorVertex(const eds_matrix &eds_segments, Vertex v,
                            int layer = 0) {
    assert(hasPredecessorVertex(v));
    if (v.index != 0) {
        return {v.segment, v.layer, v.index - 1};
    }
    return {v.segment - 1, layer,
            static_cast<int>(eds_segments[v.segment - 1][layer].size() - 1)};
}

int getPredecessorScore(const eds_matrix &eds_segments, score_matrix &scores,
                        Vertex v, bool selected, int layer = 0) {
    if (!hasPredecessorVertex(v)) {
        // BAD!
        return 0;
    }
    Vertex predecessor = getPredecessorVertex(eds_segments, v, layer);
    return scores[v.segment][v.layer][v.index][selected][layer];
}

void storePath(Vertex store_for, Vertex closest_predecessor_on_path,
               path_matrix &paths,
               belonging_path_matrix &belonging_path_starts) {
    Vertex first_on_path;
    first_on_path =
        (store_for == closest_predecessor_on_path)
            ? store_for
            : belonging_path_starts[closest_predecessor_on_path.segment]
                                  [closest_predecessor_on_path.layer]
                                  [closest_predecessor_on_path.index];
    assert(first_on_path);

    belonging_path_starts[store_for.segment][store_for.layer][store_for.index] =
        first_on_path;
    paths[first_on_path.segment][first_on_path.layer][first_on_path.index]
        .emplace_back(store_for);
}

// vector<vector<vector<vector<vector<int>>>>>
// [segment][layer][index][{SELECTED, !SELECTED}][num of layers]
score_matrix initScoreMatrix(const vector<vector<vector<int>>> &weights) {
    score_matrix scores;
    scores.resize(weights.size());
    for(int segment = 0; segment < weights.size(); segment++) {
        scores[segment].resize(weights[segment].size());
        for(int layer = 0; layer < weights[segment].size(); layer++) {
            scores[segment][layer].resize(weights[segment][layer].size());
            for(int index = 0; index < weights[segment][layer].size(); index++) {
                scores[segment][layer][index].resize(2);
                scores[segment][layer][index][0].resize(weights[segment].size());
                scores[segment][layer][index][1].resize(weights[segment].size());
            }
        }
    }
    return scores;
}

//typedef vector<vector<vector<vector<Vertex>>>> path_matrix;
// no paths in the beginning
path_matrix initPathMatrix(const vector<vector<vector<int>>> &weights) {
    path_matrix paths;
    paths.resize(weights.size());
    for(int segment = 0; segment < weights.size(); segment++) {
        paths[segment].resize(weights[segment].size());
        for(int layer = 0; layer < weights[segment].size(); layer++) {
            paths[segment][layer].resize(weights[segment][layer].size());
        }
    }
    return paths;
}

// typedef vector<vector<vector<Vertex>>> belonging_path_matrix;
belonging_path_matrix initBelongingPathMatrix(
    const vector<vector<vector<int>>> &weights) {
    belonging_path_matrix belonging_path_starts;
    belonging_path_starts.resize(weights.size());
    for (int segment = 0; segment < weights.size(); segment++) {
        belonging_path_starts[segment].resize(weights[segment].size());
        for (int layer = 0; layer < weights[segment].size(); layer++) {
            belonging_path_starts[segment][layer].resize(
                weights[segment][layer].size());
        }
    }
    return belonging_path_starts;
}

void findMaxScoringPaths(const eds_matrix &eds_segments,
                         const vector<vector<vector<int>>> &weights,
                         score_matrix &scores, path_matrix &paths,
                         belonging_path_matrix &belonging_path_starts,
                         int penalty = -1) {
    for (int segment = 0; segment < eds_segments.size(); segment++) {
        for (int layer = 0; layer < eds_segments[segment].size(); layer++) {
            for (int index = 0; index < eds_segments[segment][layer].size();
                 index++) {
                Vertex a{segment, layer, index};
                // N vertex.
                if (isNVertex(a, eds_segments)) {
                    assert(layer == 0);

                    // W(a, 1, 0) = w(a) + max(W(p, 0, 0) - x, W(p, 1, 0))
                    scores[segment][layer][index][SELECTED][layer] =
                        weights[segment][layer][index];
                    int not_selected_pred_score =
                        getPredecessorScore(eds_segments, scores, a, !SELECTED,
                                            layer) -
                        penalty;
                    int selected_pred_score = getPredecessorScore(
                        eds_segments, scores, a, SELECTED, layer);

                    if (selected_pred_score >= not_selected_pred_score) {
                        scores[segment][layer][index][SELECTED][layer] +=
                            selected_pred_score;
                        storePath(
                            a,
                            hasPredecessorVertex(a)
                                ? getPredecessorVertex(eds_segments, a, layer)
                                : a,
                            paths, belonging_path_starts);
                    }
                    else {
                        scores[segment][layer][index][SELECTED][layer] += not_selected_pred_score;
                    }

                    // W(a, 0, 0) = max(W(p, 0, 0) - x, W(a, 1, 0))


                }
            }
        }
    }
}