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

bool isFirstLayerVertex(Vertex v, const eds_matrix &eds_segments) {
    return isLayerVertex(v, eds_segments) && v.layer == 0;
}

bool isVertexFirstOnLayer(Vertex v, const eds_matrix &eds_segments) {
    return isLayerVertex(v, eds_segments) && v.index == 0;
}

bool isJVertex(Vertex v, const eds_matrix &eds_segments) {
    return !isLayerVertex(v, eds_segments) &&
           (v.index == 0 && v.segment > 0 &&
            eds_segments[v.segment - 1].size() > 1);
}

bool isNVertex(Vertex v, const eds_matrix &eds_segments) {
    return !isLayerVertex(v, eds_segments) && !isJVertex(v, eds_segments);
}

bool hasPredecessorVertex(Vertex v) { return v.segment != 0 && v.index != 0; }

Vertex getPredecessorVertex(const eds_matrix &eds_segments, Vertex v,
                            int predecessor_layer = 0) {
    assert(hasPredecessorVertex(v));
    if (v.index != 0) {
        return {v.segment, v.layer, v.index - 1};
    }
    return {v.segment - 1, predecessor_layer,
            static_cast<int>(
                eds_segments[v.segment - 1][predecessor_layer].size() - 1)};
}

int getScore(score_matrix &scores, Vertex v, bool selected,
             path_continuation layer = I) {
    return scores[v.segment][v.layer][v.index][selected][layer];
}

int getWeight(const weight_matrix &weights, Vertex v) {
    return weights[v.segment][v.layer][v.index];
}

void setScore(score_matrix &scores, int score, Vertex v, bool selected,
              path_continuation layer = I) {
    scores[v.segment][v.layer][v.index][selected][layer] = score;
}

// vector<vector<vector<vector<vector<int>>>>>
// [segment][layer][index][{SELECTED, !SELECTED}][I, E]
score_matrix initScoreMatrix(const weight_matrix &weights) {
    score_matrix scores;
    scores.resize(weights.size());
    for (int segment = 0; segment < weights.size(); segment++) {
        scores[segment].resize(weights[segment].size());
        for (int layer = 0; layer < weights[segment].size(); layer++) {
            scores[segment][layer].resize(weights[segment][layer].size());
            for (int index = 0; index < weights[segment][layer].size();
                 index++) {
                scores[segment][layer][index].resize(2);
                scores[segment][layer][index][0].resize(2);
                scores[segment][layer][index][1].resize(2);
            }
        }
    }
    return scores;
}

void findMaxScoringPaths(const eds_matrix &eds_segments,
                         const weight_matrix &weights, score_matrix &scores,
                         path_matrix &paths,
                         belonging_path_matrix &belonging_path_starts,
                         int penalty = -1) {
    for (int segment = 0; segment < eds_segments.size(); segment++) {
        for (int layer = 0; layer < eds_segments[segment].size(); layer++) {
            for (int index = 0; index < eds_segments[segment][layer].size();
                 index++) {
                Vertex a{segment, layer, index};
                int weight_a = getWeight(weights, a);

                // First vertex of the graph.
                if (!hasPredecessorVertex(a)) {
                    // W(a, 1) = w(a) - x
                    int score_a_1 = weight_a - penalty;
                    setScore(scores, score_a_1, a, SELECTED);
                    // W(a, 0) = max{0, W(a, 1)}
                    setScore(scores, max(0, score_a_1), a, !SELECTED);
                    continue;
                }

                // N vertex.
                if (isNVertex(a, eds_segments)) {
                    assert(layer == 0);
                    // W(a, 1) = w(a) + max{W(p, 0) - x, W(p, 1)}
                    Vertex p = getPredecessorVertex(eds_segments, a);
                    int score_p_0 = getScore(scores, p, !SELECTED);
                    int score_p_1 = getScore(scores, p, SELECTED);
                    int score_a_1 =
                        weight_a + max(score_p_0 - penalty, score_p_1);
                    setScore(scores, score_a_1, a, SELECTED);

                    // W(a, 0) = max{W(p, 0), W(a, 1)}
                    setScore(scores, max(score_p_0, score_a_1), a, !SELECTED);
                }

                // 1st (lowest) layer vertex.
                else if (isFirstLayerVertex(a, eds_segments)) {
                    // 1_first vertex.
                    if (isVertexFirstOnLayer(a, eds_segments)) {
                        Vertex p = getPredecessorVertex(eds_segments, a);
                        int score_p_0 = getScore(scores, p, !SELECTED);
                        int score_p_1 = getScore(scores, p, SELECTED);
                        // W(a, 1, I) = w(a) + max{W(p, 0) - x, W(p, 1)}
                        int score_a_1_I =
                            weight_a + max(score_p_0 - penalty, score_p_1);
                        setScore(scores, score_a_1_I, a, SELECTED, I);

                        // W(a, 0, I) = max{W(p, 0), W(a, 1, I)}
                        setScore(scores, max(score_p_0, score_a_1_I), a,
                                 !SELECTED, I);

                        // W(a, 1, E) = w(a) + W(p, 1) - x
                        int score_a_1_E = weight_a + score_p_1 - penalty;
                        setScore(scores, score_a_1_E, a, SELECTED, E);

                        // W(a, 0, E) = max{W(p, 1), W(a, 1, E)}
                        setScore(scores, max(score_p_1, score_a_1_E), a,
                                 !SELECTED, E);
                    }
                    // 1_later vertex.
                    else {
                        Vertex p = getPredecessorVertex(eds_segments, a);
                        int score_p_0_I = getScore(scores, p, !SELECTED, I);
                        // W(a, 1, I) = w(a) + max{W(p, 0, I) - x, W(p, 1, I)}
                        int score_a_1_I =
                            weight_a + max(score_p_0_I - penalty,
                                           getScore(scores, p, SELECTED, I));
                        setScore(scores, score_a_1_I, a, SELECTED, I);

                        // W(a, 0, I) = max{W(p, 0, I), W(a, 1, I)}
                        setScore(scores, max(score_p_0_I, score_a_1_I), a,
                                 !SELECTED, I);

                        // W(a, 1, E) = w(a) + max{W(p, 0, E) - x, W(p, 1, E)}
                        int score_p_0_E = getScore(scores, p, !SELECTED, E);
                        int score_a_1_E =
                            weight_a + max(score_p_0_E - penalty,
                                           getScore(scores, p, SELECTED, E));
                        setScore(scores, score_a_1_E, a, SELECTED, E);

                        // W(a, 0, E) = max{W(p, 0, E), W(a, 1, E)}
                        setScore(scores, max(score_p_0_E, score_a_1_E), a,
                                 !SELECTED, E);
                    }
                }
                // L in {2, ..., n} vertex.
                else if (isLayerVertex(a, eds_segments)) {
                    // L_first vertex.
                    if (isVertexFirstOnLayer(a, eds_segments)) {
                        // W(a, 1, I) = w(a)
                        int score_a_1_I = weight_a;
                        setScore(scores, score_a_1_I, a, SELECTED, I);

                        // W(a, 0, I) = W(a, 1, I)
                        setScore(scores, score_a_1_I, a, !SELECTED, I);

                        // W(a, 1, E) = w(a) - x
                        int score_a_1_E = weight_a - penalty;
                        setScore(scores, score_a_1_E, a, SELECTED, E);

                        // W(a, 0, E) = max{0, W(a, 1, E)}
                        setScore(scores, max(0, score_a_1_E), a, !SELECTED, E);
                    }
                    // L_later vertex.
                    else {
                        Vertex p = getPredecessorVertex(eds_segments, a);
                        int score_p_0_I = getScore(scores, p, !SELECTED, I);
                        // W(a, 1, I) = w(a) + max{W(p, 0, I) - x, W(p, 1, I)}
                        int score_a_1_I =
                            weight_a + max(score_p_0_I - penalty,
                                           getScore(scores, p, SELECTED, I));
                        setScore(scores, score_a_1_I, a, SELECTED, I);

                        // W(a, 0, I) = max{W(p, 0, I), W(a, 1, I)}
                        setScore(scores, max(score_p_0_I, score_a_1_I), a,
                                 !SELECTED, I);

                        // W(a, 1, E) = w(a) + max{W(p, 0, E) - x, W(p, 1, E)}
                        int score_p_0_E = getScore(scores, p, !SELECTED, E);
                        int score_a_1_E =
                            weight_a + max(score_p_0_E - penalty,
                                           getScore(scores, p, SELECTED, E));
                        setScore(scores, score_a_1_E, a, SELECTED, E);

                        // W(a, 0, E) = max{W(p, 0, E), W(a, 1, E)}
                        setScore(scores, max(score_p_0_E, score_a_1_E), a,
                                 !SELECTED, E);
                    }
                }
                // J vertex.
                else if (isJVertex(a, eds_segments)) {
                    // W(a, 1) = max{ W(p1, 0, I) + W(p2, 0, E) +...+ W(pb, 0,
                    // E) - x,
                    //                                  ...
                    //                W(p1, 0, E) + W(p2, 0, E) +...+ W(pb, 0,
                    //                I) - x, W(p1, 1, I) + W(p2, 0, E) +...+
                    //                W(pb, 0, E),
                    //                                  ...
                    //                W(p1, 1, E) + W(p2, 0, E) +...+ W(pb, 0,
                    //                I),
                    //                               .  .  .  .
                    //                W(p1, 0, I) + W(p2, 0, E) +...+ W(pb, 1,
                    //                E),
                    //                                  ...
                    //                W(p1, 0, E) + W(p2, 0, E) +...+ W(pb, 1,
                    //                I) }

                    int score_a_1 = INT_MIN;
                    // To avoid n^3 runtime complexity of W(a, 1), calculate a
                    // base score that is going to be modified: 
                    // W(p1, 0, I) +  W(p2, 0, E) + ... + W(pb, 0, E)
                    int base_score = 0;
                    // Add W(p1, 0, I).
                    base_score += getScore(
                        scores, getPredecessorVertex(eds_segments, a, 0),
                        !SELECTED, I);
                    // Add scores from remaining predecessors.
                    for (int i = 1; i < eds_segments[segment].size(); i++) {
                        Vertex p_i = getPredecessorVertex(eds_segments, a, i);
                        base_score += getScore(scores, p_i, !SELECTED, E);
                    }
                }
            }
        }
    }
}