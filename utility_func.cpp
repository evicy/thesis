#include "utility_func.hpp"

#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

// Returns the one line elastic-degenerate string file as a string.
// Adds a letter N to the start and also to the end of the string to make sure
// the elastic-degenerate string equals with the n-layered bubble graph
// definition. These two letters should get score 0.
string readEDSFile(const string &file_path) {
    ifstream input_stream(file_path);

    if (!input_stream.good()) {
        cout << "Reading of file " << file_path << " failed." << endl;
        return "";
    }
    assert(input_stream.good());
    return EMPTY_STR +
           static_cast<stringstream const &>(stringstream()
                                             << input_stream.rdbuf())
               .str() +
           EMPTY_STR;
}

// Store the elastic degenerate string `text` in a 2D matrix.
// eds_matrix[s][l]
eds_matrix EDSToMatrix(const string &EDS) {
    eds_matrix eds_segments;
    vector<string> current_segment;
    string current_string;
    bool in_nondet_segment = false;

    for (int i = 0; i < EDS.length(); i++) {
        char c = EDS[i];
        // Inside a string or segment.
        if (c != '{' && c != '}') {
            // Commas can be only inside segments.
            assert(in_nondet_segment || c != ',');

            if (c == ',') {
                // Empty layers are are denoted by a vertex with value "_".
                if (current_string.empty()) {
                    current_string += EMPTY_STR;
                }
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
            // Starts after another non-deterministic
            // segment: separate them with with an empty deterministic segment.
            if (i > 0 && EDS[i - 1] == '}') {
                eds_segments.emplace_back(vector<string>{string(1, EMPTY_STR)});
            }
            // Starts after a deterministic segment.
            if (!current_string.empty()) {
                eds_segments.emplace_back(vector<string>{move(current_string)});
                current_string.clear();
            }
        }
        // End of a segment.
        else {
            assert(c == '}' && in_nondet_segment && !current_segment.empty());
            // Empty layers are are denoted by a vertex with value "_".
            if (current_string.empty()) {
                current_string += EMPTY_STR;
            }
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
vector<vector<vector<int>>> getGCContentWeights(const eds_matrix &eds_segments,
                                                int match, int non_match) {
    vector<vector<vector<int>>> weights;
    weights.reserve(eds_segments.size());
    for (const auto &segment : eds_segments) {
        vector<vector<int>> w_segment;
        w_segment.reserve(segment.size());
        for (const auto &str : segment) {
            vector<int> w_str(str.size(), 0);
            for (int i = 0; i < str.size(); i++) {
                if (str[i] == EMPTY_STR) {
                    w_str[i] = 0;
                } else {
                    w_str[i] =
                        (str[i] == 'G' || str[i] == 'C') ? match : non_match;
                }
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

std::ostream &operator<<(std::ostream &os, Vertex const &v) { 
    return os << "(" << v.segment << "," << v.layer << "," << v.index << ")";
}

// Helper functions for Vertex type.
Vertex getLastVertex(const eds_matrix &eds_segments) {
    int last_segment = eds_segments.size() - 1;
    int last_index = eds_segments[last_segment][0].size() - 1;
    return Vertex(last_segment, 0, last_index);
}

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

bool hasPredecessorVertex(Vertex v) {
    return !(v.segment == 0 && v.index == 0);
}

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
             path_continuation path_goes = I) {
    return scores[v.segment][v.layer][v.index][selected][path_goes];
}

int getChoice(score_matrix &choices, Vertex v, bool selected,
              path_continuation path_goes = I) {
    return choices[v.segment][v.layer][v.index][selected][path_goes];
}

int getWeight(const weight_matrix &weights, Vertex v) {
    return weights[v.segment][v.layer][v.index];
}

pair<int, int> max_score(int first_score, int second_score) {
    return first_score > second_score ? make_pair(first_score, FIRST)
                                      : make_pair(second_score, SECOND);
}

void setScoreAndChoice(score_matrix &scores, score_matrix &choices,
                       pair<int, int> score_choice, Vertex v, bool selected,
                       path_continuation layer = I) {
    scores[v.segment][v.layer][v.index][selected][layer] = score_choice.first;
    choices[v.segment][v.layer][v.index][selected][layer] = score_choice.second;
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

int findMaxScoringPaths(const eds_matrix &eds_segments,
                        const weight_matrix &weights, score_matrix &scores,
                        score_matrix &choices, int penalty = -1) {
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
                    setScoreAndChoice(scores, choices,
                                      make_pair(score_a_1, FIRST), a, SELECTED);
                    // W(a, 0) = max{0, W(a, 1)}
                    setScoreAndChoice(scores, choices, max_score(0, score_a_1),
                                      a, !SELECTED);
                }
                // N vertex.
                else if (isNVertex(a, eds_segments)) {
                    assert(layer == 0);
                    // W(a, 1) = w(a) + max{W(p, 0) - x, W(p, 1)}
                    Vertex p = getPredecessorVertex(eds_segments, a);
                    int score_p_0 = getScore(scores, p, !SELECTED);
                    int score_p_1 = getScore(scores, p, SELECTED);
                    setScoreAndChoice(scores, choices,
                                      max_score(weight_a + score_p_0 - penalty,
                                                weight_a + score_p_1),
                                      a, SELECTED);

                    // W(a, 0) = max{W(p, 0), W(a, 1)}
                    setScoreAndChoice(
                        scores, choices,
                        max_score(score_p_0, getScore(scores, a, SELECTED)), a,
                        !SELECTED);
                }
                // 1st (lowest) layer vertex.
                else if (isFirstLayerVertex(a, eds_segments)) {
                    // 1_first vertex.
                    if (isVertexFirstOnLayer(a, eds_segments)) {
                        Vertex p = getPredecessorVertex(eds_segments, a);
                        int score_p_0 = getScore(scores, p, !SELECTED);
                        int score_p_1 = getScore(scores, p, SELECTED);
                        // W(a, 1, I) = w(a) + max{W(p, 0) - x, W(p, 1)}
                        setScoreAndChoice(
                            scores, choices,
                            max_score(weight_a + score_p_0 - penalty,
                                      weight_a + score_p_1),
                            a, SELECTED, I);

                        // W(a, 0, I) = max{W(p, 0), W(a, 1, I)}
                        setScoreAndChoice(
                            scores, choices,
                            max_score(score_p_0,
                                      getScore(scores, a, SELECTED, I)),
                            a, !SELECTED, I);

                        // W(a, 1, E) = w(a) + W(p, 1) - x
                        int score_a_1_E = weight_a + score_p_1 - penalty;
                        setScoreAndChoice(scores, choices,
                                          make_pair(score_a_1_E, FIRST), a,
                                          SELECTED, E);

                        // W(a, 0, E) = max{W(p, 1), W(a, 1, E)}
                        setScoreAndChoice(scores, choices,
                                          max_score(score_p_1, score_a_1_E), a,
                                          !SELECTED, E);
                    }
                    // 1_later vertex.
                    else {
                        Vertex p = getPredecessorVertex(eds_segments, a);
                        int score_p_0_I = getScore(scores, p, !SELECTED, I);
                        int score_p_1_I = getScore(scores, p, SELECTED, I);
                        // W(a, 1, I) = w(a) + max{W(p, 0, I) - x, W(p, 1, I)}
                        setScoreAndChoice(
                            scores, choices,
                            max_score(weight_a + score_p_0_I - penalty,
                                      weight_a + score_p_1_I),
                            a, SELECTED, I);

                        // W(a, 0, I) = max{W(p, 0, I), W(a, 1, I)}
                        setScoreAndChoice(
                            scores, choices,
                            max_score(score_p_0_I,
                                      getScore(scores, a, SELECTED, I)),
                            a, !SELECTED, I);

                        // W(a, 1, E) = w(a) + max{W(p, 0, E) - x, W(p, 1, E)}
                        int score_p_0_E = getScore(scores, p, !SELECTED, E);
                        int score_p_1_E = getScore(scores, p, SELECTED, E);
                        setScoreAndChoice(
                            scores, choices,
                            max_score(weight_a + score_p_0_E - penalty,
                                      weight_a + score_p_1_E),
                            a, SELECTED, E);

                        // W(a, 0, E) = max{W(p, 0, E), W(a, 1, E)}
                        setScoreAndChoice(
                            scores, choices,
                            max_score(score_p_0_E,
                                      getScore(scores, a, SELECTED, E)),
                            a, !SELECTED, E);
                    }
                }
                // L in {2, ..., n} vertex.
                else if (isLayerVertex(a, eds_segments)) {
                    // L_first vertex.
                    if (isVertexFirstOnLayer(a, eds_segments)) {
                        // W(a, 1, I) = w(a)
                        int score_a_1_I = weight_a;
                        setScoreAndChoice(scores, choices,
                                          make_pair(score_a_1_I, FIRST), a,
                                          SELECTED, I);

                        // W(a, 0, I) = W(a, 1, I)
                        setScoreAndChoice(scores, choices,
                                          make_pair(score_a_1_I, FIRST), a,
                                          !SELECTED, I);

                        // W(a, 1, E) = w(a) - x
                        int score_a_1_E = weight_a - penalty;
                        setScoreAndChoice(scores, choices,
                                          make_pair(score_a_1_E, FIRST), a,
                                          SELECTED, E);

                        // W(a, 0, E) = max{0, W(a, 1, E)}
                        setScoreAndChoice(scores, choices,
                                          max_score(0, score_a_1_E), a,
                                          !SELECTED, E);
                    }
                    // L_later vertex.
                    else {
                        Vertex p = getPredecessorVertex(eds_segments, a);
                        int score_p_0_I = getScore(scores, p, !SELECTED, I);
                        int score_p_1_I = getScore(scores, p, SELECTED, I);
                        // W(a, 1, I) = w(a) + max{W(p, 0, I) - x, W(p, 1, I)}
                        setScoreAndChoice(
                            scores, choices,
                            max_score(weight_a + score_p_0_I - penalty,
                                      weight_a + score_p_1_I),
                            a, SELECTED, I);

                        // W(a, 0, I) = max{W(p, 0, I), W(a, 1, I)}
                        setScoreAndChoice(
                            scores, choices,
                            max_score(score_p_0_I,
                                      getScore(scores, a, SELECTED, I)),
                            a, !SELECTED, I);

                        // W(a, 1, E) = w(a) + max{W(p, 0, E) - x, W(p, 1, E)}
                        int score_p_0_E = getScore(scores, p, !SELECTED, E);
                        int score_p_1_E = getScore(scores, p, SELECTED, E);
                        setScoreAndChoice(
                            scores, choices,
                            max_score(weight_a + score_p_0_E - penalty,
                                      weight_a + score_p_1_E),
                            a, SELECTED, E);

                        // W(a, 0, E) = max{W(p, 0, E), W(a, 1, E)}
                        setScoreAndChoice(
                            scores, choices,
                            max_score(score_p_0_E,
                                      getScore(scores, a, SELECTED, E)),
                            a, !SELECTED, E);
                    }
                }
                // J vertex.
                else if (isJVertex(a, eds_segments)) {
                    // W(a, 1) = w(a) + max{
                    //      W(p1, 0, I) + W(p2, 0, E) +...+ W(pb, 0, E) - x,
                    //                         ...
                    //______W(p1, 0, E) + W(p2, 0, E) +...+ W(pb, 0, I) - x,
                    //      W(p1, 1, I) + W(p2, 0, E) +...+ W(pb, 0, E),
                    //                         ...
                    //______W(p1, 1, E) + W(p2, 0, E) +...+ W(pb, 0, I),
                    //______               .  .  .  .
                    //      W(p1, 0, I) + W(p2, 0, E) +...+ W(pb, 1, E),
                    //                         ...
                    //      W(p1, 0, E) + W(p2, 0, E) +...+ W(pb, 1, I) }
                    // where b is the number of layers in the current bubble.
                    int num_preds = eds_segments[segment - 1].size();
                    // Get all predecessors.
                    vector<Vertex> a_preds(num_preds);
                    for (int i = 0; i < num_preds; i++) {
                        a_preds[i] = getPredecessorVertex(eds_segments, a, i);
                    }
                    // To avoid n^3 runtime complexity of W(a, 1), calculate a
                    // base score that is going to be modified:
                    // w(a) + W(p1, 0, E) +  W(p2, 0, E) + ... + W(pb, 0, E)
                    int base_score = 0;
                    for (const Vertex &p : a_preds) {
                        base_score += getScore(scores, p, !SELECTED, E);
                    }
                    // Choose max sum for W(a, 1) based on the rules.
                    // 1st group in max:
                    int choice_a_1 = -1;
                    int score_a_1 = INT_MIN;
                    for (int path_I = 0; path_I < num_preds; path_I++) {
                        int current_score =
                            base_score -
                            getScore(scores, a_preds[path_I], !SELECTED, E) +
                            getScore(scores, a_preds[path_I], !SELECTED, I) -
                            penalty;
                        if (score_a_1 < current_score) {
                            score_a_1 = current_score;
                            choice_a_1 = path_I;
                        }
                    }
                    // Remaining groups:
                    for (int selected_pred = 0; selected_pred < num_preds;
                         selected_pred++) {
                        int temp_base_score =
                            base_score -
                            getScore(scores, a_preds[selected_pred], !SELECTED,
                                     E) +
                            getScore(scores, a_preds[selected_pred], SELECTED,
                                     E);
                        for (int path_I = 0; path_I < num_preds; path_I++) {
                            int current_score =
                                temp_base_score -
                                getScore(scores, a_preds[path_I],
                                         selected_pred == path_I ? SELECTED
                                                                 : !SELECTED,
                                         E) +
                                getScore(scores, a_preds[path_I],
                                         selected_pred == path_I ? SELECTED
                                                                 : !SELECTED,
                                         I);
                            if (score_a_1 < current_score) {
                                score_a_1 = current_score;
                                choice_a_1 =
                                    (selected_pred + 1) * num_preds + path_I;
                            }
                        }
                    }
                    setScoreAndChoice(
                        scores, choices,
                        make_pair(weight_a + score_a_1, choice_a_1), a,
                        SELECTED);

                    // W(a, 0) = max{
                    //      W(p1, 0, I) + W(p2, 0, E) +...+ W(pb, 0, E),
                    //                         ...
                    //______W(p1, 0, E) + W(p2, 0, E) +...+ W(pb, 0, I),
                    //      W(a, 1)
                    // where b is the number of layers in the current bubble.
                    int score_a_0 = INT_MIN;
                    int choice_a_0 = 0;
                    for (int path_I = 0; path_I < num_preds; path_I++) {
                        int current_score =
                            base_score -
                            getScore(scores, a_preds[path_I], !SELECTED, E) +
                            getScore(scores, a_preds[path_I], !SELECTED, I);
                        if (score_a_0 < current_score) {
                            score_a_0 = current_score;
                            choice_a_0 = path_I;
                        }
                    }
                    if (score_a_0 < score_a_1) {
                        score_a_0 = score_a_1;
                        // For non-ambiguous decoding to get the paths, store
                        // choice `choice_a_1` here "moved" by `num_preds`: i.e.
                        // if choice_a_0 >= num_preds, then decode choice_a_0 -
                        // num_preds based on W(a, 1) rules.
                        choice_a_0 = num_preds + choice_a_1;
                    }
                    setScoreAndChoice(scores, choices,
                                      make_pair(score_a_0, choice_a_0), a,
                                      !SELECTED);
                } else {
                    assert(false);
                }
            }
        }
    }

    // Get the max score from the last vertex of the graph. The last vertex has
    // value 0, therefore, it is unnecessary to select it
    Vertex last = getLastVertex(eds_segments);
    assert(isNVertex(last, eds_segments) || isJVertex(last, eds_segments));
    const auto &last_data = scores[last.segment][last.layer][last.index];
    return max(last_data[!SELECTED][I], last_data[!SELECTED][E]);
}

vector<vector<Vertex>> getPaths(const eds_matrix &eds_segments,
                                score_matrix &scores, score_matrix &choices) {
    vector<vector<Vertex>> paths;
    // The last vertex was synthetically added to the pangenome-graph and has
    // value 0. Therefore, it is unnecessary to select it.
    Vertex a = getLastVertex(eds_segments);
    bool was_a_selected = !SELECTED;
    vector<Vertex> current_path;
    while (hasPredecessorVertex(a)) {
        // N vertex.
        if (isNVertex(a, eds_segments)) {
            // W(a, 1) = w(a) + max{W(p, 0) - x, W(p, 1)}
            // W(a, 0) = max{W(p, 0), W(a, 1)}
            // Regardless whether `a` was selected: the first choice corresponds
            // to not selecting the previos vertex `p`.
            Vertex p = getPredecessorVertex(eds_segments, a);
            int choice = getChoice(choices, a, was_a_selected);
            if (choice == FIRST) {
                was_a_selected = !SELECTED;
                if (!current_path.empty()) {
                    paths.emplace_back(current_path);
                    current_path.clear();
                }
            } else {
                current_path.emplace_back(p);
                was_a_selected = SELECTED;
            }
            a = p;
        }
        // J vertex.
        // Handle the full bubble, `a` is the start vertex of the bubble after
        // this code part.
        else if (isJVertex(a, eds_segments)) {
            Vertex j = a;
            int choice = getChoice(choices, a, was_a_selected);
            // Get all predecessors.
            int num_preds = eds_segments[a.segment - 1].size();
            vector<Vertex> j_preds(num_preds);
            for (int layer = 0; layer < num_preds; layer++) {
                j_preds[layer] = getPredecessorVertex(eds_segments, a, layer);
            }

            // W(a, 1) = w(a) + max{
            //      W(p1, 0, I) + W(p2, 0, E) +...+ W(pb, 0, E) - x,
            //                         ...
            //______W(p1, 0, E) + W(p2, 0, E) +...+ W(pb, 0, I) - x,
            //      W(p1, 1, I) + W(p2, 0, E) +...+ W(pb, 0, E),
            //                         ...
            //______W(p1, 1, E) + W(p2, 0, E) +...+ W(pb, 0, I),
            //______               .  .  .  .
            //      W(p1, 0, I) + W(p2, 0, E) +...+ W(pb, 1, E),
            //                         ...
            //      W(p1, 0, E) + W(p2, 0, E) +...+ W(pb, 1, I) }
            // where b is the number of layers in the current bubble.
            //______________________________________________________
            // W(a, 0) = max{
            //      W(p1, 0, I) + W(p2, 0, E) +...+ W(pb, 0, E),
            //                         ...
            //______W(p1, 0, E) + W(p2, 0, E) +...+ W(pb, 0, I),
            //      W(a, 1)
            // where b is the number of layers in the current bubble.

            // From the saved choice, we can get the line of the rule which
            // encodes the layer, L, where the path continues from the start
            // vertex. Note that on every layer for every predecessor p we can
            // determine whether it was selected or not by comparing W(p, 0, _)
            // and W(p, 1, _).
            int rule_line = choice % num_preds;
            bool current_path_was_continued = false;
            vector<Vertex> path_through_start_vertex;
            // Iterate through the layers and store the paths on each layer,
            // handle the first layer last. Continue current path with one of
            // the predecessors.
            for (int layer = num_preds - 1; layer >= 0; layer--) {
                vector<Vertex> layer_path;
                int path_cont_layer = rule_line == path_cont_layer ? I : E;
                bool is_j_selected =
                    getScore(scores, j_preds[layer], SELECTED,
                             path_cont_layer) < getScore(scores, j_preds[layer],
                                                         !SELECTED,
                                                         path_cont_layer)
                        ? !SELECTED
                        : SELECTED;

                if (is_j_selected) {
                    layer_path.emplace_back(j_preds[layer]);
                }

                Vertex v = j_preds[layer];
                bool is_v_selected = is_j_selected;
                while (isLayerVertex(v, eds_segments)) {
                    int choice_v =
                        getChoice(choices, v, is_v_selected, path_cont_layer);
                    Vertex v_pred = getPredecessorVertex(eds_segments, v);
                    
                    // 1_first vertex: this is the last and vertex to be
                    // processed in the J vertex block.
                    if (isFirstLayerVertex(v, eds_segments) &&
                        isVertexFirstOnLayer(v, eds_segments)) {
                        assert(layer == 0);

                        // Set `a` to point to the start vertex of the bubble;
                        a = v_pred;
                        if (path_cont_layer == I) {
                            was_a_selected =
                                choice_v == FIRST ? !SELECTED : SELECTED;
                        } else {
                            was_a_selected = SELECTED;
                        }
                    }

                    // First vertex on any layer.
                    if (isVertexFirstOnLayer(v, eds_segments)) {
                        if (path_cont_layer == I &&
                            ((isFirstLayerVertex(v, eds_segments) &&
                              was_a_selected) ||
                             !isFirstLayerVertex(v, eds_segments))) {
                            layer_path.emplace_back(v_pred);
                            if (!current_path.empty()) {
                                // Maybe `layer_path` should continue
                                // `current_path`. If the bubble's end vertex's
                                // (the processed J vertex's) predecessor
                                if (current_path.back() == j &&
                                    layer_path[0] == j_preds[layer]) {
                                    current_path.insert(current_path.end(),
                                                        layer_path.begin(),
                                                        layer_path.end());
                                } else {
                                    paths.emplace_back(current_path);
                                    current_path = layer_path;
                                }
                            } else {
                                current_path = layer_path;
                            }
                        } else if (!layer_path.empty()) {
                            paths.emplace_back(layer_path);
                            layer_path.clear();
                        }
                        v = v_pred;
                    }
                    // Later vertex on any layer.
                    else {
                        if (choice == FIRST) {
                            is_v_selected = !SELECTED;
                            if (!layer_path.empty()) {
                                // Try to continue current_path if possible.
                                if (!current_path.empty() &&
                                    current_path.back() == j &&
                                    layer_path[0] == j_preds[layer]) {
                                    current_path.insert(current_path.end(),
                                                        layer_path.begin(),
                                                        layer_path.end());
                                    paths.emplace_back(current_path);
                                    current_path.clear();
                                } else {
                                    paths.emplace_back(layer_path);
                                    layer_path.clear();
                                }
                            }
                        } else {
                            layer_path.emplace_back(v_pred);
                            is_v_selected = SELECTED;
                        }
                        v = v_pred;
                    }
                }
            }
        } 
        else {
            assert(false);
        }
    }
    if (!current_path.empty()) {
        paths.emplace_back(current_path);
    }
    return paths;
}

void printPaths(vector<vector<Vertex>> paths) {
    for (const auto &path : paths) {
        for (const auto &vertex : path) {
            cout << vertex; 
        }
        cout << endl;
    }
}