#ifndef MAXSCOREPATH_UTILITY_FUNC_HEADER
#define MAXSCOREPATH_UTILITY_FUNC_HEADER

#include <string>
#include <vector>

using namespace std;

// This file contains a dynamic programming algorithm for the maximum-score
// disjoint paths problem in n-layered bubble graphs and utility functions.
//
// Elastic-degenerate strings, EDS for short, can be used to represent
// pan-genomes. See more: https://doi.org/10.48550/arXiv.1610.08111.
//
// EDS text is a sequence of segments with the following format:
// `{A,CCC,GA,}ATA{A,}{GG,CC}A`. Braces determine the start and end of
// non-deterministic segments. Inside the braces we can see comma-separated
// segment variants. In non-deterministic segments, a trailing comma or a comma
// without a preceding segment variant encodes an empty segment variant.
// See SOPanG, a tool used for pattern matching in EDS texts:
// https://doi.org/10.1093/bioinformatics/bty506.
//
// To make EDS texts fully correspond to n-layered bubble graphs, we need to
// make small modifications. Each character of the EDS text corresponds to a
// vertex. The non-deterministic segments inside braces correspond to bubble
// layers. In an n-layered bubble graph, all bubbles have a start and end
// vertex. To fulfill this condition, a special separation character is inserted
// in between adjacent non-deterministic segments and also at the beginning and
// end of the text. Furthermore, the start and end vertex of a bubble are not
// adjacent, therefore, empty segment variants are also replaced by a the
// special separation character. The separation character has no biological
// significance, therefore, it has weight 0.

// Store the EDS text inside a 2D matrix of strings. The character on
// `eds_matrix[segment][layer][index]` corresponds to the character in segment
// `segment` on layer `layer` on position `index`.
typedef vector<vector<string>> eds_matrix;
// Each character from the `eds_matrix` has a weight, given with the following
// matrix.
typedef vector<vector<vector<int>>> weight_matrix;

// Separation character between adjacent non-deterministic segments, start and
// end of the EDS text, and for empty segment variants in non-deterministic
// segments.
#define EMPTY_STR '_'

// Reads a line containing EDS text from `file_path` and returns it as a string.
string readEDSFile(const string &file_path);

// Store the EDS text in an `eds_matrix`.
eds_matrix EDSToMatrix(const string &EDS);

// Return a `weight matrix` of weights to the given `eds_matrix` based on the GC
// content. Assign scores:
// - bases G and C get score `match`;
// - bases A and T (or N) get score `non_match`;
// - the separation character `EMPTY_STR` gets score 0.
vector<vector<vector<int>>> getGCContentWeights(const eds_matrix &eds_segments,
                                                int match = 1,
                                                int non_match = -1);

// Each character in the EDS text represents a vertex.
struct Vertex {
    Vertex() : segment(-1), layer(-1), index(-1) {}
    Vertex(int segment, int layer, int index)
        : segment(segment), layer(layer), index(index) {}

    explicit operator bool() {
        return !(segment == -1 && layer == -1 && index == -1);
    }

    int segment;
    int layer;
    int index;
};
bool operator==(const Vertex &a, const Vertex &b);
bool operator!=(const Vertex &a, const Vertex &b);

// DP score calculation.
// For each vertex we consider the score when it is selected or not selected.
// Index `SURELY_SELECTED = 1` corresponds to the score when the vertex is
// selected. Its negation `!SURELY_SELECTED = 0` index contains the maximum
// score of selecting or not selecting the vertex.
#define SURELY_SELECTED true
// The path from the start vertex of a bubble can continue on different layers.
// For each layer vertex the path continuation on its, I(dentical), layer or a
// different, E(lse), layer is considered.
typedef bool path_continuation;
#define I 0
#define E 1
// The DP algorithm uses the following score vector to store the score for each
// vertex if the vertex is/is not selected and the path continues/does not
// continue on it's layer (in case of layer vertices). Each dimension has the
// following meaning:
// score_matrix[segment][layer][index][(!)SURELY_SELECTED][path_continuation].
typedef vector<vector<vector<vector<vector<int>>>>> score_matrix;

// Helper functions for the `Vertex`.
// Return the last vertex of the n-layered bubble graph.
Vertex getLastVertex(const eds_matrix &eds_segments);

// Return true if the vertex is on one of the layers.
bool isLayerVertex(Vertex v, const eds_matrix &eds_segments);

// Returns true if the vertex is on the first layer.
bool isFirstLayerVertex(Vertex v, const eds_matrix &eds_segments);

// Returns true if the vertex is the first on any layer.
bool isVertexFirstOnLayer(Vertex v, const eds_matrix &eds_segments);

// Returns whether a vertex is a J vertex, i.e. end vertex of a bubble.
bool isJVertex(Vertex v, const eds_matrix &eds_segments);

// Returns whether a vertex is a "normal" vertex, i.e. not J or layer vertex.
bool isNVertex(Vertex v, const eds_matrix &eds_segments);

// Returns true if the vertex has at least one predecessor vertex.
bool hasPredecessorVertex(Vertex v);

// Returns the predecessor vertex of vertex `v`. For J vertices, the layer can
// be specified in `predecessor_layer`.
Vertex getPredecessorVertex(const eds_matrix &eds_segments, Vertex v,
                            int predecessor_layer = 0);

int getScore(score_matrix &scores, Vertex v, bool surely_selected,
             path_continuation path_goes = I);

// The majority of rules in the DP algorithm have only two choices for score
// calculation.
#define FIRST 0
#define SECOND 1

int getChoice(score_matrix &choices, Vertex v, bool surely_selected,
              path_continuation path_goes = I);

// Returns a tuple containing the maximum score and wether this was the first or
// second parameter.
pair<int, int> max_score(int first_score, int second_score);

void setScoreAndChoice(score_matrix &scores, score_matrix &choices,
                       pair<int, int> score_choice, Vertex v, bool selected,
                       path_continuation layer);

int getWeight(const weight_matrix &weights, Vertex v);

// Initializes the score matrix to its known size.
score_matrix initScoreMatrix(const weight_matrix &weights);

// Dynamic programming algorithm to maximize the score of selected disjoint
// paths while each path incurs a penalty. Fills the `scores` score matrix and
// returns the best score from the last vertex which is the maximal score of
// selecting disjoint paths. Fills also the `choices` matrix which contains
// which previous score was used when calculating the current score.
int findMaxScoringPaths(const eds_matrix &eds_segments,
                        const weight_matrix &weights, score_matrix &scores,
                        score_matrix &choices, int penalty);

// Based on the `choices` that were filled by `findMaxScoringPaths()`, return
// all the selected paths.
vector<vector<Vertex>> getPaths(const eds_matrix &eds_segments,
                                score_matrix &scores, score_matrix &choices);

// Prints out the paths that were found by `getPaths()`.
void printPaths(vector<vector<Vertex>> paths);

#endif