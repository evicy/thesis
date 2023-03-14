#ifndef MAXSCOREPATH_UTILITY_FUNC_HEADER
#define MAXSCOREPATH_UTILITY_FUNC_HEADER

#include <boost/algorithm/string/trim.hpp>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

typedef vector<vector<string>> eds_matrix;

// Returns the one line elastic degenerate string file as a string.
string readEDSFile(const string &file_path);

// Store the elastic degenerate string `text` in a 2D matrix.
eds_matrix EDSToMatrix(const string &EDS);

// Assign score based on GC content:
// - bases G and C get score 1
// - empty non-deterministic segment parts get score 0
// - every other character gets score 0.
vector<vector<vector<int>>> getGCContentWeights(const eds_matrix &eds_segments);

// Find maximum-scoring paths.

struct Vertex {
    int segment;
    int layer;
    int index;

    Vertex() : segment(-1), layer(-1), index(-1) {}
    Vertex(int segment, int layer, int index)
        : segment(segment), layer(layer), index(index) {}

    explicit operator bool() {
        return !(segment == -1 && layer == -1 && index == -1);
    }
};

bool operator==(const Vertex &a, const Vertex &b);

#define SELECTED true
typedef vector<vector<vector<vector<vector<int>>>>> score_matrix;
typedef vector<vector<vector<vector<Vertex>>>> path_matrix;
typedef vector<vector<vector<Vertex>>> belonging_path_matrix;
typedef vector<vector<vector<int>>> weight_matrix;

typedef int path_continuation;
#define I 0
#define E 1

// Helper functions for vertex type.
bool isLayerVertex(Vertex v, const eds_matrix &eds_segments);

bool isJVertex(Vertex v, const eds_matrix &eds_segments);

bool isNVertex(Vertex v, const eds_matrix &eds_segments);

bool isFirstLayerVertex(Vertex v, const eds_matrix &eds_segments);

bool hasPredecessorVertex(Vertex v);

Vertex getPredecessorVertex(const eds_matrix &eds_segments, Vertex v,
                            int predecessor_layer);

int getScore(score_matrix &scores, Vertex v, bool selected,
             path_continuation layer);

void setScore(score_matrix &scores, int score, Vertex v, bool selected,
              path_continuation layer);

int getWeight(const weight_matrix &weights, Vertex v);

score_matrix initScoreMatrix(const weight_matrix &weights);

void findMaxScoringPaths(const eds_matrix &eds_segments,
                         const weight_matrix &weights,
                         score_matrix &scores, path_matrix &paths,
                         belonging_path_matrix &belonging_path_start,
                         int penalty);

#endif