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
vector<vector<vector<int>>> getGCContentWeights(
    const eds_matrix &eds_segments);

// Find maximum-scoring paths.

struct Vertex {
    int segment;
    int layer;
    int index;

    Vertex() : segment(-1), layer(-1), index(-1) {}
    Vertex(int segment, int layer, int index) : segment(segment), layer(layer), index(index) {}

    explicit operator bool() { return !(segment == -1 && layer == -1 && index == -1); }
};

bool operator==(const Vertex &a, const Vertex &b);

#define SELECTED true
typedef vector<vector<vector<vector<vector<int>>>>> score_matrix;
typedef vector<vector<vector<vector<Vertex>>>> path_matrix;
typedef vector<vector<vector<Vertex>>> belonging_path_matrix;

// Helper functions for vertex type.
bool isLayerVertex(Vertex v, const eds_matrix &eds_segments);

bool isJVertex(Vertex v, const eds_matrix &eds_segments);

bool isNVertex(Vertex v, const eds_matrix &eds_segments);

bool isBaseLayerVertex(Vertex v, const eds_matrix &eds_segments);

bool hasPredecessorVertex(Vertex v);

Vertex getPredecessorVertex(const eds_matrix &eds_segments,
                            Vertex v, int layer);

int getPredecessorScore(const eds_matrix &eds_segments,
                        score_matrix &scores,
                        Vertex v, bool selected, int layer);

void storePath(Vertex store_for, Vertex closest_predecessor_on_path,
               path_matrix &paths,
               belonging_path_matrix &belonging_path_start);

score_matrix initScoreMatrix(const vector<vector<vector<int>>> &weights);
path_matrix initPathMatrix(const vector<vector<vector<int>>> &weights);
belonging_path_matrix initBelongingPathMatrix(
    const vector<vector<vector<int>>> &weights);

void findMaxScoringPaths(const eds_matrix &eds_segments,
                         const vector<vector<vector<int>>> &weights,
                         score_matrix &scores,
                         path_matrix &paths,
                         belonging_path_matrix &belonging_path_start,
                         int penalty);

#endif