#ifndef MAXSCOREPATH_UTILITY_FUNC_HEADER
#define MAXSCOREPATH_UTILITY_FUNC_HEADER

#include <boost/algorithm/string/trim.hpp>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

using namespace std;

typedef vector<vector<string>> eds_matrix;
typedef vector<vector<vector<vector<vector<int>>>>> score_matrix;

// Returns the one line elastic degenerate string file as a string.
string readEDSFile(const string &file_path);

// Store the elastic degenerate string `text` in a 2D matrix.
eds_matrix EDSToMatrix(const string &EDS);

// Assign score based on GC content:
// - bases G and C get score 1
// - empty non-deterministic segment parts get score 0
// - every other character gets score 0.
vector<vector<vector<int>>> getGCContentWeights(
    const vector<vector<string>> &eds_segments);

// Find maximum-scoring paths.
#define SELECTED true

struct vertex {
    int segment;
    int layer;
    int index;
};

bool operator==(const vertex &a, const vertex &b);

// Helper functions for vertex type.
bool isLayerVertex(vertex v, const vector<vector<string>> &eds_segments);

bool isJVertex(vertex v, const vector<vector<string>> &eds_segments);

bool isNVertex(vertex v, const vector<vector<string>> &eds_segments);

bool isBaseLayerVertex(vertex v, const vector<vector<string>> &eds_segments);

bool hasPredecessorVertex(vertex v);

vertex getPredecessorVertex(const vector<vector<string>> &eds_segments,
                            vertex v, int layer);

int getPredecessorScore(const vector<vector<string>> &eds_segments,
                        vector<vector<vector<vector<vector<int>>>>> &scores,
                        vertex v, bool selected, int layer);

void storePath(vertex store_for, vertex closest_predecessor_on_path,
               unordered_map<vertex, vector<vertex>> &paths,
               unordered_map<vertex, vertex> &belonging_path_start);

void findMaxScoringPaths(const vector<vector<string>> &eds_segments,
                         const vector<vector<vector<int>>> &weights,
                         score_matrix &scores,
                         unordered_map<vertex, vector<vertex>> &paths,
                         unordered_map<vertex, vertex> &belonging_path_start,
                         int penalty);

#endif