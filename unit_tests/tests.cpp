#include <gtest/gtest.h>

#include <iostream>
#include <set>

#include "../utility_func.hpp"

using namespace std;

TEST(InputProcessing, readEDSFileTest) {
    EXPECT_EQ(readEDSFile("../unit_tests/test_inputs/input_01.txt"),
              "_GG{AGAA,GGGA,,ACCCCC}{AG,G}AGG{A,G}{C,}{A,AG}G{A,GA,CCC}{,A}_");
}

TEST(InputProcessing, EDSToMatrixTest) {
    string EDS =
        "_GG{AGAA,GGGA,,ACCCCC}{AG,G}AGG{A,G}{C,}{A,AG}G{A,GA,CCC}{,A}_";
    eds_matrix eds_segments = EDSToMatrix(EDS);

    eds_matrix expected;
    expected.emplace_back(vector<string>{"_GG"});
    expected.emplace_back(vector<string>{"AGAA", "GGGA", "_", "ACCCCC"});
    expected.emplace_back(vector<string>{"_"});
    expected.emplace_back(vector<string>{"AG", "G"});
    expected.emplace_back(vector<string>{"AGG"});
    expected.emplace_back(vector<string>{"A", "G"});
    expected.emplace_back(vector<string>{"_"});
    expected.emplace_back(vector<string>{"C", "_"});
    expected.emplace_back(vector<string>{"_"});
    expected.emplace_back(vector<string>{"A", "AG"});
    expected.emplace_back(vector<string>{"G"});
    expected.emplace_back(vector<string>{"A", "GA", "CCC"});
    expected.emplace_back(vector<string>{"_"});
    expected.emplace_back(vector<string>{"_", "A"});
    expected.emplace_back(vector<string>{"_"});

    EXPECT_EQ(eds_segments, expected);
}

class PathFindingTest : public testing::Test {
   public:
    PathFindingTest() {}
    void GetPaths(const string& EDS, int penalty = 2, int match = 1,
                  int non_match = -1) {
        eds_matrix eds_segments = EDSToMatrix(EDS);
        weight_matrix weights =
            getGCContentWeights(eds_segments, match, non_match);

        score_matrix scores = initScoreMatrix(weights);
        score_matrix choices = initScoreMatrix(weights);

        score = findMaxScoringPaths(eds_segments, weights, scores, choices,
                                    penalty);
        paths = getPaths(eds_segments, scores, choices);
    }

    int score;
    vector<vector<Vertex>> paths;
};

TEST_F(PathFindingTest, LinearGraphTest01) {
    GetPaths("_GGCAGGGAAGAAGGA_");
    EXPECT_EQ(this->score, 3);
    set<vector<Vertex>> expected;
    expected.insert({Vertex(0, 0, 1), Vertex(0, 0, 2), Vertex(0, 0, 3),
                     Vertex(0, 0, 4), Vertex(0, 0, 5), Vertex(0, 0, 6),
                     Vertex(0, 0, 7)});
    expected.insert({Vertex(0, 0, 13), Vertex(0, 0, 14)});
    EXPECT_EQ(set<vector<Vertex>>(this->paths.begin(), this->paths.end()),
              expected);
}

TEST_F(PathFindingTest, LinearGraphTest02) {
    GetPaths("_GGCGAAAGGGA_");
    EXPECT_EQ(this->score, 3);
    set<vector<Vertex>> expected;
    expected.insert(
        {Vertex(0, 0, 1), Vertex(0, 0, 2), Vertex(0, 0, 3), Vertex(0, 0, 4)});
    expected.insert({Vertex(0, 0, 8), Vertex(0, 0, 9), Vertex(0, 0, 10)});
    EXPECT_EQ(set<vector<Vertex>>(this->paths.begin(), this->paths.end()),
              expected);
}

TEST_F(PathFindingTest, BubbleContinuationUp01) {
    GetPaths("_GG{GA,AAGG,TTAACAG,,ACTCCTT}_");
    EXPECT_EQ(this->score, 1);
    set<vector<Vertex>> expected;
    expected.insert({Vertex(0, 0, 1), Vertex(0, 0, 2), Vertex(1, 0, 0)});
    expected.insert({Vertex(1, 1, 2), Vertex(1, 1, 3)});
    expected.insert(
        {Vertex(1, 4, 1), Vertex(1, 4, 2), Vertex(1, 4, 3), Vertex(1, 4, 4)});
    EXPECT_EQ(set<vector<Vertex>>(this->paths.begin(), this->paths.end()),
              expected);
}

TEST_F(PathFindingTest, BubbleContinuationUp02) {
    GetPaths("_GG{GGCA,AAGG,CATT,TATTTA,,GATTGTTG}_");
    EXPECT_EQ(this->score, 3);
    set<vector<Vertex>> expected;
    expected.insert({Vertex(0, 0, 1), Vertex(0, 0, 2), Vertex(1, 0, 0),
                     Vertex(1, 0, 1), Vertex(1, 0, 2)});
    expected.insert({Vertex(1, 1, 2), Vertex(1, 1, 3)});
    EXPECT_EQ(set<vector<Vertex>>(this->paths.begin(), this->paths.end()),
              expected);
}

TEST_F(PathFindingTest, BubbleContinuationUpAndMultipleLayerPaths01) {
    GetPaths("_GG{GGCAAAGGGAGG,AAGG,TTACTTACTTCT,}_");
    EXPECT_EQ(this->score, 5);
    set<vector<Vertex>> expected;
    expected.insert({Vertex(0, 0, 1), Vertex(0, 0, 2), Vertex(1, 0, 0),
                     Vertex(1, 0, 1), Vertex(1, 0, 2)});
    expected.insert({Vertex(1, 0, 6), Vertex(1, 0, 7), Vertex(1, 0, 8),
                     Vertex(1, 0, 9), Vertex(1, 0, 10), Vertex(1, 0, 11)});
    expected.insert({Vertex(1, 1, 2), Vertex(1, 1, 3)});
    EXPECT_EQ(set<vector<Vertex>>(this->paths.begin(), this->paths.end()),
              expected);
}

TEST_F(PathFindingTest, BubbleContinuationUpAndMultipleLayerPaths02) {
    GetPaths("_GG{GGCAAAGGGAAACGCAAAA,AAGG,GTTGTGT,TGGTGGCAATG,,ATGTGGC}_");
    EXPECT_EQ(this->score, 8);
    set<vector<Vertex>> expected;
    expected.insert({Vertex(0, 0, 1), Vertex(0, 0, 2), Vertex(1, 0, 0),
                     Vertex(1, 0, 1), Vertex(1, 0, 2)});
    expected.insert({Vertex(1, 0, 6), Vertex(1, 0, 7), Vertex(1, 0, 8)});
    expected.insert({Vertex(1, 0, 12), Vertex(1, 0, 13), Vertex(1, 0, 14)});
    expected.insert({Vertex(1, 1, 2), Vertex(1, 1, 3)});
    expected.insert({Vertex(1, 3, 1), Vertex(1, 3, 2), Vertex(1, 3, 3),
                     Vertex(1, 3, 4), Vertex(1, 3, 5), Vertex(1, 3, 6)});
    expected.insert({Vertex(1, 5, 2), Vertex(1, 5, 3), Vertex(1, 5, 4),
                     Vertex(1, 5, 5), Vertex(1, 5, 6)});
    EXPECT_EQ(set<vector<Vertex>>(this->paths.begin(), this->paths.end()),
              expected);
}

TEST_F(PathFindingTest, BubbleContinuationUpAndMultipleLayerPaths03) {
    GetPaths(
        "_GG{GGCAAAGGGAAACGCAAAA,AAGGAGGG,GTTGTGT,TGGTGGCAATG,,ATGTGGC,AAGG}_");
    EXPECT_EQ(this->score, 10);
    set<vector<Vertex>> expected;
    expected.insert({Vertex(0, 0, 1), Vertex(0, 0, 2), Vertex(1, 0, 0),
                     Vertex(1, 0, 1), Vertex(1, 0, 2)});
    expected.insert({Vertex(1, 0, 6), Vertex(1, 0, 7), Vertex(1, 0, 8)});
    expected.insert({Vertex(1, 0, 12), Vertex(1, 0, 13), Vertex(1, 0, 14)});
    expected.insert({Vertex(1, 1, 2), Vertex(1, 1, 3), Vertex(1, 1, 4),
                     Vertex(1, 1, 5), Vertex(1, 1, 6), Vertex(1, 1, 7)});
    expected.insert({Vertex(1, 3, 1), Vertex(1, 3, 2), Vertex(1, 3, 3),
                     Vertex(1, 3, 4), Vertex(1, 3, 5), Vertex(1, 3, 6)});
    expected.insert({Vertex(1, 5, 2), Vertex(1, 5, 3), Vertex(1, 5, 4),
                     Vertex(1, 5, 5), Vertex(1, 5, 6)});
    expected.insert({Vertex(1, 6, 2), Vertex(1, 6, 3)});
    EXPECT_EQ(set<vector<Vertex>>(this->paths.begin(), this->paths.end()),
              expected);
}

TEST_F(PathFindingTest, BubbleContinuationUpAndMultipleLayerPaths04) {
    GetPaths(
        "_GG{GGCAAAGGGAAACGCAAAA,AAGGAGGGAAAGGG,GTTGTGT,TGGTGGCAATGCGT,,"
        "ATGTGGC,AAGGAGCAATGGGGG}_");
    EXPECT_EQ(this->score, 16);
    set<vector<Vertex>> expected;
    expected.insert({Vertex(0, 0, 1), Vertex(0, 0, 2), Vertex(1, 0, 0),
                     Vertex(1, 0, 1), Vertex(1, 0, 2)});
    expected.insert({Vertex(1, 0, 6), Vertex(1, 0, 7), Vertex(1, 0, 8)});
    expected.insert({Vertex(1, 0, 12), Vertex(1, 0, 13), Vertex(1, 0, 14)});
    expected.insert({Vertex(1, 1, 2), Vertex(1, 1, 3), Vertex(1, 1, 4),
                     Vertex(1, 1, 5), Vertex(1, 1, 6), Vertex(1, 1, 7)});
    expected.insert({Vertex(1, 1, 11), Vertex(1, 1, 12), Vertex(1, 1, 13)});
    expected.insert({Vertex(1, 3, 1), Vertex(1, 3, 2), Vertex(1, 3, 3),
                     Vertex(1, 3, 4), Vertex(1, 3, 5), Vertex(1, 3, 6)});
    expected.insert({Vertex(1, 3, 10), Vertex(1, 3, 11), Vertex(1, 3, 12)});
    expected.insert({Vertex(1, 5, 2), Vertex(1, 5, 3), Vertex(1, 5, 4),
                     Vertex(1, 5, 5), Vertex(1, 5, 6)});
    expected.insert({Vertex(1, 6, 2), Vertex(1, 6, 3), Vertex(1, 6, 4),
                     Vertex(1, 6, 5), Vertex(1, 6, 6)});
    expected.insert({Vertex(1, 6, 10), Vertex(1, 6, 11), Vertex(1, 6, 12),
                     Vertex(1, 6, 13), Vertex(1, 6, 14)});
    EXPECT_EQ(set<vector<Vertex>>(this->paths.begin(), this->paths.end()),
              expected);
}

TEST_F(PathFindingTest, BubbleContinuationNonUp01) {
    GetPaths("_GG{GA,GGAA,TTAACAG,,ACTCCTT}_");
    EXPECT_EQ(this->score, 2);
    set<vector<Vertex>> expected;
    expected.insert(
        {Vertex(0, 0, 1), Vertex(0, 0, 2), Vertex(1, 1, 0), Vertex(1, 1, 1)});
    expected.insert(
        {Vertex(1, 4, 1), Vertex(1, 4, 2), Vertex(1, 4, 3), Vertex(1, 4, 4)});
    EXPECT_EQ(set<vector<Vertex>>(this->paths.begin(), this->paths.end()),
              expected);
}

TEST_F(PathFindingTest, BubbleContinuationNonUp02) {
    GetPaths("_GG{GA,AAGG,TTAACAG,,ACTCCTT,GGAG}_");
    EXPECT_EQ(this->score, 2);
    set<vector<Vertex>> expected;
    expected.insert({Vertex(0, 0, 1), Vertex(0, 0, 2), Vertex(1, 5, 0),
                     Vertex(1, 5, 1), Vertex(1, 5, 2), Vertex(1, 5, 3)});
    expected.insert({Vertex(1, 1, 2), Vertex(1, 1, 3)});
    expected.insert(
        {Vertex(1, 4, 1), Vertex(1, 4, 2), Vertex(1, 4, 3), Vertex(1, 4, 4)});
    EXPECT_EQ(set<vector<Vertex>>(this->paths.begin(), this->paths.end()),
              expected);
}

TEST_F(PathFindingTest, BubbleContinuationNotUpAndMultipleLayerPaths01) {
    GetPaths(
        "_GG{TCTAAAGGGAAACGCAAAA,AAGGAGGGAAAGGG,GTTGTGT,TGGTGGCAATGCGT,"
        "GGCAAAGGGAAACGCAAAA,ATGTGGC,AAGGAGCAATGGGGG,}_");
    EXPECT_EQ(this->score, 18);
    set<vector<Vertex>> expected;
    expected.insert({Vertex(1, 0, 6), Vertex(1, 0, 7), Vertex(1, 0, 8)});
    expected.insert({Vertex(1, 0, 12), Vertex(1, 0, 13), Vertex(1, 0, 14)});
    expected.insert({Vertex(1, 1, 2), Vertex(1, 1, 3), Vertex(1, 1, 4),
                     Vertex(1, 1, 5), Vertex(1, 1, 6), Vertex(1, 1, 7)});
    expected.insert({Vertex(1, 1, 11), Vertex(1, 1, 12), Vertex(1, 1, 13)});
    expected.insert({Vertex(1, 3, 1), Vertex(1, 3, 2), Vertex(1, 3, 3),
                     Vertex(1, 3, 4), Vertex(1, 3, 5), Vertex(1, 3, 6)});
    expected.insert({Vertex(1, 3, 10), Vertex(1, 3, 11), Vertex(1, 3, 12)});
    expected.insert({Vertex(0, 0, 1), Vertex(0, 0, 2), Vertex(1, 4, 0),
                     Vertex(1, 4, 1), Vertex(1, 4, 2)});
    expected.insert({Vertex(1, 4, 6), Vertex(1, 4, 7), Vertex(1, 4, 8)});
    expected.insert({Vertex(1, 4, 12), Vertex(1, 4, 13), Vertex(1, 4, 14)});
    expected.insert({Vertex(1, 5, 2), Vertex(1, 5, 3), Vertex(1, 5, 4),
                     Vertex(1, 5, 5), Vertex(1, 5, 6)});
    expected.insert({Vertex(1, 6, 2), Vertex(1, 6, 3), Vertex(1, 6, 4),
                     Vertex(1, 6, 5), Vertex(1, 6, 6)});
    expected.insert({Vertex(1, 6, 10), Vertex(1, 6, 11), Vertex(1, 6, 12),
                     Vertex(1, 6, 13), Vertex(1, 6, 14)});
    EXPECT_EQ(set<vector<Vertex>>(this->paths.begin(), this->paths.end()),
              expected);
}

TEST_F(PathFindingTest, BubbleContinuationNotUpAndMultipleLayerPaths02) {
    GetPaths(
        "_GG{GTGAAAGGGAAACGCAAAA,GGGCT,AAGG,GTTGTGT,TGGTGGCAATG,,ATGTGGC}_");
    EXPECT_EQ(this->score, 9);
    set<vector<Vertex>> expected;
    expected.insert({Vertex(0, 0, 1), Vertex(0, 0, 2), Vertex(1, 1, 0),
                     Vertex(1, 1, 1), Vertex(1, 1, 2), Vertex(1, 1, 3)});
    expected.insert({Vertex(1, 0, 6), Vertex(1, 0, 7), Vertex(1, 0, 8)});
    expected.insert({Vertex(1, 0, 12), Vertex(1, 0, 13), Vertex(1, 0, 14)});
    expected.insert({Vertex(1, 2, 2), Vertex(1, 2, 3)});
    expected.insert({Vertex(1, 4, 1), Vertex(1, 4, 2), Vertex(1, 4, 3),
                     Vertex(1, 4, 4), Vertex(1, 4, 5), Vertex(1, 4, 6)});
    expected.insert({Vertex(1, 6, 2), Vertex(1, 6, 3), Vertex(1, 6, 4),
                     Vertex(1, 6, 5), Vertex(1, 6, 6)});
    EXPECT_EQ(set<vector<Vertex>>(this->paths.begin(), this->paths.end()),
              expected);
}

TEST_F(PathFindingTest, BubbleContinuationToJVertex01) {
    GetPaths("_TGGTG{GA,TACTCT,TC,GAG,TCTT}CG");
    EXPECT_EQ(this->score, 3);
    set<vector<Vertex>> expected;
    expected.insert({Vertex(0, 0, 2), Vertex(0, 0, 3), Vertex(0, 0, 4),
                     Vertex(0, 0, 5), Vertex(1, 3, 0), Vertex(1, 3, 1),
                     Vertex(1, 3, 2), Vertex(2, 0, 0), Vertex(2, 0, 1)});
    EXPECT_EQ(set<vector<Vertex>>(this->paths.begin(), this->paths.end()),
              expected);
}

TEST_F(PathFindingTest, BubbleContinuationToJVertex02) {
    GetPaths("_TGGTG{GGCATT,TACTCT,TC,GAAGG,TCTT}CG");
    EXPECT_EQ(this->score, 5);
    set<vector<Vertex>> expected;
    expected.insert({Vertex(0, 0, 2), Vertex(0, 0, 3), Vertex(0, 0, 4),
                     Vertex(0, 0, 5), Vertex(1, 0, 0), Vertex(1, 0, 1),
                     Vertex(1, 0, 2)});
    expected.insert(
        {Vertex(1, 3, 3), Vertex(1, 3, 4), Vertex(2, 0, 0), Vertex(2, 0, 1)});
    EXPECT_EQ(set<vector<Vertex>>(this->paths.begin(), this->paths.end()),
              expected);
}

TEST_F(PathFindingTest, BubbleContinuationToJVertex03) {
    GetPaths("_TGATA{GGCATT,TACTCTC,TC,GAAGG,TCTC,}CG");
    EXPECT_EQ(this->score, 3);
    set<vector<Vertex>> expected;
    expected.insert({Vertex(1, 0, 0), Vertex(1, 0, 1), Vertex(1, 0, 2)});
    expected.insert(
        {Vertex(1, 3, 3), Vertex(1, 3, 4), Vertex(2, 0, 0), Vertex(2, 0, 1)});
    EXPECT_EQ(set<vector<Vertex>>(this->paths.begin(), this->paths.end()),
              expected);
}

TEST_F(PathFindingTest, Mix01) {
    GetPaths("_GG{AAGG,GGA}_");
    EXPECT_EQ(this->score, 2);
    set<vector<Vertex>> expected;
    expected.insert(
        {Vertex(0, 0, 1), Vertex(0, 0, 2), Vertex(1, 1, 0), Vertex(1, 1, 1)});
    expected.insert({Vertex(1, 0, 2), Vertex(1, 0, 3)});
    EXPECT_EQ(set<vector<Vertex>>(this->paths.begin(), this->paths.end()),
              expected);
}

TEST_F(PathFindingTest, Mix02) {
    GetPaths("_GG{AAGG,GGA}_");
    EXPECT_EQ(this->score, 2);
    set<vector<Vertex>> expected;
    expected.insert(
        {Vertex(0, 0, 1), Vertex(0, 0, 2), Vertex(1, 1, 0), Vertex(1, 1, 1)});
    expected.insert({Vertex(1, 0, 2), Vertex(1, 0, 3)});
    EXPECT_EQ(set<vector<Vertex>>(this->paths.begin(), this->paths.end()),
              expected);
}

TEST_F(PathFindingTest, Mix03) {
    GetPaths("_GGG{CCG,AGGGA}A_");
    EXPECT_EQ(this->score, 5);
    set<vector<Vertex>> expected;
    expected.insert({Vertex(0, 0, 1), Vertex(0, 0, 2), Vertex(0, 0, 3),
                     Vertex(1, 0, 0), Vertex(1, 0, 1), Vertex(1, 0, 2)});
    expected.insert({Vertex(1, 1, 1), Vertex(1, 1, 2), Vertex(1, 1, 3)});
    EXPECT_EQ(set<vector<Vertex>>(this->paths.begin(), this->paths.end()),
              expected);
}

TEST_F(PathFindingTest, Mix04) {
    GetPaths("_GGG{ATT,ATA}AGCGC_");
    EXPECT_EQ(this->score, 3);
    set<vector<Vertex>> expected;
    expected.insert({Vertex(0, 0, 1), Vertex(0, 0, 2), Vertex(0, 0, 3)});
    expected.insert({Vertex(2, 0, 1), Vertex(2, 0, 2), Vertex(2, 0, 3),
                     Vertex(2, 0, 4), Vertex(2, 0, 5)});
    EXPECT_EQ(set<vector<Vertex>>(this->paths.begin(), this->paths.end()),
              expected);
}

TEST_F(PathFindingTest, Mix05) {
    GetPaths("_AG{GGG,,CCC}{AG,GCGG,AA}A{A,G}{G,CC}{AAAA,}_");
    EXPECT_EQ(this->score, 9);
    set<vector<Vertex>> expected;
    expected.insert({Vertex(0, 0, 2), Vertex(1, 0, 0), Vertex(1, 0, 1),
                     Vertex(1, 0, 2), Vertex(2, 0, 0), Vertex(3, 1, 0),
                     Vertex(3, 1, 1), Vertex(3, 1, 2), Vertex(3, 1, 3),
                     Vertex(4, 0, 0), Vertex(5, 1, 0), Vertex(6, 0, 0),
                     Vertex(7, 1, 0), Vertex(7, 1, 1)});
    expected.insert({Vertex(1, 2, 0), Vertex(1, 2, 1), Vertex(1, 2, 2)});
    EXPECT_EQ(set<vector<Vertex>>(this->paths.begin(), this->paths.end()),
              expected);
}

TEST_F(PathFindingTest, Mix06) {
    GetPaths(
        "_GG{AGAA,GGGA,,ACCCCC}{AG,G}AGG{C,A,GT}{C,}{A,AG,TA}G{A,GA,CCC,"
        "TTTAGTG}{,A}_");
    EXPECT_EQ(this->score, 14);
    set<vector<Vertex>> expected;
    expected.insert({Vertex(0, 0, 1), Vertex(0, 0, 2), Vertex(1, 1, 0),
                     Vertex(1, 1, 1), Vertex(1, 1, 2)});
    expected.insert(
        {Vertex(1, 3, 1),  Vertex(1, 3, 2),  Vertex(1, 3, 3),  Vertex(1, 3, 4),
         Vertex(1, 3, 5),  Vertex(2, 0, 0),  Vertex(3, 1, 0),  Vertex(4, 0, 0),
         Vertex(4, 0, 1),  Vertex(4, 0, 2),  Vertex(5, 0, 0),  Vertex(6, 0, 0),
         Vertex(7, 0, 0),  Vertex(8, 0, 0),  Vertex(9, 1, 0),  Vertex(9, 1, 1),
         Vertex(10, 0, 0), Vertex(11, 2, 0), Vertex(11, 2, 1), Vertex(11, 2, 2),
         Vertex(12, 0, 0), Vertex(13, 0, 0)});
    EXPECT_EQ(set<vector<Vertex>>(this->paths.begin(), this->paths.end()),
              expected);
}

TEST_F(PathFindingTest, Mix07) {
    GetPaths("_GGAG{CCAAA,AGGGA,TATGC,ATGAC,GAC}CAGTTT{TC,GC,TT,GAC}CAGG_");
    EXPECT_EQ(this->score, 6);
    set<vector<Vertex>> expected;
    printPaths(paths);
    expected.insert({Vertex(0, 0, 1), Vertex(0, 0, 2), Vertex(0, 0, 3),
                     Vertex(0, 0, 4), Vertex(1, 0, 0), Vertex(1, 0, 1)});
    expected.insert({Vertex(1, 1, 1), Vertex(1, 1, 2), Vertex(1, 1, 3)});
    expected.insert({Vertex(1, 2, 3), Vertex(1, 2, 4), Vertex(2, 0, 0),
                     Vertex(2, 0, 1), Vertex(2, 0, 2)});
    expected.insert({Vertex(3, 1, 0), Vertex(3, 1, 1), Vertex(4, 0, 0),
                     Vertex(4, 0, 1), Vertex(4, 0, 2), Vertex(4, 0, 3),
                     Vertex(4, 0, 4)});
    EXPECT_EQ(set<vector<Vertex>>(this->paths.begin(), this->paths.end()),
              expected);
}