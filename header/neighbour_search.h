#ifndef NEIGHBOUR_SEARCH_H
#define NEIGHBOUR_SEARCH_H

#include <vector>
#include <list>
#include <chrono>
#include "Point.h"
#include "Dataset.h"

using namespace std;


#define USE_QUERY_TRICK true


struct Neighbour {
    string neighbourID;
    long double dist;
    Neighbour(string neighbourID, long double dist)
        : neighbourID(neighbourID), dist(dist) {}
};

struct Result {
    string queryID;
    vector<Neighbour> *NNs;
    long long int dt;   // in ms
};


class NearestNeighboursSearch {
    Dataset &index_dataset;
    Dataset &query_dataset;
    long double **exact_distances;                      // not used
    long double **lsh_distances;                        // not used
    long double **cube_distances;                       // not used
public:
    NearestNeighboursSearch(Dataset &index_dataset, Dataset &query_dataset);
    ~NearestNeighboursSearch();
    // Note: ended up not using these two
    long long int calculate_exact_distances();
    long long int calculate_only_lsh_distances();      // only calculate distances for same buckets
    vector<Result *> *calculate_exact_NN(int k);
    vector<Result *> *calculate_approximate_NN(int k, bool use_query_trick=USE_QUERY_TRICK);
    vector<vector<string> *> *range_search(long double r, bool use_query_trick=USE_QUERY_TRICK);
    list<Point*> LSH_search(long double r, Point *q, bool use_query_trick=USE_QUERY_TRICK);
     list<Point *> cube_search(unsigned int M, int probes, int cur_ver, int hamming, long double r, Point *q);
    vector<Result *> *calculate_approximate_NN_using_cube(int k, int probes, int M);
    vector<vector<string> *> *range_search_using_cube(long double r, int probes, int M);
};

// for sorting
struct sort_pred {
    bool operator()(const pair<string, long double> &left, const pair<string, long double> &right) {
        return left.second < right.second;
    }
};

struct sort_pred2 {
    bool operator()(const pair<Point *, long double> &left, const pair<Point *, long double> &right) {
        return left.second < right.second;
    }
};

#endif
