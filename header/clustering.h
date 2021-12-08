#ifndef CLUSTERING_H
#define CLUSTERING_H

//#include "Dataset.h"
//#include "Point.h"
#include <unordered_map>
#include "neighbour_search.h"


#define CLASSIC 1
#define USE_RANGE_LSH 2
#define USE_HYPERCUBE 3

#define DIVIDE_DATASET_FOR_HASHTABLE_SIZE 8
#define MINIMUM_NUM_CLUSTERS_FOR_LSH 32
#define CENTROID_CHANGE_THRESHOLD 5.0


struct Centroid : public Point {
    long double silhouette;
    Centroid(vector<double> *coordinates, string id = "", int pos = -1)
        : Point(coordinates, id, pos), silhouette(0.0) {}
};


struct Cluster {
    vector<Point *> points;
    unordered_map<string, Point*> point_map;
    vector<long double> points_silhouette;
    Centroid *centroid;
    Cluster(Centroid *c) : centroid(c) {}
    void add(Point *p);
    void clear_cluster();
    pair<Centroid *, long double> recalculate_centroid(unsigned int dim);
    ~Cluster() {
        delete centroid;
    }
};


class Clustering {
    Dataset *data;
    vector<Centroid *> centroids;
    vector<Cluster *> clusters;
    void initialize_clusters(unsigned int k);
    vector<long double> silhouette;
    const int method;
    NearestNeighboursSearch search_functions;
public:
    Clustering(Dataset *data, int method = CLASSIC) : data(data), method(method), search_functions(*data,*data) {}
    ~Clustering();
    double perform_kMeans(unsigned int k, unsigned int M, unsigned int probes);   // number of clusters
    const vector<Cluster *> &get_clusters() const { return clusters; }
    vector<long double> get_silhouette() const { return this->silhouette; }
    double initialize_radius();
    void calculate_silhouette();
};


#endif
