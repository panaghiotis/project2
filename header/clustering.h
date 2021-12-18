#ifndef CLUSTERING_H
#define CLUSTERING_H

//#include "Dataset.h"
//#include "Point.h"
#include <unordered_map>
#include "neighbour_search.h"


#define CLASSIC 1
#define USE_RANGE_LSH 2
#define USE_HYPERCUBE 3
#define USE_LSH_FRECHET 4
#define MEAN_VECTOR 5
#define MEAN_FRECHET 6

#define DIVIDE_DATASET_FOR_HASHTABLE_SIZE 8
#define CENTROID_CHANGE_THRESHOLD 5.0
#define MAX_LENGTH 1000

struct Centroid : public Point {
    long double silhouette;
    Centroid(vector<double> *coordinates, string id = "", int pos = -1)
        : Point(coordinates, id, pos), silhouette(0.0) {}
};

struct Cluster {
    vector<Point *> points;
    vector<Curve *> curves;
    unordered_map<string, Point *> point_map;
    unordered_map<string, Curve *> curve_map;
    vector<long double> points_silhouette;
    Centroid *centroid;
    Curve *meanCurve;
    Cluster(Centroid *c) : centroid(c) {meanCurve = NULL;}
    Cluster(Curve *c) : meanCurve(c) {centroid = NULL;}
    void add(Point *p);
    void addCurve(Curve *c);
    void clear_cluster();
    pair<Centroid *, long double> recalculate_centroid(unsigned int dim);
    void curve_map_to_vec();        // curves vector gets the same curves as curve map (for safety reasons)

     //optimal traversal computation of two curves
    vector<pair< pair<double,double>, pair<double,double> >> OptimalTraversal(vector<pair<double,double>> c, vector<pair<double,double>> centroid);

    //BST traversal that creates coordinates for new mean curve in root.(Starting from leaves-cluster_curves to root)
    int PostOrderTraversal(int node, int leaves_places, vector<vector<pair<double,double>>> &tree);

    //initializes tree with leaves the curves of the cluster.Calls post order BST traversal, gets the root and creates the new mean curve
    pair<Curve *, long double> recalculate_meanCurve();

    ~Cluster();
};


class Clustering {
    Dataset *data;
    vector<Centroid *> centroids;
    vector<Curve *> meanCurves;
    vector<Cluster *> clusters;
    void initialize_clusters(unsigned int k);
    void initialize_R2clusters(unsigned int k);
    vector<long double> silhouette;
    const int method;
    const int update;
    NearestNeighboursSearch search_functions;
public:
    Clustering(Dataset *data, int method = CLASSIC, int update = MEAN_VECTOR) : data(data), method(method), update(update), search_functions(*data,*data) {}
    ~Clustering();
    double perform_kMeans(unsigned int k, unsigned int M, unsigned int probes);   // k is number of clusters
    double perform_R2kMeans(unsigned int k);                                      // k is number of clusters, kMeans for curves in R2
    const vector<Cluster *> &get_clusters() const { return clusters; }
    vector<long double> get_silhouette() const { return this->silhouette; }
    double initialize_radius();
    double initialize_radiusR2();           //radius for curves
    void calculate_silhouette();
    void calculate_silhouetteR2();          //silhouette for curves
};


#endif
