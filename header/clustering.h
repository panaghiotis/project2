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

//struct Mean_Curve : public Curve {
//    int max_length;
//    unordered_map<string,long double **> backtrack_map;
//    //Mean_Curve(Curve &curve): Curve(curve) {}
//    Mean_Curve(vector<pair<double,double>> *curve_coords ,int dim = 0, string name = "")
//        : Curve(dim,name,curve_coords), max_length(MAX_LENGTH) {}
//    ~Mean_Curve();
//    void set_backtrack(Curve *curve, int mean_size);    //set the array for later optimal traversal
//};

struct Cluster {
    vector<Point *> points;
    vector<Curve *> curves;
    unordered_map<string, Point *> point_map;
    unordered_map<string, Curve *> curve_map;
    vector<long double> points_silhouette;
    Centroid *centroid;
    /*Mean_Curve*/Curve *meanCurve;
    unordered_map<string, long double **> backtrack_map;
    vector<Curve *> *meanTree;
    Cluster(Centroid *c) : centroid(c) {}
    Cluster(/*Mean_*/Curve *c) : meanCurve(c) {}
    void add(Point *p);
    void clear_cluster();
    pair<Centroid *, long double> recalculate_centroid(unsigned int dim);

    //set the array for later optimal traversal
    void set_backtrack(Curve *curve, int mean_size);

    //optimal traversal computation
    vector<pair< pair<double,double>, pair<double,double> >> *OptimalTraversal(Curve *c, Curve *centroid);

    //size is cluster size for tree height computation
    void initialize_tree();
    pair</*Mean_*/Curve *, long double> recalculate_meanCurve();

    //Tree traversal
    /*Mean_Curve **/int PostOrderTraversal(int node);

    ~Cluster();
};


class Clustering {
    Dataset *data;
    vector<Centroid *> centroids;
    vector</*Mean_*/Curve *> meanCurves;
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
    //double perform_R2kMeans(unsigned int k);                                      // k is number of clusters
    const vector<Cluster *> &get_clusters() const { return clusters; }
    vector<long double> get_silhouette() const { return this->silhouette; }
    double initialize_radius();     //TODO: radius for curves
    void calculate_silhouette();    //TODO: silhouette for curves
};


#endif
