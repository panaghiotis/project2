#ifndef VECTOR_OPERATIONS_H
#define VECTOR_OPERATIONS_H


#include <utility>
#include "Point.h"
#include "clustering.h"


/***
 * Wrapper to use for distances elsewhere.
 * If we want to change the distance used we can do this here.
 */
long double used_distance(const Point &p1, const Point &p2, bool isFrechet = false);

long double L2_distance(const Point &p1, const Point &p2);
long double euclidean2d_distance(pair<int,double> p, pair<int,double> q);

//Curve *curve(this curve is a mean curve) is only used for saving c array for clustering
long double discreteFrechet_distance(vector<pair<double,double>> p, vector<pair<double,double>> q, Curve *curve = NULL);
double** get_C(vector<pair<double,double>> p, vector<pair<double,double>> q);

double continuousFrechet_distanceByFred(Curve &c1, Curve &c2);
long double dot_product(const Point &p1, const Point &p2);
long unsigned int int_dot_product(const vector<unsigned int> &v1, const vector<unsigned int> &v2);
Point *vec_add(const Point &p1, const Point &p2);
Centroid *vec_avg(vector<Point *> &points, unsigned int dim);

//filtering mean curve from root of BST if root has more coordinates than MAX_LENGTH number
vector<pair<double,double>> filtering(vector<pair<double,double>> R2_coords, double epsilon);

//getting reverse coordinates from optimal traversal of two curves. Sets coordinates for new node until we reach root
vector<pair<double,double>> vec_avg_curve(vector<pair< pair<double,double>, pair<double,double> >> reverse_traversal, double epsilon);

//MAF functions
double approx_factor(double exact, double approx);
double max_approx_factor(vector<double> factor_arr);

#endif
