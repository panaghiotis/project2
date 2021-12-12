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
long double discreteFrechet_distance(vector<pair<int,double>> p, vector<pair<int,double>> q);
double continuousFrechet_distanceByFred(Curve &c1, Curve &c2);
long double dot_product(const Point &p1, const Point &p2);
long unsigned int int_dot_product(const vector<unsigned int> &v1, const vector<unsigned int> &v2);
Point *vec_add(const Point &p1, const Point &p2);
Centroid *vec_avg(vector<Point *> &points, unsigned int dim);
double approx_factor(double exact, double approx);
double max_approx_factor(vector<double> factor_arr);

#endif
