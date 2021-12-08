#include <cmath>
#include <vector>
#include<algorithm>
#include "../header/vector_operations.h"

using namespace std;


/***
 * Wrapper to use for distances elsewhere.
 * If we want to change the distance used we can do this here.
 */
long double used_distance(const Point &p1, const Point &p2) {
    return L2_distance(p1, p2);   // use Euclidean distance
}

long double L2_distance(const Point &p1, const Point &p2) {
    long double sum = 0.0;
    for (int i = 0 ; i < p1.coords->size() && i < p2.coords->size(); i++) {
        sum += pow(p1.coords->at(i) - p2.coords->at(i), 2);
    }
    return sqrt(sum);
}

//euclidian distance for Discrete Frechet Distance implementation
long double euclidean2d_distance(pair<int,double> p, pair<int,double> q) {
    long double sum = 0.0;
    long double x = p.first - q.first;
    long double y = p.second - q.second;
    sum = pow(x,2) + pow(y,2);
    return sqrt(sum);
}

long double discreteFrechet_distance(vector<pair<int,double>> p, vector<pair<int,double>> q) {
    long double result = 0.0;

    //c(i,j) initialization
    long double **c = new long double*[p.size()];
    for(int i = 0; i < p.size(); ++i)
        c[i] = new long double[p.size()];
    //computing the DFD

    // fill the first value with the distance between the first two points in p and q
    c[0][0] = euclidean2d_distance(p[0], q[0]);

    // load the first column and first row with distances (memorize)
    for (int i=1; i < p.size(); i++)
        c[i][0] = max(c[i-1][0], euclidean2d_distance(p[i], q[0]));
    for (int j=1; j < q.size(); j++)
        c[0][j] = max(c[0][j-1], euclidean2d_distance(p[0], q[j]));

    // load random column and row with distances
    for (int i=1; i < p.size(); i++)
        for (int j=1; j < q.size(); j++)
            c[i][j] = max(min(min(c[i-1][j], c[i][j-1]), c[i-1][j-1]), euclidean2d_distance(p[i], q[j]));

    //get discrete frechet distance
    result= c[p.size()-1][q.size()-1];

    for(int i = 0; i < p.size();++i)
        delete[] c[i];
    delete[] c;

    return result;
}

long double dot_product(const Point &p1, const Point &p2) {
    long double sum = 0.0;
    for (int i = 0 ; i < p1.coords->size() && i < p2.coords->size(); i++) {
        sum += p1.coords->at(i) * p2.coords->at(i);
    }
    return sum;
}

long unsigned int int_dot_product(const vector<unsigned int> &v1, const vector<unsigned int> &v2) {
    long unsigned int sum = 0;
    for (int i = 0 ; i < v1.size() && i < v2.size(); i++) {
        sum += v1.at(i) * v2.at(i);
    }
    return sum;
}

Point *vec_add(const Point &p1, const Point &p2) {
    vector<double> *coords = new vector<double>();
    for (int i = 0 ; i < p1.coords->size() && i < p2.coords->size(); i++) {
        coords->push_back(p1.coords->at(i) + p2.coords->at(i));
    }
    return new Point(coords);
}

Centroid *vec_avg(vector<Point *> &points, unsigned int dim) {
    vector<double> *sums = new vector<double>(dim);
    for (int j = 0 ; j < points[0]->coords->size() ; j++) {
        (*sums)[j] = 0.0;
    }
    if (points.empty()){          // empty cluster
        cerr << "Warning: cluster with no points assigned" << endl;
        return NULL;
    }
    for (int i = 0 ; i < points.size() ; i++) {
        for (int j = 0 ; j < dim ; j++) {
            (*sums)[j] += points[i]->coords->at(j);
        }
    }
    for (int j = 0 ; j < dim ; j++) {
        (*sums)[j] /= ((double) points.size());
    }
    return new Centroid(sums);
}

ostream& operator<<(ostream& os, const Point& p) {
    os << "[ ";
    for (int i = 0 ; i < p.coords->size() ; i++) {
        os << p.coords->at(i) << ", ";
    }
    os << "]";
    return os;
}
