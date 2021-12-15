#include <cmath>
#include <vector>
#include<algorithm>
#include <cfloat>
#include "../header/vector_operations.h"
#include "../FredFilesNeeded/include/frechet.hpp"       //Fred lib


using namespace std;


/***
 * Wrapper to use for distances elsewhere.
 * If we want to change the distance used we can do this here.
 */
long double used_distance(const Point &p1, const Point &p2, bool isFrechet) {
    if(isFrechet) {
        //cout << p1.curve->get_curve_coords()->at(0).first << " " << p1.curve->get_curve_coords()->at(0).second << endl;
        return discreteFrechet_distance(*p1.curve->get_curve_coords(), *p2.curve->get_curve_coords());
    }
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

long double discreteFrechet_distance(vector<pair<double,double>> p, vector<pair<double,double>> q, Curve *curve) {
    long double result = 0.0;

    //c(i,j) initialization
    long double **c = new long double*[p.size()];
    for(int i = 0; i < p.size(); ++i)
        //c[i] = new long double[p.size()];
        c[i] = new long double[q.size()];
    //computing the DFD

    // fill the first value with the distance between the first two points in p and q
    c[0][0] = euclidean2d_distance(p[0], q[0]);

    // load the first column and first row with distances (memorize)
    for (int i=1; i < p.size(); i++)
        c[i][0] = max(c[i-1][0], euclidean2d_distance(p[i], q[0]));
    for (int j=1; j < q.size(); j++)
        c[0][j] = max(c[0][j-1], euclidean2d_distance(p[0], q[j]));

    // load every other column and row with distances
    for (int i=1; i < p.size(); i++)
        for (int j=1; j < q.size(); j++)
            c[i][j] = max(min(min(c[i-1][j], c[i][j-1]), c[i-1][j-1]), euclidean2d_distance(p[i], q[j]));

    //get discrete frechet distance
    result= c[p.size()-1][q.size()-1];

    //keep c(i,j) for optimal traversal if we are doing clustering
    if(curve != NULL) {
        //initialize c(i,j) for optimal traversal
       // Mean_c = new long double*[p.size()];
       // for(int i = 0; i < p.size(); ++i)
       //     Mean_c[i] = new long double[q.size()];

        //load c(i,j) for optimal traversal
        for(int i=0; i < p.size(); i++)
            for(int j= 0; j < q.size(); j++) {
                curve->C[i][j] = c[i][j];
            }
    }

    for(int i = 0; i < p.size();++i)
        delete[] c[i];
    delete[] c;

    return result;
}

double continuousFrechet_distanceByFred(Curve &c1, Curve &c2) {
    //set R coordinates of each curve in Fred Curves
    Fred_Points f_p1_arr(c1.get_filtered_coords()->size(),1);
    for(int i=0; i < c1.get_filtered_coords()->size(); i++){
        Fred_Point f_p1 = Fred_Point(1);
        f_p1.set(0,c1.get_filtered_coords()->at(i));
        f_p1_arr.push_back(f_p1);
    }

    Fred_Points f_p2_arr(c2.get_filtered_coords()->size(),1);
    for(int i=0; i < c2.get_filtered_coords()->size(); i++){
        Fred_Point f_p2 = Fred_Point(1);
        f_p2.set(0,c2.get_filtered_coords()->at(i));
        f_p2_arr.push_back(f_p2);
    }

    // Construct Fred Curves
    Fred_Curve f_c1 = Fred_Curve(f_p1_arr,c1.get_id());
    Fred_Curve f_c2 = Fred_Curve(f_p2_arr,c2.get_id());

    //calculate Continuous Frechet Distance using Fred's Algorithm
    return Frechet::Continuous::distance(f_c1, f_c2).value;
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

double approx_factor(double exact, double approx) {
    if(exact == 0 && approx == 0)
        return 1;
    if(exact == 0)
        return DBL_MAX;

    return approx / exact;
}

double max_approx_factor(vector<double> factor_arr) {
    sort(factor_arr.begin(), factor_arr.end());
    return factor_arr.at(factor_arr.size()-1);
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

Curve *vec_avg_curve(vector<pair< pair<double,double>, pair<double,double> >> &reverse_traversal) {
    vector<pair<double,double>> *new_coords = new vector <pair<double,double>>(reverse_traversal.size());
    pair<double,double> new_coord;
    double sum_x, sum_y;

    for(int i = (reverse_traversal.size() - 1); i >= 0; i--) {
        //curve on x
        sum_x = reverse_traversal.at(i).first.first + reverse_traversal.at(i).second.first;
        sum_x = sum_x / 2 ;

        //curve on y
        sum_y = reverse_traversal.at(i).first.second + reverse_traversal.at(i).second.second;
        sum_y = sum_y / 2 ;

        //set and load new mean curve coordinate
        new_coord.first = sum_x;
        new_coord.second = sum_y;
        new_coords->push_back(new_coord);
    }

    //return new Mean Curve for tree traversal
    return new /*Mean_*/Curve(new_coords);
}

ostream& operator<<(ostream& os, const Point& p) {
    os << "[ ";
    for (int i = 0 ; i < p.coords->size() ; i++) {
        os << p.coords->at(i) << ", ";
    }
    os << "]";
    return os;
}
