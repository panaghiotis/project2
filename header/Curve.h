#ifndef PROJECT2_CURVE_H
#define PROJECT2_CURVE_H

#include <iostream>

using namespace std;


class Curve {
private:
    int y_dim;                                  // dimension of coordinates in y
    string id;                                  // name
    double e;                                   // ε for Continuous Frechet
    vector<pair<int,double>> *R2_coords;        // coordinates in R2 for Discrete Frechet
    vector<double> *R_coords;                   // coordinates in R for Continuous Frechet
    vector<vector<double>> *coords_arr;         // coordinates for saving curves in LSH
public:
    Curve(int dim, string name, vector<pair<int,double>> *curve_coords, vector<double> *ccurve_coords = NULL);
    ~Curve();
    void Grid_hash(double **t);                 //Grid hash function for Discrete Frechet
    void R_Grid_hash();                         //Grid hash function for Continuous Frechet
    vector<vector<double>> *get_grid_coords();
    vector<pair<int,double>> *get_curve_coords();
    string get_id();
    void print(bool isCont = false);
};


#endif
