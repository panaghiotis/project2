#ifndef PROJECT2_CURVE_H
#define PROJECT2_CURVE_H

#include <iostream>

using namespace std;


class Curve {
private:
    int y_dim;
    string id;
    vector<pair<int,double>> *R2_coords;
    vector<vector<double>> *coords_arr;
public:
    Curve(int dim, string name, vector<pair<int,double>> *curve_coords);
    ~Curve();
    void Grid_hash(double **t);
    vector<vector<double>> *get_grid_coords();
    vector<pair<int,double>> *get_curve_coords();
    string get_id();
    void print();
};


#endif
