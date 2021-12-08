#ifndef PROJECT2_CURVE_H
#define PROJECT2_CURVE_H

#include <iostream>

using namespace std;


class Curve {
private:
    int y_dim;
    string id;
    vector<pair<int,double>> *R2_coords;
    vector<double> *coords;
public:
    Curve(int dim, string name, vector<pair<int,double>> *curve_coords);
    ~Curve();
    void Grid_hash();
    vector<double> *get_grid_coords();
    void print();
};


#endif
