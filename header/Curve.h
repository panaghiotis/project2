#ifndef PROJECT2_CURVE_H
#define PROJECT2_CURVE_H

#include <iostream>

using namespace std;


class Curve {
private:
    int y_dim;                                  // dimension of coordinates in y
    string id;                                  // name
    vector<pair<double,double>> *R2_coords;     // coordinates in R2 for Discrete Frechet
    vector<double> *R_coords;                   // coordinates in R for Continuous Frechet
    vector<double> *filtered_coords;            // filtered coordinates in R for Continuous Frechet
    vector<vector<double>> *coords_arr;         // coordinates for saving curves in LSH

public:
    Curve(vector<pair<double,double>> *curve_coords, int dim=0, string name= "",  vector<double> *ccurve_coords = NULL);
    ~Curve();
    void Grid_hash(double **t, int max_len=0);  //Grid hash function for Discrete Frechet. Max len is for padding and saving of mean curves
    void R_Grid_hash(double &t);                //Grid hash function for Continuous Frechet
    void R_Filtering();                         //Filtering curves for Continuous Frechet
    vector<vector<double>> *get_grid_coords();
    vector<pair<double,double>> *get_curve_coords();
    vector<double> *get_Rcoords();
    vector<double> *get_filtered_coords();
    string get_id();
    int get_dimension();
    void print(bool isCont = false);
};


#endif
