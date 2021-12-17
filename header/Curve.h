#ifndef PROJECT2_CURVE_H
#define PROJECT2_CURVE_H

#include <iostream>

using namespace std;


class Curve {
private:
    int y_dim;                                  // dimension of coordinates in y
    string id;                                  // name
    double e;                                   // ε for Continuous Frechet
    vector<pair<double,double>> *R2_coords;     // coordinates in R2 for Discrete Frechet
    vector<double> *R_coords;                   // coordinates in R for Continuous Frechet
    vector<double> *filtered_coords;            // filtered coordinates in R for Continuous Frechet
    vector<vector<double>> *coords_arr;         // coordinates for saving curves in LSH

public:
    //long double **C;
    //vector<vector<long double>> *C;
    Curve(vector<pair<double,double>> *curve_coords, int dim=0, string name= "",  vector<double> *ccurve_coords = NULL);
    ~Curve();
    void Grid_hash(double **t);                 //Grid hash function for Discrete Frechet
    void R_Grid_hash(double &t);                //Grid hash function for Continuous Frechet
    void R_Filtering();                         //Filtering curves for Continuous Frechet
    void R2_Filtering(unsigned int endfilt);    //Filtering Mean Curves for Clustering
    void Remove_duplicates();                   //only removing duplicates just for Clustering. For search duplicates are removed by Grid_hash
    vector<vector<double>> *get_grid_coords();
    vector<pair<double,double>> *get_curve_coords();
    vector<double> *get_Rcoords();
    vector<double> *get_filtered_coords();
    string get_id();
    int get_dimension();
    //set_dimension();
    //void set_C(int mean_size);        //for clustering
    void set_epsilon(int &e);         //for clustering
    void print(bool isCont = false);
};


#endif
