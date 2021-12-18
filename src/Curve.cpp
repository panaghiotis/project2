#include <vector>
#include <algorithm>
#include <random>
#include <ctime>
#include <cmath>
#include <limits>
#include <cfloat>
#include <cstdlib>
#include "../header/Curve.h"

extern double delta;
extern double epsilon; // ε for Continuous Frechet

Curve::Curve(vector<pair<double, double>> *curve_coords,int dim, string name, vector<double> *ccurve_coords) {
    this->y_dim = dim;
    this->id = name;
    this->R2_coords = curve_coords;
    this->R_coords = ccurve_coords;
    this->coords_arr = new vector<vector<double>>();
    this->filtered_coords = new vector<double>();

}

Curve::~Curve() {
    delete this->coords_arr;
    if(this->R2_coords != NULL)
        delete this->R2_coords;
    if(this->R_coords != NULL)
        delete this->R_coords;
    delete this->filtered_coords;
}

void Curve::Grid_hash(double **t, int max_len) {
    // randomly generate all t from [0, delta)
    uniform_real_distribution<double> dist(0.0, delta);

    // snap to grid
    vector<pair<double,double>> snapping_arr;
    pair<double,double> xy;
    for(int i = 0 ; i < this->R2_coords->size() ; i++) {
        double x = this->R2_coords->at(i).first;
        double y = this->R2_coords->at(i).second;

        // snap x and y
        x = (floor(abs(x-(*t)[0])/delta + 0.5)*delta + (*t)[0]);
        y = (floor(abs(y-(*t)[1])/delta + 0.5)*delta + (*t)[1]);
        xy.first = x;
        xy.second = y;

        //get snapped curve
        snapping_arr.push_back(xy);
    }

    // remove consecutive duplicates and do padding for each duplicate you remove
    int snapped_size = snapping_arr.size();
    pair<double,double> padding_num;
    padding_num.first = DBL_MAX;
    padding_num.second = DBL_MAX;
    for(int i=0; i < snapped_size; i++) {
        if(snapping_arr.at(i) == padding_num)
            break;
        else {
            if(i != (snapped_size-1) && snapping_arr.at(i) == snapping_arr.at(i+1)) {
                snapping_arr.erase(snapping_arr.begin()+i);

                //padding
                snapping_arr.push_back(padding_num);

                //check xy in i position again since we erased previous duplicate in position i
                i = i-1;
            }
        }
    }

    // If we are doing clustering do padding until max length of mean curve for saving in LSH
    if(max_len !=0 && max_len != snapping_arr.size()) {
        if(max_len > snapping_arr.size()) {
            for(int i = (snapping_arr.size() - 1); i < max_len; i++)
                snapping_arr.push_back(padding_num);
        } else if(max_len < snapping_arr.size()) {
            for(int i = snapping_arr.size() - 1; i < max_len; i--)
                snapping_arr.erase(snapping_arr.begin()+i);
        }
    }

    //get 1d coordinates of curve for saving in LSH
    vector<double> R_coords;
    for(int i=0; i < snapping_arr.size(); i++) {
        R_coords.push_back(snapping_arr.at(i).first);
        R_coords.push_back(snapping_arr.at(i).second);
    }
    this->coords_arr->push_back(R_coords);
}

void Curve::R_Filtering() {
    // get R coordinates
    vector<double> temp_coords;
    for(int i=0; i < this->R_coords->size(); i++)
        temp_coords.push_back(this->R_coords->at(i));

    // filter coordinates
    for(int i=0; i < temp_coords.size(); i++) {
        if( i < temp_coords.size() - 2) {
            // |a-b| <= ε and |b-c| <= ε ,remove b
            if(abs(temp_coords.at(i) - temp_coords.at(i+1)) <= epsilon && abs(temp_coords.at(i+1) - temp_coords.at(i+2)) <= epsilon) {
                temp_coords.erase(temp_coords.begin()+(i+1));

                //check i again with i+2 next to it this time
                i = i - 1;
            }
        }
    }
    //save filtered coordinates
    for(int i=0; i < temp_coords.size(); i++) {
        this->filtered_coords->push_back(temp_coords.at(i));
    }
}

void Curve::R_Grid_hash(double &t) {
    // copy filtered coordinates in coordinates to be saved later as Grid keys in LSH
    vector<double> temp_coords;
    for(int i=0; i < this->filtered_coords->size(); i++)
        temp_coords.push_back(this->filtered_coords->at(i));

    // snap to grid
    for(int i=0; i < temp_coords.size(); i++) {
        // snap each x
        double x = floor((temp_coords.at(i) + t) / delta) * delta;
        temp_coords.at(i) = x;
    }

    // key implementation of minima maxima
    vector<double> key;
    for(int i=1; i < temp_coords.size() - 1; i++) {   //start from the second num and finish loop with the second from last
        // { minima , maxima } of vi
        key.push_back(min(min(temp_coords.at(i-1), temp_coords.at(i)), temp_coords.at(i+1)));
        key.push_back(max(max(temp_coords.at(i-1), temp_coords.at(i)), temp_coords.at(i+1)));
    }

    // padding key
    int k_last = key.size() - 1;
    for(int i = k_last; i < ((2*y_dim)-4) ; i++)        // minima maxima gives from a num 2 except from the first and the last one
        key.push_back(DBL_MAX);

    //save key for LSH
    this->coords_arr->push_back(key);
}

vector<vector<double>> *Curve::get_grid_coords() {
    return this->coords_arr;
}

vector<pair<double, double>> *Curve::get_curve_coords() {
    return this->R2_coords;
}

vector<double> *Curve::get_Rcoords() {
    return this->R_coords;
}

vector<double> *Curve::get_filtered_coords() {
    return this->filtered_coords;
}

string Curve::get_id() {
    return this->id;
}

int Curve::get_dimension() {
    return this->y_dim;
}

void Curve::print(bool isCont) {
    if(!isCont)
        for(int i=0; i < 5 ; i++)
            cout << "id: " << id << i <<" -coords: " << R2_coords->at(i).first <<" "<< R2_coords->at(i).second << endl;
    else {
        cout << "id: " << id << " - 1st coord: " << R_coords->at(0) << " - second coord: " << R_coords->at(1) << endl;
    }
}


