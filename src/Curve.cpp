#include <vector>
#include <algorithm>
#include <random>
#include <ctime>
#include <cmath>
#include <limits>
#include <float.h>
#include "../header/Curve.h"

extern double delta; // 0.1 or 0.01

// random seeded generator for discrete curves
default_random_engine generator_2d = default_random_engine(time(NULL));


Curve::Curve(int dim, string name, vector<pair<int, double>> *curve_coords) {
    this->y_dim = dim;
    this->id = name;
    this->R2_coords = curve_coords;
    this->coords = new vector<double>();
}

Curve::~Curve() {
    delete this->coords;
    delete this->R2_coords;
}

void Curve::Grid_hash() {
    // randomly generate all t from [0, delta)
    uniform_real_distribution<double> dist(0.0, delta);
    vector<double> *t = new vector<double>();
    int t_size = 2*(this->y_dim);
    for (int i = 0 ; i < t_size ; i++) {
        t->push_back(dist(generator_2d));
    }

    // snap to grid
    int t_placement=0;
    vector<pair<double,double>> snapping_arr;
    pair<double,double> xy;
    for(int i = 0 ; i < this->y_dim ; i++) {
        double x = (double) this->R2_coords->at(i).first;
        double y = this->R2_coords->at(i).second;

        // snap x and y
        x = (floor((x-(t->at(t_placement)))/delta + 0.5)*delta + t->at(t_placement));
        t_placement++;
        y = (floor((y-(t->at(t_placement)))/delta + 0.5)*delta + t->at(t_placement));
        t_placement++;
        xy.first = x;
        xy.second = y;

        //get snapped curve
        snapping_arr.push_back(xy);
    }

    // remove consecutive duplicates and do padding for each duplicate you remove
    int snapped_size = snapping_arr.size();
    pair<double,double> padding_num;
    padding_num.first = DBL_MAX;                //LDBL_MAX
    padding_num.second = DBL_MAX;
    for(int i=0; i < snapped_size; i++) {
        if(snapping_arr.at(i) == padding_num)
            break;
        else {
            if(i != (snapped_size-1) && snapping_arr.at(i) == snapping_arr.at(i+1)) {
                snapping_arr.erase(snapping_arr.begin()+i);

                //padding
                snapping_arr.push_back(padding_num);

                //check xy in i position again since we erased previous duplicate in postition i
                i = i-1;
            }
        }
    }

    //get 1d coordinates of curve for saving in LSH
    for(int i=0; i < snapping_arr.size(); i++) {
        this->coords->push_back(snapping_arr.at(i).first);
        this->coords->push_back(snapping_arr.at(i).second);
    }
}

vector<double> *Curve::get_grid_coords() {
    return this->coords;
}

void Curve::print() {
    cout << "id: " << id << " - 1st coord: " << R2_coords->at(0).first <<" "<< R2_coords->at(0).second <<" "<< coords->at(0)  << endl;
}


