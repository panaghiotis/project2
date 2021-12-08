#include <iostream>
#include "../header/Dataset.h"


extern unsigned int k;
extern unsigned int L;


Dataset::Dataset(ifstream &file, bool isFrechet) : size(0), dim(0) {
    string token;
    int i = 0;
    while (file >> token) {
        // read line tokens one-by-one
        string ID = token;
        vector<double> *coords = new vector<double>();
        while ((file.peek() != '\n' && file.peek() != '\r') && (file >> token)) {   // '\r' for windows, '\n' for unix
            coords->push_back((double) atoi(token.c_str()));   // convert to real numbers
            if (file.peek() == '\t' || file.peek() == ' ') {   // ignore these
                file.ignore(1);
            }
        }
        if (isFrechet) {
            vector<pair<int, double>> *curve_coords = new vector<pair<int, double>>;
            pair<int, double> coord_2d;
            for (int i = 0; i < coords->size(); i++) {
                coord_2d.first = i + 1;
                coord_2d.second = coords->at(i);
                curve_coords->push_back(coord_2d);
            }
            curves.push_back(new Curve(coords->size(), ID, curve_coords));

            //Grid hash and save in 1d LSH as new points only used for saving curves
            curves.at(curves.size() - 1)->Grid_hash();
            dim = 2 * (coords->size());
            points.push_back(
                    new Point(curves.at(curves.size() - 1)->get_grid_coords(), ID, i, curves.at(curves.size() - 1)));
        } else {
            dim = coords->size();   // should be the same for every point
            points.push_back(new Point(coords, ID, i));        // save pos i for later
        }
        i++;
        // ignore anything until '\n'
        file.ignore(1024, '\n');
    }
    size = i;
    cout << "Loaded " << size << " points." << endl;
    file.close();
}


void Dataset::print(bool isFrechet) const {
    if(isFrechet) {
        for(int i=0; i < curves.size(); i++) {
            curves[i]->print();
        }
    }
    else {
        for (int i = 0; i < points.size(); i++){
            points[i]->print();
        }
    }
}

Dataset::~Dataset() {
    for (int i = 0; i < points.size(); i++){
        delete points[i];
    }
    for(int i = 0; i < curves.size(); i++){
        delete curves[i];
    }
    for (int i = 0; i < hashTables.size(); i++){
        delete hashTables[i];
    }
}

void Dataset::index_LSH(unsigned int hashtable_size) {
    // init hashtables
    for (int i = 0 ; i < L ; i++) {
        hashTables.push_back(new HashTable(hashtable_size, k, dim));
    }
    // hash all the points into the hash tables
    for (int i = 0 ; i < points.size() ; i++) {
        for (int j = 0 ; j < L ; j++) {
            hashTables[j]->hash_and_add(points[i]);
        }
    }
}

Bucket **Dataset::get_buckets_for_point(Point *p) {
    Bucket **result = new Bucket*[hashTables.size()];
    for (int i = 0 ; i < hashTables.size() ; i++) {
        unsigned int bucket_num = hashTables[i]->g(*p);
        result[i] = hashTables[i]->get_bucket(bucket_num);
    }
    return result;
}

void Dataset::index_HyperCube(unsigned int k) {
    this->hypercube = new HyperCube(k,dim);
    // hash all the points into hypercube
    for (int i = 0 ; i < points.size() ; i++){
        hypercube->hyper_hash_and_add(points[i]);
    }
}

Bucket *Dataset::get_vertex_for_point(Point *p) {
    unsigned int vertex_num = hypercube->hyper_hash(*p);
    return hypercube->get_vertex(vertex_num);
}

Bucket *Dataset::get_vertex_by_number(int num) {
    return hypercube->get_vertex(num);
}

unsigned int Dataset::get_vertex_number(Point *p) {
    return hypercube->hyper_hash(*p);
}

int Dataset::get_dimension_from_cube() {
    return hypercube->get_hyper_dimension();
}