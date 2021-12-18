#include <iostream>
#include "../header/Dataset.h"


extern unsigned int k;
extern unsigned int L;


Dataset::Dataset(ifstream &file, bool isFrechet, bool isCont) : size(0), dim(0) {
    string token;
    int i = 0;
    while (file >> token) {
        // read line tokens one-by-one
        string ID = token;
        vector<double> *coords = new vector<double>();
        while ((file.peek() != '\n' && file.peek() != '\r') && (file >> token)) {   // '\r' for windows, '\n' for unix
            coords->push_back((double) atof(token.c_str()));   // convert to real numbers
            if (file.peek() == '\t' || file.peek() == ' ') {   // ignore these
                file.ignore(1);
            }
        }
        if(isFrechet) {
            if (!isCont) {
                vector<pair<double, double>> *curve_coords = new vector<pair<double, double>>;
                pair<double, double> coord_2d;
                for (int i = 0; i < coords->size(); i++) {
                    coord_2d.first = i + 1;
                    coord_2d.second = coords->at(i);
                    curve_coords->push_back(coord_2d);
                }
                curves.push_back(new Curve(curve_coords,coords->size(), ID));

                //Grid hash and save in 1d LSH as new points only used for saving curves
                dim = 2 * (coords->size());     //TODO: Is that needed?
                //free memory
                delete coords;
            } else {
                dim = coords->size();   // should be the same for every curve
                curves.push_back(new Curve(NULL,dim,ID,coords));
                curves.at(curves.size()-1)->R_Filtering();
            }
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


void Dataset::print(bool isFrechet, bool isCont) const {
    if(isFrechet) {
        if(isCont) {
            for(int i=0; i < curves.size(); i++) {
                curves[i]->print(true);
            }
        }
        else {
            for(int i=0; i < curves.size(); i++) {
                curves[i]->print();
            }
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

void Dataset::index_LSH(unsigned int hashtable_size, bool isFrechet, bool isCont, int max_len) { //max len is the max length of mean curves
    // init hashtables
    if(isFrechet){
        for (int i = 0 ; i < L ; i++) {
            hashTables.push_back(new HashTable(hashtable_size, k, dim));

            // Grid hash for Discrete or Continuous Frechet
            for(int j=0; j < curves.size(); j++) {
                if(!isCont){// R2 Grid
                    if(max_len !=0) // we are doing clustering
                        curves.at(j)->Grid_hash((hashTables.at(i)->get_Frechet_t()), max_len);
                    else
                        curves.at(j)->Grid_hash((hashTables.at(i)->get_Frechet_t()));
                }
                else       // R Grid
                    curves.at(j)->R_Grid_hash(*hashTables.at(i)->get_Frechet_t()[0]);
            }
        }

        // Save into L LSH tables as new points only used for saving curves
        for(int i=0; i < curves.size(); i++) {
            for(int j = 0 ; j < L ; j++){
                points.push_back(new Point(&curves.at(i)->get_grid_coords()->at(j), curves.at(i)->get_id(), i, curves.at(i)));
                hashTables[j]->hash_and_add(points.at(points.size()-1));
            }
        }
    }
    else{ //curves as points
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
}

Bucket **Dataset::get_buckets_for_point(Point *p) {
    Bucket **result = new Bucket*[hashTables.size()];
    for (int i = 0 ; i < hashTables.size() ; i++) {
        unsigned int bucket_num = hashTables[i]->g(*p);
        result[i] = hashTables[i]->get_bucket(bucket_num);
    }
    return result;
}

Bucket *Dataset::get_bucket_for_curve(Point *p, int num) {
    unsigned int bucket_num = hashTables[num]->g(*p);
    return hashTables[num]->get_bucket(bucket_num);
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

double **Dataset::get_Fr_t_of_htable(int num) {
    return this->hashTables.at(num)->get_Frechet_t();
}

void Dataset::print_LSH() {
    for(int i=0; i < hashTables.size(); i++) {
        cout << "Hash Table numbber " << i << " :" << endl;
        hashTables[i]->print();
    }
}