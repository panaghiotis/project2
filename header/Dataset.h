#ifndef DATASET_H
#define DATASET_H

#include <vector>
#include <fstream>
#include "Point.h"
//#include "HashTable.h"
#include "Hypercube.h"


using namespace std;


class Dataset{
public:
    vector<Point *> points;
    vector<Curve *> curves;
private:
    unsigned int dim, size;
    vector<HashTable *> hashTables;   // empty if not called index
    HyperCube *hypercube;
public:
    Dataset() {};
    Dataset(ifstream &file, bool isFrechet = false, bool isCont = false);
    ~Dataset();
    void print(bool isFrechet = false, bool isCont = false) const;
    void index_LSH(unsigned int hashtable_size, bool isFrechet = false, bool isCont = false, int max_len = 0);
    Bucket **get_buckets_for_point(Point *p);
    Bucket *get_bucket_for_curve(Point *p, int num);
    unsigned int get_dim() const { return dim; }
    void set_dim(unsigned int d) { dim = d; }
    unsigned int get_size() const { return size; }
    void set_size(unsigned int sz) { size = sz; }
    unsigned int get_hashtables_count() const { return hashTables.size(); }
    void index_HyperCube(unsigned int k);
    Bucket *get_vertex_for_point(Point *p);
    Bucket *get_vertex_by_number(int num);
    unsigned int get_vertex_number(Point *p);
    int get_dimension_from_cube();
    double **get_Fr_t_of_htable(int num);
    void print_LSH();
};


#endif
