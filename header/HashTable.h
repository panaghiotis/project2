#ifndef HASHTABLE_H
#define HASHTABLE_H

#include <vector>
#include "Point.h"


#define DEFAULT_W (3*500)  // in [2, 6]?
#define m 4294967291       // 2^32 - 5

struct Bucket {
    vector<Point *> points;
    void add(Point *p) {
        points.push_back(p);
    }
    unsigned int get_size() const {
        return points.size();
    }
};


class HashFunction {
    Point *v;
    double t, w;
public:
    HashFunction(unsigned int dim, double w = DEFAULT_W);
    ~HashFunction();
    unsigned int hash(const Point &p) const; // maybe unsigned long??
};


class HashTable {
    unsigned int hashtable_size;
    Bucket *buckets;
    vector<unsigned int> r;
    vector<HashFunction *> h;     // h.size() == r.size() == # num of buckets
public:
    HashTable(unsigned int size, unsigned int k, unsigned int dim);
    ~HashTable();
    unsigned int g(Point &p);     // hash into ints in [0, hashtable_size - 1]
    void hash_and_add(Point *p);
    Bucket *get_bucket(unsigned int num);
    void print() const;
};

#endif
