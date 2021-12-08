#include <random>
#include <ctime>
#include <cmath>
#include "../header/HashTable.h"
#include "../header/vector_operations.h"

using namespace std;


// random seeded generator seeded once at start of program
default_random_engine generator = default_random_engine(time(NULL));


HashFunction::HashFunction(unsigned int dim, double w) : w(w) {
    // randomly generate vector v from U(0, 1)
    uniform_real_distribution<double> dist(0.0, 1.0);
    vector<double> *coords = new vector<double>();
    for (int i = 0 ; i < dim ; i++) {
        coords->push_back(dist(generator));
    }
    v = new Point(coords);
    // pick t from [0, w)
    uniform_real_distribution<double> dist2(0.0, w);
    t = dist2(generator);
}

HashFunction::~HashFunction() { delete v; }

unsigned int HashFunction::hash(const Point &p) const {
    return floor((dot_product(p, *v) + t) / w);
}


HashTable::HashTable(unsigned int size, unsigned int k, unsigned int dim) : hashtable_size(size) {
    buckets = new Bucket[size];
    // generate r vector
    uniform_int_distribution<unsigned int> dist(0, UINT32_MAX );
    for (int i = 0 ; i < k ; i++) {
        r.push_back(dist(generator));
    }
    // create hash functions
    for (int i = 0 ; i < k ; i++) {
        h.push_back(new HashFunction(dim));
    }
}

HashTable::~HashTable() {
    delete[] buckets;
    for (int i = 0; i < h.size(); i++) {
        delete h[i];
    }
}

unsigned int HashTable::g(Point &p) {
    vector<unsigned int> hashes;
    for (int i = 0; i < h.size(); i++) {
        hashes.push_back(h[i]->hash(p));  // h_i(p)
    }
    unsigned int hashedID = int_dot_product(r, hashes) % m;
    // store these for later trick
    p.hashed_ID = hashedID;
    p.hashed = true;
    // return bucked number
    return hashedID % hashtable_size;
}

void HashTable::hash_and_add(Point *p) {
    unsigned int bucket_num = g(*p);
    if (bucket_num >= hashtable_size) cerr << "ERROR in g()'s calculation" << endl;
    buckets[bucket_num].add(p);
}

Bucket *HashTable::get_bucket(unsigned int num) {
    if (num >= hashtable_size) { cerr << "Warning: invalid bucket num arg" << endl; return NULL; }
    return &buckets[num];
}

void HashTable::print() const {
    int count = 0;
    for (int i = 0 ; i < hashtable_size ; i++) {
        cout << "Bucket[" << i << "] size:" << buckets[i].get_size() << endl;
        count += buckets[i].get_size();
    }
    cout << "Total count:" << count << endl;
}
