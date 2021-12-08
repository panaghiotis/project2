#include <random>
#include <ctime>
#include <cmath>
#include <bits/stdc++.h>
//#include <string>
#include "../header/Hypercube.h"
//#include "../header/vector_operations.h"

using namespace std;

// random seeded generator seeded once at start of program
default_random_engine coin_generator = default_random_engine(time(NULL));



HyperCube::HyperCube(unsigned int k, unsigned int dim) : k(k) , dim(dim){
    //calculate hypercube size
    this->hyper_size = pow(2,k);
    this->vertices = new Bucket[this->hyper_size];

    // create hash functions
    for (int i = 0 ; i < k ; i++) {
        this->h.push_back(new HashFunction(dim));
    }

    // Initialize hash maps for f functions
    this->f = new unordered_map<unsigned int, int>[k];
}

HyperCube::~HyperCube() {
    delete[] this->vertices;
    delete[] this->f;
    for (int i = 0; i < h.size(); i++) {
        delete h[i];
    }
}

unsigned int HyperCube::hyper_hash(Point &p) {
    uniform_int_distribution<int> dist(0, 1);

    vector<unsigned int> inverted_hashes;
    for (int i = (h.size()-1) ; i >= 0; i--) {
        inverted_hashes.push_back(h[i]->hash(p));  // h_i(p)
    }

    unsigned int ver=0;
    for(int i=0 ; i < this->k ; i++)
    {
        // Check if we have already generated a bit for the d' - i th bucket
        if(this->f[this->k -1 -i].find(inverted_hashes[i]) == this->f[this->k -1 -i].end()){
            //if not generate it now
            this->f[this->k -1 -i][inverted_hashes[i]] = dist(coin_generator);
        }
        ver |= (this->f[this->k -1 -i][inverted_hashes[i]] << i);
    }

    p.hashed = true;

    return ver;
}

void HyperCube::hyper_hash_and_add(Point *p) {
    unsigned int vertex_num = hyper_hash(*p);
    if (vertex_num >= this->hyper_size) cerr << "ERROR in hyper_hash()'s calculation" << endl;
    this->vertices[vertex_num].add(p);
}

Bucket *HyperCube::get_vertex(unsigned int num) {
    if(num >= this->hyper_size) { cerr << "Warning: invalid vertex num arg" << endl; return NULL; }
    return &vertices[num];
}

int HyperCube::get_hyper_dimension() {
    return this->k;
}

void HyperCube::hyper_print() const {
    int count = 0;
    for (int i = 0 ; i < this->hyper_size ; i++) {
        cout << "Vertex[" << i << "] size:" << vertices[i].get_size() << endl;
        count += vertices[i].get_size();
    }
    cout << "Total count:" << count << endl;
}