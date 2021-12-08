#include <unordered_map>
#include <list>
#include <string>
#include "HashTable.h"
//#include "../headers/dataset.h"

class HyperCube
{
    private:
        int k,w,dim;              // k is dimension d' , dim is coords dimension of a point
        vector<HashFunction *> h;
        std::unordered_map<unsigned int, int> *f;
        unsigned int hyper_size;
        Bucket *vertices;
    public:
        HyperCube(unsigned int k, unsigned int dim);
        unsigned int hyper_hash(Point &p);
        void hyper_hash_and_add(Point *p);
        Bucket *get_vertex(unsigned int num);
        int get_hyper_dimension();
        void hyper_print() const;
        ~HyperCube();
};