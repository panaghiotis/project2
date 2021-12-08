#include <vector>
#include <chrono>
#include <algorithm>
#include <unordered_set>
#include "../header/neighbour_search.h"
#include "../header/vector_operations.h"
#define INF 1.0/0.0

using namespace std;


NearestNeighboursSearch::NearestNeighboursSearch(Dataset &index_dataset, Dataset &query_dataset)
        : index_dataset(index_dataset), query_dataset(query_dataset), exact_distances(NULL), lsh_distances(NULL), cube_distances(NULL) {

}

NearestNeighboursSearch::~NearestNeighboursSearch() {
    if (exact_distances != NULL) {
        for (int i = 0; i < query_dataset.points.size(); i++) {
            delete[] exact_distances[i];
        }
        delete[] exact_distances;
    }
    if (lsh_distances != NULL) {
        for (int i = 0; i < query_dataset.points.size(); i++) {
            delete[] lsh_distances[i];
        }
        delete[] lsh_distances;
    }
    if (cube_distances != NULL) {
        for (int i = 0; i < query_dataset.points.size(); i++) {
            delete[] cube_distances[i];
        }
    }
}

long long int NearestNeighboursSearch::calculate_exact_distances() {
    unsigned int q_size = query_dataset.points.size();
    unsigned int ind_size = index_dataset.points.size();
    exact_distances = new long double*[q_size];
    for (int i = 0 ; i < q_size ; i++) {
        exact_distances[i] = new long double[ind_size];
    }
    chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    for (int i = 0 ; i < q_size ; i++) {
        for (int j = 0 ; j < ind_size ; j++) {
            exact_distances[i][j] = used_distance(*query_dataset.points.at(i), *index_dataset.points.at(j));
        }
    }
    chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    return chrono::duration_cast<chrono::microseconds>(end - begin).count();
}

long long int NearestNeighboursSearch::calculate_only_lsh_distances() {
    unsigned int q_size = query_dataset.points.size();
    unsigned int ind_size = index_dataset.points.size();
    lsh_distances = new long double*[q_size];
    for (int i = 0 ; i < q_size ; i++) {
        lsh_distances[i] = new long double[ind_size];
        for (int j = 0 ; j < ind_size ; j++) {
            lsh_distances[i][j] = -1.0;   // not calculated
        }
    }
    // calculate and fill only distances between points that hash into the same bucket for the L hash tables
    chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    for (int i = 0 ; i < q_size ; i++) {
        Point *p1 = query_dataset.points.at(i);
        Bucket **buckets = index_dataset.get_buckets_for_point(p1);           // L buckets
        for (int j = 0 ; j < index_dataset.get_hashtables_count(); j++) {
            for (int k = 0 ; k < buckets[j]->points.size() ; k++) {           // might be empty
                const Point *p2 = buckets[j]->points[k];
                if (p1->hashed_ID == p2->hashed_ID) {                         // query trick
                    // update appropriate distance at [i, p2->pos]
                    lsh_distances[i][p2->pos] = used_distance(*p1, *p2);
                }
            }
        }
        delete[] buckets;
    }
    chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    return chrono::duration_cast<chrono::microseconds>(end - begin).count();
}

vector<Result *> *NearestNeighboursSearch::calculate_exact_NN(int k) {
    vector<Result *> *results = new vector<Result *>();
    unsigned int q_size = query_dataset.points.size();
    unsigned int ind_size = index_dataset.points.size();
    for (int i = 0 ; i < q_size ; i++) {
        // create new result entry
        Result *res = new Result();

        // time starts
        chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

        // calculate all distances
        vector<pair<string, long double>> dists;
        dists.reserve(ind_size);
        for (int j = 0 ; j < ind_size ; j++) {
            dists.emplace_back(index_dataset.points.at(j)->id, used_distance(*query_dataset.points.at(i), *index_dataset.points.at(j)));
        }

        // sort distances with their ids
        std::sort(dists.begin(), dists.end(), sort_pred());

        // time ends
        chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        long long int dt = chrono::duration_cast<chrono::microseconds>(end - begin).count();

        // add kNNs to the result
        vector<Neighbour> *NNs = new vector<Neighbour>();
        for (int j = 0 ; j < k && j < dists.size() ; j++) {
            NNs->emplace_back(dists[j].first, dists[j].second);
        }

        // set result fields
        res->queryID = query_dataset.points[i]->id;
        res->NNs = NNs;
        res->dt = dt;
        results->push_back(res);
    }
    return results;
}

vector<Result *> *NearestNeighboursSearch::calculate_approximate_NN(int k,  bool use_query_trick) {
    vector<Result *> *results = new vector<Result *>();
    unsigned int q_size = query_dataset.points.size();
    for (int i = 0 ; i < q_size ; i++) {
        // create new result entry
        Result *res = new Result();

        // time starts
        chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

        // calculate only specific distances according to LSH
        vector<pair<string, long double>> dists;
        Point *p1 = query_dataset.points.at(i);
        Bucket **buckets = index_dataset.get_buckets_for_point(p1);           // L buckets
        unordered_set<string> seen;                                           // multiple hash tables -> may find an element multiple times
        for (int j = 0 ; j < index_dataset.get_hashtables_count(); j++) {
            for (int z = 0 ; z < buckets[j]->points.size() ; z++) {           // might be empty
                const Point *p2 = buckets[j]->points[z];
                if (!use_query_trick || (p1->hashed && p2->hashed && p1->hashed_ID == p2->hashed_ID)) {     // query trick
                    if (seen.find(p2->id) == seen.end()) {                    // if not already seen this point
                        dists.emplace_back(p2->id, used_distance(*p1, *p2));
                        seen.insert(p2->id);
                    }
                }
            }
        }
        delete[] buckets;

        // sort distances with their ids
        std::sort(dists.begin(), dists.end(), sort_pred());

        // time ends
        chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        long long int dt = chrono::duration_cast<chrono::microseconds>(end - begin).count();

        // add kNNs to the result
        vector<Neighbour> *NNs = new vector<Neighbour>();
        for (int j = 0 ; j < k && j < dists.size() ; j++) {     // k might be > dists.size() if no neighbors
            NNs->emplace_back(dists[j].first, dists[j].second);
        }

        // set result fields
        res->queryID = query_dataset.points[i]->id;
        res->NNs = NNs;
        res->dt = dt;
        results->push_back(res);
    }
    return results;
}

vector<vector<string> *> *NearestNeighboursSearch::range_search(long double r, bool use_query_trick) {
    vector<vector<string> *> *combined_result = new vector<vector<string> *>();
    unsigned int q_size = query_dataset.points.size();
    for (int i = 0 ; i < q_size ; i++) {
        vector<string> *result = new vector<string>();
        // calculate only specific distances according to LSH
        Point *p1 = query_dataset.points.at(i);
        Bucket **buckets = index_dataset.get_buckets_for_point(p1);           // L buckets
        unordered_set<string> seen;                                           // multiple hash tables -> may find an element multiple times
        for (int j = 0 ; j < index_dataset.get_hashtables_count(); j++) {
            for (int z = 0 ; z < buckets[j]->points.size() ; z++) {           // might be empty
                const Point *p2 = buckets[j]->points[z];
                if (!use_query_trick || (p1->hashed && p2->hashed && p1->hashed_ID == p2->hashed_ID)) {     // query trick
                    if (seen.find(p2->id) == seen.end()) {                    // if not already seen this point
                        if (used_distance(*p1, *p2) < r) {                    // if in range
                            result->push_back(p2->id);
                        }
                        seen.insert(p2->id);
                    }
                }
            }
        }
        delete[] buckets;
        combined_result->push_back(result);
    }
    return combined_result;
}

//range search for clustering
list<Point*> NearestNeighboursSearch::LSH_search(long double r, Point *q, bool use_query_trick) {
    list<Point *> result;
    // calculate only specific distances according to LSH
    Bucket **buckets = index_dataset.get_buckets_for_point(q);           // L buckets
    unordered_set<string> seen;                                           // multiple hash tables -> may find an element multiple times
    for (int j = 0 ; j < index_dataset.get_hashtables_count(); j++) {
        for (int z = 0 ; z < buckets[j]->points.size() ; z++) {           // might be empty
            Point *p = buckets[j]->points[z];
            if (!use_query_trick || (q->hashed && p->hashed && q->hashed_ID == p->hashed_ID)) {     // query trick
                if (seen.find(p->id) == seen.end()) {                                                // if not already seen this point
                    if (used_distance(*q, *p) < r) {                                                 // if in range
                        result.push_back(p);
                    }
                    seen.insert(p->id);
                }
            }
        }
    }
    delete[] buckets;
    return result;
}


list<Point*> NearestNeighboursSearch::cube_search(unsigned int M, int probes, int cur_ver, int hamming, long double r, Point *q) {
    // Base case
    if (probes <= 0) {
        list<Point*> emptyList;
        return emptyList;
    }
    // Search current vertex
   list<Point*> result;
    for (int i = 0; i < index_dataset.get_vertex_by_number(cur_ver)->get_size(); i++) {
        // Check if current point distance to q lies in range R
        Point *p = index_dataset.get_vertex_by_number(cur_ver)->points[i];
        if (used_distance(*q,*p) <= r) {
            result.push_back(p);
        }
    }

    if (result.size() > M) {
        result.resize(M);
    }
    M -= result.size();
    probes--;

    // Recursively search non searched neighbour vertices with hamming distance
    for (int i = 0; M > 0 && probes > 0 && i < index_dataset.get_dimension_from_cube() ; i++) {
        // Check if ith bit was not previously hammed
        if (!(hamming & (1 << i))) {
            list<Point*> rec = cube_search(M, probes - index_dataset.get_dimension_from_cube(), cur_ver ^ (1 << i), hamming | (1 << i), r, q);
            M -= rec.size();
            probes--;
            result.splice(result.end(), rec);
        }
    }

    return result;
}

vector<Result *> *NearestNeighboursSearch::calculate_approximate_NN_using_cube(int k, int probes, int M) {
    vector<Result *> *results = new vector<Result *>();
    unsigned int q_size = query_dataset.points.size();
    for (int i = 0 ; i < q_size ; i++) {
        // create new result entry
        Result *res = new Result();

        // time starts
        chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

        vector<pair<string, long double>> dists;
        Point *q = query_dataset.points.at(i);
        list<Point*> p_to_search = this->cube_search(M,probes, index_dataset.get_vertex_number(q),0,INF,q);
        for(list<Point*>::iterator it = p_to_search.begin(); it != p_to_search.end(); it++)
            dists.emplace_back((*it)->id, used_distance(*q,*(*it)));


        // sort distances with their ids
        std::sort(dists.begin(), dists.end(), sort_pred());

        // time ends
        chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        long long int dt = chrono::duration_cast<chrono::microseconds>(end - begin).count();

        // add kNNs to the result
        vector<Neighbour> *NNs = new vector<Neighbour>();
        for (int j = 0 ; j < k && j < dists.size() ; j++) {
            NNs->emplace_back(dists[j].first, dists[j].second);
        }

        // set result fields
        res->queryID = query_dataset.points[i]->id;
        res->NNs = NNs;
        res->dt = dt;
        results->push_back(res);
    }
    return results;
}

vector<vector<string> *> *NearestNeighboursSearch::range_search_using_cube(long double r, int probes, int M) {
    vector<vector<string> *> *combined_result = new vector<vector<string> *>();
    unsigned int q_size = query_dataset.points.size();
    for (int i = 0 ; i < q_size ; i++) {
        vector<string> *result = new vector<string>();
        Point *q = query_dataset.points.at(i);
        list<Point*> p_to_search = this->cube_search(M,probes, index_dataset.get_vertex_number(q),0,r,q);
        for(list<Point*>::iterator it = p_to_search.begin(); it != p_to_search.end(); it++) {
            result->push_back((*it)->id);
        }
        combined_result->push_back(result);
    }
    return  combined_result;
}
