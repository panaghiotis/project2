#include <iostream>
#include <random>
#include <cfloat>
//#include <unordered_map>
#include <algorithm>
#include <unordered_set>
#include "../header/clustering.h"
//#include "../header/neighbour_search.h"
#include "../header/vector_operations.h"

using namespace std;


extern default_random_engine generator;   // initialized from another .cpp


template<typename T>
vector<T> *create_copy(vector<T> const &vec){         // to copy coords and not have double deletes
    vector<T> *v = new vector<T>(vec.size());
    copy(vec.begin(), vec.end(), v->begin());
    return v;
}

Clustering::~Clustering() {
    for (int i = 0 ; i < clusters.size() ; i++) {
        delete clusters[i];
    }
}

void Clustering::initialize_clusters(unsigned int k) {    // initialization++
    // Step 1: pick uniformly random point in data:
    srand(time(NULL));  // seed
    int indx = rand() % data->get_size();
    Centroid *c = new Centroid(create_copy(*data->points[indx]->coords));
    centroids.push_back(c);
    clusters.push_back(new Cluster(c));

    unordered_set<unsigned int> ignore_indexes;
    ignore_indexes.insert(indx);

    // Step 2-k: choose the other centers by chance depending on their distance from already selected centers
    for (int j = 1 ; j < k ; j++) {
        vector<long double> D;
        D.reserve(data->points.size());
        // calculate D(i) for each unassigned point i
        for (int i = 0 ; i < data->points.size() ; i++) {
            if (ignore_indexes.find(i) != ignore_indexes.end()) {   // if index already picked
                D.push_back(0.0);                                   // assign zero probabilty
            } else {
                long double minDist = FLT_MAX;
                for (int z = 0 ; z < j ; z++) {    // for already picked centroids
                    long double dist = used_distance(*data->points[i], *centroids[z]);
                    if (dist < minDist) {
                        minDist = dist;
                    }
                }
                D.push_back(pow(minDist, 2));    // no need to normalize, dist does it for us, only put D(i)^2
            }
        }
        // create discrete dist based on D(i)
        discrete_distribution<size_t> dist{D.begin(), D.end()};
        // assign new centroid
        unsigned int selected_index = dist(generator);
        ignore_indexes.insert(selected_index);                 // ignore this in the future
        c = new Centroid(create_copy(*data->points[selected_index]->coords));
        // add to clusters
        centroids.push_back(c);
        clusters.push_back(new Cluster(c));
    }
}


double Clustering::perform_kMeans(unsigned int k, unsigned int M, unsigned int probes) {
    chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    // initialize clusters and centroids
    cout << "Initializing clusters..." << endl;
    this->initialize_clusters(k);
    cout << "Done!" << endl;

    double R;
    if(method == USE_RANGE_LSH || method == USE_HYPERCUBE)
        R = initialize_radius();

    int count = 0;
    while (true) {
        count++;

        // clear any previous clusters  (!) important
        if( method == CLASSIC) {
            for (int i = 0; i < clusters.size(); i++) {
                clusters[i]->clear_cluster();
            }
        }

        cout << "E-step" << endl;
        // E-step: Assign each point to nearest centroid's cluster
        unsigned int *assign_counts = new unsigned int[k];
        for (int i = 0 ; i < k ; i++) {
            assign_counts[i] = 0;
        }
        switch (method){
            case CLASSIC:
                // exact NNs
                for (int i = 0 ; i < data->get_size() ; i++){    // for each point
                    long double minDist = FLT_MAX;
                    int minCluster = -1;
                    for (int j = 0 ; j < k ; j++) {     // calculate all distances and keep the minimum
                        long double dist = used_distance(*data->points[i], *centroids[j]);
                        if (dist < minDist) {
                            minDist = dist;
                            minCluster = j;
                        }
                    }
                    clusters[minCluster]->add(data->points[i]);   // assign to nearest cluster
                    assign_counts[minCluster]++;
                }
                break;
            case USE_RANGE_LSH:
            {
                for(int i=0 ; i < k; i++) {     // for every cluster get points within range
                    list<Point*> tempPoints;
                    tempPoints = search_functions.LSH_search(R,centroids[i]);
                    if(!tempPoints.empty()) {   //avoid buckets with very few points
                        for (list<Point *>::iterator it = tempPoints.begin(); it != tempPoints.end(); it++) {
                            //if we just started the algorithm
                            if(clusters[i]->point_map.empty()) {
                                if(i == 0) {                            //add point to the first cluster
                                    clusters[i]->add((*it));
                                    assign_counts[i]++;
                                }else {                                 //else check if point is already assigned in another cluster
                                    int mark_assign = 0;                //flag already assigned point
                                    for( int j=0; j < k; j++) {
                                        //if cluster isn't empty check if point is assigned in this cluster
                                        if(!clusters[j]->point_map.empty()) {
                                            if (clusters[j]->point_map.find((*it)->id) != clusters[j]->point_map.end()) {
                                                //if it is assigned mark it and check distances to choose if it changes clusters
                                                mark_assign = 1;
                                                long double dist_i = used_distance(*(*it), *centroids[i]);
                                                long double dist_j = used_distance(*(*it), *centroids[j]);
                                                if (dist_i < dist_j) {
                                                    clusters[i]->add((*it));
                                                    clusters[j]->point_map.erase((*it)->id);
                                                    assign_counts[i]++;
                                                }
                                            }
                                        }
                                    }
                                    if(!mark_assign) {
                                        //unassigned point by all clusters, assign it
                                        clusters[i]->add((*it));
                                        assign_counts[i]++;
                                    }
                                }
                            }else {
                                //check if current point is assigned in cluster i
                                if( clusters[i]->point_map.find((*it)->id) == clusters[i]->point_map.end()) {
                                    int mark_assign = 0;
                                    for ( int j=0; j < k; j++) {
                                        //if it is unassigned check if it is assigned in other clusters
                                        if(j != i) {
                                            if( clusters[j]->point_map.find((*it)->id) != clusters[j]->point_map.end()) {
                                                // if it is assigned in another cluster check its distance between each centroid
                                                // and assign point in cluster with the closest centroid
                                                mark_assign = 1;
                                                long double dist_i = used_distance(*(*it), *centroids[i]);
                                                long double dist_j = used_distance(*(*it), *centroids[j]);
                                                if (dist_i < dist_j) {
                                                    clusters[i]->add((*it));
                                                    clusters[j]->point_map.erase((*it)->id);
                                                    assign_counts[i]++;
                                                }
                                            }
                                        }
                                    }
                                    if(!mark_assign) {
                                        //unassigned point by all clusters, assign it
                                        clusters[i]->add((*it));
                                        assign_counts[i]++;
                                    }
                                }
                            }
                        }
                    }
                }
            }
                break;
            case USE_HYPERCUBE:
            {
                for(int i=0 ; i < k; i++) {     // for every cluster get points within range
                    list<Point*> tempPoints;
                    tempPoints = search_functions.cube_search(M,probes,data->get_vertex_number(centroids[i]),0,R,centroids[i]);
                    if(!tempPoints.empty()) {   //avoid buckets with very few points
                        for (list<Point *>::iterator it = tempPoints.begin(); it != tempPoints.end(); it++) {
                            //if we just started the algorithm
                            if(clusters[i]->point_map.empty()) {
                                if(i == 0) {                            //add point to the first cluster
                                    clusters[i]->add((*it));
                                    assign_counts[i]++;
                                }else {                                 //else check if point is already assigned in another cluster
                                    int mark_assign = 0;                //flag already assigned point
                                    for( int j=0; j < k; j++) {
                                        //if cluster isn't empty check if point is assigned in this cluster
                                        if(!clusters[j]->point_map.empty()) {
                                            if (clusters[j]->point_map.find((*it)->id) != clusters[j]->point_map.end()) {
                                                //if it is assigned mark it and check distances to choose if it changes clusters
                                                mark_assign = 1;
                                                long double dist_i = used_distance(*(*it), *centroids[i]);
                                                long double dist_j = used_distance(*(*it), *centroids[j]);
                                                if (dist_i < dist_j) {
                                                    clusters[i]->add((*it));
                                                    clusters[j]->point_map.erase((*it)->id);
                                                    assign_counts[i]++;
                                                }
                                            }
                                        }
                                    }
                                    if(!mark_assign) {
                                        //unassigned point by all clusters, assign it
                                        clusters[i]->add((*it));
                                        assign_counts[i]++;
                                    }
                                }
                            }else {
                                //check if current point is assigned in cluster i
                                if( clusters[i]->point_map.find((*it)->id) == clusters[i]->point_map.end()) {
                                    int mark_assign = 0;
                                    for ( int j=0; j < k; j++) {
                                        //if it is unassigned check if it is assigned in other clusters
                                        if(j != i) {
                                            if( clusters[j]->point_map.find((*it)->id) != clusters[j]->point_map.end()) {
                                                // if it is assigned in another cluster check its distance between each centroid
                                                // and assign point in cluster with the closest centroid
                                                mark_assign = 1;
                                                long double dist_i = used_distance(*(*it), *centroids[i]);
                                                long double dist_j = used_distance(*(*it), *centroids[j]);
                                                if (dist_i < dist_j) {
                                                    clusters[i]->add((*it));
                                                    clusters[j]->point_map.erase((*it)->id);
                                                    assign_counts[i]++;
                                                }
                                            }
                                        }
                                    }
                                    if(!mark_assign) {
                                        //unassigned point by all clusters, assign it
                                        clusters[i]->add((*it));
                                        assign_counts[i]++;
                                    }
                                }
                            }
                        }
                    }
                }
            }
                break;
            default:
                cerr << "Error: invalid method selected" << endl;
                return 0.0;
        }
        for (int i = 0 ; i < k ; i++) {
            cout << "Cluster " << (i+1) << " has " << assign_counts[i] << " points" << endl;
        }
        delete[] assign_counts;

        // M-step: Recalculate centroids
        cout << "M-step" << endl;
        long double total_change = 0.0;
        for (int i = 0 ; i < clusters.size() ; i++) {
            pair<Centroid *, long double> temp = clusters[i]->recalculate_centroid(data->get_dim());
            total_change += temp.second;
            if (temp.first != NULL){   // if centroid changed
                // update centroid in both locations
                delete centroids[i];
                centroids[i] = temp.first;
                clusters[i]->centroid = temp.first;
            }
        }

        // stop if centroids didn't change much
        cout << "Total change: " << total_change << endl;
        if (total_change < CENTROID_CHANGE_THRESHOLD) {    // should be decreasing
            break;
        }
        if( method == USE_RANGE_LSH || method == USE_HYPERCUBE)
            R *= 2.0;
    }

    // assign every left unassigned pointer to clusters by comparing distances to all centroids
    //like classic Lloyds algorithm
    if(method == USE_RANGE_LSH || method == USE_HYPERCUBE) {
        for( int i=0; i < data->get_size(); i++) {
            long double minDist = FLT_MAX;
            int minCluster = -1;
            for( int j=0; j<k; j++) {
                if(clusters[j]->point_map.find(data->points[i]->id) == clusters[j]->point_map.end()) {
                    if( j == k-1) {
                        for (int n = 0; n < k; n++) {
                            long double dist = used_distance(*data->points[i], *centroids[n]);
                            if (dist < minDist) {
                                minDist = dist;
                                minCluster = n;
                            }
                        }
                        clusters[minCluster]->add(data->points[i]);
                    }
                }
            }
        }
    }

    cout << "DONE" << endl;

    chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    cout << "kMeans iterations: " << count << endl;

    return chrono::duration_cast<chrono::seconds>(end - begin).count();
}

void Cluster::add(Point *p) {
    points.push_back(p);
    //point_map[p->id] = p;
    point_map.insert({p->id, p});
}

void Cluster::clear_cluster() {
    points.clear();
    point_map.clear();
}

pair<Centroid *, long double> Cluster::recalculate_centroid(unsigned int dim) {
    // calculate new centroid
    Centroid *new_centroid = vec_avg(points, dim);
    // calculate change
    long double change = 0.0;
    if (new_centroid != NULL) {
        change = used_distance(*centroid, *new_centroid);
    }
    return make_pair(new_centroid, change);
}

double Clustering::initialize_radius() {
    double minDist=1.0/0.0;
    for(int i=0; i < centroids.size(); i++) {
        for(int j = i+1; j < centroids.size(); j++) {
            double temp_dist = used_distance(*centroids[i],*centroids[j]);
            if(temp_dist < minDist)
                minDist = temp_dist;
        }
    }
    return minDist/2.0;
}

void Clustering::calculate_silhouette() {
    if(method == CLASSIC) {
        int loop_count = 0;
        int loop2_count=0;
        int loop3_count=0;
        for (int i = 0; i < clusters.size(); i++) {
            //for each point in cluster calculate silhouette
            for( int j =0; j < clusters[i]->points.size();j++) {
                long double a = 0.0;
                for( int z = 0; z < clusters[i]->points.size();z++)
                {
                    if( clusters[i]->points[z]->id != clusters[i]->points[j]->id ) {
                        a += used_distance(*clusters[i]->points[j], *clusters[i]->points[z]);
                    }
                    if(loop2_count == 200) {
                        loop2_count=0;
                        break;
                    }
                    loop2_count++;
                }
                a = a / clusters[i]->points.size();
                long double b = 0.0;
                long double minDist = 1.0/0.0;
                int second_cluster;
                for(int n = 0; n < clusters.size();n++) {
                    if (n != i) {
                        long double dist = used_distance(*clusters[i]->points[j],*centroids[n]);
                        if( dist < minDist) {
                            minDist = dist;
                            second_cluster = n;
                        }
                    }
                }
                for ( int l = 0; l < clusters[second_cluster]->points.size(); l++) {
                    b += used_distance(*clusters[i]->points[j],*clusters[second_cluster]->points[l]);
                    if(loop3_count == 200) {
                        loop3_count = 0;
                        break;
                    }
                    loop3_count++;
                }
                b = b / clusters[second_cluster]->points.size();
                long double s = (b-a)/ max(a,b);
                clusters[i]->points_silhouette.push_back(s);
                if(loop_count == 200) {
                    loop_count=0;
                    break;
                }
                loop_count++;
            }
            //for each cluster calculate average silhouette
            long double avg_cl_s = 0.0;
            for ( int j = 0; j < clusters[i]->points_silhouette.size();j++) {
                avg_cl_s += clusters[i]->points_silhouette[j];
            }
            avg_cl_s = avg_cl_s /clusters[i]->points_silhouette.size();
            this->silhouette.push_back(avg_cl_s);
        }
        long double total_silhouette = 0.0;
        for(int i = 0; i < this->silhouette.size(); i++) {
            total_silhouette += this->silhouette[i];
        }
        total_silhouette = total_silhouette / this->silhouette.size();
        this->silhouette.push_back(total_silhouette);
    }
    else {
        int loop_count = 0;
        int loop2_count=0;
        int loop3_count=0;
        for (int i = 0; i < clusters.size(); i++) {
            //for each point in cluster calculate silhouette
            for( auto j = clusters[i]->point_map.begin(); j != clusters[i]->point_map.end(); j++) {
                long double a = 0.0;
                for( auto z = clusters[i]->point_map.begin(); z != clusters[i]->point_map.end();z++)
                {
                    if( z->first != j->first) {
                        a += used_distance(*j->second, *z->second);
                    }
                    if(loop2_count == 200) {
                        loop2_count=0;
                        break;
                    }
                    loop2_count++;
                }
                a = a / clusters[i]->point_map.size();
                long double b = 0.0;
                long double minDist = 1.0/0.0;
                int second_cluster;
                for(int n = 0; n < clusters.size();n++) {
                    if (n != i) {
                        long double dist = used_distance(*j->second,*centroids[n]);
                        if( dist < minDist) {
                            minDist = dist;
                            second_cluster = n;
                        }
                    }
                }
                for ( auto l = clusters[second_cluster]->point_map.begin(); l != clusters[second_cluster]->point_map.end(); l++) {
                    b += used_distance(*j->second,*l->second);
                    if(loop3_count == 200) {
                        loop3_count = 0;
                        break;
                    }
                    loop3_count++;
                }
                b = b / clusters[second_cluster]->point_map.size();
                long double s = (b-a)/ max(a,b);
                clusters[i]->points_silhouette.push_back(s);
                if(loop_count == 200) {
                    loop_count=0;
                    break;
                }
                loop_count++;
            }
            //for each cluster calculate average silhouette
            long double avg_cl_s = 0.0;
            for ( int j = 0; j < clusters[i]->points_silhouette.size();j++) {
                avg_cl_s += clusters[i]->points_silhouette[j];
            }
            avg_cl_s = avg_cl_s /clusters[i]->points_silhouette.size();
            this->silhouette.push_back(avg_cl_s);
        }
        long double total_silhouette = 0.0;
        for(int i = 0; i < this->silhouette.size(); i++) {
            total_silhouette += this->silhouette[i];
        }
        total_silhouette = total_silhouette / this->silhouette.size();
        this->silhouette.push_back(total_silhouette);
    }
}
