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

extern double epsilon;                    //for filtering mean curves
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

//Same initialization of points centroids but now for mean curves-centroids in R2
void Clustering::initialize_R2clusters(unsigned int k) {
    // Step 1: pick uniformly random curve in data:
    srand(time(NULL));  // seed
    int indx = rand() % data->get_size();
    Curve *c = new Curve(create_copy(*data->curves[indx]->get_curve_coords()));
    meanCurves.push_back(c);
    clusters.push_back(new Cluster(c));

    unordered_set<unsigned int> ignore_indexes;
    ignore_indexes.insert(indx);

    // Step 2-k: choose the other centers by chance depending on their distance from already selected centers
    for (int j = 1 ; j < k ; j++) {
        vector<long double> D;
        D.reserve(data->curves.size());
        // calculate D(i) for each unassigned point i
        for (int i = 0 ; i < data->curves.size() ; i++) {
            if (ignore_indexes.find(i) != ignore_indexes.end()) {   // if index already picked
                D.push_back(0.0);                                   // assign zero probabilty
            } else {
                long double minDist = FLT_MAX;
                for (int z = 0 ; z < j ; z++) {    // for already picked mean curves
                    long double dist = discreteFrechet_distance(*data->curves[i]->get_curve_coords(), *meanCurves[z]->get_curve_coords());
                    if (dist < minDist) {
                        minDist = dist;
                    }
                }
                D.push_back(pow(minDist, 2));    // no need to normalize, dist does it for us, only put D(i)^2
            }
        }
        // create discrete dist based on D(i)
        discrete_distribution<size_t> dist{D.begin(), D.end()};
        // assign new mean curve
        unsigned int selected_index = dist(generator);
        ignore_indexes.insert(selected_index);                 // ignore this in the future
        c = new Curve(create_copy(*data->curves[selected_index]->get_curve_coords()));
        // add to clusters
        meanCurves.push_back(c);
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

// K Means for mean curves. Same implementation only for Classic and LSH Frechet methods
double Clustering::perform_R2kMeans(unsigned int k) {
    chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    // initialize clusters and centroids
    cout << "Initializing clusters..." << endl;
    this->initialize_R2clusters(k);
    cout << "Done!" << endl;

    double R;
    if(method == USE_LSH_FRECHET)
        R = initialize_radiusR2();

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
        // E-step: Assign each curve to the nearest mean curve's cluster
        unsigned int *assign_counts = new unsigned int[k];
        for (int i = 0 ; i < k ; i++) {
            assign_counts[i] = 0;
        }
        switch (method){
            case CLASSIC:
                // exact NNs
                for (int i = 0 ; i < data->get_size() ; i++){    // for each curve
                    long double minDist = FLT_MAX;
                    int minCluster = -1;
                    for (int j = 0 ; j < k ; j++) {                       // calculate all distances and keep the minimum
                        long double dist = discreteFrechet_distance(*data->curves[i]->get_curve_coords(), *meanCurves[j]->get_curve_coords());
                        if (dist < minDist) {
                            minDist = dist;
                            minCluster = j;
                        }
                    }
                    clusters[minCluster]->addCurve(data->curves[i]);   // assign to the nearest cluster
                    assign_counts[minCluster]++;
                }
                break;
            case USE_LSH_FRECHET:
            {
                for(int i=0 ; i < k; i++) {     // for every cluster get curves within range
                    list<Curve*> tempCurves;
                    tempCurves = search_functions.LSH_searchR2(R,meanCurves[i],MAX_LENGTH);
                    if(!tempCurves.empty()) {   //avoid buckets with very few points
                        for (list<Curve *>::iterator it = tempCurves.begin(); it != tempCurves.end(); it++) {
                            //if we just started the algorithm
                            if(clusters[i]->curve_map.empty()) {
                                if(i == 0) {                            //add curve to the first cluster
                                    clusters[i]->addCurve((*it));
                                    assign_counts[i]++;
                                }else {                                 //else check if curve is already assigned in another cluster
                                    int mark_assign = 0;                //flag already assigned curve
                                    for( int j=0; j < k; j++) {
                                        //if cluster isn't empty check if curve is assigned in this cluster
                                        if(!clusters[j]->curve_map.empty()) {
                                            if (clusters[j]->curve_map.find((*it)->get_id()) != clusters[j]->curve_map.end()) {
                                                //if it is assigned mark it and check distances to choose if it changes clusters
                                                mark_assign = 1;
                                                long double dist_i = discreteFrechet_distance(*(*it)->get_curve_coords(), *meanCurves[i]->get_curve_coords());
                                                long double dist_j = discreteFrechet_distance(*(*it)->get_curve_coords(), *meanCurves[j]->get_curve_coords());
                                                if (dist_i < dist_j) {
                                                    clusters[i]->addCurve((*it));
                                                    clusters[j]->curve_map.erase((*it)->get_id());
                                                    clusters[j]->curves.erase(remove(clusters[j]->curves.begin(), clusters[j]->curves.end(), (*it)), clusters[j]->curves.end());
                                                    assign_counts[i]++;
                                                }
                                            }
                                        }
                                    }
                                    if(!mark_assign) {
                                        //unassigned curve by all clusters, assign it
                                        clusters[i]->addCurve((*it));
                                        assign_counts[i]++;
                                    }
                                }
                            }else {
                                //check if current curve is assigned in cluster i
                                if( clusters[i]->curve_map.find((*it)->get_id()) == clusters[i]->curve_map.end()) {
                                    int mark_assign = 0;
                                    for ( int j=0; j < k; j++) {
                                        //if it is unassigned check if it is assigned in other clusters
                                        if(j != i) {
                                            if( clusters[j]->curve_map.find((*it)->get_id()) != clusters[j]->curve_map.end()) {
                                                // if it is assigned in another cluster check its distance between each mean curve
                                                // and assign curve in cluster with the closest mean curve
                                                mark_assign = 1;
                                                long double dist_i = discreteFrechet_distance(*(*it)->get_curve_coords(), *meanCurves[i]->get_curve_coords());
                                                long double dist_j = discreteFrechet_distance(*(*it)->get_curve_coords(), *meanCurves[j]->get_curve_coords());
                                                if (dist_i < dist_j) {
                                                    clusters[i]->addCurve((*it));
                                                    clusters[j]->curve_map.erase((*it)->get_id());
                                                    clusters[j]->curves.erase(remove(clusters[j]->curves.begin(), clusters[j]->curves.end(), (*it)), clusters[j]->curves.end());
                                                    assign_counts[i]++;
                                                }
                                            }
                                        }
                                    }
                                    if(!mark_assign) {
                                        //unassigned curve by all clusters, assign it
                                        clusters[i]->addCurve((*it));
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
            if(method == USE_LSH_FRECHET) {
                clusters[i]->curve_map_to_vec();
                assign_counts[i] = clusters[i]->curve_map.size();
            }
            cout << "Cluster " << (i+1) << " has " << assign_counts[i] << " curves" << endl;
        }
        delete[] assign_counts;

        // M-step: Recalculate mean curves
        cout << "M-step" << endl;
        long double total_change = 0.0;
        for (int i = 0 ; i < clusters.size() ; i++) {
            pair<Curve *, long double> temp = clusters[i]->recalculate_meanCurve();
            total_change += temp.second;
            //clusters[i]->clear_tree();             //Tree must be cleared after every loop
            if (temp.first != NULL){   // if mean curve changed
                // update mean curve in both locations
                delete meanCurves[i];
                meanCurves[i] = temp.first;
                clusters[i]->meanCurve = temp.first;
            }
            else {
                total_change = 0.0;
                cout <<"Found an empty cluster. End Kmeans." << endl;
                break;
            }
        }

        // stop if mean curves didn't change much
        cout << "Total change: " << total_change << endl;
        if (total_change < CENTROID_CHANGE_THRESHOLD) {
            break;
        }
        else if ( count > 1 && total_change >= 5.0 ){
            int break_flag_3 = 0;
            int break_flag_2 = 0;
            for(int i=0; i < clusters.size(); i++) {       // in Frechet mean curves tend to be a lot different after every loop but clusters usually
                if (clusters[i]->curves.size() >= 10)      // have almost the same number of curves. So no need to go through a lot of loops.
                    break_flag_3++;                        // If we find at least 3 clusters with 10 or more curves inside them, stop the loop,
                if (clusters[i]->curves.size() >= 20)      // else if we find at least 2 clusters with 20 or more curves, stop the loop.
                    break_flag_2++;
            }
            if(break_flag_3 > 2 || break_flag_2 > 1){
                cout << "Found good clusters!" << endl;
                break;
            }
        }
        if (count == 10) { //end loop
            cout<<"Taking too much time. End it" << endl;
            break;
        }

        if( method == USE_LSH_FRECHET)
            R *= 2.0;
    }
    // assign every left unassigned pointer to clusters by comparing distances to all mean curves
    //like classic Lloyd's algorithm
    if(method == USE_LSH_FRECHET) {
        for( int i=0; i < data->get_size(); i++) {
            long double minDist = FLT_MAX;
            int minCluster = -1;
            for( int j=0; j<k; j++) {
                if(clusters[j]->curve_map.find(data->curves[i]->get_id()) == clusters[j]->curve_map.end()) {
                    if( j == k-1) {
                        for (int n = 0; n < k; n++) {
                            long double dist = discreteFrechet_distance(*data->curves[i]->get_curve_coords(), *meanCurves[n]->get_curve_coords());
                            if (dist < minDist) {
                                minDist = dist;
                                minCluster = n;
                            }
                        }
                        clusters[minCluster]->addCurve(data->curves[i]);
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

Cluster::~Cluster() {
    if(centroid != NULL)
        delete centroid;
    if(meanCurve != NULL)
        delete meanCurve;
}

void Cluster::add(Point *p) {
    points.push_back(p);
    point_map.insert({p->id, p});
}

void Cluster::addCurve(Curve *c) {
    curves.push_back(c);
    curve_map.insert({c->get_id(),c});
}

void Cluster::clear_cluster() {
    if(!points.empty())
        points.clear();
    if(!curves.empty())
        curves.clear();
    if(!point_map.empty())
        point_map.clear();
    if(!curve_map.empty())
        curve_map.clear();
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

void Cluster::curve_map_to_vec() {
    this->curves.clear();
    for( auto i = curve_map.begin(); i != curve_map.end();i++) {
        curves.push_back( i->second );
    }
}

pair</*Mean_*/Curve *, long double> Cluster::recalculate_meanCurve() {
    //initialize binary tree and do PostTraversal order in order to get new mean curve(centroid)
    if (curves.empty()){          // empty cluster
        cerr << "Warning: cluster with no curves assigned!" << endl;
        pair <Curve *, long double> null_pair;
        null_pair.first = NULL;
        null_pair.second = 0.0;
        return null_pair;
    }

    srand(time(NULL));  // seed

    //compute height of tree ceiling of lg(size of cluster)
    unsigned int size;
    size = curves.size();  //leaves size
    unsigned int h = ceil(log2(size));

    //find size of tree
    int tree_size = pow(2,(h+1)) -1;
    int leaves_sum = pow(2,h);          //how many leaves does the BST have

    //initialize tree
    vector<vector<pair<double,double>>> meanTree(tree_size);
    //initialize tree leaves randomly
    vector<int> seen;
    int last_leaf = leaves_sum + size - 1; //last leaf placement

    for(int i=0; i < (leaves_sum-1); i++)
        meanTree.at(i).push_back({0,0});
    for(int i=last_leaf; i < meanTree.size(); i++)
        meanTree.at(i).push_back({0,0});

    for(int i = (leaves_sum-1); i < last_leaf ; i++){
        //different curve inserted in tree's leaves every time
        int indx = rand() % size;
        for(int j=0; j < seen.size(); j++) {
            if(seen.at(j) == indx) {
                indx = rand() % size;
                j=0;
            }
        }
        seen.push_back(indx);
        for(int j=0; j < curves[indx]->get_curve_coords()->size(); j++)
            meanTree.at(i).push_back(curves[indx]->get_curve_coords()->at(j));
    }

    cout<< "Post Order Traversal tree method for each cluster..." << endl;
    vector<pair<double, double>> coords = meanTree.at(PostOrderTraversal(0,last_leaf,meanTree)); //get root mean curve
    cout << "Done!" << endl;
    vector<pair<double, double>> *curve_coords = new vector<pair<double, double>>();
    for(int i=0; i < coords.size(); i++)
        curve_coords->push_back(coords.at(i));
    Curve *new_centroid = new Curve(curve_coords);

    //no need for BST anymore
    meanTree.clear();

    // calculate change
    long double change = 0.0;
    if (new_centroid != NULL) {
        change = discreteFrechet_distance(*meanCurve->get_curve_coords(), *new_centroid->get_curve_coords());
    }
    return make_pair(new_centroid, change);
}

vector<pair< pair<double,double>, pair<double,double> >> Cluster::OptimalTraversal(vector<pair<double,double>> c, vector<pair<double,double>> centroid) {
    //traversal is an empty list of pairs
    vector<pair< pair<double,double>, pair<double,double> >> traversal;
    double **C = get_C(c,centroid);

    //initialization for traversal
    int P_i, Q_i;
    P_i = c.size() -1;
    Q_i = centroid.size() -1;
    pair <double, double> p, q;
    pair<pair <double, double>, pair <double, double> > p_q;
    p = c.at(P_i);
    q = centroid.at(Q_i);
    p_q.first = p;
    p_q.second = q;

    //first value set
    traversal.push_back(p_q);

    //loop until one of those curves end their traversal
    while(P_i && Q_i) {
        if(traversal.size() >= 1000)
            break;
        //minIdx = index of min([Pi − 1, Qi], C[Pi, Qi − 1], C[Pi − 1, Qi − 1])
        int minInd;
        if(C[P_i - 1][Q_i] <= C[P_i][Q_i - 1]) {
            if(C[P_i - 1][Q_i] < C[P_i - 1][Q_i - 1])
                minInd = 0;
            else
                minInd = 2;
        } else {
            if(C[P_i][Q_i - 1] < C[P_i - 1][Q_i - 1])
                minInd = 1;
            else
                minInd = 2;
        }
        if(!minInd) {
            --P_i;
            p = c.at(P_i);
            q = centroid.at(Q_i);
            p_q.first = p;
            p_q.second = q;
            traversal.push_back(p_q);
        } else if(minInd == 1) {
            --Q_i;
            p = c.at(P_i);
            q = centroid.at(Q_i);
            p_q.first = p;
            p_q.second = q;
            traversal.push_back(p_q);
        } else {
            --P_i;
            --Q_i;
            p = c.at(P_i);
            q = centroid.at(Q_i);
            p_q.first = p;
            p_q.second = q;
            traversal.push_back(p_q);
        }
    }

    for(int i = 0; i < c.size();++i)
        delete[] C[i];
    delete[] C;

    //traversal is reversed. Mean curve will not be!
    return traversal;
}

//implementation of BS tree is structured on array of Curves
int Cluster::PostOrderTraversal(int node, int leaves_places, vector<vector<pair<double,double>>> &meanTree) {
    //if node is leaf
    if(node < leaves_places && node >= (leaves_places - curves.size())) {
        return node;
    }else {
        int leftNode = (2*node)+1;
        //corner case
        if(leftNode >= leaves_places) { //if left node is NULL then right will also be NULL so just get left of parent
            meanTree.at(node) = meanTree.at(node-1);
            return node;
        }
        leftNode = PostOrderTraversal(leftNode, leaves_places, meanTree);
        int rightNode = leftNode+1;
        if(rightNode < leaves_places) { //if right node isn't empty
            rightNode = PostOrderTraversal(rightNode, leaves_places, meanTree);

            //corner case we explained in left node check! If pointers have the same address just give curve to their parent node
            if(meanTree.at(leftNode) == meanTree.at(rightNode)) {
                meanTree.at(node) = meanTree.at(leftNode);
                return node;
            }
            vector<pair<double,double>> nodeCurve = vec_avg_curve(OptimalTraversal(meanTree.at(leftNode),meanTree.at(rightNode)), epsilon);
            meanTree.at(node) = nodeCurve;                 //put new mean curve in current tree node
            return node;
        }
        else {
            //Mean Curve only for the existing node so just go back to parent with this value
            meanTree.at(node) = meanTree.at(leftNode);
            return node;
        }
    }
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

double Clustering::initialize_radiusR2() {
    double minDist=1.0/0.0;
    for(int i=0; i < meanCurves.size(); i++) {
        for(int j = i+1; j < meanCurves.size(); j++) {
            double temp_dist = discreteFrechet_distance(*meanCurves[i]->get_curve_coords(),*meanCurves[j]->get_curve_coords());
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

void Clustering::calculate_silhouetteR2() {
    if(method == CLASSIC) {
        int empty_flag=0;
        int loop_count = 0;
        int loop2_count=0;
        int loop3_count=0;
        for (int i = 0; i < clusters.size(); i++) {
            while(clusters[i]->curves.empty()) { //empty cluster
                i++;
                if(i == clusters.size()) {
                    empty_flag = 1;
                    break;
                }
            }
            if(empty_flag)
                break;
            //for each curve in cluster calculate silhouette
            for( int j =0; j < clusters[i]->curves.size();j++) {
                long double a = 0.0;
                for( int z = 0; z < clusters[i]->curves.size();z++)
                {
                    if( clusters[i]->curves[z]->get_id() != clusters[i]->curves[j]->get_id() ) {
                        a += discreteFrechet_distance(*clusters[i]->curves[j]->get_curve_coords(), *clusters[i]->curves[z]->get_curve_coords());
                    }
                    if(loop2_count == 20) {
                        loop2_count=0;
                        break;
                    }
                    loop2_count++;
                }
                a = a / clusters[i]->curves.size();
                long double b = 0.0;
                long double minDist = 1.0/0.0;
                int second_cluster;
                int second_flag = 0;
                for(int n = 0; n < clusters.size();n++) {
                    while(clusters[n]->curves.empty()) { //empty cluster
                        n++;
                        if(n == clusters.size()) {
                            second_flag = 1;
                            break;
                        }
                    }
                    if(second_flag)
                        break;
                    if (n != i) {
                        long double dist = discreteFrechet_distance(*clusters[i]->curves[j]->get_curve_coords(),*meanCurves[n]->get_curve_coords());
                        if( dist < minDist) {
                            minDist = dist;
                            second_cluster = n;
                        }
                    }
                }
                for ( int l = 0; l < clusters[second_cluster]->curves.size(); l++) {
                    b += discreteFrechet_distance(*clusters[i]->curves[j]->get_curve_coords(),*clusters[second_cluster]->curves[l]->get_curve_coords());
                    if(loop3_count == 20) {
                        loop3_count = 0;
                        break;
                    }
                    loop3_count++;
                }
                b = b / clusters[second_cluster]->curves.size();
                long double s = (b-a)/ max(a,b);
                clusters[i]->points_silhouette.push_back(s);
                if(loop_count == 20) {
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
        int empty_flag=0;
        int loop_count = 0;
        int loop2_count=0;
        int loop3_count=0;
        for (int i = 0; i < clusters.size(); i++) {
            while(clusters[i]->curve_map.empty()) { //empty cluster
                i++;
                if(i == clusters.size()) {
                    empty_flag = 1;
                    break;
                }
            }
            if(empty_flag)
                break;
            //for each point in cluster calculate silhouette
            for( auto j = clusters[i]->curve_map.begin(); j != clusters[i]->curve_map.end(); j++) {
                long double a = 0.0;
                for( auto z = clusters[i]->curve_map.begin(); z != clusters[i]->curve_map.end();z++)
                {
                    if( z->first != j->first) {
                        a += discreteFrechet_distance(*j->second->get_curve_coords(), *z->second->get_curve_coords());
                    }
                    if(loop2_count == 20) {
                        loop2_count=0;
                        break;
                    }
                    loop2_count++;
                }
                a = a / clusters[i]->curve_map.size();
                long double b = 0.0;
                long double minDist = 1.0/0.0;
                int second_cluster;
                int second_flag = 0;
                for(int n = 0; n < clusters.size();n++) {
                    while(clusters[n]->curves.empty()) { //empty cluster
                        n++;
                        if(n == clusters.size()) {
                            second_flag = 1;
                            break;
                        }
                    }
                    if(second_flag)
                        break;
                    if (n != i) {
                        long double dist = discreteFrechet_distance(*j->second->get_curve_coords(),*meanCurves[n]->get_curve_coords());
                        if( dist < minDist) {
                            minDist = dist;
                            second_cluster = n;
                        }
                    }
                }
                for ( auto l = clusters[second_cluster]->curve_map.begin(); l != clusters[second_cluster]->curve_map.end(); l++) {
                    b += discreteFrechet_distance(*j->second->get_curve_coords(),*l->second->get_curve_coords());
                    if(loop3_count == 20) {
                        loop3_count = 0;
                        break;
                    }
                    loop3_count++;
                }
                b = b / clusters[second_cluster]->curve_map.size();
                long double s = (b-a)/ max(a,b);
                clusters[i]->points_silhouette.push_back(s);
                if(loop_count == 20) {
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