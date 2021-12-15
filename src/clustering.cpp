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

void Clustering::initialize_R2clusters(unsigned int k) {
    // Step 1: pick uniformly random point in data:
    srand(time(NULL));  // seed
    int indx = rand() % data->get_size();
    /*Mean_*/Curve *c = new /*Mean_*/Curve(create_copy(*data->curves[indx]->get_curve_coords()));
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
                for (int z = 0 ; z < j ; z++) {    // for already picked centroids
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
        // assign new centroid
        unsigned int selected_index = dist(generator);
        ignore_indexes.insert(selected_index);                 // ignore this in the future
        c = new /*Mean_*/Curve(create_copy(*data->curves[selected_index]->get_curve_coords()));
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
        if (total_change < /*CENTROID_CHANGE_THRESHOLD*/5.0) {    // should be decreasing
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
        // E-step: Assign each point to nearest centroid's cluster
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
                        //data->curves[i]->set_C(meanCurves[j]->get_curve_coords()->size()); //set C(i,j) for optimal traversal
                        long double dist = discreteFrechet_distance(*data->curves[i]->get_curve_coords(), *meanCurves[j]->get_curve_coords());
                        //clusters[j]->set_backtrack(data->curves[i]);  //needed for optimal traversal
                        if (dist < minDist) {
                            minDist = dist;
                            minCluster = j;
                        }
                    }
                    clusters[minCluster]->addCurve(data->curves[i]);   // assign to nearest cluster
                    assign_counts[minCluster]++;
                }
                break;
            case USE_LSH_FRECHET:
            {
                for(int i=0 ; i < k; i++) {     // for every cluster get curves within range
                    list<Curve*> tempCurves;
                    tempCurves = search_functions.LSH_searchR2(R,meanCurves[i]);
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
                                                //(*it)->set_C(meanCurves[i]->get_curve_coords()->size()); //set C(i,j) for optimal traversal
                                                long double dist_i = discreteFrechet_distance(*(*it)->get_curve_coords(), *meanCurves[i]->get_curve_coords());
                                                //clusters[i]->set_backtrack(*it); //needed for optimal traversal

                                                //(*it)->set_C(meanCurves[j]->get_curve_coords()->size()); //set C(i,j) for optimal traversal
                                                long double dist_j = discreteFrechet_distance(*(*it)->get_curve_coords(), *meanCurves[j]->get_curve_coords());
                                                //clusters[j]->set_backtrack(*it); //needed for optimal traversal
                                                if (dist_i < dist_j) {
                                                    clusters[i]->addCurve((*it));
                                                    clusters[j]->curve_map.erase((*it)->get_id());
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

                                                //(*it)->set_C(meanCurves[i]->get_curve_coords()->size()); //set C(i,j) for optimal traversal
                                                long double dist_i = discreteFrechet_distance(*(*it)->get_curve_coords(), *meanCurves[i]->get_curve_coords());
                                                //clusters[i]->set_backtrack(*it);

                                                //(*it)->set_C(meanCurves[j]->get_curve_coords()->size()); //set C(i,j) for optimal traversal
                                                long double dist_j = discreteFrechet_distance(*(*it)->get_curve_coords(), *meanCurves[j]->get_curve_coords());
                                                //clusters[j]->set_backtrack(*it);
                                                if (dist_i < dist_j) {
                                                    clusters[i]->addCurve((*it));
                                                    clusters[j]->curve_map.erase((*it)->get_id());
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
            cout << "Cluster " << (i+1) << " has " << assign_counts[i] << " points" << endl;
        }
        delete[] assign_counts;

        // M-step: Recalculate mean curves
        cout << "M-step" << endl;
        long double total_change = 0.0;
        for (int i = 0 ; i < clusters.size() ; i++) {
            pair<Curve *, long double> temp = clusters[i]->recalculate_meanCurve();
            total_change += temp.second;
            if (temp.first != NULL){   // if mean curve changed
                // update mean curve in both locations
                delete meanCurves[i];
                meanCurves[i] = temp.first;
                clusters[i]->meanCurve = temp.first;
            }
            else { //bad cluster.conf given in such small inputs!
                if(method == CLASSIC) {
                total_change = 0.0;
                break;
                }
            }
        }

        // stop if centroids didn't change much
        cout << "Total change: " << total_change << endl;
        if (total_change < /*CENTROID_CHANGE_THRESHOLD*/5.0) {    // should be decreasing
            break;
        }
        if( method == USE_LSH_FRECHET)
            R *= 2.0;
    }

    // assign every left unassigned pointer to clusters by comparing distances to all mean curves
    //like classic Lloyds algorithm
    if(method == USE_LSH_FRECHET) {
        for( int i=0; i < data->get_size(); i++) {
            long double minDist = DBL_MAX;
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
    if(meanTree != NULL)
        delete meanTree;
//    if(!backtrack_map.empty()) {
//        for(auto j = backtrack_map.begin(); j != backtrack_map.end(); j++){
//            for(int i = 0; i < (sizeof(j->second) / sizeof(j->second[0])) ;++i)
//                delete[] j->second[i];
//            delete[] j->second;
//        }
//    }
}

void Cluster::add(Point *p) {
    points.push_back(p);
    //point_map[p->id] = p;
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

//void Cluster::set_backtrack(Curve *curve/*, int mean_size*/) {
//    //initialize c(i,j) for optimal traversal
//    long double **Mean_c = new long double*[curve->get_curve_coords()->size()];
//    for(int i = 0; i < curve->get_curve_coords()->size(); ++i)
//        Mean_c[i] = new long double[mean_size];

    //set it to backtracking map with its curve
    //backtrack_map[curve->get_id()] = /*Mean_c*/meanCurve->C;
//    backtrack_map.insert({curve->get_id(), /*meanCurve->C*/curve->C});
//}

vector<pair< pair<double,double>, pair<double,double> >> *Cluster::OptimalTraversal(/*Mean_*/Curve *c, /*Mean_*/Curve *centroid) {
    //traversal is an empty list of pairs
    vector<pair< pair<double,double>, pair<double,double> >> *traversal = new vector<pair< pair<double,double>, pair<double,double> >>();
    //cout << "start traversal!" << endl;

    //initialization for traversal
    int P_i, Q_i;
    P_i = c->get_curve_coords()->size() -1;
    Q_i = centroid->get_curve_coords()->size() -1;
    pair <double, double> p, q;
    pair<pair <double, double>, pair <double, double> > p_q;
    p = c->get_curve_coords()->at(P_i);
    q = centroid->get_curve_coords()->at(Q_i);
    p_q.first = p;
    p_q.second = q;

    //first value set
    traversal->push_back(p_q);

    //loop until one of those curves end their traversal
    while(P_i && Q_i) {
        //minIdx = index of min([Pi − 1, Qi], C[Pi, Qi − 1], C[Pi − 1, Qi − 1])
        int minInd;
        if(c->C[P_i - 1][Q_i] <= c->C[P_i][Q_i - 1]) {
            if(c->C[P_i - 1][Q_i] < c->C[P_i - 1][Q_i - 1])
                minInd = 0;
            else
                minInd = 2;
        } else {
            if(c->C[P_i][Q_i - 1] < c->C[P_i - 1][Q_i - 1])
                minInd = 1;
            else
                minInd = 2;
        }
        if(!minInd) {
            --P_i;
            p = c->get_curve_coords()->at(P_i);
            q = centroid->get_curve_coords()->at(Q_i);
            p_q.first = p;
            p_q.second = q;
            traversal->push_back(p_q);
            //cout << "(" << traversal->at(0).first.first << " , " << traversal->at(0).first.second << "), (" << traversal->at(0).second.first << " , " << traversal->at(0).second.second <<")" <<endl;
        } else if(minInd == 1) {
            --Q_i;
            p = c->get_curve_coords()->at(P_i);
            q = centroid->get_curve_coords()->at(Q_i);
            p_q.first = p;
            p_q.second = q;
            traversal->push_back(p_q);
        } else {
            --P_i;
            --Q_i;
            p = c->get_curve_coords()->at(P_i);
            q = centroid->get_curve_coords()->at(Q_i);
            p_q.first = p;
            p_q.second = q;
            traversal->push_back(p_q);
        }
    }

    //traversal is reversed. Mean curve will not be!
    //cout << "not killed in traversal!" << endl;
    return traversal;
}

int Cluster::initialize_tree() {
    srand(time(NULL));  // seed
    cout<< "initialize tree" << endl;

    //compute height of tree ceiling of lg(size of cluster)
    unsigned int size = curves.size();  //leaves size
   // if(size % 2 != 0)       //sum of leaves can't be an odd number in Binary Search Tree
   //     size = size+1;
    unsigned int h = ceil(log2(size));
    cout << "height: " << h << endl;

    //find size of tree
    int tree_size = pow(2,(h+1)) -1;
    cout << "power: " << pow(2,(h+1)) << ", tree_size:" << tree_size << endl;
    int leaves_sum = pow(2,h);          //how many leaves does the bst have
    cout << "leaves sum: " << leaves_sum << endl;
    //cout << "tree size= " << tree_size << "curves size= " << size << endl;

    //initialize tree
    //cout << "vector tree initialization" << endl;
    meanTree = new vector<Curve *>(tree_size);
//    for(int i=0; i < meanTree->size(); i++)
//        meanTree->push_back(NULL);

    //initialize tree leaves randomly
    vector<int> seen;
    int last_leaf = leaves_sum + size - 1; //last leaf placement
    for(int i = (leaves_sum-1); i < last_leaf ; i++){

        //different curve inserted in tree every time
        int indx = rand() % size;
        for(int j=0; j < seen.size(); j++) {
            if(seen.at(j) == indx) {
                indx = rand() % size;
                j=0;
            }
        }
        seen.push_back(indx);
        meanTree->at(i) = curves[indx];

        //Curve *leaf = vec_avg_curve(*OptimalTraversal(curves[indx],meanCurve));
        //leaf->Remove_duplicates();
        //leaf->R2_Filtering(MAX_LENGTH);
    }
    return last_leaf;
}

//implementation of BS tree is structured on array of Curves
int Cluster::PostOrderTraversal(int node, int leaves_places) {
    //if node is leaf
    if(node < leaves_places && node >= (leaves_places - curves.size())) {
        //cout << "leaves: " << node << endl;
        //meanTree->at(node)->print();
        return node;
    }else {
        int leftNode = (2*node)+1;
        //corner case
        if(leftNode >= leaves_places) { //if left node is NULL then right will also be NULL so just get left of parent
            meanTree->at(node) = meanTree->at(node-1);
            return node;
        }
        leftNode = PostOrderTraversal(leftNode, leaves_places);
        int rightNode = leftNode+1;
        if(rightNode < leaves_places) { //if right node isn't empty
            rightNode = PostOrderTraversal(rightNode, leaves_places);
            //cout << "segm is coming" << endl;

            //corner case we explained in left node check! If pointers have the same address just give curve to their parent node
            if(meanTree->at(leftNode) == meanTree->at(rightNode)) {
                meanTree->at(node) = meanTree->at(leftNode);
                //cout <<"got in?" << endl;
                return node;
            }
            //cout <<"segm came!" << endl;

            //cout << "left - right: " << leftNode << " - " << rightNode << endl;
            //cout << "node: " << node << endl;
//            cout << "right node size: " << meanTree->at(rightNode)->get_curve_coords()->size() << endl;

            meanTree->at(leftNode)->set_C((meanTree->at(rightNode)->get_curve_coords()->size()));
//            cout << "after set" << endl;


            double distance = discreteFrechet_distance(*meanTree->at(leftNode)->get_curve_coords(),*meanTree->at(rightNode)->get_curve_coords(),meanTree->at(leftNode));
            Curve *nodeCurve = vec_avg_curve(*OptimalTraversal(meanTree->at(leftNode),meanTree->at(rightNode)));
            nodeCurve->Remove_duplicates();                    //remove duplicates
            nodeCurve->R2_Filtering(MAX_LENGTH);        //filter mean curve
            meanTree->at(node) = nodeCurve;                 //put new mean curve in current tree node
            return node;
        }
        else {
            //Mean Curve only for the existing node so just go back to parent with this value
            //cout << "else right null: " << node << endl;
            meanTree->at(node) = meanTree->at(leftNode);
            return node;
        }
    }
}

pair</*Mean_*/Curve *, long double> Cluster::recalculate_meanCurve() {
    //initialize binary tree and do PostTraversal order to get new mean curve(centroid)
   // cout << "initialize tree..." << endl;

    if (curves.empty()){          // empty cluster
        cerr << "Warning: cluster with no curves assigned!" << endl;
        pair <Curve *, long double> null_pair;
        null_pair.first = NULL;
        null_pair.second = 0.0;
        return null_pair;
    }

    int place = initialize_tree();
  //  cout << "done!" << endl;
    Curve *new_centroid = meanTree->at(PostOrderTraversal(0,place)); //0 because our root is the new centroid

    // calculate change
    long double change = 0.0;
    if (new_centroid != NULL) {
        change = discreteFrechet_distance(*meanCurve->get_curve_coords(), *new_centroid->get_curve_coords());
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
        int loop_count = 0;
        int loop2_count=0;
        int loop3_count=0;
        for (int i = 0; i < clusters.size(); i++) {
            //for each curve in cluster calculate silhouette
            for( int j =0; j < clusters[i]->curves.size();j++) {
                long double a = 0.0;
                for( int z = 0; z < clusters[i]->curves.size();z++)
                {
                    if( clusters[i]->curves[z]->get_id() != clusters[i]->curves[j]->get_id() ) {
                        a += discreteFrechet_distance(*clusters[i]->curves[j]->get_curve_coords(), *clusters[i]->curves[z]->get_curve_coords());
                    }
                    if(loop2_count == 200) {
                        loop2_count=0;
                        break;
                    }
                    loop2_count++;
                }
                a = a / clusters[i]->curves.size();
                long double b = 0.0;
                long double minDist = 1.0/0.0;
                int second_cluster;
                for(int n = 0; n < clusters.size();n++) {
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
                    if(loop3_count == 200) {
                        loop3_count = 0;
                        break;
                    }
                    loop3_count++;
                }
                b = b / clusters[second_cluster]->curves.size();
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
            for( auto j = clusters[i]->curve_map.begin(); j != clusters[i]->curve_map.end(); j++) {
                long double a = 0.0;
                for( auto z = clusters[i]->curve_map.begin(); z != clusters[i]->curve_map.end();z++)
                {
                    if( z->first != j->first) {
                        a += discreteFrechet_distance(*j->second->get_curve_coords(), *z->second->get_curve_coords());
                    }
                    if(loop2_count == 200) {
                        loop2_count=0;
                        break;
                    }
                    loop2_count++;
                }
                a = a / clusters[i]->curve_map.size();
                long double b = 0.0;
                long double minDist = 1.0/0.0;
                int second_cluster;
                for(int n = 0; n < clusters.size();n++) {
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
                    if(loop3_count == 200) {
                        loop3_count = 0;
                        break;
                    }
                    loop3_count++;
                }
                b = b / clusters[second_cluster]->curve_map.size();
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


