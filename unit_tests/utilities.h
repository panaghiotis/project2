#include <iostream>
#include <vector>
#include <utility>
#include <algorithm>
#include <cmath>


using namespace std;

//euclidian(L2) distance for Discrete Frechet Distance implementation
long double euclidean2d_distance(pair<int,double> p, pair<int,double> q) {
    long double sum = 0.0;
    long double x = p.first - q.first;
    long double y = p.second - q.second;
    sum = pow(x,2) + pow(y,2);
    return sqrt(sum);
}

long double discreteFrechet_distance(vector<pair<double,double>> p, vector<pair<double,double>> q) {
    long double result = 0.0;

    //c(i,j) initialization
    long double **c = new long double*[p.size()];
    for(int i = 0; i < p.size(); ++i)
        c[i] = new long double[q.size()];

    // fill the first value with the distance between the first two points in p and q
    c[0][0] = euclidean2d_distance(p[0], q[0]);

    // load the first column and first row with distances (memorize)
    for (int i=1; i < p.size(); i++)
        c[i][0] = max(c[i-1][0], euclidean2d_distance(p[i], q[0]));
    for (int j=1; j < q.size(); j++)
        c[0][j] = max(c[0][j-1], euclidean2d_distance(p[0], q[j]));

    // load every other column and row with distances
    for (int i=1; i < p.size(); i++)
        for (int j=1; j < q.size(); j++)
            c[i][j] = max(min(min(c[i-1][j], c[i][j-1]), c[i-1][j-1]), euclidean2d_distance(p[i], q[j]));

    //get discrete frechet distance
    result= c[p.size()-1][q.size()-1];

    //free memory
    for(int i = 0; i < p.size();++i)
        delete[] c[i];
    delete[] c;

    //return distance
    return result;
}

//implementation of BS tree is structured on array of Curves
int PostOrderTraversal(int node, int leaves_places, vector<int> &meanTree) {
    //if node is leaf
    if(node < leaves_places && node >= (leaves_places - 2)) {
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
            int nodeCurve = 0;
            meanTree[node] = nodeCurve;                 //put new mean curve in current tree node
            return node;
        }
        else {
            //Mean Curve only for the existing node so just go back to parent with this value
            meanTree.at(node) = meanTree.at(leftNode);
            return node;
        }
    }
}

//return Discrete Frechet C(i,j) array for clustering
double** get_C(vector<pair<double,double>> p, vector<pair<double,double>> q){

    //c(i,j) initialization
    double **c = new double*[p.size()];
    for(int i = 0; i < p.size(); ++i)
        c[i] = new double[q.size()];

    // fill the first value with the distance between the first two points in p and q
    c[0][0] = (double)euclidean2d_distance(p[0], q[0]);

    // load the first column and first row with distances (memorize)
    for (int i=1; i < p.size(); i++)
        c[i][0] = max(c[i-1][0], (double)euclidean2d_distance(p[i], q[0]));
    for (int j=1; j < q.size(); j++)
        c[0][j] = max(c[0][j-1], (double)euclidean2d_distance(p[0], q[j]));

    // load every other column and row with distances
    for (int i=1; i < p.size(); i++)
        for (int j=1; j < q.size(); j++)
            c[i][j] = max(min(min(c[i-1][j], c[i][j-1]), c[i-1][j-1]), (double)euclidean2d_distance(p[i], q[j]));

    return c;
}

//filtering mean curve for clustering
vector<pair<double,double>> filtering(vector<pair<double,double>> R2_coords, double epsilon) {
    //if mean curve is smaller than max mean curve length don't filter
    if(R2_coords.size() <=10)
        return R2_coords;

    // filter coordinates
    for(int i=0; i < R2_coords.size(); i++) {
        if(R2_coords.size() <= 10)       //max mean curve length
            break;
        if( i < R2_coords.size() - 2) {
            // |a-b| <= ε and |b-c| <= ε ,remove b
            if(abs(R2_coords.at(i).second - R2_coords.at(i+1).second) <= epsilon && abs(R2_coords.at(i+1).second - R2_coords.at(i+2).second) <= epsilon) {
                R2_coords.erase(R2_coords.begin()+(i+1));
                //check i again with i+2 next to it this time
                i = i - 1;
            }
        }
    }
    //get filtered coordinates
    return R2_coords;
}