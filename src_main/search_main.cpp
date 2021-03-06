#include <iostream>
#include <fstream>
#include <cfloat>
#include "../header/neighbour_search.h"
#include "../header/InputParser.h"

using namespace std;

// globals
unsigned int k = 0;           // number of hi hash functions, also dimension d' for hypercube
unsigned int L = 0;           // number of g()s , hash tables and grids
unsigned int N = 1;           // number of neighbours to search
unsigned int M = 0;           // number of points to be checked
unsigned int probes = 0;      // number of vertices to check
double delta = 0.0;           //delta for Frechet
double epsilon = 1;        //epsilon for curve filtering

#define LSH 1
#define HYPERCUBE 2
#define FRECHET 3
#define DISCRETE 4
#define CONTINUOUS 5

#define DEFAULT_k 4           // number of hash functions per hash table
#define DEFAULT_L 5           // number of hash tables
#define DEFAULT_CONTI_L 1     // L must be 1 for continuous Frechet
#define DEFAULT_CUBE_k 14     // number of dimension d' (also k: number of h)
#define DEFAULT_M 10          // number of points to be checked
#define DEFAULT_PR 2          // number of vertices to be checked
#define DEFAULT_Delta 0.1     // delta number

#define DIVIDE_DATASET_FOR_HASHTABLE_SIZE 8


int main(int argc, char **argv) {
    // read input
    InputParser input(argc, argv);
    const string &input_file = input.getCmdOption("-i");
    const string &query_file = input.getCmdOption("-q");
    const string &output_file = input.getCmdOption("-o");
    const string &algorithm = input.getCmdOption("-algorithm");
    const string &metric = input.getCmdOption("-metric");

    int method, frechet_metric;
    if (algorithm == "LSH"){
        method = LSH;
    }else if(algorithm == "Hypercube"){
        method = HYPERCUBE;
    }else if (algorithm == "Frechet"){
        method = FRECHET;
        if(metric == "discrete")
            frechet_metric = DISCRETE;
        else if(metric == "continuous")
            frechet_metric = CONTINUOUS;
        else {
            string answer;
            while (1) {
                cerr << "No correct metric specified for Frechet" << endl;
                cout << "Please choose between discrete or continuous: ";
                cin >> answer;
                if(answer == "discrete"){
                    frechet_metric = DISCRETE;
                    break;
                } else if (answer == "continuous") {
                    frechet_metric = CONTINUOUS;
                    break;
                }
            }
        }
    }else {
        string answer;
        while (1) {
            cerr << "No correct method specified for search" << endl;
            cout << "Please choose from options LSH , Hypercube , Frechet: ";
            cin >> answer;
            if(answer == "LSH") {
                method = LSH;
                break;
            } else if (answer == "Hypercube") {
                method = HYPERCUBE;
                break;
            } else if (answer == "Frechet") {
                method = FRECHET;
                string answer;
                while(1) {
                    cout << "Please choose between discrete or continuous for Frechet metric: ";
                    cin >> answer;
                    if(answer == "discrete"){
                        frechet_metric = DISCRETE;
                        break;
                    } else if (answer == "continuous") {
                        frechet_metric = CONTINUOUS;
                        break;
                    } else {
                        cerr << "No correct metric specified for Frechet" << endl;
                    }
                }
                break;
            }
        }
    }

    cout << "Input files from command line:" << endl;
    cout << "input_file:" << input_file << endl;
    cout << "query_file:" << query_file << endl;
    cout << "output_file:" << output_file << endl;

    // check input files from cdl
    ifstream input_stream, query_stream;
    if (!input_file.empty()) {
        if (!input_stream.good()) {
            cerr << "Error: input file specified does not exist" << endl;
            return -1;
        } else {
            input_stream.open(input_file, std::ifstream::in);
        }
    }
    if (!query_file.empty()) {
        if (!query_stream.good()) {
            cerr << "Error: query file specified does not exist" << endl;
            return -1;
        } else {
            query_stream.open(query_file);
        }
    }

    // read numeric parameters (globals)
    k = atoi(input.getCmdOption("-k").c_str());
    L = atoi(input.getCmdOption("-L").c_str());
    M = atoi(input.getCmdOption("-M").c_str());
    probes = atoi(input.getCmdOption("-probes").c_str());
    delta = atof(input.getCmdOption("-delta").c_str());

    // fall back to defaults if not given
    if ( method == LSH || method == FRECHET && k == 0) k = DEFAULT_k;
    if ( method == HYPERCUBE && k == 0) k = DEFAULT_CUBE_k;
    if( method == FRECHET && frechet_metric == CONTINUOUS) L = DEFAULT_CONTI_L;
    else {
        if (L == 0) L = DEFAULT_L;
    }
    if (M == 0) M = DEFAULT_M;
    if (probes == 0) probes = DEFAULT_PR;
    if(delta == 0.0) delta = DEFAULT_Delta;

    // get input file if not given
    if (input_file.empty()) {
        do {
            cout << "Please specify input file path:" << endl;
            string input_file_in;
            cin >> input_file_in;
            // check
            input_stream.open(input_file_in);
            if (!input_stream.good()) {
                cout << "Invalid file path. Try again." << endl;
            } else {
                break;
            }
        } while (true);
    }

    // load input dataset
    Dataset *input_dataset;
    if(method == LSH || method == HYPERCUBE)
        input_dataset = new Dataset(input_stream);
    else if (method == FRECHET && frechet_metric == DISCRETE)
    {
        input_dataset = new Dataset(input_stream, true);
    }
    else if(method == FRECHET && frechet_metric == CONTINUOUS)
    {
        input_dataset = new Dataset(input_stream, true, true);
    }
    else {cerr<< "Something went wrong with input Dataset!" << endl; return -1;}
    cout << "Indexing input dataset ";
    if( method == LSH) {
        cout <<"using LSH vector algorithm..." << endl;
        input_dataset->index_LSH(input_dataset->get_size() / DIVIDE_DATASET_FOR_HASHTABLE_SIZE);
    }
    else if ( method == HYPERCUBE) {
        cout <<"using Hypercube algorithm..." << endl;
        input_dataset->index_HyperCube(k);
    }else if( method == FRECHET && frechet_metric == DISCRETE){
        cout << "using LSH Discrete Frechet algorithm..." << endl;
        input_dataset->index_LSH((input_dataset->get_size()) / DIVIDE_DATASET_FOR_HASHTABLE_SIZE, true);
    }
    else if( method == FRECHET && frechet_metric == CONTINUOUS){
        cout << "using LSH Continuous Frechet algorithm..." << endl;
        input_dataset->index_LSH((input_dataset->get_size()) / DIVIDE_DATASET_FOR_HASHTABLE_SIZE, true, true);
    }
    else {
        cerr << "Something went wrong with methods!" <<endl;
        return -1;
    }
    cout << "Done!" << endl;

    // get query file if not given
    if (query_file.empty()) {
        do {
            cout << "Please specify query file path:" << endl;
            string query_file_in;
            cin >> query_file_in;
            // check
            query_stream.open(query_file_in);
            if (!query_stream.good()) {
                cout << "Invalid file path. Try again." << endl;
            } else {
                break;
            }
        } while (true);
    }

    // output stream
    ofstream output_stream;
    if (!output_file.empty()) {
        output_stream.open(output_file.c_str());
    }
        // get output file if not given
    else {
        do {
            cout << "Please specify output file path:" << endl;
            string output_file_in;
            cin >> output_file_in;
            // check
            output_stream.open(output_file_in);
            if (!output_stream.good()) {
                cout << "Invalid file path. Try again." << endl;
            } else {
                break;
            }
        } while (true);
    }

    bool repeat;
    do {
        // load query dataset
        Dataset *query_dataset;
        if(method == LSH || method == HYPERCUBE) {
            query_dataset = new Dataset(query_stream);
        }
        else if (method == FRECHET) {
            if(frechet_metric == DISCRETE)
                query_dataset = new Dataset(query_stream, true);
            else
                query_dataset = new Dataset(query_stream, true, true);
        }
        else {cerr<< "Something went wrong with query Dataset!" << endl; return -1;}
        // initialize average time
        double t_approx = 0.0;
        double t_true = 0.0;
        int q_size = 0;

        /**
        * Both datasets are ready!
        */

        NearestNeighboursSearch nns(*input_dataset, *query_dataset);

        cout << "Calculating exact neighbours";
        vector<Result *> *exactResult;
        if(method == FRECHET && frechet_metric == DISCRETE) {
            cout << " with discrete frechet..." << endl;
            exactResult = nns.calculate_exact_NN_curves();
        }
        else if(method == FRECHET && frechet_metric == CONTINUOUS) {
            cout << " with continuous frechet..." << endl;
            exactResult = nns.calculate_exact_NN_curves(true);
        }
        else {
            cout << "..."<<endl;
            exactResult = nns.calculate_exact_NN(N);
        }
        cout << "Done!" << endl;

        vector<Result *> *approximateResult;
        if(method == LSH) {
            cout << "Calculating approximate neighbours using LSH..." << endl;
            approximateResult = nns.calculate_approximate_NN(N);
        }
        else if(method == HYPERCUBE) {
            cout << "Calculating approximate neighbours using Hypercube..." << endl;
            approximateResult = nns.calculate_approximate_NN_using_cube(N, probes, M);
        }
        else if(method == FRECHET && frechet_metric == DISCRETE){
            cout << "Calculating approximate neighbours using Discrete Frechet..." << endl;
            approximateResult = nns.calculate_approximate_NN_curves();
        }
        else if(method == FRECHET && frechet_metric == CONTINUOUS){
            cout << "Calculating approximate neighbours using Continuous Frechet..." << endl;
            approximateResult = nns.calculate_approximate_NN_curves(true);
        }
        else {cerr<< "Something went wrong with approximate nearest neighbor result!" << endl; return -1;}
        cout << "Done!" << endl;

        // query size found
        q_size = approximateResult->size();

        cout << "Calculating Max Approximation Factor..." << endl;
        vector<double> factors;
        for(int i=0; i < q_size; i++) {
            if(approximateResult->at(i)->NNs->empty())
                factors.push_back(nns.approx_get(exactResult->at(i)->NNs->at(0).dist,0));
            else
                factors.push_back(nns.approx_get(exactResult->at(i)->NNs->at(0).dist,approximateResult->at(i)->NNs->at(0).dist));
        }
        double maf = nns.max_approx_get(factors);
        cout << "Done!" << endl;

        for (int i = 0; i < exactResult->size() && i < approximateResult->size() ; i++) {   // for each query point
            Result *exact_res = exactResult->at(i);
            Result *apprx_res = approximateResult->at(i);
            output_stream << "Query: " << apprx_res->queryID << endl;
            if (method == LSH)
                output_stream << "Algorithm: {LSH_Vector}" << endl;
            else if(method == HYPERCUBE)
                output_stream << "Algorithm: {Hypercube}" << endl;
            else if(method == FRECHET && frechet_metric == DISCRETE)
                output_stream << "Algorithm: {LSH_Frechet_Discrete}" << endl;
            else
                output_stream << "Algorithm: {LSH_Frechet_Continuous}" << endl;
            for (int j = 0; j < apprx_res->NNs->size(); j++) {
                Neighbour &approx_nn = apprx_res->NNs->at(j);
                Neighbour &real_nn = exact_res->NNs->at(j);
                if(!approximateResult->at(i)->NNs->empty())
                    output_stream << "Approximate Nearest neighbor: "<< approx_nn.neighbourID << endl;
                output_stream << "True Nearest neighbor: "<< real_nn.neighbourID << endl;
                if(!approximateResult->at(i)->NNs->empty())
                    output_stream << "distanceApproximate: " << approx_nn.dist << endl;
                output_stream << "distanceTrue: " << real_nn.dist << endl;
            }
            t_approx += apprx_res->dt;
            t_true += exact_res->dt;

            output_stream << endl;

            // delete
            delete exact_res->NNs;
            delete apprx_res->NNs;
            delete exact_res;
            delete apprx_res;
        }
        //count and print average time
        t_approx = t_approx / q_size;
        t_true = t_true / q_size;
        output_stream << "tApproximateAverage: "<< t_approx << " seconds" << endl;
        output_stream << "tTrueAverage: "<< t_true << " seconds" << endl;
        if(maf == DBL_MAX)
            output_stream << "MAF: Infinite" << endl;
        else
            output_stream << "MAF: "<< maf << endl;

        delete exactResult;
        delete approximateResult;

        //delete query dataset
        delete query_dataset;

        string answer;
        cout << "Do you want to enter another query (Y/N)?  ";

        do {
            cin >> answer;
        } while(answer != "Y" && answer != "N" && cout << "Invalid answer! Try again." << endl);

        repeat = (answer == "Y");
        if(repeat) {
            string query_file_in;
            cout << "Query File: ";
            cin >> query_file_in;
            query_stream.open(query_file_in);
            while (!query_stream.good()) {
                cout << "Invalid file path. Try again." << endl;
                cout << "Query File: ";
                cin >> query_file_in;
                query_stream.open(query_file_in);
            }
        }
    } while (repeat);

    //delete input dataset
    delete input_dataset;

    // close output stream
    output_stream.close();

    return 0;
}
