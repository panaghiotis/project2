#include <iostream>
#include <fstream>
#include "../header/neighbour_search.h"
#include "../header/InputParser.h"
#include "../header/Dataset.h"

using namespace std;

// globals
unsigned int k = 0;           // number of hi hash functions, also dimension d'
unsigned int M = 0;           // number of points to be checked
unsigned int probes = 0;      // number of vertices to check
unsigned int N = 0;           // number of neighbours to search
long double R = 0;            // radius for range search
unsigned int L = 0;           // not really needed here (it's for Dataset.cpp)

#define DEFAULT_k 14          // number of dimension d' (also k: number of h)
#define DEFAULT_M 10          // number of points to be checked
#define DEFAULT_PR 2          // number of vertices to be checked
#define DEFAULT_N 1           // number of nearest neighbours
#define DEFAULT_R 10000       // radius for range search

#define DIVIDE_DATASET_FOR_HASHTABLE_SIZE 8


int main(int argc, char **argv) {
    // read input
    InputParser input(argc, argv);
    const string &input_file = input.getCmdOption("-i");
    const string &query_file = input.getCmdOption("-q");
    const string &output_file = input.getCmdOption("-o");

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
    M = atoi(input.getCmdOption("-M").c_str());
    probes = atoi(input.getCmdOption("-probes").c_str());
    N = atoi(input.getCmdOption("-N").c_str());
    R = atoi(input.getCmdOption("-R").c_str());

    // fall back to defaults if not given
    if (k == 0) k = DEFAULT_k;
    if (M == 0) M = DEFAULT_M;
    if (probes == 0) probes = DEFAULT_PR;
    if (N == 0) N = DEFAULT_N;
    if (R == 0) R = DEFAULT_R;

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
    Dataset input_dataset(input_stream);
    cout << "Indexing input dataset..." << endl;
    input_dataset.index_HyperCube(k);
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
        Dataset *query_dataset = new Dataset(query_stream);
        /**
        * Both datasets are ready!
        */

        NearestNeighboursSearch nns(input_dataset, *query_dataset);

        cout << "Calculating exact neighbours..." << endl;
        vector<Result *> *exactResult = nns.calculate_exact_NN(N);
        cout << "Done!" << endl;

        cout << "Calculating approximate neighbours using Hypercube..." << endl;
        vector<Result *> *approximateResult = nns.calculate_approximate_NN_using_cube(N, probes, M);
        cout << "Done!" << endl;

        cout << "Calculating range search using Hypercube..." << endl;
        vector<vector<string> *> *range_search_res = nns.range_search_using_cube(R, probes, M);
        cout << "Done!" << endl;

        for (int i = 0; i < exactResult->size() && i < approximateResult->size() && i < range_search_res->size(); i++) {   // for each query point
            Result *exact_res = exactResult->at(i);
            Result *apprx_res = approximateResult->at(i);
            output_stream << "Query: " << apprx_res->queryID << endl;
            for (int j = 0; j < apprx_res->NNs->size(); j++) {
                Neighbour &approx_nn = apprx_res->NNs->at(j);
                Neighbour &real_nn = exact_res->NNs->at(j);
                output_stream << "Nearest neighbor-" << (j + 1) << ": " << approx_nn.neighbourID << endl;
                output_stream << "distanceHypercube: " << approx_nn.dist << endl;
                output_stream << "distanceTrue: " << real_nn.dist << endl;
            }
            output_stream << "tHypercube: " << apprx_res->dt << " micro secs" << endl;
            output_stream << "tTrue: " << exact_res->dt << " micro secs" << endl;

            vector<string> *r_nns = range_search_res->at(i);
            output_stream << "R-near neighbours:" << endl;
            for (int j = 0; j < r_nns->size(); j++) {
                output_stream << r_nns->at(j) << endl;
            }

            output_stream << endl;

            // delete
            delete exact_res->NNs;
            delete apprx_res->NNs;
            delete exact_res;
            delete apprx_res;
            delete r_nns;
        }

        delete exactResult;
        delete approximateResult;
        delete range_search_res;

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

    // close output stream
    output_stream.close();

    return 0;
}