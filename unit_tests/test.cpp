#include "../include/acutest.h"
#include "../include/utilities.h"
#include <exception>
#include <stdexcept>
#include <string>

void test_discrete_frechet(void) {
    char expected[] = "Discrete Frechet returns positive number (greater than 0)";
    char error[] = "Discrete Frechet returns negative number (or 0)";
    char produced[100];

    //initialize values for our DiscreteFrechet function
    vector<pair<double,double>> p , q;
    for(int i=0; i < 10; i++) {
        p.push_back({i,2*i});
        q.push_back({i+4, 10*(i+6)});
    }
    //Test if it returns a positive (greater than 0) number as for distance
    long double dist = discreteFrechet_distance(p,q);
    if(dist <= 0)
        strcpy(produced,error);
    else
        strcpy(produced,expected);

    //Do the test!
    TEST_CHECK((strcmp(produced,expected)==0));
    cout << endl;
    cout << "Expected: " << expected << endl;
    cout << "Produced: " << produced << endl;

    // my compiler couldn't use TEST_MSG
    //TEST_MSG("Expected: %s", expected);
    //TEST_MSG("Produced: %s", produced);
}

void test_tree_recursion(void) {
    int produced; // we except value 0 which is the root's placement in tree

    //initialize dump values for a simple tree
    vector<int> tree(7);

    //Test if our BST function returns the root (node num = 0)
    produced = PostOrderTraversal(0,5,tree);

    //Do the test!
    TEST_CHECK(produced == 0);
    cout << endl;
    cout << "Expected: 0" << endl;
    cout << "Produced: " << produced << endl;

    //TEST_MSG("Expected: 0");
    //TEST_MSG("Produced: %d", produced);
}

void test_c_array(void) {
    double **produced;

    //initialize values for our get_C function
    vector<pair<double,double>> p , q;
    for(int i=0; i < 10; i++) {
        p.push_back({i,2*i});
        q.push_back({i+4, 10*(i+6)});
    }

    //test if array is allocated correctly
    produced = get_C(p,q);

    //Do the test!
    TEST_CHECK(produced != NULL);
    cout << endl;
    cout <<"Expected a float number." << endl;
    cout<< "Produced: " << produced[0][0] << endl;

    //TEST_MSG("Produced: %f", produced[0][0]);
}

void test_mean_curve_filtering(void) {
    //initialize values for filtering
    vector<pair<double,double>> produced;
    for(int i=0; i < 20; i++) // we will commit 20 values with epsilon=1 so as to cut some of those
        produced.push_back({i,i});

    //call the function
    produced = filtering(produced,1);

    //Do the test!
    TEST_CHECK(produced.size() < 20);
    cout << endl;
    cout << "Expected less than 20 coordinates: " << endl;
    cout << "Produced " << produced.size() << " coordinates" << endl;
}

TEST_LIST = {
        { "test_discrete_frechet", test_discrete_frechet },
        { "test_tree_recursion", test_tree_recursion },
        { "test_c_array", test_c_array },
        { "test_mean_curve_filtering", test_mean_curve_filtering },
        { NULL, NULL }
};

















