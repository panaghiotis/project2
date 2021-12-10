#ifndef PROJECTEMIRI_POINT_H
#define PROJECTEMIRI_POINT_H

#include <iostream>
#include "Curve.h"

using namespace std;


struct Point{
    unsigned int hashed_ID;
    bool hashed;
    int pos;
    string id;
    vector<double> *coords;
    Curve *curve;
    Point(vector<double> *coordinates, string id = "", int pos = -1, Curve *curve = NULL) :
        id(id), coords(coordinates), hashed_ID(-1), hashed(false), pos(pos), curve(curve) {}
    ~Point() {
        //delete coords;
    }
    void print() const {
        cout << "id:" << id <<endl;
        for(int i=0; i < 5;i++)
            cout << coords->at(i) << " ";
        cout << endl;
    }
    friend ostream& operator<<(ostream& os, const Point& p);
};

#endif
