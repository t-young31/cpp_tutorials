//
// Created by tom on 09/01/2021.
//

#include <sstream>
#include "utils.h"


vector<string> utils::split(const string &s, char delim) {
    /*********************************************************
     *  Split a string by a delimiter into a vector of strings
     ********************************************************/
    stringstream stream(s);
    string item;
    vector<string> elems;
    while (getline(stream, item, delim)) {
        if (!item.empty()){
            elems.push_back(item);
        }
    }
    return elems;
}

