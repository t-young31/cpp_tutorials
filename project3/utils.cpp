//
// Created by tom on 09/01/2021.
//

#include <sstream>
#include <fstream>
#include <iostream>

#include <Eigen/Dense>
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

Eigen::MatrixXd utils::symmetric_matrix_from_file(const string &filename){
    /*******************************************************
     *  Extract a symmetric matrix from a file, where the
     *  file is formatted as:
     *
     *  idx1  idx2 value
     ******************************************************/
    string line, item;
    ifstream input_file(filename);

    std::vector<long> idxs1;
    std::vector<long> idxs2;
    std::vector<double> values;

    // Iterate through the file line by line
    while (getline(input_file, line, '\n')) {

        // Ignore any blank lines
        if (line.empty()) {
            continue;
        }
        vector<string> items = utils::split(line, ' ');
        if (items.size() != 3){
            throw runtime_error("Overlap matrix file not correctly formatted."
                                "Expecting lines in the format: idx1  idx2 value");
        }
        idxs1.push_back(stol(items[0]));
        idxs2.push_back(stol(items[1]));
        values.push_back(stod(items[2]));

    }
    long max_idx1 = *max_element(begin(idxs1), end(idxs1));
    long min_idx1 = *min_element(begin(idxs1), end(idxs1));
    long size = max_idx1 - min_idx1 + 1;     // Allows for indexing from 0 or 1 in the file

    Eigen::MatrixXd matrix(size, size);
    for (int n = 0; n < values.size(); n++){

        long i = idxs1[n] - min_idx1;
        long j = idxs2[n] - min_idx1;

        matrix(i, j) = values[n];
        // And because the matrix is symmetric the flipped element can also be set
        matrix(j, i) = values[n];
    }

    return matrix;
}


Eigen::VectorXd utils::ee_from_file(const string &filename, const int &length){
    /*******************************************************
     *  Extract a 1D vector of values for 4 index electron
     *  repulsion integrals from a file
     *
     *  idx1  idx2 idx3 idx4  value
     ******************************************************/
    // Set up the vector we're going to populate
    int M = (length*(length+1))/2;
    int max_ijkl = (M*(M+1))/2;
    Eigen::VectorXd ee_vector(max_ijkl+1);

    string line;
    int i, j, k, l;
    double value;

    ifstream input_file(filename);

    // Iterate through the file line by line
    while (getline(input_file, line, '\n')) {

        // Ignore any blank lines
        if (line.empty()) {
            continue;
        }
        vector<string> items = utils::split(line, ' ');
        if (items.size() != 5){
            throw runtime_error("Unexpected line length");
        }

        // Indexes from 1
        i = stoi(items[0]) - 1;
        j = stoi(items[1]) - 1;
        k = stoi(items[2]) - 1;
        l = stoi(items[3]) - 1;

        value = stod(items[4]);

        if (i < j){
            swap(i, j);
        }
        if (k < l){
            swap(k, l);
        }

        int ij = i * (i + 1) / 2 + j;
        int kl = k * (k + 1) / 2 + l;

        if (ij < kl){
            swap(ij, kl);
        }
        int idx = ij*(ij+1)/2+kl;

        if (idx >= ee_vector.size()){
            throw runtime_error("Index not settable, exceeded the "
                                "size of the array");
        }
        ee_vector(idx) = value;
    }
    return ee_vector;
}


