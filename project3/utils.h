//
// Created by tom on 09/01/2021.
//

#ifndef CPP_TUTORIAL_UTILS_H
#define CPP_TUTORIAL_UTILS_H
#include <vector>
#include <string>
#include <Eigen/Dense>

using namespace std;

namespace utils {

    vector<string> split(const string &s, char delim);

    Eigen::MatrixXd symmetric_matrix_from_file(const string &filename);

    Eigen::VectorXd ee_from_file(const string &filename, const int &length);

}


#endif //CPP_TUTORIAL_UTILS_H
