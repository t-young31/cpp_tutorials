//
// Created by tom on 09/01/2021.
//

#ifndef CPP_TUTORIAL_MOLECULE_H
#define CPP_TUTORIAL_MOLECULE_H
#include <Eigen/Dense>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include "utils.h"


using namespace std;


namespace base_mol {

    class Molecule {

    public:
        vector<int> atomic_numbers;
        Eigen::MatrixXd coords;
        int n_atoms = 0;
        int n_elecs = 0;

        explicit Molecule(string xyz_filename) {                         // Constructor
            extract_from_xyz_file(xyz_filename);

            // Set the number of electrons as the number of protons, assuming neutral
            for (auto z : atomic_numbers){
                n_elecs += z;
            }
        }

        void extract_from_xyz_file(string &xyz_filename) {
            /*******************************************************
             *  Set coordinates and atomic numbers from a .xyz file
             *  where the units are a0 *not* angstroms
             ******************************************************/

            string line, item;
            ifstream xyz_file(xyz_filename);

            bool assigned_n_atoms = false;

            std::vector<float> coord_list;

            // Iterate through the xyz file
            while (getline(xyz_file, line, '\n')) {

                // Ignore any blank lines etc.
                if (line.empty()) {
                    continue;
                }

                // Assign the number of atoms
                if (!assigned_n_atoms) {
                    n_atoms = stoi(line);
                    assigned_n_atoms = true;
                    continue;
                }

                vector<string> xyz_items = utils::split(line, ' ');

                atomic_numbers.push_back(stoi(xyz_items[0]));
                for (int i = 1; i < 4; i++) {
                    coord_list.push_back(stod(xyz_items[i]));
                }

            }
            xyz_file.close();

            if (atomic_numbers.size() != n_atoms) {
                cout << atomic_numbers.size() << " " << n_atoms << endl;
                throw runtime_error("Number of atoms not equal to the number declared");
            }

            // Now we can assign the coordinate as a Nx3 matrix, after defining it's shape
            coords.resize(n_atoms, 3);

            for (int i = 0; i < coord_list.size(); i++) {
                coords(i / 3, i % 3) = coord_list[i];
            }
        }

        double distance_ij(int i, int j){
            /******************************************************
             *  Calculate the distance between to atoms i and j as
             *  √((x_i - x_j)^2 + (y_i - y_j)^2 + (z_i - z_j)^2)
             *****************************************************/

            return sqrt(pow((coords(i, 0) - coords(j, 0)), 2)
                        + pow((coords(i, 1) - coords(j, 1)), 2)
                        + pow((coords(i, 2) - coords(j, 2)), 2));
        }

        Eigen::MatrixXd distance_matrix(){
            /************************************************
             *  Calculate the distance matrix (N x N)
             ***********************************************/
            double dist;

            // N x N matrix of zeros
            Eigen::MatrixXd dist_mat = Eigen::MatrixXd::Zero(n_atoms, n_atoms);

            for (int i=0; i < n_atoms; i++){
                for (int j=i+1; j < n_atoms; j++){

                    dist = distance_ij(i, j);

                    // And set the values of the symmetric matrix
                    dist_mat(i, j) = dist;
                    dist_mat(j, i) = dist;
                }
            }
            return dist_mat;
        }
    };
}

#endif //CPP_TUTORIAL_MOLECULE_H
