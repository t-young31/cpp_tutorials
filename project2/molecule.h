//
// Created by tom on 09/01/2021.
//

#ifndef CPP_TUTORIAL_MOLECULE_H
#define CPP_TUTORIAL_MOLECULE_H
#include <Eigen/Dense>
#include <cmath>
#include <fstream>
#include <iostream>
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

        explicit Molecule(string xyz_filename) {                         // Constructor
            extract_from_xyz_file(xyz_filename);
        }

        void extract_from_xyz_file(string &xyz_filename) {
            /*******************************************************
             *  Set coordinates and atomic numbers from a .xyz file
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
    };
}

#endif //CPP_TUTORIAL_MOLECULE_H
