//
// Created by tom on 09/01/2021.
//
#include <utility>
#include "molecule.h"


Eigen::MatrixXd symmetric_matrix_from_file(const string &filename){
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


class Molecule : public base_mol::Molecule{
public:

    double E_nuc;         // Nuclear repulsion energy
    Eigen::MatrixXd S;    // Overlap matrix with elements

    explicit Molecule(string xyz_filename, const string &overlap_filename):
        base_mol::Molecule(move(xyz_filename)){

        set_e_nuc();
        S = symmetric_matrix_from_file(overlap_filename);
    }

    void set_e_nuc(){
        /*******************************************************
         *  Calculate the nuclear repulsion energy from the
         *  atomic coordinates in atomic units
         *
         *  E_nuc = Σ_i Σ_j>i Z_i Z_j/R_ij
         ******************************************************/

        double total_E_nuc = 0;
        const double ang_to_a0 = 1.8897259886;              // a_0 A^-1

        // Enumerate over the unique pairs of atoms (i, j)
        for (int i = 0; i < n_atoms; i++){
            for (int j = i+1; j < n_atoms; j++){

                cout << i << ' ' << j << endl;
                total_E_nuc += (atomic_numbers[i] * atomic_numbers[j]
                                / distance_ij(i, j));

            }
        }
        // Finally set the total value
        E_nuc = total_E_nuc;
    }

};



int main(){

    Molecule mol = Molecule("/home/tom/repos/cpp_tutorial/project3/h2o/geom.xyz",
                            "/home/tom/repos/cpp_tutorial/project3/h2o/STO-3G/s.dat");

    cout << mol.E_nuc << endl;

    return 0;
}

