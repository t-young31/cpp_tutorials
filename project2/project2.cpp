//
// Created by tom on 09/01/2021.
//
#include <utility>
#include "molecule.h"
#include "atoms.h"

class Molecule : public base_mol::Molecule{         // New molecule class inherited from base class
public:

    Eigen::MatrixXd hessian;

    explicit Molecule(string xyz_filename, string hess_filename):
        base_mol::Molecule(move(xyz_filename)){     // Call the constructor of the base class

        // Populate the mass-weighted hessian matrix
        extract_hessian_from_file(hess_filename);
        mass_weight_hessian();
    }

    void extract_hessian_from_file(string &hess_filename){
        /*******************************************************
         *  Set the components of the hessian matrix from a file
         ******************************************************/
        string line, item;
        ifstream hess_file(hess_filename);

        std::vector<double> hess_list;
        int line_n = 0;

        // Iterate through the xyz file
        while (getline(hess_file, line, '\n')) {

            // Ignore any blank lines etc.
            if (line.empty()) {
                continue;
            }

            // Check the number of atoms matches
            if (line_n == 0) {
                if (stoi(line) != n_atoms){
                    throw runtime_error("Number of atoms in the hessian doesn't match the number"
                                        "declared in the xyz file");
                }
                line_n += 1;
                continue;
            }

            // Split the line on spaces and add the floats to the hessian list
            vector<string> hess_items = utils::split(line, ' ');
            for (int i = 0; i < 3; i++) {
                hess_list.push_back(stod(hess_items[i]));
            }

            // Increment the line number
            line_n += 1;
        }
        hess_file.close();

        // From the 1D vector of Hessian elements ordered F_x1x1, F_x1y1, F_x1z1, F_x1x2 ...
        // generate a 3Nx3N matrix
        int size = 3 * n_atoms;
        hessian.resize(size, size);

        for (int i = 0; i < hess_list.size(); i++) {
            hessian(i / size, i % size) = hess_list[i];
        }
    }

    void mass_weight_hessian(){
        /*****************************************************
         * Divide the hessian matrix by the product of atomic
         * weights for each atom as
         *
         *                F_ij
         * F_ij^M  =  ------------
         *             √(m_i x m_j)
         ****************************************************/
         for (int i=0; i < n_atoms; i ++){
             double m_i = atomic_weight(atomic_numbers[i]);

             // Small matrix so set all the values
             for (int j=0; j < n_atoms; j++){

                 double m_j = atomic_weight(atomic_numbers[j]);
                 double weight = 1.0 / sqrt(m_i * m_j);

                 for (int k_i=0; k_i < 3; k_i ++){
                     for (int k_j=0; k_j < 3; k_j ++){
                         hessian(3*i+k_i, 3*j+k_j) *= weight;
                     }
                 }

             }
         }
     // now mass weighted!
    }

    vector<double> harmonic_freqs(){
        /*****************************************************
         * Compute the harmonic normal modes by diagonalising
         * the mass weighted Hessian matrix. Only the
         * rotational and vibrational modes are returned
         ****************************************************/
        vector<double> freqs;

        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(hessian);
        Eigen::VectorXd lambda = solver.eigenvalues().real();

        double amu_to_kg = 1.66053906660e-27;                            // kg Da^-1
        double eh_to_j = 4.3597447222071e-18;                            // J E_h^-1
        double a0_to_m = 5.291772109E-11;                                // m a_0^-1
        double c_cm = 2.998e10;                                          // cm s^-1
        double two_pi = 2.0 * 3.14159265359;                             // 2π

        double conversion = sqrt(eh_to_j / amu_to_kg) / (two_pi * c_cm * a0_to_m);

        // Miss the fist 3 translational modes, which are zero
        for (int i=2; i<lambda.size(); i++){
            freqs.push_back(conversion * sqrt(lambda[i]));
        }
        return freqs;
    }

};



int main(){

    Molecule mol = Molecule("/home/tom/repos/cpp_tutorial/project2/h2o.xyz",
                            "/home/tom/repos/cpp_tutorial/project2/h2o_hessian.txt");

    for (auto &freq : mol.harmonic_freqs()){
        cout << freq  << " cm^-1" << endl;
    }

    return 0;
}

