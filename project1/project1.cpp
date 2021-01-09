#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include <stdexcept>
#include "atoms.h"


using namespace std;
using namespace Eigen;

vector<string> split(const string &s, char delim) {
    stringstream ss(s);
    string item;
    vector<string> elems;
    while (getline(ss, item, delim)) {
        if (!item.empty()){
            elems.push_back(item);
        }
    }
    return elems;
}


class Molecule{

    public:
        vector<int> atomic_numbers;
        MatrixXd coords;
        int n_atoms = 0;

        explicit Molecule(string xyz_filename){                         // Constructor
                extract_from_xyz_file(xyz_filename);
        }

    void extract_from_xyz_file(string &xyz_filename){
        /************************************************
         *  Set coordinates and atomic numbers from a .xyz file
         ***********************************************/

        string line, item;
        ifstream xyz_file (xyz_filename);

        bool assigned_n_atoms = false;

        vector<float> coord_list;

        // Iterate through the xyz file
        while (getline(xyz_file, line, '\n')){

            // Ignore any blank lines etc.
            if (line.empty()){
                continue;
            }

            // Assign the number of atoms
            if (!assigned_n_atoms){
                n_atoms = stoi(line);
                assigned_n_atoms = true;
                continue;
            }

            vector<string> xyz_items = split(line, ' ');

            atomic_numbers.push_back(stoi(xyz_items[0]));
            for (int i=1; i < 4; i++){
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

        for (int i=0; i < coord_list.size(); i++){
            coords(i/3, i%3) = coord_list[i];
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

    double angle_ijk(int i, int j, int k){
        /******************************************************
        *  Calculate the angle in radians between three atoms
        *  raises a runtime error if any indices are the same
        *
        *         i    k
        *          \  /            j is the mid atom
        *           j
        *
        *                    ( v_ij . v_jk  )
        *          θ = arccos(--------------)
          *                  ( |v_ij||v_jk| )
        *
        * where v_ij and v_jk are unit vectors
        *****************************************************/
        if (i==j | i==k | j==k){
            throw runtime_error("Angle must be calcd. with three different "
                                "indices");
            }

        // Calculate the dot product over the three components of the vector
        double dot_product = 0;

        for (int c=0; c < 3; c++){
            dot_product += (coords(i, c) - coords(j, c)) * (coords(k, c) - coords(j, c));
        }

        return acos(dot_product / (distance_ij(i, j) * distance_ij(j, k)));

        }

    MatrixXd distance_matrix(){
        /************************************************
         *  Calculate the distance matrix (N x N)
         ***********************************************/
        double dist;

        // N x N matrix of zeros
        MatrixXd dist_mat = MatrixXd::Zero(n_atoms, n_atoms);

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

    void shift_to_com(){
        /************************************************
        *  Shift the molecule so that the center of mass_i
        *  (COM) is centered at the origin
        ***********************************************/
        // Center of mass as a column vector
        Vector3d com = Vector3d::Zero();

        double mass_i;
        double total_mass = 0;

        // COM = Σ_i m_i v_i /  Σ_i m_i
        for (int i=0; i<n_atoms; i++){

            mass_i = atomic_weight(atomic_numbers[i]);
            total_mass += mass_i;

            com += mass_i * coords.row(i);
        }

        // Subtract the COM from each row of the coordinates
        // COM needs transposing into a row vector_
        coords.rowwise() -= com.transpose() / total_mass;
    }

    Vector3d moments_of_inertia(){
        /************************************************
        *  Calculate the moments of inertia as the
        *  eigenvalues of the moment of inertia tensor
        ***********************************************/
        shift_to_com();

        VectorXd atomic_weights(n_atoms);

        for (int i=0; i<n_atoms; i++) atomic_weights(i) = atomic_weight(atomic_numbers[i]);

        auto x = coords.col(0);
        auto y = coords.col(1);
        auto z = coords.col(2);

        VectorXd x_sq = x.array().pow(2);
        VectorXd y_sq = y.array().pow(2);
        VectorXd z_sq = z.array().pow(2);

        Matrix3d I_mat = Matrix3d::Zero();
        // Calculate the diagonal elements
        I_mat(0, 0) = atomic_weights.adjoint() * (y_sq + z_sq);
        I_mat(1, 1) = atomic_weights.adjoint() * (x_sq + z_sq);
        I_mat(2, 2) = atomic_weights.adjoint() * (x_sq + y_sq);

        // Calculate the off diagonals
        I_mat(0, 1 ) = -atomic_weights.adjoint() * (x.cwiseProduct(y));
        I_mat(0, 2 ) = -atomic_weights.adjoint() * (x.cwiseProduct(z));
        I_mat(1, 2 ) = -atomic_weights.adjoint() * (y.cwiseProduct(z));

        // Assign the rest of the bottom left of the symmetric matrix
        I_mat(1, 0) = I_mat(0, 1);
        I_mat(2, 0) = I_mat(0, 2);
        I_mat(2, 1) = I_mat(1, 2);

        EigenSolver<MatrixXd> solver(I_mat);
        return solver.eigenvalues().real();
    }

};


int main(){

    Molecule mol = Molecule("../project1.xyz");
    cout << mol.moments_of_inertia().transpose() << endl;

    return 0;
}
