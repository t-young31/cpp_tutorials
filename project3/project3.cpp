//
// Created by tom on 09/01/2021. following https://github.com/CrawfordGroup/ProgrammingProjects/tree/master/Project%2303
//
#include <utility>
#include "molecule.h"
#include "utils.h"



Eigen::MatrixXd matrix_minus_half(Eigen::MatrixXd const &matrix){
    /*******************************************************
     *  Compute the minus half power of a matrix as
     *
     *  M^(-1/2) = A D^(-1/2) A^T
     *
     *  where M is the matrix, A the matrix of eigenvectors
     *  and D a diagonal matrix of eigenvalues
     ******************************************************/
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(matrix);
    Eigen::VectorXd lambda_mhalf = solver.eigenvalues().array().inverse().sqrt();

    return solver.eigenvectors() * lambda_mhalf.asDiagonal() * solver.eigenvectors().transpose();
}


Eigen::MatrixXd eigenvectors(Eigen::MatrixXd const&matrix){
    /*******************************************************
     *  Calculate the matrix of eigenvectors (columns)
     *
     *  M = A D A^T
     ******************************************************/
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(matrix);
    return solver.eigenvectors();
}


class Molecule : public base_mol::Molecule{
public:

    int n_occ;                 // Number of occupied orbitals
    int n_ao;                 // Number of atomic orbitals

    double E_nuc;              // Nuclear repulsion energy
    Eigen::MatrixXd S;         // Overlap matrix
    Eigen::MatrixXd S_orthog;  // S^(-1/2) orthogonalising matrix
    Eigen::MatrixXd T;         // Kinetic energy matrix
    Eigen::MatrixXd V_en;      // Electron-nuclear attraction matrix
    Eigen::VectorXd V_ee_vec;  // Flat (1D) array of electron-electron repulsion integrals

    Eigen::MatrixXd H_core;     // Core Hamiltonian matrix
    Eigen::MatrixXd F;          // Fock matrix

    Eigen::MatrixXd C_prime;    // Matrix of molecular coefficients in the orthogonal AO basis
    Eigen::MatrixXd C;          // .                                in the non-orthogonal AO basis

    Eigen::MatrixXd D;          // Density matrix


    explicit Molecule(string xyz_filename,
                      const string &overlap_filename,
                      const string &kinetic_filename,
                      const string &elec_nuc_filename,
                      const string &elec_elec_filename):
        base_mol::Molecule(move(xyz_filename)){

        n_occ = n_elecs / 2;   // Number of occupied spatial orbitals is half the number of electrons
        set_e_nuc();           // Calculate and set the nuclear repulsion (E_h)

        S = utils::symmetric_matrix_from_file(overlap_filename);
        S_orthog = matrix_minus_half(S);
        n_ao = S.rows();      // n_rows = n_cols = number of atomic orbitals

        T = utils::symmetric_matrix_from_file(kinetic_filename);
        V_en = utils::symmetric_matrix_from_file(elec_nuc_filename);
        V_ee_vec = utils::ee_from_file(elec_elec_filename, n_ao);

        H_core = T + V_en;

        // Construct an initial guess Fock matrix
        F = H_core;

        // And build an initial guess of the density matrix
        C_prime = eigenvectors(S_orthog.transpose() * F * S_orthog);
        // back transform
        C = S_orthog * C_prime;

        D.resize(n_ao, n_ao);
        set_density_matrix();
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

                total_E_nuc += (atomic_numbers[i] * atomic_numbers[j]
                                / distance_ij(i, j));

            }
        }
        // Finally set the total value
        E_nuc = total_E_nuc;
    }

    void set_density_matrix(){
        /*******************************************************
         *  Set the elements of the density matrix given a set
         *  of MO coefficients
         ******************************************************/
        double d_mu_nu;

         for (int mu = 0; mu < n_ao; mu ++){
             for (int nu = 0; nu < n_ao; nu++){
                 d_mu_nu = 0;

                 // Sum over occupied orbitals
                 for (int m = 0; m < n_occ; m++){
                     d_mu_nu += C(mu, m) * C(nu, m);
                 }

                 D(mu, nu) = d_mu_nu;
             }
         }

    }

    void set_fock_matrix(){
        /*******************************************************
        *  Update the elements of the Fock matrix as
        *
        * F_ij = H^core_ij + Σ_nm D_nm(2xV_ijnm - V_injm)
        *
        * where D is the current desnity matrix and V the
        * tensor of two index integrals
        ******************************************************/
        double sum;

        for (int mu = 0; mu < n_ao; mu++) {
            for (int nu = mu; nu < n_ao; nu++) {
                sum = 0;

                for (int lambda = 0; lambda < n_ao; lambda++){
                    for (int sigma = 0; sigma < n_ao; sigma++){
                        sum += D(lambda, sigma) * (2*V_ee(mu, nu, lambda, sigma) - V_ee(mu, lambda, nu, sigma));
                    } // sigma
                } // lambda

                F(mu, nu) = H_core(mu, nu) + sum;
                // set the equiv. element of the symmetric matrix
                F(nu, mu) = H_core(mu, nu) + sum;

            } // nu
        } // mu
    }

    double V_ee(int mu, int nu, int lambda, int sigma){
        /*******************************************************
        *  Return an element of the V_ee tensor stored as a
        *  1D vector of only the unique elements
        ******************************************************/
        if (mu < nu){
            swap(mu, nu);
        }
        if (lambda < sigma){
            swap(lambda, sigma);
        }

        int mu_nu = mu * (mu + 1)/2 + nu;
        int lambda_sigma = lambda * (lambda + 1)/2 + sigma;

        if (mu_nu < lambda_sigma){
            swap(mu_nu, lambda_sigma);
        }
        // Return the element in the vector with the compound index
        return V_ee_vec(mu_nu*(mu_nu+1)/2+lambda_sigma);
    }

    double E(){
        /*******************************************************
        *  Compute the total electronic energy, including the
        *  nuclear repulsion
        ******************************************************/
        double E_elec = 0;

        for (int mu = 0; mu < n_ao; mu++){
            for (int nu = 0; nu < n_ao; nu++){
                E_elec += D(mu, nu) * (H_core(mu, nu) + F(mu, nu));
            }
        }
        return E_elec + E_nuc;
    }

    void do_scf(){
        /*******************************************************
        *  Do the self-consistent field (SCF) cycles
        ******************************************************/
        const double delta_E_tol = 1E-6;      // Tolerance on the energy convergence
        double prev_E = 0;
        double curr_E = E();
        double delta_E = fabs(curr_E - prev_E);

        cout << "Energy_i       Energy_i-1    Delta(Energy)"  << endl;

        while (delta_E > delta_E_tol){
            prev_E = curr_E;

            set_fock_matrix();
            C_prime = eigenvectors(S_orthog.transpose() * F * S_orthog);
            C = S_orthog * C_prime;
            set_density_matrix();

            curr_E = E();
            cout << setprecision(9) << curr_E << "   " << prev_E << "   " << delta_E << endl;
            delta_E = fabs(curr_E - prev_E);
        }
    }

};



int main(){

    Molecule mol = Molecule("/home/tom/repos/cpp_tutorial/project3/h2o/geom.xyz",
                            "/home/tom/repos/cpp_tutorial/project3/h2o/STO-3G/s.dat",
                            "/home/tom/repos/cpp_tutorial/project3/h2o/STO-3G/t.dat",
                            "/home/tom/repos/cpp_tutorial/project3/h2o/STO-3G/v.dat",
                            "/home/tom/repos/cpp_tutorial/project3/h2o/STO-3G/eri.dat");

    // cout << mol.V_ee(5, 0, 2, 1) << endl;

    // cout << "should be the same.." << endl;
    // cout << mol.V_ee(0, 5, 2, 1) << endl;
    // cout << mol.V_ee(5, 0, 1, 2) << endl;
    // cout << mol.V_ee(0, 5, 2, 1) << endl;
    // cout << mol.V_ee(2, 1, 5, 0) << endl;
    // cout << mol.V_ee(1, 2, 0, 5) << endl;
    // cout << mol.V_ee(2, 1, 0, 5) << endl;
    // cout << "should not be the same.." << endl;
    // cout << mol.V_ee(5, 2, 0, 1) << endl;

    //for (int i=0; i < 7; i ++){
    //    cout << mol.D.row(i) << endl;
    //}

    //cout << mol.E() << endl;

    mol.do_scf();

    return 0;
}

