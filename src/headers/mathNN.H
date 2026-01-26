/*
 * mathNN.H (clean)
 *
 * Linear-algebra / ODE utilities based on Eigen, meant to replace math33.H
 * for the NN/3x3 routines used across projects.
 */

#ifndef MATHNN_H
#define MATHNN_H

#include <Eigen/Dense>
#include <complex>
#include <iostream>
#include <cstdlib>

inline std::complex<double> determinant(std::complex<double> **a, int N = 3) {
    Eigen::MatrixXcd Am(N, N);

    // Copy data from input array a into Eigen matrix Am
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            Am(i, j) = a[i][j];

    std::complex<double> dtr = Am.determinant();
    return dtr;
}

inline std::complex<double>* eigenvalues(std::complex<double> **a, int N = 3) {
    std::complex<double>* eva = new std::complex<double>[N];
    Eigen::MatrixXcd Am(N, N);

    // Copy input matrix a into Eigen matrix Am
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            Am(i, j) = a[i][j];

    // Compute eigenvalues
    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> ces;
    ces.compute(Am);
    if (ces.info() != Eigen::Success) {
        std::cerr << "Error: Failed to calculate eigenvalues." << std::endl;
        std::exit(EXIT_FAILURE);
    }

    // Copy eigenvalues to output array eva
    for (int i = 0; i < N; i++) {
        eva[i] = ces.eigenvalues()[i];
    }

    return eva;
}

inline std::complex<double>** eigenvectors(std::complex<double> **a, int N = 3) {
    Eigen::MatrixXcd Am(N, N);

    // Copy data from input array a into Eigen matrix Am
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            Am(i, j) = a[i][j];
        }
    }

    // Calculate the eigenvalues and eigenvectors of Am
    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> ces;
    ces.compute(Am);
    if (ces.info() != Eigen::Success) {
        std::cerr << "Error: Failed to calculate eigenvectors." << std::endl;
        std::exit(EXIT_FAILURE);  // here returning a null pointer
    }

    // Allocate memory for the output eigenvectors (std::complex<double>**)
    std::complex<double>** eve = new std::complex<double>*[N];
    for (int i = 0; i < N; i++) {
        eve[i] = new std::complex<double>[N];
    }

    // Copy data from the Eigen matrix to the dynamically allocated array eve
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            eve[i][j] = ces.eigenvectors()(i, j);
        }
    }

    return eve;
}

inline std::complex<double>* find_constants(std::complex<double> **EVE, std::complex<double> *qss, std::complex<double> *q0, int N = 3) {
    std::complex<double>* c = new std::complex<double>[N];

    Eigen::VectorXcd q_stst(N);
    Eigen::VectorXcd q_init(N);
    Eigen::VectorXcd noto(N);
    Eigen::VectorXcd cnstnt(N);
    Eigen::MatrixXcd eigvec(N, N);

    // Fill the vectors and matrix from input
    for (int i = 0; i < N; i++) q_init(i) = q0[i];
    for (int i = 0; i < N; i++) q_stst(i) = qss[i];
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++) eigvec(i, j) = EVE[i][j];

    if (norm(eigvec.determinant()) == 0) { std::cout << "EVE is singular!\n"; std::exit(EXIT_FAILURE); }
    // Calculate the deviation from steady state
    noto = q_init - q_stst;

    // Solve the linear system to find the constants
    cnstnt = eigvec.partialPivLu().solve(noto);

    // Convert to std::complex array for return
    for (int i = 0; i < N; i++) c[i] = cnstnt(i);

    return c;
}

inline std::complex<double>* steady_state_solution(std::complex<double> **a, std::complex<double> *b, int N = 3) {
    std::complex<double>* qss = new std::complex<double>[N];
    Eigen::VectorXcd q_stst(N);
    Eigen::VectorXcd b_dsh(N);
    Eigen::MatrixXcd Am(N, N);

    for (int i = 0; i < N; i++) b_dsh(i) = -b[i];

    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++) Am(i, j) = a[i][j];

    // Solve the linear system for steady state solution
    q_stst = Am.partialPivLu().solve(b_dsh);

    for (int i = 0; i < N; i++) qss[i] = q_stst(i);
    return qss;
}

// Compatibility wrapper (math33.H-style):
// c = V^{-1} (q0 - qss), where qss solves A qss + b = 0 and V are eigenvectors of A.
inline std::complex<double>* constants(std::complex<double> **a, std::complex<double> *b, std::complex<double> *q0, int N = 3) {
    std::complex<double>* qss = steady_state_solution(a, b, N);
    std::complex<double>** eve = eigenvectors(a, N);
    std::complex<double>* c   = find_constants(eve, qss, q0, N);
    delete[] qss;
    for (int i = 0; i < N; ++i) delete[] eve[i];
    delete[] eve;
    return c;
}

/*
 Coupled RK4: To solve coupled equations of the form:
          
          df/dt=A*f+C,
          
          where f=(f0, f1, f2),
          
              |a00  a01 a02|
          A = |a10  a11 a12|  and C=(c0, c1, c2)
              |a20  a21 a22|
*/    

inline std::complex<double> * Runge_Kutta_4 (std::complex<double>* f, std::complex<double> **a, std::complex<double>* b, double dt, int N = 3){
    std::complex<double> k[4][N];
    std::complex<double> sum[N] = {0};
    double d;
    
    for (int i=0; i<=3; i++){
        if (i==0) d=0.;
            else if (i==3) d=1.;
            else d=0.5;
        for (int j = 0 ; j < N; j++) sum[j] = std::complex<double> (0.,0.);
        for (int j = 0 ; j < N; j++){
            for (int jj = 0 ; jj < N; jj++) sum[j]=sum[j]+a[j][jj]*(f[jj]+d*k[i-1][jj]);
            k[i][j]=dt*(sum[j]+b[j]);
            }
                }
    for (int j = 0; j < N; j++)
        f[j]=f[j]+(1/6.)*(k[0][j]+2.*k[1][j]+2.*k[2][j]+k[3][j]);

    return f;
}

/*
  Single complex RK4: To solve equations of the form:
        df/dt=a*f+b.
*/
std::complex<double> Runge_Kutta_moco_4 (std::complex<double> f, 
                     std::complex<double> a, std::complex<double> b, double dt){
    
    std::complex<double> k0=dt*(a*f+b);
    std::complex<double> k1=dt*(a*(f+0.5*k0)+b);
    std::complex<double> k2=dt*(a*(f+0.5*k1)+b);
    std::complex<double> k3=dt*(a*(f+k2)+b);
    
    f=f+(1/6.)*(k0+2.*k1+2.*k2+k3);
    
    return f;
}

/*
  Single real RK4: To solve equations of the form:
        df/dt=a*f+b.
*/
double Runge_Kutta_mono_4 (double f, double a, double b, double dt){
    
    double k0=dt*(a*f+b);
    double k1=dt*(a*(f+0.5*k0)+b);
    double k2=dt*(a*(f+0.5*k1)+b);
    double k3=dt*(a*(f+k2)+b);
    
    f=f+(1/6.)*(k0+2.*k1+2.*k2+k3);
    
    return f;
    }

#endif // MATHNN_H