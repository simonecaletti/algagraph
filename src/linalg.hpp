// Small dense linear algebra for the 2x2 ... 4x4 symmetric matrices that
// describe conics and quadrics. Replaces the ROOT classes (TMatrixDSym,
// TMatrixDSymEigen) used by the 2015 version.
#pragma once

#include <cstddef>
#include <vector>

namespace alga {

class Matrix {
public:
    Matrix() = default;
    Matrix(std::size_t rows, std::size_t cols, double value = 0.0);

    static Matrix identity(std::size_t n);

    std::size_t rows() const { return m_rows; }
    std::size_t cols() const { return m_cols; }

    double& operator()(std::size_t i, std::size_t j) { return m_data[i * m_cols + j]; }
    double operator()(std::size_t i, std::size_t j) const { return m_data[i * m_cols + j]; }

    Matrix transposed() const;
    // Largest absolute value of the entries (used to scale tolerances).
    double maxAbs() const;

private:
    std::size_t m_rows = 0, m_cols = 0;
    std::vector<double> m_data;
};

Matrix operator*(const Matrix& a, const Matrix& b);

// Determinant by Gaussian elimination with partial pivoting.
double determinant(Matrix m);

// Sum of all principal minors of order k (k = 1 is the trace, k = n the determinant).
double principalMinorSum(const Matrix& m, std::size_t k);

struct EigenSystem {
    std::vector<double> values; // sorted in decreasing order
    Matrix vectors;             // eigenvectors in columns, det = +1
};

// Eigenvalues and orthonormal eigenvectors of a symmetric matrix (cyclic Jacobi).
EigenSystem symmetricEigen(const Matrix& m);

} // namespace alga
