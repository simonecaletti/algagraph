// Conics and quadrics given by a second-degree equation in x, y (and z).
#pragma once

#include "linalg.hpp"

#include <array>
#include <map>
#include <string>
#include <vector>

namespace alga {

// Exponents of x, y, z in a monomial.
using Monomial = std::array<int, 3>;

// A polynomial of degree <= 2 in x, y, z, stored as monomial -> coefficient.
struct Polynomial {
    std::map<Monomial, double> terms;
    bool mentionsZ = false; // true if the variable z appears in the input text

    double coefficient(int ex, int ey, int ez) const;
    double evaluate(double x, double y, double z = 0.0) const;
};

// Parses equations such as "x^2+3*x*y-2y^2+5 = 0", "2x^2 + y^2 = 1" or
// "x^2+y^2-z". A missing right-hand side means "= 0". Throws
// std::invalid_argument with a readable message on syntax errors.
Polynomial parseEquation(const std::string& text);

// Canonical (normal) form reached by a rotation and a translation:
//   sum_i squares[i] * X_i^2 + 2 * linear * X_k = 0    (if hasLinear)
//   sum_i squares[i] * X_i^2 + constant = 0            (otherwise)
// where X_k is the first variable not used by the square terms.
struct NormalForm {
    std::vector<double> squares; // non-zero eigenvalues, decreasing
    bool hasLinear = false;
    double linear = 0.0;
    double constant = 0.0;
};

class Alga {
public:
    explicit Alga(const std::string& equation);

    const std::string& equation() const { return m_equation; }
    const Polynomial& polynomial() const { return m_poly; }
    bool isQuadric() const { return m_dim == 3; }
    // Number of variables: 2 for a conic, 3 for a quadric.
    int dim() const { return m_dim; }

    // Complete (n+1)x(n+1) matrix A' and the nxn matrix A of the quadratic part.
    const Matrix& matrix() const { return m_full; }
    const Matrix& quadraticMatrix() const { return m_quad; }

    // Orthogonal invariants in decreasing order: {I3, I2, I1} for a conic,
    // {I4, I3, I2, I1} for a quadric.
    std::vector<double> invariants() const;

    // Eigenvalues (decreasing) and the rotation matrix P of eigenvectors in columns.
    const EigenSystem& eigen() const { return m_eigen; }

    NormalForm normalForm() const;
    // e.g. "2*x^2+1.80278*y^2-1.80278*z^2-6.30769 = 0"
    std::string normalFormString() const;
    // Polynomial of the normal form, ready to be drawn.
    Polynomial normalFormPolynomial() const;

    // Human readable classification, e.g. "It's a hyperboloid of one sheet."
    std::string type() const;

    // Values smaller than this are treated as zero.
    double tolerance() const { return m_tol; }
    // Rounds numerical noise (relative to the matrix entries) to exactly zero.
    double chop(double value, int degree = 1) const;

private:
    std::string m_equation;
    Polynomial m_poly;
    int m_dim = 2;
    Matrix m_full, m_quad;
    EigenSystem m_eigen;
    double m_tol = 0.0;
};

// Formats sum_i coeff_i * name_i plus a constant, e.g. "x^2-2*y^2+3".
std::string formatPolynomial(const std::vector<std::pair<double, std::string>>& terms,
                             double constant);
// Formats a number with 6 significant digits, without "-0".
std::string formatNumber(double value);

} // namespace alga
