#include "linalg.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <utility>

namespace alga {

Matrix::Matrix(std::size_t rows, std::size_t cols, double value)
    : m_rows(rows), m_cols(cols), m_data(rows * cols, value) {}

Matrix Matrix::identity(std::size_t n) {
    Matrix id(n, n);
    for (std::size_t i = 0; i < n; ++i) id(i, i) = 1.0;
    return id;
}

Matrix Matrix::transposed() const {
    Matrix t(m_cols, m_rows);
    for (std::size_t i = 0; i < m_rows; ++i)
        for (std::size_t j = 0; j < m_cols; ++j) t(j, i) = (*this)(i, j);
    return t;
}

double Matrix::maxAbs() const {
    double m = 0.0;
    for (double v : m_data) m = std::max(m, std::abs(v));
    return m;
}

Matrix operator*(const Matrix& a, const Matrix& b) {
    if (a.cols() != b.rows()) throw std::invalid_argument("matrix size mismatch");
    Matrix c(a.rows(), b.cols());
    for (std::size_t i = 0; i < a.rows(); ++i)
        for (std::size_t j = 0; j < b.cols(); ++j) {
            double s = 0.0;
            for (std::size_t k = 0; k < a.cols(); ++k) s += a(i, k) * b(k, j);
            c(i, j) = s;
        }
    return c;
}

double determinant(Matrix m) {
    const std::size_t n = m.rows();
    double det = 1.0;
    for (std::size_t col = 0; col < n; ++col) {
        std::size_t pivot = col;
        for (std::size_t r = col + 1; r < n; ++r)
            if (std::abs(m(r, col)) > std::abs(m(pivot, col))) pivot = r;
        if (m(pivot, col) == 0.0) return 0.0;
        if (pivot != col) {
            for (std::size_t j = 0; j < n; ++j) std::swap(m(pivot, j), m(col, j));
            det = -det;
        }
        det *= m(col, col);
        for (std::size_t r = col + 1; r < n; ++r) {
            const double f = m(r, col) / m(col, col);
            for (std::size_t j = col; j < n; ++j) m(r, j) -= f * m(col, j);
        }
    }
    return det;
}

double principalMinorSum(const Matrix& m, std::size_t k) {
    const std::size_t n = m.rows();
    // Enumerate the k-subsets of {0..n-1} through a selection mask.
    std::vector<bool> mask(n, false);
    std::fill(mask.begin(), mask.begin() + static_cast<long>(k), true);
    double sum = 0.0;
    do {
        std::vector<std::size_t> idx;
        for (std::size_t i = 0; i < n; ++i)
            if (mask[i]) idx.push_back(i);
        Matrix sub(k, k);
        for (std::size_t i = 0; i < k; ++i)
            for (std::size_t j = 0; j < k; ++j) sub(i, j) = m(idx[i], idx[j]);
        sum += determinant(sub);
    } while (std::prev_permutation(mask.begin(), mask.end()));
    return sum;
}

EigenSystem symmetricEigen(const Matrix& input) {
    const std::size_t n = input.rows();
    Matrix a = input;
    Matrix v = Matrix::identity(n);

    for (int sweep = 0; sweep < 100; ++sweep) {
        double off = 0.0;
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = i + 1; j < n; ++j) off += a(i, j) * a(i, j);
        if (off < 1e-30 * std::max(1.0, a.maxAbs() * a.maxAbs())) break;

        for (std::size_t p = 0; p < n; ++p)
            for (std::size_t q = p + 1; q < n; ++q) {
                if (a(p, q) == 0.0) continue;
                const double theta = (a(q, q) - a(p, p)) / (2.0 * a(p, q));
                const double t = (theta >= 0 ? 1.0 : -1.0) /
                                 (std::abs(theta) + std::sqrt(theta * theta + 1.0));
                const double c = 1.0 / std::sqrt(t * t + 1.0);
                const double s = t * c;
                for (std::size_t k = 0; k < n; ++k) { // A <- A J
                    const double akp = a(k, p), akq = a(k, q);
                    a(k, p) = c * akp - s * akq;
                    a(k, q) = s * akp + c * akq;
                }
                for (std::size_t k = 0; k < n; ++k) { // A <- J^T A
                    const double apk = a(p, k), aqk = a(q, k);
                    a(p, k) = c * apk - s * aqk;
                    a(q, k) = s * apk + c * aqk;
                }
                for (std::size_t k = 0; k < n; ++k) { // V <- V J
                    const double vkp = v(k, p), vkq = v(k, q);
                    v(k, p) = c * vkp - s * vkq;
                    v(k, q) = s * vkp + c * vkq;
                }
            }
    }

    std::vector<std::size_t> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(),
                     [&](std::size_t i, std::size_t j) { return a(i, i) > a(j, j); });

    EigenSystem es;
    es.vectors = Matrix(n, n);
    for (std::size_t c = 0; c < n; ++c) {
        es.values.push_back(a(order[c], order[c]));
        // Deterministic orientation: the largest component of each vector is positive.
        std::size_t big = 0;
        for (std::size_t r = 1; r < n; ++r)
            if (std::abs(v(r, order[c])) > std::abs(v(big, order[c])) + 1e-12) big = r;
        const double sign = v(big, order[c]) < 0 ? -1.0 : 1.0;
        for (std::size_t r = 0; r < n; ++r) es.vectors(r, c) = sign * v(r, order[c]);
    }
    // P must be a rotation (right-handed frame), as stated in the manual.
    if (n > 0 && determinant(es.vectors) < 0)
        for (std::size_t r = 0; r < n; ++r) es.vectors(r, n - 1) = -es.vectors(r, n - 1);
    return es;
}

} // namespace alga
