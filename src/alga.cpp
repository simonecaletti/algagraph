#include "alga.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <sstream>
#include <stdexcept>

namespace alga {

// ---------------------------------------------------------------- Polynomial

double Polynomial::coefficient(int ex, int ey, int ez) const {
    auto it = terms.find({ex, ey, ez});
    return it == terms.end() ? 0.0 : it->second;
}

double Polynomial::evaluate(double x, double y, double z) const {
    double sum = 0.0;
    for (const auto& [m, c] : terms) {
        double v = c;
        for (int i = 0; i < m[0]; ++i) v *= x;
        for (int i = 0; i < m[1]; ++i) v *= y;
        for (int i = 0; i < m[2]; ++i) v *= z;
        sum += v;
    }
    return sum;
}

// ---------------------------------------------------------------- Parser

namespace {

class Parser {
public:
    explicit Parser(std::string s) : m_s(std::move(s)) {}

    // side := term (('+'|'-') term)*   with an optional leading sign
    void parseSide(Polynomial& out, double sideSign) {
        if (atEnd()) fail("missing expression");
        bool first = true;
        while (!atEnd()) {
            double sign = 1.0;
            bool hasSign = false;
            while (!atEnd() && (peek() == '+' || peek() == '-')) {
                if (get() == '-') sign = -sign;
                hasSign = true;
            }
            if (!first && !hasSign) fail("expected '+' or '-'");
            first = false;
            parseTerm(out, sign * sideSign);
        }
    }

private:
    // term := factor (['*'] factor | '/' number)*
    void parseTerm(Polynomial& out, double sign) {
        if (atEnd()) fail("missing term after sign");
        double coeff = sign;
        Monomial exps{0, 0, 0};
        bool any = false;
        while (!atEnd() && peek() != '+' && peek() != '-') {
            if (any) {
                if (peek() == '*') {
                    get();
                    if (atEnd()) fail("missing factor after '*'");
                } else if (peek() == '/') {
                    get();
                    if (atEnd() || !startsNumber()) fail("expected a number after '/'");
                    const double d = parseNumber();
                    if (d == 0.0) fail("division by zero");
                    coeff /= d;
                    continue;
                }
            }
            if (startsNumber()) {
                coeff *= parseNumber();
            } else if (peek() == 'x' || peek() == 'y' || peek() == 'z') {
                const char var = get();
                if (var == 'z') out.mentionsZ = true;
                int power = 1;
                if (!atEnd() && (peek() == '^' || m_s.compare(m_pos, 2, "**") == 0)) {
                    m_pos += peek() == '^' ? 1 : 2;
                    if (atEnd() || !std::isdigit(static_cast<unsigned char>(peek())))
                        fail("expected an integer exponent");
                    power = 0;
                    while (!atEnd() && std::isdigit(static_cast<unsigned char>(peek())))
                        power = power * 10 + (get() - '0');
                }
                exps[static_cast<std::size_t>(var - 'x')] += power;
            } else {
                fail(std::string("unexpected character '") + peek() + "'");
            }
            any = true;
        }
        if (!any) fail("missing term");
        if (exps[0] + exps[1] + exps[2] > 2)
            fail("terms of degree higher than 2 are not allowed");
        out.terms[exps] += coeff;
    }

    bool startsNumber() const {
        return std::isdigit(static_cast<unsigned char>(peek())) ||
               (peek() == '.' && m_pos + 1 < m_s.size() &&
                std::isdigit(static_cast<unsigned char>(m_s[m_pos + 1])));
    }

    double parseNumber() {
        const std::size_t start = m_pos;
        while (!atEnd() && (std::isdigit(static_cast<unsigned char>(peek())) || peek() == '.'))
            ++m_pos;
        // Scientific notation: only when 'e' is followed by an exponent.
        if (!atEnd() && (peek() == 'e' || peek() == 'E')) {
            std::size_t p = m_pos + 1;
            if (p < m_s.size() && (m_s[p] == '+' || m_s[p] == '-')) ++p;
            if (p < m_s.size() && std::isdigit(static_cast<unsigned char>(m_s[p]))) {
                m_pos = p;
                while (!atEnd() && std::isdigit(static_cast<unsigned char>(peek()))) ++m_pos;
            }
        }
        const std::string tok = m_s.substr(start, m_pos - start);
        std::size_t used = 0;
        double v = 0.0;
        try {
            v = std::stod(tok, &used);
        } catch (const std::exception&) {
            used = 0;
        }
        if (used != tok.size()) {
            m_pos = start;
            fail("malformed number '" + tok + "'");
        }
        return v;
    }

    bool atEnd() const { return m_pos >= m_s.size(); }
    char peek() const { return atEnd() ? '\0' : m_s[m_pos]; }
    char get() { return m_s[m_pos++]; }

    [[noreturn]] void fail(const std::string& what) const {
        throw std::invalid_argument(what + " (at position " + std::to_string(m_pos + 1) +
                                    " of \"" + m_s + "\")");
    }

    std::string m_s;
    std::size_t m_pos = 0;
};

} // namespace

Polynomial parseEquation(const std::string& text) {
    std::string s;
    for (char c : text)
        if (!std::isspace(static_cast<unsigned char>(c)))
            s += static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    if (s.empty()) throw std::invalid_argument("empty equation");

    const std::size_t eq = s.find('=');
    if (eq != std::string::npos && s.find('=', eq + 1) != std::string::npos)
        throw std::invalid_argument("more than one '=' in the equation");

    Polynomial p;
    Parser(s.substr(0, eq)).parseSide(p, 1.0);
    if (eq != std::string::npos) Parser(s.substr(eq + 1)).parseSide(p, -1.0);
    return p;
}

// ---------------------------------------------------------------- Formatting

std::string formatNumber(double value) {
    std::ostringstream os;
    os.precision(6);
    os << value;
    std::string s = os.str();
    if (s == "-0") s = "0";
    return s;
}

std::string formatPolynomial(const std::vector<std::pair<double, std::string>>& terms,
                             double constant) {
    std::string out;
    auto append = [&out](double c, const std::string& name) {
        std::string num = formatNumber(c);
        if (num == "0") return;
        const bool negative = num[0] == '-';
        if (negative) num.erase(0, 1);
        if (negative)
            out += '-';
        else if (!out.empty())
            out += '+';
        if (name.empty())
            out += num;
        else if (num == "1")
            out += name;
        else
            out += num + "*" + name;
    };
    for (const auto& [c, name] : terms) append(c, name);
    append(constant, "");
    return out.empty() ? "0" : out;
}

// ---------------------------------------------------------------- Alga

Alga::Alga(const std::string& equation) : m_equation(equation), m_poly(parseEquation(equation)) {
    // As in the original program: the equation is a quadric as soon as z appears.
    m_dim = m_poly.mentionsZ ? 3 : 2;
    const std::size_t n = static_cast<std::size_t>(m_dim);

    // f = X^T A' X with X = (x, y[, z], 1)
    m_full = Matrix(n + 1, n + 1);
    for (std::size_t i = 0; i < n; ++i) {
        Monomial sq{0, 0, 0};
        sq[i] = 2;
        m_full(i, i) = m_poly.coefficient(sq[0], sq[1], sq[2]);
        for (std::size_t j = i + 1; j < n; ++j) {
            Monomial mixed{0, 0, 0};
            mixed[i] = mixed[j] = 1;
            m_full(i, j) = m_full(j, i) = m_poly.coefficient(mixed[0], mixed[1], mixed[2]) / 2.0;
        }
        Monomial lin{0, 0, 0};
        lin[i] = 1;
        m_full(i, n) = m_full(n, i) = m_poly.coefficient(lin[0], lin[1], lin[2]) / 2.0;
    }
    m_full(n, n) = m_poly.coefficient(0, 0, 0);

    m_quad = Matrix(n, n);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) m_quad(i, j) = m_full(i, j);

    if (m_quad.maxAbs() == 0.0)
        throw std::invalid_argument("\"" + equation +
                                    "\" is not a second-degree equation (no quadratic terms)");

    m_tol = 1e-9 * m_full.maxAbs();
    m_eigen = symmetricEigen(m_quad);
    for (double& v : m_eigen.values) v = chop(v);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j)
            if (std::abs(m_eigen.vectors(i, j)) < 1e-12) m_eigen.vectors(i, j) = 0.0;
}

double Alga::chop(double value, int degree) const {
    const double scale = std::pow(m_full.maxAbs(), degree);
    return std::abs(value) <= 1e-9 * scale ? 0.0 : value;
}

std::vector<double> Alga::invariants() const {
    std::vector<double> inv;
    inv.push_back(chop(determinant(m_full), m_dim + 1));
    for (int k = m_dim; k >= 1; --k)
        inv.push_back(chop(principalMinorSum(m_quad, static_cast<std::size_t>(k)), k));
    return inv;
}

NormalForm Alga::normalForm() const {
    const std::size_t n = static_cast<std::size_t>(m_dim);
    const Matrix& P = m_eigen.vectors;

    // Rotation x = P x': f = sum l_i x'_i^2 + 2 b'.x' + c with b' = P^T b.
    std::vector<double> b(n), bp(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) b[i] = m_full(i, n);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t k = 0; k < n; ++k) bp[i] += P(k, i) * b[k];

    NormalForm nf;
    double constant = m_full(n, n);
    double constScale = std::abs(constant);
    double linear2 = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const double l = m_eigen.values[i];
        if (l != 0.0) {
            // Translation completes the square along x'_i.
            nf.squares.push_back(l);
            constant -= bp[i] * bp[i] / l;
            constScale = std::max(constScale, std::abs(bp[i] * bp[i] / l));
        } else {
            linear2 += bp[i] * bp[i];
        }
    }
    const double linear = std::sqrt(linear2);
    if (linear > m_tol) {
        // A rotation inside the null space of A collects the linear terms on
        // one axis, and a translation along it absorbs the constant.
        nf.hasLinear = true;
        nf.linear = linear;
    } else {
        constScale = std::max(constScale, m_full.maxAbs());
        nf.constant = std::abs(constant) <= 1e-9 * constScale ? 0.0 : constant;
    }
    return nf;
}

namespace {
const char* const kSquareNames[] = {"x^2", "y^2", "z^2"};
const char* const kLinearNames[] = {"x", "y", "z"};
} // namespace

std::string Alga::normalFormString() const {
    const NormalForm nf = normalForm();
    std::vector<std::pair<double, std::string>> terms;
    for (std::size_t i = 0; i < nf.squares.size(); ++i)
        terms.emplace_back(nf.squares[i], kSquareNames[i]);
    if (nf.hasLinear) terms.emplace_back(2.0 * nf.linear, kLinearNames[nf.squares.size()]);
    return formatPolynomial(terms, nf.hasLinear ? 0.0 : nf.constant) + " = 0";
}

Polynomial Alga::normalFormPolynomial() const {
    const NormalForm nf = normalForm();
    Polynomial p;
    p.mentionsZ = isQuadric();
    for (std::size_t i = 0; i < nf.squares.size(); ++i) {
        Monomial m{0, 0, 0};
        m[i] = 2;
        p.terms[m] = nf.squares[i];
    }
    if (nf.hasLinear) {
        Monomial m{0, 0, 0};
        m[nf.squares.size()] = 1;
        p.terms[m] = 2.0 * nf.linear;
    } else {
        p.terms[{0, 0, 0}] = nf.constant;
    }
    return p;
}

std::string Alga::type() const {
    const NormalForm nf = normalForm();
    const std::size_t rank = nf.squares.size();
    int positive = 0;
    for (double l : nf.squares)
        if (l > 0) ++positive;
    const bool sameSign = positive == 0 || positive == static_cast<int>(rank);
    // Number of square coefficients with the sign of -constant, i.e. on the
    // "right" side of  sum l_i X_i^2 = -constant.
    int s = 0;
    for (double l : nf.squares)
        if ((l > 0) == (nf.constant < 0)) ++s;
    const bool allEqual =
        std::all_of(nf.squares.begin(), nf.squares.end(),
                    [&](double l) { return std::abs(l - nf.squares[0]) <= 1e-9 * std::abs(l); });

    if (!isQuadric()) {
        if (rank == 2) {
            if (nf.constant != 0.0) {
                if (!sameSign)
                    return std::abs(nf.squares[0] + nf.squares[1]) <= m_tol
                               ? "It is an equilateral hyperbola."
                               : "It is a generic hyperbola.";
                if (s == 2)
                    return allEqual ? "It is an ellipse in real points (a circle)."
                                    : "It is an ellipse in real points.";
                return "It is an ellipse in complex points (no real points).";
            }
            return sameSign ? "Degenerate conic: pair of complex conjugate lines meeting in a "
                              "real point."
                            : "Degenerate conic: pair of distinct real lines meeting in a point.";
        }
        if (nf.hasLinear) return "It is a generic parabola.";
        if (nf.constant == 0.0) return "Degenerate conic: pair of coincident lines (double line).";
        return s == 1 ? "Degenerate conic: pair of distinct real parallel lines."
                      : "Degenerate conic: pair of complex conjugate parallel lines (no real "
                        "points).";
    }

    if (rank == 3) {
        if (nf.constant == 0.0)
            return sameSign ? "Degenerate quadric: imaginary cone (the only real point is the "
                              "vertex)."
                            : "Degenerate quadric: real quadric cone.";
        switch (s) {
        case 3: return allEqual ? "It's an ellipsoid in real points (a sphere)."
                                : "It's an ellipsoid in real points.";
        case 0: return "It's an ellipsoid in complex points (no real points).";
        case 2: return "It's a hyperboloid of one sheet (hyperbolic hyperboloid).";
        default: return "It's a hyperboloid of two sheets (elliptic hyperboloid).";
        }
    }
    if (rank == 2) {
        if (nf.hasLinear)
            return sameSign ? "It's an elliptic paraboloid."
                            : "It's a hyperbolic paraboloid (saddle paraboloid).";
        if (nf.constant == 0.0)
            return sameSign ? "Degenerate quadric: pair of complex conjugate planes meeting in "
                              "a real line."
                            : "Degenerate quadric: pair of distinct real planes meeting in a "
                              "line.";
        if (!sameSign) return "Degenerate quadric: hyperbolic cylinder.";
        return s == 2 ? "Degenerate quadric: elliptic cylinder in real points."
                      : "Degenerate quadric: elliptic cylinder in complex points (no real "
                        "points).";
    }
    if (nf.hasLinear) return "Degenerate quadric: parabolic cylinder.";
    if (nf.constant == 0.0) return "Degenerate quadric: pair of coincident planes (double plane).";
    return s == 1 ? "Degenerate quadric: pair of distinct real parallel planes."
                  : "Degenerate quadric: pair of complex conjugate parallel planes (no real "
                    "points).";
}

} // namespace alga
