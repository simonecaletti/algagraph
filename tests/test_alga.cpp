// Unit tests for the ALGAgraph core. Reference values for the quadric
// x^2+2y^2-z^2-3xz+2x-6 are the ones worked out by hand in the 2015 manual.
#include "alga.hpp"
#include "plot.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace alga;

namespace {

int failures = 0;

void check(bool ok, const std::string& what) {
    if (!ok) {
        ++failures;
        std::cerr << "FAIL: " << what << "\n";
    }
}

void near(double got, double want, const std::string& what, double tol = 1e-5) {
    check(std::abs(got - want) <= tol * std::max(1.0, std::abs(want)),
          what + ": got " + std::to_string(got) + ", want " + std::to_string(want));
}

void equal(const std::string& got, const std::string& want, const std::string& what) {
    check(got == want, what + ": got \"" + got + "\", want \"" + want + "\"");
}

void typeIs(const std::string& eq, const std::string& want) {
    equal(Alga(eq).type(), want, "type of " + eq);
}

void throws(const std::string& eq) {
    try {
        Alga a(eq);
        check(false, "expected a parse error for \"" + eq + "\"");
    } catch (const std::invalid_argument&) {
    }
}

void testManualExample() {
    const Alga q("x^2+2*y^2-3*x*z-z^2+2*x-6");
    check(q.isQuadric(), "manual example is a quadric");
    const auto inv = q.invariants();
    near(inv[0], 41, "I4");
    near(inv[1], -6.5, "I3");
    near(inv[2], -3.25, "I2");
    near(inv[3], 2, "I1");
    const auto& ev = q.eigen().values;
    near(ev[0], 2, "lambda1");
    near(ev[1], std::sqrt(13.0) / 2, "lambda2");
    near(ev[2], -std::sqrt(13.0) / 2, "lambda3");
    near(determinant(q.eigen().vectors), 1.0, "det P");
    // P^T A P must be diagonal with the eigenvalues on the diagonal.
    const Matrix d = q.eigen().vectors.transposed() * q.quadraticMatrix() * q.eigen().vectors;
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j)
            near(d(i, j), i == j ? ev[i] : 0.0, "P^T A P");
    equal(q.type(), "It's a hyperboloid of one sheet (hyperbolic hyperboloid).", "manual type");
    equal(q.normalFormString(), "2*x^2+1.80278*y^2-1.80278*z^2-6.30769 = 0", "manual NF");
}

void testParser() {
    const Polynomial p = parseEquation(" 2x^2 + 3 y*x - y**2 = 4 - x + 1.5e1 ");
    near(p.coefficient(2, 0, 0), 2, "2x^2");
    near(p.coefficient(1, 1, 0), 3, "y*x -> x*y");
    near(p.coefficient(0, 2, 0), -1, "y**2");
    near(p.coefficient(1, 0, 0), 1, "-x moved to the left");
    near(p.coefficient(0, 0, 0), -19, "constant");
    near(parseEquation("x^2+x^2-x*x/4+3/2*y").coefficient(2, 0, 0), 1.75, "like terms");
    near(parseEquation("x^2+3/2*y").coefficient(0, 1, 0), 1.5, "fraction");
    check(!parseEquation("x^2+y^2-1").mentionsZ, "conic has no z");
    check(parseEquation("x^2+y^2+0*z").mentionsZ, "0*z forces a quadric");
    throws("");
    throws("x^3+y^2");
    throws("x^2*y");
    throws("x^2+w");
    throws("x^2+y^2=1=2");
    throws("x^2+");
    throws("x^2 y^2");  // whitespace is ignored, so this is x^2*y^2: degree 4
    throws("3x+2y-1");  // no quadratic part
    throws("x^2+y^2/0");
}

void testConics() {
    // From the manual and the FAQ.
    typeIs("x^2+3*y^2-x*y-5", "It is an ellipse in real points.");
    typeIs("x^2-2*y^2-3*y-4", "It is a generic hyperbola.");
    typeIs("x^2+3*x*y+2*y^2-5*x+2*y-1", "It is a generic hyperbola.");
    const Alga h("x^2-2*y^2-3*y-4");
    near(h.eigen().values[0], 1, "hyperbola lambda1");
    near(h.eigen().values[1], -2, "hyperbola lambda2");

    typeIs("x^2+y^2-1", "It is an ellipse in real points (a circle).");
    typeIs("x^2+y^2+1", "It is an ellipse in complex points (no real points).");
    typeIs("x*y-1", "It is an equilateral hyperbola.");
    typeIs("x^2-y", "It is a generic parabola.");
    typeIs("0.3x^2+0.6x*y+0.3y^2+0.1x-0.7", "It is a generic parabola.");
    typeIs("x^2-y^2", "Degenerate conic: pair of distinct real lines meeting in a point.");
    typeIs("x^2+y^2", "Degenerate conic: pair of complex conjugate lines meeting in a real point.");
    typeIs("x^2-1", "Degenerate conic: pair of distinct real parallel lines.");
    typeIs("x^2+2x*y+y^2+1",
           "Degenerate conic: pair of complex conjugate parallel lines (no real points).");
    typeIs("x^2-2x+1", "Degenerate conic: pair of coincident lines (double line).");

    // Normal forms (central: constant = I3/I2; parabola: 2*sqrt(-I3/l)).
    equal(Alga("x^2+y^2-1").normalFormString(), "x^2+y^2-1 = 0", "circle NF");
    equal(Alga("y^2-4x").normalFormString(), "x^2+4*y = 0", "parabola NF");
    const Alga e("x^2+3*y^2-x*y-5");
    const auto inv = e.invariants();
    near(e.normalForm().constant, inv[0] / inv[1], "ellipse constant = I3/I2");
    const Alga par("x^2+2x*y+y^2+3x-y");
    near(par.normalForm().linear, std::sqrt(-par.invariants()[0] / par.normalForm().squares[0]),
         "parabola linear coefficient");
}

void testQuadrics() {
    typeIs("x^2+y^2+z^2-1", "It's an ellipsoid in real points (a sphere).");
    typeIs("x^2+2y^2+3z^2-1", "It's an ellipsoid in real points.");
    typeIs("x^2+y^2+z^2+1", "It's an ellipsoid in complex points (no real points).");
    typeIs("x^2+y^2-z^2-1", "It's a hyperboloid of one sheet (hyperbolic hyperboloid).");
    typeIs("x^2+y^2-z^2+1", "It's a hyperboloid of two sheets (elliptic hyperboloid).");
    typeIs("x^2+y^2-z", "It's an elliptic paraboloid.");
    typeIs("x^2-y^2-z", "It's a hyperbolic paraboloid (saddle paraboloid).");
    typeIs("x*y-z", "It's a hyperbolic paraboloid (saddle paraboloid).");
    typeIs("x^2+y^2-z^2", "Degenerate quadric: real quadric cone.");
    typeIs("x^2+y^2+z^2", "Degenerate quadric: imaginary cone (the only real point is the vertex).");
    typeIs("x^2+y^2-1+0*z", "Degenerate quadric: elliptic cylinder in real points.");
    typeIs("x^2+y^2+1+0z", "Degenerate quadric: elliptic cylinder in complex points (no real points).");
    typeIs("x^2-z^2-1", "Degenerate quadric: hyperbolic cylinder.");
    typeIs("x^2-y+0z", "Degenerate quadric: parabolic cylinder.");
    typeIs("x^2-z^2", "Degenerate quadric: pair of distinct real planes meeting in a line.");
    typeIs("x^2+z^2",
           "Degenerate quadric: pair of complex conjugate planes meeting in a real line.");
    typeIs("z^2-1", "Degenerate quadric: pair of distinct real parallel planes.");
    typeIs("z^2+1", "Degenerate quadric: pair of complex conjugate parallel planes (no real points).");
    typeIs("x^2+y^2+z^2+2xy+2xz+2yz","Degenerate quadric: pair of coincident planes (double plane).");

    // The paraboloid normal form is linear in z (the 2015 drawNF used z^2).
    equal(Alga("x^2+y^2-z").normalFormString(), "x^2+y^2+z = 0", "paraboloid NF");
    const Alga s("x^2+2y^2-z^2-3x*z+2x+y-z-6");
    const auto inv = s.invariants();
    near(s.normalForm().constant, inv[0] / inv[1], "quadric constant = I4/I3");
}

void testDrawing() {
    // Circle of radius 1: every segment endpoint lies (almost) on the circle.
    const auto segs = conicCurve(parseEquation("x^2+y^2-1"), {-2, 2}, {-2, 2}, 200);
    check(segs.size() > 100, "circle has segments");
    for (const Segment& sg : segs)
        near(std::hypot(sg.a[0], sg.a[1]), 1.0, "point on circle", 1e-3);
    check(conicCurve(parseEquation("x^2+y^2+1"), {-2, 2}, {-2, 2}, 50).empty(),
          "imaginary ellipse draws nothing");
    const auto sphere = quadricWireframe(parseEquation("x^2+y^2+z^2-4"), {-3, 3}, {-3, 3},
                                         {-3, 3}, 10, 60);
    check(!sphere.empty(), "sphere has segments");
    for (const Segment& sg : sphere)
        near(std::sqrt(sg.a[0] * sg.a[0] + sg.a[1] * sg.a[1] + sg.a[2] * sg.a[2]), 2.0,
             "point on sphere", 1e-2);
}

} // namespace

int main() {
    testManualExample();
    testParser();
    testConics();
    testQuadrics();
    testDrawing();
    if (failures) {
        std::cerr << failures << " check(s) failed\n";
        return 1;
    }
    std::cout << "All tests passed\n";
    return 0;
}
