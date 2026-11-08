// Drawing of conics and quadrics with gnuplot (replaces ROOT's TF2/TF3).
//
// The zero set of the polynomial is extracted with marching squares: a
// conic becomes a set of line segments, a quadric a wireframe made of its
// sections with planes parallel to the coordinate planes.
#pragma once

#include "alga.hpp"

#include <array>
#include <cstddef>
#include <string>
#include <vector>

namespace alga {

struct Range {
    double min = -10.0, max = 10.0;
};

struct Segment {
    std::array<double, 3> a, b;
};

// Segments approximating f(x, y) = 0 on an n x n grid (z = 0).
std::vector<Segment> conicCurve(const Polynomial& p, Range x, Range y, int n = 800);

// Sections of f(x, y, z) = 0 with `slices` planes per axis, each traced on a
// grid x grid mesh.
std::vector<Segment> quadricWireframe(const Polynomial& p, Range x, Range y, Range z,
                                      int slices = 24, int grid = 100);

// Draws the zero set of p. `ranges` holds x, y (and z for a quadric).
// An empty `output` opens an interactive gnuplot window (the program does not
// wait for it); otherwise the plot is written to a .png, .svg or .pdf file.
// Returns the number of drawn segments (0 means: no real points in range,
// nothing was drawn). Throws std::runtime_error if gnuplot fails.
std::size_t draw(const Polynomial& p, bool quadric, const std::string& title,
                 const std::array<Range, 3>& ranges, const std::string& output);

// Throws std::runtime_error if `output` is not empty and not a .png, .svg or .pdf file.
void checkOutputFormat(const std::string& output);

// Removes the temporary gnuplot scripts.
void cleanupPlots();

} // namespace alga
