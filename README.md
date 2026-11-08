# ALGAgraph

Interactive program for the study of **conics** and **quadrics**.

You enter a second-degree equation such as

```
x^2+3*x*y+2*y^2-5*x+2*y-1          (conic, in x and y)
x^2+2*y^2-3*x*z-z^2+2*x-6          (quadric, in x, y and z)
```

and ALGAgraph computes the invariants, eigenvalues and eigenvectors, tells you
what kind of curve or surface it is, gives a normal (canonical) form, and draws
it both as entered and in normal form.

ALGAgraph was written in 2015 as the final project of the *Laboratorio di
Calcolo* course (a.a. 2014/2015), on top of CERN ROOT 5. Version 2.0 is a
rewrite in plain C++17. It does the same things, has no ROOT dependency and
fixes the bugs of the original (see [What changed](#what-changed-in-version-20)).
The original source is kept in [`legacy/ALGAgraph.cpp`](legacy/ALGAgraph.cpp).

## Building

Requirements: a C++17 compiler, CMake ≥ 3.14 and, for drawing,
[gnuplot](http://www.gnuplot.info/) (`sudo apt install gnuplot`,
`brew install gnuplot`, ...).

```bash
cmake -S . -B build
cmake --build build
(cd build && ctest)           # optional: run the unit tests
./build/algagraph
```

You can also pass the first equation on the command line:
`./build/algagraph "x^2+y^2-z^2-1"`.

## Writing equations

* Variables are `x`, `y` and `z`. The equation is treated as a quadric as soon
  as `z` appears (write `+0*z` to study e.g. a cylinder `x^2+y^2-1+0*z`).
* Monomials may come in any order and may repeat (like terms are summed).
  `y*x` is the same as `x*y`.
* The `*` is optional (`2xy`, `3x^2`), powers can be written `^2` or `**2`,
  and coefficients can be decimals, fractions or in scientific notation
  (`0.5x^2`, `3/2*y`, `1e-3*z`).
* Without a right-hand side the equation is `... = 0`, but you can also write
  `x^2 + y^2 = 4`.

## Commands

| Command          | Description |
|------------------|-------------|
| `help`           | List of the commands. |
| `type`           | Type of the conic/quadric (ellipse, hyperbola, parabola, ellipsoid, hyperboloids, paraboloids, and also the degenerate cases: cones, cylinders, pairs of lines/planes). |
| `invar`          | Orthogonal invariants in decreasing order (I4, I3, I2, I1). |
| `eigenVal`       | Eigenvalues of the matrix of the quadratic part. |
| `eigenVect`      | Matrix P of the normalized eigenvectors in columns (det P = +1, so P is a rotation). |
| `matrix`         | Complete matrix A' associated to the conic/quadric. |
| `normalForm`     | A normal form of the conic/quadric. |
| `draw [file]`    | Draws the conic/quadric in a range chosen by you (press Enter for the default range [-10, 10]). The gnuplot window can be zoomed, and rotated for quadrics, with the mouse. With a file name (`.png`, `.svg`, `.pdf`) the drawing is saved instead. |
| `drawNF [file]`  | Same as `draw`, for the normal form. |
| `new [eq]`       | Enter a new equation; the previous ones stay in memory. |
| `list`           | List of all the equations entered (`*` marks the current one). |
| `select [n]`     | Select equation number `n` of the list. |
| `clear`          | Delete all the equations from memory. |
| `.q`             | Quit (also `quit`, `exit` or Ctrl-D). |

Command names are case-insensitive.

## Example

This is the quadric worked out by hand in the original manual:

```
ALGAgraph[0]->Enter the eq. :x^2+2*y^2-3*x*z-z^2+2*x-6
ALGAgraph[1]->invar
Invariants:
I4: 41
I3: -6.5
I2: -3.25
I1: 2
ALGAgraph[2]->eigenVal
EigenValues:
2
1.80278
-1.80278
ALGAgraph[3]->type
It's a hyperboloid of one sheet (hyperbolic hyperboloid).
ALGAgraph[4]->normalForm
The equation in normal form is:
2*x^2+1.80278*y^2-1.80278*z^2-6.30769 = 0
ALGAgraph[5]->draw
Define the range in which to graph the quadric: (xmin xmax ymin ymax zmin zmax) [default -10 10 ...]: -20 20 -20 20 -20 20
```

As the manual's FAQ already said, a plot in a small range can be misleading:
`x^2+3*x*y+2*y^2-5*x+2*y-1` looks like a parabola between -5 and 5, but `type`
tells you it is a hyperbola, which you can see by drawing it between -200 and 200.

## How it works

* The equation is parsed into its coefficients and written as
  `X^T A' X = 0` with `X = (x, y[, z], 1)`. `A` is the block of `A'` that
  belongs to the quadratic part.
* The invariants are `det A'` and the sums of the principal minors of `A`
  (`I1 = tr A`, ..., `det A`).
* Eigenvalues and eigenvectors of `A` come from the Jacobi method. The
  eigenvalues are sorted in decreasing order and P is oriented so that it is
  a rotation.
* The normal form follows the rotation `x = P x'`. Completing the squares
  gives `Σ λᵢ Xᵢ² + c = 0` (with `c = I4/I3`, or `I3/I2` for conics, in the
  non-degenerate case) or `Σ λᵢ Xᵢ² + 2m Y = 0` for parabolas, paraboloids and
  parabolic cylinders. The type is read off the signs of the λᵢ and of `c`.
* Drawing extracts the zero set with *marching squares*. A conic is a set of
  short segments; a quadric is drawn as a wireframe of its sections with planes
  parallel to the coordinate planes. gnuplot then renders the result.

## What changed in version 2.0

* **No ROOT.** The linear algebra (`TMatrixDSym`, `TMatrixDSymEigen`) is
  replaced by a small built-in module, and drawing (`TF2`/`TF3`) by gnuplot.
* **Parser rewritten.** It now accepts spaces, `y*x`, `2x`, `x**2`, repeated
  monomials, decimals and fractions, and an `=` sign. Invalid input gets an
  error message; before, it was silently misread.
* **Drawing fixed.** The contour level passed to `SetContour` was never set
  (uninitialized array). The normal form of a paraboloid was drawn as
  `...+2*z^2` instead of `...+2*z`.
* **Normal form fixed.** For parabolas with the zero eigenvalue in the first
  position, two branches ran and the second divided by I2 = 0. Every case,
  degenerate ones included, now has a normal form.
* **Classification.** Zero tests use a tolerance instead of exact `== 0`
  comparisons on floating-point numbers, so non-integer coefficients work.
  Degenerate quadrics are classified too; before, `type` printed nothing.
* **Session commands fixed.** `select` used the wrong index and accepted
  invalid numbers, `clear` led to out-of-range accesses, and `new` leaked
  memory.
* Drawing no longer blocks the prompt. Plots can be saved to file. New
  `matrix` command. End of input (Ctrl-D) is handled.

## Project layout

```
src/linalg.*     small matrices, determinant, Jacobi eigenvalues
src/alga.*       equation parser, invariants, normal form, classification
src/plot.*       marching squares and gnuplot output
src/main.cpp     interactive command loop
tests/           unit tests (run with ctest)
legacy/          the original 2015 ROOT-based source
```

## License

MIT, see [LICENSE](LICENSE).
