// ALGAgraph: interactive study of conics and quadrics.
//
// Originally written in 2015 for the "Laboratorio di Calcolo" course on top of
// CERN ROOT; this version is plain C++17 and draws with gnuplot.
#include "alga.hpp"
#include "plot.hpp"

#include <algorithm>
#include <cctype>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace alga;

namespace {

const char* const kVersion = "2.0";
const char* const kDefaultEquation = "x^2+y^2-1";

std::string trim(const std::string& s) {
    const auto b = s.find_first_not_of(" \t\r\n");
    if (b == std::string::npos) return "";
    return s.substr(b, s.find_last_not_of(" \t\r\n") - b + 1);
}

std::string lower(std::string s) {
    for (char& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return s;
}

bool readLine(const std::string& prompt, std::string& line) {
    std::cout << prompt << std::flush;
    if (!std::getline(std::cin, line)) {
        std::cout << '\n';
        return false;
    }
    line = trim(line);
    return true;
}

void printHelp() {
    std::cout <<
        "Enter a second-degree equation in x, y (conic) or x, y, z (quadric), e.g.\n"
        "  x^2+3*x*y+2*y^2-5*x+2*y-1      or      x^2 + 2y^2 - z^2 - 3xz + 2x = 6\n"
        "Monomials can be written in any order (y*x is the same as x*y); a missing\n"
        "right-hand side means '= 0'. The equation is a quadric as soon as z appears.\n"
        "\n"
        "->type:\n  Returns the type of the conic/quadric (ellipse, hyperbola, elliptic\n"
        "  paraboloid, cone, cylinder, ...), degenerate cases included.\n"
        "->eigenVal:\n  Returns the eigenvalues of the matrix of the quadratic part.\n"
        "->eigenVect:\n  Returns the matrix P of the normalized eigenvectors put in columns\n"
        "  (det P = +1, so P is a rotation).\n"
        "->invar:\n  Returns the orthogonal invariants in decreasing order (I4, I3, I2, I1).\n"
        "->matrix:\n  Returns the complete matrix A' associated to the conic/quadric.\n"
        "->normalForm:\n  Returns a normal (canonical) form of the conic/quadric.\n"
        "->draw [file]:\n  Draws the conic/quadric in a range chosen by the user. The window can\n"
        "  be zoomed (and rotated, for quadrics) with the mouse. With a file name\n"
        "  (.png, .svg or .pdf) the drawing is saved instead.\n"
        "->drawNF [file]:\n  Like draw, but for the conic/quadric in normal form.\n"
        "->new [eq]:\n  Enter a new equation to study; the previous ones stay in memory.\n"
        "->list:\n  Returns the ordered list of all the equations entered so far.\n"
        "->select [n]:\n  Selects equation number n of the list as the current one.\n"
        "->clear:\n  Deletes all the equations from memory (the current one becomes "
     << kDefaultEquation << ").\n"
        "->.q:\n  Quits ALGAgraph.\n";
}

void printMatrix(const Matrix& m) {
    std::vector<std::string> cells;
    std::size_t width = 0;
    for (std::size_t i = 0; i < m.rows(); ++i)
        for (std::size_t j = 0; j < m.cols(); ++j) {
            cells.push_back(formatNumber(m(i, j)));
            width = std::max(width, cells.back().size());
        }
    for (std::size_t i = 0; i < m.rows(); ++i) {
        for (std::size_t j = 0; j < m.cols(); ++j) {
            const std::string& c = cells[i * m.cols() + j];
            std::cout << std::string(width - c.size() + (j ? 2 : 0), ' ') << c;
        }
        std::cout << '\n';
    }
}

// Asks for the drawing range; an empty answer keeps the default [-10, 10].
bool askRanges(bool quadric, std::array<Range, 3>& ranges) {
    const int n = quadric ? 6 : 4;
    std::string line;
    if (!readLine(quadric ? "Define the range in which to graph the quadric: "
                            "(xmin xmax ymin ymax zmin zmax) [default -10 10 ...]: "
                          : "Define the range in which to graph the conic: "
                            "(xmin xmax ymin ymax) [default -10 10 -10 10]: ",
                  line))
        return false;
    for (char& c : line)
        if (c == ',') c = ' ';
    std::istringstream is(line);
    std::vector<double> v;
    double d;
    while (is >> d) v.push_back(d);
    if (line.empty()) return true;
    if (!is.eof() || static_cast<int>(v.size()) != n) {
        std::cout << "Please enter " << n << " numbers.\n";
        return false;
    }
    for (int i = 0; i < n / 2; ++i) {
        if (!(v[2 * i] < v[2 * i + 1])) {
            std::cout << "Each range must satisfy min < max.\n";
            return false;
        }
        ranges[i] = {v[2 * i], v[2 * i + 1]};
    }
    return true;
}

void drawCommand(const Alga& q, bool normalForm, const std::string& output) {
    try {
        checkOutputFormat(output);
    } catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << "\n";
        return;
    }
    std::array<Range, 3> ranges;
    if (!askRanges(q.isQuadric(), ranges)) return;
    const std::string title = normalForm ? "Normal form: " + q.normalFormString()
                                         : q.equation();
    const Polynomial p = normalForm ? q.normalFormPolynomial() : q.polynomial();
    try {
        if (draw(p, q.isQuadric(), title, ranges, output) == 0) {
            std::cout << "Nothing to draw: there are no real points in the selected range.\n"
                         "Try a larger range, and check the 'type': the figure may have no "
                         "real points\nor reduce to a single point, a double line or a double "
                         "plane.\n";
            return;
        }
        if (!output.empty()) std::cout << "Saved to " << output << "\n";
    } catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << "\n";
    }
}

struct Session {
    std::vector<Alga> equations;
    std::size_t current = 0;

    const Alga& q() const { return equations[current]; }

    void add(const Alga& a) {
        equations.push_back(a);
        current = equations.size() - 1;
    }
};

// Reads equations until one parses; false on end of input.
bool askEquation(const std::string& prompt, Session& s, std::string text = "") {
    while (true) {
        if (text.empty() && !readLine(prompt, text)) return false;
        try {
            s.add(Alga(text));
            return true;
        } catch (const std::invalid_argument& e) {
            std::cout << "Invalid equation: " << e.what() << "\n";
            text.clear();
        }
    }
}

} // namespace

int main(int argc, char** argv) {
    std::string first;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "-h" || a == "--help") {
            std::cout << "Usage: " << argv[0] << " [equation]\n\n";
            printHelp();
            return 0;
        }
        if (a == "-v" || a == "--version") {
            std::cout << "ALGAgraph " << kVersion << "\n";
            return 0;
        }
        first += (first.empty() ? "" : " ") + a;
    }

    std::cout << "\n\t~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
              << "\t\t\tALGAgraph ver " << kVersion << ":\n"
              << "\twrite 'help' after entering the first equation for more info.\n"
              << "\t~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n";

    Session s;
    if (!askEquation("ALGAgraph[0]->Enter the eq. :", s, first)) return 0;

    for (int times = 1;; ++times) {
        std::string line;
        if (!readLine("ALGAgraph[" + std::to_string(times) + "]->", line)) break;
        if (line.empty()) {
            --times;
            continue;
        }
        std::istringstream is(line);
        std::string command;
        is >> command;
        const std::string arg = trim(line.substr(command.size()));
        command = lower(command);
        const Alga& q = s.q();

        if (command == "eigenval") {
            std::cout << "EigenValues:\n";
            for (double v : q.eigen().values) std::cout << formatNumber(v) << "\n";
        } else if (command == "eigenvect") {
            std::cout << "EigenVectors:\n";
            printMatrix(q.eigen().vectors);
        } else if (command == "invar") {
            const std::vector<double> inv = q.invariants();
            std::cout << "Invariants:\n";
            for (std::size_t i = 0; i < inv.size(); ++i)
                std::cout << "I" << inv.size() - i << ": " << formatNumber(inv[i]) << "\n";
        } else if (command == "matrix") {
            std::cout << "Associated matrix A':\n";
            printMatrix(q.matrix());
        } else if (command == "normalform") {
            std::cout << "The equation in normal form is:\n" << q.normalFormString() << "\n";
        } else if (command == "draw") {
            drawCommand(q, false, arg);
        } else if (command == "drawnf") {
            drawCommand(q, true, arg);
        } else if (command == "type") {
            std::cout << q.type() << "\n";
        } else if (command == "new") {
            if (!askEquation("Enter the eq. : ", s, arg)) break;
        } else if (command == "list") {
            for (std::size_t i = 0; i < s.equations.size(); ++i)
                std::cout << i + 1 << (i == s.current ? "*-> " : " -> ")
                          << s.equations[i].equation() << "\n";
        } else if (command == "select") {
            std::string which = arg;
            if (which.empty() && !readLine("Equation number: ", which)) break;
            std::size_t pos = 0;
            long n = 0;
            try {
                n = std::stol(which, &pos);
            } catch (const std::exception&) {
                pos = 0;
            }
            if (pos == 0 || pos != which.size() || n < 1 ||
                n > static_cast<long>(s.equations.size())) {
                std::cout << "The selected equation does not exist (choose 1-"
                          << s.equations.size() << ").\n";
            } else {
                s.current = static_cast<std::size_t>(n - 1);
                std::cout << "Selected: " << s.q().equation() << "\n";
            }
        } else if (command == "clear") {
            s.equations.clear();
            s.add(Alga(kDefaultEquation));
            std::cout << "Memory cleared. The current equation is now " << kDefaultEquation
                      << " (use 'new' to enter another one).\n";
        } else if (command == "help") {
            printHelp();
        } else if (command == ".q" || command == "quit" || command == "exit") {
            break;
        } else {
            std::cout << "Command not recognized (write 'help' for the list of commands).\n";
        }
    }
    std::cout << "Thanks for using ALGAgraph.\n";
    cleanupPlots();
    return 0;
}
