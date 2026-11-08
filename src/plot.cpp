#include "plot.hpp"

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <unistd.h>

namespace alga {

namespace {

// Marching squares on [u0,u1]x[v0,v1]; emits segments in (u, v) coordinates.
template <class F, class Emit>
void marchingSquares(F f, Range u, Range v, int n, Emit emit) {
    const double du = (u.max - u.min) / n, dv = (v.max - v.min) / n;
    std::vector<double> prev(static_cast<std::size_t>(n) + 1), cur(prev.size());
    for (int i = 0; i <= n; ++i) prev[static_cast<std::size_t>(i)] = f(u.min + i * du, v.min);

    for (int j = 0; j < n; ++j) {
        const double va = v.min + j * dv, vb = va + dv;
        for (int i = 0; i <= n; ++i) cur[static_cast<std::size_t>(i)] = f(u.min + i * du, vb);
        for (int i = 0; i < n; ++i) {
            const double ua = u.min + i * du, ub = ua + du;
            // Corners counter-clockwise from the bottom-left one.
            const double c[4] = {prev[static_cast<std::size_t>(i)],
                                 prev[static_cast<std::size_t>(i) + 1],
                                 cur[static_cast<std::size_t>(i) + 1],
                                 cur[static_cast<std::size_t>(i)]};
            const double pu[4] = {ua, ub, ub, ua}, pv[4] = {va, va, vb, vb};
            const bool in[4] = {c[0] > 0, c[1] > 0, c[2] > 0, c[3] > 0};

            std::array<double, 2> cross[4];
            bool has[4] = {false, false, false, false};
            int count = 0;
            for (int e = 0; e < 4; ++e) { // edge e joins corner e and e+1
                const int k = (e + 1) % 4;
                if (in[e] == in[k]) continue;
                const double t = c[e] / (c[e] - c[k]);
                cross[e] = {pu[e] + t * (pu[k] - pu[e]), pv[e] + t * (pv[k] - pv[e])};
                has[e] = true;
                ++count;
            }
            if (count == 2) {
                int e1 = -1, e2 = -1;
                for (int e = 0; e < 4; ++e)
                    if (has[e]) (e1 < 0 ? e1 : e2) = e;
                emit(cross[e1], cross[e2]);
            } else if (count == 4) { // saddle: decide with the centre value
                const bool centre = f((ua + ub) / 2, (va + vb) / 2) > 0;
                if (centre == in[0]) { // corners 1 and 3 are isolated
                    emit(cross[0], cross[1]);
                    emit(cross[2], cross[3]);
                } else { // corners 0 and 2 are isolated
                    emit(cross[3], cross[0]);
                    emit(cross[1], cross[2]);
                }
            }
        }
        prev.swap(cur);
    }
}

std::filesystem::path tempDir() {
    static const std::filesystem::path dir =
        std::filesystem::temp_directory_path() / ("algagraph-" + std::to_string(getpid()));
    std::filesystem::create_directories(dir);
    return dir;
}

std::string shellQuote(const std::string& s) {
    std::string q = "'";
    for (char c : s) q += c == '\'' ? std::string("'\\''") : std::string(1, c);
    return q + "'";
}

std::string gnuplotQuote(const std::string& s) {
    std::string q = "'";
    for (char c : s) q += c == '\'' ? std::string("''") : std::string(1, c);
    return q + "'";
}

std::string terminalFor(const std::string& output) {
    const std::string ext = std::filesystem::path(output).extension().string();
    if (ext == ".png") return "pngcairo enhanced size 1000,800";
    if (ext == ".svg") return "svg enhanced size 1000,800";
    if (ext == ".pdf") return "pdfcairo enhanced size 7in,5.6in";
    throw std::runtime_error("unsupported output format '" + ext + "' (use .png, .svg or .pdf)");
}

} // namespace

std::vector<Segment> conicCurve(const Polynomial& p, Range x, Range y, int n) {
    std::vector<Segment> out;
    marchingSquares([&](double u, double v) { return p.evaluate(u, v); }, x, y, n,
                    [&](const std::array<double, 2>& a, const std::array<double, 2>& b) {
                        out.push_back({{a[0], a[1], 0.0}, {b[0], b[1], 0.0}});
                    });
    return out;
}

std::vector<Segment> quadricWireframe(const Polynomial& p, Range x, Range y, Range z, int slices,
                                      int grid) {
    std::vector<Segment> out;
    const Range r[3] = {x, y, z};
    for (int axis = 0; axis < 3; ++axis) {
        // The section plane is orthogonal to `axis`; (u, v) are the other two.
        const int ua = (axis + 1) % 3, va = (axis + 2) % 3;
        for (int k = 0; k <= slices; ++k) {
            const double w = r[axis].min + k * (r[axis].max - r[axis].min) / slices;
            auto point = [&](double u, double v) {
                std::array<double, 3> pt{};
                pt[static_cast<std::size_t>(axis)] = w;
                pt[static_cast<std::size_t>(ua)] = u;
                pt[static_cast<std::size_t>(va)] = v;
                return pt;
            };
            marchingSquares(
                [&](double u, double v) {
                    const auto pt = point(u, v);
                    return p.evaluate(pt[0], pt[1], pt[2]);
                },
                r[ua], r[va], grid,
                [&](const std::array<double, 2>& a, const std::array<double, 2>& b) {
                    out.push_back({point(a[0], a[1]), point(b[0], b[1])});
                });
        }
    }
    return out;
}

std::size_t draw(const Polynomial& p, bool quadric, const std::string& title,
                 const std::array<Range, 3>& ranges, const std::string& output) {
    const std::vector<Segment> segments =
        quadric ? quadricWireframe(p, ranges[0], ranges[1], ranges[2])
                : conicCurve(p, ranges[0], ranges[1]);
    if (segments.empty()) return 0;

    const bool interactive = output.empty();
    if (interactive && !std::getenv("DISPLAY") && !std::getenv("WAYLAND_DISPLAY"))
        throw std::runtime_error("no graphical display available; save the plot to a file "
                                 "instead, e.g. 'draw plot.png'");
    if (std::system("command -v gnuplot > /dev/null 2>&1") != 0)
        throw std::runtime_error("gnuplot not found; please install it (e.g. 'sudo apt install "
                                 "gnuplot')");

    static int counter = 0;
    const std::filesystem::path script = tempDir() / ("plot" + std::to_string(++counter) + ".gp");
    std::ofstream gp(script);
    gp.precision(10);
    if (!interactive) {
        gp << "set terminal " << terminalFor(output) << "\n";
        gp << "set output " << gnuplotQuote(output) << "\n";
    }
    gp << "set title " << gnuplotQuote(title) << "\n"
       << "unset key\nset grid\n"
       << "set xrange [" << ranges[0].min << ":" << ranges[0].max << "]\n"
       << "set yrange [" << ranges[1].min << ":" << ranges[1].max << "]\n"
       << "set xlabel 'x'\nset ylabel 'y'\n";
    if (quadric) {
        gp << "set zrange [" << ranges[2].min << ":" << ranges[2].max << "]\n"
           << "set zlabel 'z'\nset xyplane relative 0\nset view 60,30\n";
    } else {
        gp << "set xzeroaxis lt -1\nset yzeroaxis lt -1\n";
        // Same unit on both axes, unless that would squash the picture.
        const double aspect =
            (ranges[1].max - ranges[1].min) / (ranges[0].max - ranges[0].min);
        if (aspect > 1.0 / 3 && aspect < 3.0) gp << "set size ratio -1\n";
    }
    gp << "$zeroset << EOD\n";
    for (const Segment& s : segments) {
        // splot joins blocks separated by a single blank line into a grid,
        // so every 3D segment gets its own block (two blank lines).
        if (quadric)
            gp << s.a[0] << ' ' << s.a[1] << ' ' << s.a[2] << '\n'
               << s.b[0] << ' ' << s.b[1] << ' ' << s.b[2] << "\n\n\n";
        else
            gp << s.a[0] << ' ' << s.a[1] << '\n' << s.b[0] << ' ' << s.b[1] << "\n\n";
    }
    gp << "EOD\n"
       << (quadric ? "splot" : "plot") << " $zeroset with lines lc rgb '#d62728' lw "
       << (quadric ? 1 : 2) << "\n";
    if (interactive)
        gp << "pause mouse close\n";
    else
        gp << "unset output\n";
    gp.close();

    if (interactive) {
        // Runs in background: the window stays open (and can be rotated or
        // zoomed with the mouse) while ALGAgraph keeps accepting commands.
        // gnuplot messages go to a log file next to the script.
        std::filesystem::path log = script;
        log.replace_extension(".log");
        if (std::system(("gnuplot " + shellQuote(script.string()) + " > " +
                         shellQuote(log.string()) + " 2>&1 &")
                            .c_str()) != 0)
            throw std::runtime_error("could not start gnuplot");
    } else if (std::system(("gnuplot " + shellQuote(script.string())).c_str()) != 0) {
        throw std::runtime_error("gnuplot failed to write '" + output + "'");
    }
    return segments.size();
}

void checkOutputFormat(const std::string& output) {
    if (!output.empty()) terminalFor(output);
}

void cleanupPlots() {
    std::error_code ec;
    std::filesystem::remove_all(std::filesystem::temp_directory_path() /
                                    ("algagraph-" + std::to_string(getpid())),
                                ec);
}

} // namespace alga
