#include <iostream>
#include <json/json.h>
#include <functional>
#include <fstream>
#include <filesystem>
#include <cmath>

#define out(mes, val) std::cout << "\t" << mes << ": " << val << '\n';

#define INDEX_TYPE int32_t

#define IS_DOUBLE

#ifndef IS_FLOAT
#ifndef IS_DOUBLE
#error "No type specified"
#else
typedef long long INTEGER_TYPE;
typedef double FLOAT_TYPE;
typedef double CONTAINER_TYPE;
constexpr FLOAT_TYPE eps = DBL_EPSILON;
constexpr INTEGER_TYPE precision = DBL_DIG;
#endif
#else
typedef long long INTEGER_TYPE;
typedef float FLOAT_TYPE;
typedef float CONTAINER_TYPE;
constexpr FLOAT_TYPE eps = FLT_EPSILON;
constexpr INTEGER_TYPE precision = FLT_DIG;
#endif

const FLOAT_TYPE one = FLOAT_TYPE(1);

#define to_f(x) FLOAT_TYPE(x)

#include <vector>
#include "methods.hpp"
#include "dynamic_grid.hpp"

double __fastcall two_sum (double &t, double a, double b) {
    double s = a+b;
    double bs = s-a;
    double as = s-bs;
    t = (b-bs) + (a-as);
    return s;
}

int main(int argc, char **argv) {
    std::cout << "Starting program!\n";
    if (argc != 2) {
        std::cout << "Usage: porous_media [configuration_file.json]\n";
        return -1;
    }
    if (!std::filesystem::is_directory("results") || !std::filesystem::exists("results")) {
        std::filesystem::create_directory("results");
    }

    INDEX_TYPE i, j;
    std::cout << "Configuration: " << argv[1] << '\n';
    std::ifstream configuration_json(argv[1], std::ifstream::binary);
    Json::Reader reader;
    Json::Value root;
    reader.parse(configuration_json, root);

    const FLOAT_TYPE dt = (FLOAT_TYPE) root.get("dt", FLOAT_TYPE(1)).as<FLOAT_TYPE>();
    const FLOAT_TYPE max_t = (FLOAT_TYPE) root.get("max_t", FLOAT_TYPE(1)).as<FLOAT_TYPE>();
    const FLOAT_TYPE delta_0 = (FLOAT_TYPE) root.get("delta_0", FLOAT_TYPE(1)).as<FLOAT_TYPE>();
    const FLOAT_TYPE F0 = (FLOAT_TYPE) root.get("F0", FLOAT_TYPE(0)).as<FLOAT_TYPE>();
    const FLOAT_TYPE v_0 = (FLOAT_TYPE) root.get("v0", FLOAT_TYPE(1)).as<FLOAT_TYPE>();
    const FLOAT_TYPE s_m = (FLOAT_TYPE) root.get("s_m", FLOAT_TYPE(1)).as<FLOAT_TYPE>();
    const FLOAT_TYPE c_x_0 = (FLOAT_TYPE) root.get("C_x_0", FLOAT_TYPE(0)).as<FLOAT_TYPE>();
    const FLOAT_TYPE c_0_t = (FLOAT_TYPE) root.get("C_0_t", FLOAT_TYPE(1)).as<FLOAT_TYPE>();
    const FLOAT_TYPE s_x_0 = (FLOAT_TYPE) root.get("S_x_0", FLOAT_TYPE(0)).as<FLOAT_TYPE>();

    std::cout << "Configuration:\n";
    out("dt", dt)
    out("max_t", max_t)
    out("s_m", s_m)
    out("delta_0", delta_0)
    out("c_x_0", c_x_0)
    out("c_0_t", c_0_t)
    out("s_x_0", s_x_0)
    std::cout << '\n';

    const FLOAT_TYPE A = (delta_0 / s_m) + F0;

    const std::function<FLOAT_TYPE (const FLOAT_TYPE&)> F = [F0](const FLOAT_TYPE& S) {
        return F0 * S;
    };

    const std::function<FLOAT_TYPE(FLOAT_TYPE)> S_0 = [delta_0, s_m, F0, A](FLOAT_TYPE t) {
        return -delta_0*s_m*(exp(-(delta_0+F0)*t)-1)/(delta_0+F0);
    };

    const std::function<FLOAT_TYPE(FLOAT_TYPE)> flow_rate = [delta_0, A, F0, v_0](FLOAT_TYPE t) {
//        return (v_0 / A) * (F0 + delta_0 * exp(-A * t));
        return v_0;
    };
    const std::function<FLOAT_TYPE(FLOAT_TYPE)> Delta = [delta_0, s_m](FLOAT_TYPE S) {
        return delta_0 * (s_m - S);
    };

    auto dots_G = std::vector<FLOAT_TYPE>();
    auto deltas = std::vector<FLOAT_TYPE>();

    dots_G.push_back(FLOAT_TYPE(0.));

    auto x = FLOAT_TYPE(0.);
    auto t = FLOAT_TYPE(0.);

    const FLOAT_TYPE half_dt = dt / FLOAT_TYPE(2.);

    j = 0;

    const std::function<void(const std::string&, std::ofstream&)> mapKey = [&deltas, dt](const std::string& key, std::ofstream& output) {
        std::istringstream is(key);
        INTEGER_TYPE ix, iy;
        char del;
        is >> ix >> del >> iy;
        FLOAT_TYPE x = 0, y = FLOAT_TYPE (iy) * dt;
        INTEGER_TYPE bix = ix;
        INTEGER_TYPE biy = iy;
        while (bix != 0) {
            x += deltas[biy];
            bix--;
            biy--;
        }
        output.precision(precision);
        output << ix << del << iy << del << x << del << y << del;
    };

    DynamicGrid<FLOAT_TYPE, INTEGER_TYPE , INTEGER_TYPE, precision> C = DynamicGrid<FLOAT_TYPE, INTEGER_TYPE, INTEGER_TYPE, precision>(mapKey, "conc.dat", "./results");
    DynamicGrid<FLOAT_TYPE, INTEGER_TYPE, INTEGER_TYPE, precision> S = DynamicGrid<FLOAT_TYPE, INTEGER_TYPE, INTEGER_TYPE, precision>(mapKey, "rest.dat", "./results");

    S.put(0., 0., 0.);

    std::cout << "Created grids!\n";

    // i is for time
    while (t <= max_t) {
        auto dx = integrate(flow_rate, t, t + dt);
        C.put(0, j, c_0_t);
        const auto s_0_t = S.get(0, j);
        t += dt;
        j++;
        S.put(0., j, (Delta(s_0_t) - F(s_0_t)) * dt + s_0_t);
        x += dx;
        dots_G.push_back(x);
        deltas.push_back(dx);
//        std::cout << "t: " << t << " x: " << x << " dx: " << dx << '\n';
    }

    const auto max_j = j;

    std::cout << "Max j = " << max_j << '\n';

    t = FLOAT_TYPE(0.);

    std::cout.precision(precision);

    double e = 0.;

    for (j = 0; j + 1 < max_j; j++) {
        auto next_t_s = two_sum(e, t, dt);
        const auto next_t = next_t_s + e;
        if (j % 100 == 0) { std::cout << t << '\n'; }
        for (i = 0; i <= j; i++) {
            /*
             * Iterating over each x in given time
             */

            const auto c_x_t = C.get(i, j);
            const auto s_x_t = S.get(i, j);
            auto s_x_next_t_next = FLOAT_TYPE (0.);
            if (i < j) {
                const auto x_i_j = get_x(i, j, deltas);
                const auto x_i_next_j = get_x(i + 1, j, deltas);
                const auto x_i_next_j_next = get_x(i + 1, j + 1, deltas);
                const auto s_l = s_x_t + dt * (Delta(s_x_t) * c_x_t - F(s_x_t));
                const auto s_x_next_t = S.get(i + 1, j);
                const auto c_x_next_t = C.get(i + 1, j);
                const auto s_r = s_x_next_t + dt * (Delta(s_x_next_t) * c_x_next_t - F(s_x_next_t));
                s_x_next_t_next = interpolate(s_l, s_r, x_i_j, x_i_next_j, x_i_next_j_next);
            }
            S.put(i + 1, j + 1, s_x_next_t_next);
            const auto m1 = one - half_dt * (Delta(s_x_t) - F(s_x_next_t_next) - F(s_x_t));
            const auto m2 = one + half_dt * Delta(s_x_next_t_next);
            C.put(i + 1, j + 1, c_x_t * (m1 / m2));
        }
        t = next_t;
        C.flush();
        S.flush();
        for (i = 0; i <= j; i++) {
            C.remove(i, j);
            S.remove(i, j);
        }
    }

    std::cout << "Ready!\n";

    return 0;

}