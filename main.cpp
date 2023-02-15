#include <iostream>
#include <json/json.h>
#include <fstream>
#include <filesystem>
#include <cmath>
#include <Eigen/Sparse>

#define out(mes, val) std::cout << "\t" << mes << ": " << val << '\n';

#define INDEX_TYPE long long

#define IS_FLOAT

#ifndef IS_FLOAT
#ifndef IS_DOUBLE
#error "No type specified"
#else
typedef long long INTEGER_TYPE;
typedef double FLOAT_TYPE;
typedef double CONTAINER_TYPE;
const FLOAT_TYPE eps = DBL_EPSILON;
#endif
#else
typedef long long INTEGER_TYPE;
typedef float FLOAT_TYPE;
typedef float CONTAINER_TYPE;
const FLOAT_TYPE eps = FLT_EPSILON;
#endif

bool check_differs(FLOAT_TYPE a, FLOAT_TYPE b) {
    return fabs(a - b) > eps;
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

    std::cerr.precision(20);
    std::cout.precision(20);

    INDEX_TYPE i, j;
    std::cout << "Configuration: " << argv[1] << '\n';
    std::ifstream configuration_json(argv[1], std::ifstream::binary);
    Json::Reader reader;
    Json::Value root;
    reader.parse(configuration_json, root);
    const INDEX_TYPE M = root.get("M", (INDEX_TYPE) 0).asInt64();
    const INDEX_TYPE N = root.get("N", (INDEX_TYPE) 0).asInt64();

    const FLOAT_TYPE dx = (FLOAT_TYPE) root.get("dx", double(1)).as<double>();
    const FLOAT_TYPE dt = (FLOAT_TYPE) root.get("dt", double(1)).as<double>();
    const FLOAT_TYPE delta_0 = (FLOAT_TYPE) root.get("delta_0", double(1)).as<double>();
    const FLOAT_TYPE s_m = (FLOAT_TYPE) root.get("s_m", double(1)).as<double>();
    const FLOAT_TYPE c_x_0 = (FLOAT_TYPE) root.get("C_x_0", double(0)).as<double>();
    const FLOAT_TYPE c_0_t = (FLOAT_TYPE) root.get("C_0_t", double(1)).as<double>();
    const FLOAT_TYPE s_x_0 = (FLOAT_TYPE) root.get("S_x_0", double(0)).as<double>();

    std::cout << "Configuration:\n";
    out("M", M)
    out("N", N)
    out("dx", dx)
    out("dt", dt)
    out("s_m", s_m)
    out("delta_0", delta_0)
    out("c_x_0", c_x_0)
    out("c_0_t", c_0_t)
    out("s_x_0", s_x_0)
    std::cout << '\n';

    const std::function<FLOAT_TYPE(FLOAT_TYPE)> concentration_front = [delta_0, s_m](FLOAT_TYPE t) {
        return (s_m / delta_0) * (FLOAT_TYPE(1) - exp(-delta_0 * t / s_m));
    };

    const std::function<FLOAT_TYPE(FLOAT_TYPE)> S_0 = [delta_0, s_m](FLOAT_TYPE t) {
        return s_m * (FLOAT_TYPE(1) - exp(-delta_0 * t / s_m));
    };

    const std::function<FLOAT_TYPE(FLOAT_TYPE)> flow_rate = [delta_0, S_0, s_m](FLOAT_TYPE t) {
        return exp(-delta_0 * t / s_m);
    };
    const std::function<FLOAT_TYPE(FLOAT_TYPE)> Delta = [delta_0, s_m](FLOAT_TYPE S) {
        return delta_0 * (FLOAT_TYPE(1) - S / s_m);
    };


    Eigen::SparseMatrix<FLOAT_TYPE> c(N, M);
    Eigen::SparseMatrix<FLOAT_TYPE> s(N, M);

    for (i = 0; i < N; i++) {
        // C(x, 0) = C|t=0
        c.coeffRef(i, 0) = c_x_0;
    }

    for (i = 0; i < N; i++) {
        // S(x, 0) = S|t=0
        s.coeffRef(i, 0) = s_x_0;
    }

    for (j = 0; j < M; j++) {
        // C(0, t) = C|x=0
        c.coeffRef(0, j) = c_0_t;
    }

    std::vector<FLOAT_TYPE> concentration_front_at_t;
    std::vector<FLOAT_TYPE> flow_rate_at_t;

    for (j = 0; j < M; j++) {
        // Вычисление точек фронта концентрации
        const auto concentration_front_at_t_j = concentration_front(dt * FLOAT_TYPE (j));
        concentration_front_at_t.push_back(concentration_front_at_t_j);
        // Проверка с точным решением
        // Вычисление точек скорости распространения суспензии
        const auto flow_rate_at_time_j = flow_rate(FLOAT_TYPE(j) * dt);
        flow_rate_at_t.push_back(flow_rate_at_time_j);
    }

    // 'i' is for x measure, 'j' is for time measure
    for (j = 0; j < M - 1; j++) {
        const auto flow_rate_j = flow_rate_at_t[j];
        const auto concentration_front_v = concentration_front_at_t[j];
        const auto concentration_front_next = concentration_front_at_t[j + 1];
        for (i = 0; i < N - 1; i++) {
            const auto x_i = FLOAT_TYPE(i) * dx;
            const auto x_i_next = x_i + dx;
            const auto c_i_j = c.coeff(i, j);
            const auto s_i_j = s.coeff(i, j);
            const auto s_i_next_j = s.coeff(i + 1, j);
            const auto ds_dt = Delta(s_i_j) * c_i_j;
            s.coeffRef(i, j + 1) = s_i_j + dt * ds_dt;
            const auto ds_next_dt = Delta(s_i_next_j) * c.coeff(i + 1, j);
            s.coeffRef(i + 1, j + 1) = ds_next_dt * dt + s_i_next_j;
            const bool c2 = check_differs(s_i_next_j, s.coeff(i + 1, j + 1));
            if (((concentration_front_next - x_i_next) > eps) && ((concentration_front_v - x_i) > eps)) {
                const auto c_i_j_next = c.coeff(i, j + 1);
                const auto dc_dt = (c_i_j_next - c_i_j) / dt;
                const bool c1 = check_differs(c.coeff(i + 1, j), c.coeff(i + 1, j - 1));

                const auto dc_dx = (c1 || c2) ? (-ds_dt - dc_dt) / flow_rate_j : FLOAT_TYPE(0.);

                const auto c_i_next_j_next = c_i_j_next + dx * dc_dx;
                c.coeffRef(i + 1, j + 1) = c_i_next_j_next;

                if (fabs(((c_i_j_next - c_i_j) + (flow_rate_j * (c_i_next_j_next - c_i_j_next)) +
                          (Delta(s_i_j) * c_i_j * dt))) >= std::max(dx, dt)) {
                    std::cerr << c.coeff(i + 1, j + 1) << " " << c.coeff(i + 1, j) << '\n';
                    std::cerr << "C(" << i + 1 << ", " << j + 1 << ") = " << c_i_next_j_next << '\n';
                    std::cerr << "C(" << i << ", " << j << ") = " << c_i_j << '\n';
                    std::cerr << "C(" << i << ", " << j + 1 << ") = " << c_i_j_next << '\n';
                    std::cerr << "S(" << i << ", " << j << ") = " << s_i_j << '\n';
                    std::cerr << "/\\(S) = " << Delta(s_i_j) << '\n';
                    std::cerr << "v(t) = " << flow_rate_j << '\n';
                    std::cerr << "dc/dx = " << c_i_next_j_next - c_i_j_next << '\n';
                    std::cerr << "dc/dt = " << c_i_j_next - c_i_j << '\n';
                    std::cerr << "v * dc/dx = " << flow_rate_j * (c_i_next_j_next - c_i_j_next) << '\n';
                    std::cerr << "ds/dt = " << Delta(s_i_j) * c_i_j * dt << '\n';

                    std::cerr << "dc/dt + v*dc/dx + ds/dt = "
                              << ((c_i_j_next - c_i_j) + (flow_rate_j * (c_i_next_j_next - c_i_j_next)) +
                                  (Delta(s_i_j) * c_i_j * dt)) << " = 0?" << '\n';

                    std::cerr << "t = " << dt * FLOAT_TYPE(j) << '\n';
                    std::cerr << "x = " << x_i << " Concentration front G " << concentration_front_v << '\n';
                    std::cerr << "i = " << i << " j = " << j << '\n';
                    exit(3);
                }
            }
        }
    }

    std::cout << "Ready!" << std::endl;
    std::ofstream output("results/result.csv");
    output.precision(20);
    output << "dx;\n" << dx << "\n" << "dt;\n" << dt << "\n\n\n";
    output << "x;t;c(i, j);s;v;G;\n";
    std::cout << "Starting writing results!\n";

    for (i = 0; i < N; i++) {
        const auto x = dx * FLOAT_TYPE(i);
        for (j = 0; j < M; j++) {
            const auto t = dt * FLOAT_TYPE(j);
            auto c_i_j = FLOAT_TYPE(0);
            auto s_i_j = s.coeff(i, j);
            auto dc_dt = FLOAT_TYPE(0);
            auto dc_dx = FLOAT_TYPE(0);
            auto ds_dt = FLOAT_TYPE(0);
            if (concentration_front_at_t[j] >= x) {
                c_i_j = c.coeff(i, j);
                if (i < N - 1) {
                    const auto c_i_next_j = c.coeff(i + 1, j);
                    dc_dx = (c_i_next_j - c_i_j) / dx;
                }
                if (j < M - 1) {
                    const auto c_i_j_next = c.coeff(i, j + 1);
                    dc_dt = (c_i_j_next - c_i_j) / dt;
                    const auto s_i_j_next = s.coeff(i, j + 1);
                    ds_dt = (s_i_j_next - s_i_j) / dt;
                }
            }
            output << x << ';' << t << ';'
                   << c_i_j << ';'
                   << s_i_j << ';'
                   << flow_rate_at_t[j] << ';'
                   << concentration_front_at_t[j] << ';'
                   << '\n';
        }
    }

    std::cout << "Saved!\n";
    if (output.is_open()) { output.close(); }
    return 0;
}
