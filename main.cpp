#include <iostream>
#include <json/json.h>
#include <fstream>
#include <filesystem>
#include <sstream>
#include <cmath>
#include <Eigen/Dense>

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
constexpr FLOAT_TYPE eps = DBL_EPSILON;
#endif
#else
typedef long long INTEGER_TYPE;
typedef float FLOAT_TYPE;
typedef float CONTAINER_TYPE;
constexpr FLOAT_TYPE eps = FLT_EPSILON;
#endif

const FLOAT_TYPE one = FLOAT_TYPE(1);

#define to_f(x) FLOAT_TYPE(x)

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
        return (s_m / delta_0) * (one - exp(-delta_0 * t / s_m));
    };

    const std::function<FLOAT_TYPE(FLOAT_TYPE)> S_0 = [delta_0, s_m](FLOAT_TYPE t) {
        return s_m * (one - exp(-delta_0 * t / s_m));
    };

    const std::function<FLOAT_TYPE(FLOAT_TYPE)> flow_rate = [delta_0, S_0, s_m](FLOAT_TYPE t) {
        return exp(-delta_0 * t / s_m);
    };
    const std::function<FLOAT_TYPE(FLOAT_TYPE)> Delta = [delta_0, s_m](FLOAT_TYPE S) {
        return delta_0 * (one - (S / s_m));
    };

    const auto calc = [](const FLOAT_TYPE c_i_j_next, const FLOAT_TYPE dc_dt, const FLOAT_TYPE ds_dt,
                         const FLOAT_TYPE flow_rate) {
        return c_i_j_next - (one / flow_rate) * (dc_dt + ds_dt);
    };

    std::cout << "Creating containers!\n";

    Eigen::MatrixX<FLOAT_TYPE> c(N + 2, M + 2);
    Eigen::MatrixX<FLOAT_TYPE> s(N + 2, M + 2);

    std::cout << "Created containers!\n";

    for (i = 0; i < N + 2; i++) {
        for (j = 0; j < M + 2; j++) {
            c.coeffRef(i, j) = 0;
            s.coeffRef(i, j) = 0;
        }
    }

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
        const auto s_0_j = s.coeff(0, j);
        s.coeffRef(0, j + 1) = s_0_j + dt * Delta(s_0_j) * c_0_t;
    }

    std::vector<FLOAT_TYPE> concentration_front_at_t;
    std::vector<FLOAT_TYPE> flow_rate_at_t;

    std::cout << "Starting calculations!\n";

    // 'i' is for x measure, 'j' is for time measure
    for (j = 0; j < M; j++) {
        const auto t = dt * to_f(j);
        const auto t_next = dt * to_f(j + 1);
        // Вычисление точек фронта концентрации
        const auto concentration_front_at_t_j = concentration_front(t);
        const auto concentration_front_at_j_next = concentration_front(t_next);
        concentration_front_at_t.push_back(concentration_front_at_t_j);
        // Проверка с точным решением
        // Вычисление точек скорости распространения суспензии
        const auto flow_rate_at_time_j = flow_rate(t);
        const auto flow_rate_at_time_j_next = flow_rate(t_next);
        flow_rate_at_t.push_back(flow_rate_at_time_j);

        for (i = 0; i < N; i++) {
            const auto x_next = to_f(i + 1) * dx;
            const auto x = to_f(i) * dx;
            bool f1 = (concentration_front_at_t_j > x);
            bool f2 = (concentration_front_at_j_next > x_next);
//            std::cout << "\n\n\n\n";
//            std::cout << "i = " << i << " j = " << j << '\n';
//            std::cout << "x = " << x << " G = " << concentration_front_at_t_j << " f1 = " << f1 << '\n';
//            std::cout << "x_next = " << x_next << " G_next = " << concentration_front_at_j_next << " f2 = " << f2 << '\n';
            if (f2) {
                const auto c_i_j = c.coeff(i, j);
                const auto c_i_j_next = c.coeff(i, j + 1);
                const auto s_i_j = s.coeff(i, j);
                const auto s_i_j_next = s.coeff(i, j + 1);
                const auto dc_dt = (c_i_j_next - c_i_j) / dt;
                const auto ds_dt1 = (s_i_j_next - s_i_j) / dt;
                const auto inv_speed = one / flow_rate_at_time_j_next;
                const auto c_i_next_j_next = c_i_j_next - dx * inv_speed * (dc_dt + ds_dt1);
//                std::cout << c_i_j << ' ' << c_i_j_next << '\n';
//                std::cout << s_i_j << ' ' << s_i_j_next << '\n';
//                std::cout << "dc/dt = " << dc_dt << '\n';
//                std::cout << "ds_dt = " << ds_dt1 << '\n';
//                std::cout << "1/v = " << inv_speed << '\n';
//                std::cout << "C(" << x_next << "," << t_next << ") = " << c_i_next_j_next << '\n';
                c.coeffRef(i + 1, j + 1) = c_i_next_j_next;
                const auto s_i_next_j_next = s.coeff(i + 1, j + 1);
                const auto s_i_next_j_next_next = s_i_j_next + dt * Delta(s_i_next_j_next) * c_i_next_j_next;
                s.coeffRef(i + 1, j + 2) = s_i_next_j_next_next;
            }
        }
    }

    std::cout << "Ready!\n";
    std::stringstream ss;
    ss << "results/result" << M << "_" << N << ".csv";
    std::string fname = ss.str();
    std::ofstream output(fname);
    output.precision(4);
    output << "dx;\n" << dx << "\n" << "dt;\n" << dt << "\n\n\n";
    output << "i;j;x;t;c(i, j);s;v;G;\n";
    std::cout << "Starting writing results!\n";

    for (i = 0; i < N; i++) {
        const auto x = dx * FLOAT_TYPE(i);
        for (j = 0; j < M; j++) {
            const auto t = dt * FLOAT_TYPE(j);
            auto c_i_j = c.coeff(i, j);
            auto s_i_j = s.coeff(i, j);
//            std::cout << "i = " << i << " j = " << j << '\n';
//            std::cout << c_i_j << " " << s_i_j << '\n';
            output << i << ';' << j << ';' << x << ';' << t << ';'
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
