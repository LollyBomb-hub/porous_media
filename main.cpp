#include <iostream>
#include <json/json.h>
#include <fstream>
#include <filesystem>
#include <cmath>
#include <Eigen/Sparse>
#include <gmp.h>

#define out(mes, val) std::cout << "\t" << mes << ": " << val << '\n';

#define TESTING

#define INDEX_TYPE long long

#define IS_FLOAT

#ifndef IS_FLOAT
#ifndef IS_DOUBLE
#ifndef IS_LONG_DOUBLE
#define MPZ_T
#include <string>
typedef mpz_t INTEGER_TYPE;
typedef mpf_t FLOAT_TYPE;
typedef std::string CONTAINER_TYPE;
#else
typedef long long INTEGER_TYPE;
typedef long double FLOAT_TYPE;
typedef long double CONTAINER_TYPE;
#endif
#else
typedef long long INTEGER_TYPE;
typedef double FLOAT_TYPE;
typedef double CONTAINER_TYPE;
#endif
#else
typedef long long INTEGER_TYPE;
typedef float FLOAT_TYPE;
typedef float CONTAINER_TYPE;
#endif

#ifndef MPZ_T
const FLOAT_TYPE eps = FLOAT_TYPE(1e-5);
const FLOAT_TYPE one = FLOAT_TYPE(1);
const FLOAT_TYPE zero = FLOAT_TYPE(0);

FLOAT_TYPE
concentration_front(const std::function<FLOAT_TYPE(FLOAT_TYPE)> &f_flow_rate, FLOAT_TYPE h, INDEX_TYPE n, INDEX_TYPE multiplier = 100,
                    FLOAT_TYPE result_0 = FLOAT_TYPE(0)) {
    FLOAT_TYPE result = result_0;
    n *= multiplier;
    h /= FLOAT_TYPE(multiplier);
    for (INDEX_TYPE i = 0; i < n; i++) {
        result += h * f_flow_rate(h * FLOAT_TYPE(i));
    }
    return result;
}


FLOAT_TYPE exact_flow_rate(FLOAT_TYPE t, FLOAT_TYPE s_m, FLOAT_TYPE delta_0) {
//    return (TYPE(1) - delta_0 * t / s_m);
    return exp(- delta_0 * t / s_m);
}


FLOAT_TYPE exact_concentration_front(FLOAT_TYPE t, FLOAT_TYPE s_m, FLOAT_TYPE delta_0) {
//    if (fabs(delta_0) < eps) {
//        return t;
//    }
//    const auto T = s_m / delta_0;
////    return t - ((t * t) / 2 * T);
//    return T * (FLOAT_TYPE(1) - exp(-t / T));
    return t;
}
#else
FLOAT_TYPE eps;
FLOAT_TYPE one;
FLOAT_TYPE zero;
FLOAT_TYPE exponent_const;
#endif


int main(int argc, char **argv) {
    std::cout << "Starting program!\n";
#ifdef MPZ_T
    try {
        mpf_init2(eps, 64);
        mpf_init2(one, 64);
        mpf_init2(zero, 64);
        mpf_init2(exponent_const, 64);
        mpf_set_str(eps, "0.0001", 10);
        mpf_set_str(one, "1", 10);
        mpf_set_str(zero, "0", 10);
        mpf_set_str(exponent_const, "2.718281828459045235360287471352", 10);
    } catch (const std::exception& exception) {
        std::cout << exception.what() << '\n';
        return 1;
    }
    mp_exp_t mp_exp = 1;
    std::cout << mpf_get_str(nullptr, &mp_exp, 10, 20, exponent_const) << "E" << mp_exp << '\n';
#endif
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


#ifndef MPZ_T
    const FLOAT_TYPE dx = (FLOAT_TYPE)root.get("dx", double(1)).as<double>();
    const FLOAT_TYPE dt = (FLOAT_TYPE)root.get("dt", double(1)).as<double>();
    const FLOAT_TYPE delta_0 = (FLOAT_TYPE)root.get("delta_0", double(1)).as<double>();
    const FLOAT_TYPE s_m = (FLOAT_TYPE)root.get("s_m", double(1)).as<double>();
    const FLOAT_TYPE c_x_0 = (FLOAT_TYPE)root.get("C_x_0", double(0)).as<double>();
    const FLOAT_TYPE c_0_t = (FLOAT_TYPE)root.get("C_0_t", double(1)).as<double>();
    const FLOAT_TYPE s_x_0 = (FLOAT_TYPE)root.get("S_x_0", double(0)).as<double>();
    const FLOAT_TYPE a = (FLOAT_TYPE) root.get("a", double(1.)).as<double>();
    const FLOAT_TYPE b = (FLOAT_TYPE) root.get("b", double(2.)).as<double>();

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
#else
    const auto dx_string = root.get("dx", "0.001").as<std::string>();
    const auto dt_string = root.get("dt", "0.0005").as<std::string>();
    const auto delta_0_string = root.get("delta_0", "1").as<std::string>();
    const auto s_m_string = root.get("s_m", "1").as<std::string>();
    const auto c_x_0_string = root.get("C_x_0", "0").as<std::string>();
    const auto c_0_t_string = root.get("C_0_t", "1").as<std::string>();
    const auto s_x_0_string = root.get("S_x_0", "0").as<std::string>();

    FLOAT_TYPE dx;
    mpf_init2(dx, 64);
    mpf_set_str(dx, dx_string.c_str(), 10);
    FLOAT_TYPE dt;
    mpf_init2(dt, 64);
    mpf_set_str(dt, dt_string.c_str(), 10);
    FLOAT_TYPE delta_0;
    mpf_init2(delta_0, 64);
    mpf_set_str(delta_0, delta_0_string.c_str(), 10);
    FLOAT_TYPE s_m;
    mpf_init2(s_m, 64);
    mpf_set_str(s_m, s_m_string.c_str(), 10);
    FLOAT_TYPE c_x_0;
    mpf_init2(c_x_0, 64);
    mpf_set_str(c_x_0, c_x_0_string.c_str(), 10);
    FLOAT_TYPE c_0_t;
    mpf_init2(c_0_t, 64);
    mpf_set_str(c_0_t, c_0_t_string.c_str(), 10);
    FLOAT_TYPE s_x_0;
    mpf_init2(s_x_0, 64);
    mpf_set_str(s_x_0, s_x_0_string.c_str(), 10);
#endif

    const auto S_0 = [delta_0, s_m, a, b](FLOAT_TYPE t) {
#ifdef MPZ_T
        FLOAT_TYPE delta;
        FLOAT_TYPE exponent;
        FLOAT_TYPE exponent_argument;
        mpf_init2(delta, 64);
        mpf_init2(exponent, 64);
        mpf_init2(exponent_argument, 64);
        mpf_mul(exponent_argument, delta_0, t);
        mpf_div(exponent_argument, exponent_argument, s_m);
        mpf_neg(exponent_argument, exponent_argument);

        mpf_sub(delta, one, exponent);
        mpf_mul(delta, delta, s_m);
        return delta;
#else
        return (a/b - exp(-b * t) / b);
#endif
    };

#ifndef MPZ_T
    const std::function<FLOAT_TYPE(FLOAT_TYPE)> flow_rate = [S_0, s_m, a, b](FLOAT_TYPE t) {
//        return 1;
//        return v(t, s_0, s_m, delta_0);
//        return (FLOAT_TYPE(1) - S_0(t) / (a / b));
        return FLOAT_TYPE(0.5);
    };
    const std::function<FLOAT_TYPE(FLOAT_TYPE)> Delta = [delta_0, s_m, a, b](FLOAT_TYPE S) {
//        return (FLOAT_TYPE(1) - FLOAT_TYPE(2) * S);
//        return delta(S, delta_0, s_m);
        return (a - b * S);
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
        const auto concentration_front_at_t_j = concentration_front(flow_rate, dt, j);
        concentration_front_at_t.push_back(concentration_front_at_t_j);
        // Проверка с точным решением
#ifdef CHECK_EXACT_VALUES
        if (fabs(concentration_front_at_t_j - exact_concentration_front(FLOAT_TYPE(j) * dt, s_m, delta_0)) > eps) {
            std::cerr << j << '\n';
            std::cerr << concentration_front_at_t_j << " <> " << exact_concentration_front(FLOAT_TYPE(j) * dt, s_m, delta_0)
                      << '\n';
            throw std::runtime_error("Concentration front mismatch!");
        }
#endif
        // Вычисление точек скорости распространения суспензии
        const auto flow_rate_at_time_j = flow_rate(FLOAT_TYPE(j) * dt);
        flow_rate_at_t.push_back(flow_rate_at_time_j);
        // Сравнение
#ifdef CHECK_EXACT_VALUES
        if (fabs(flow_rate_at_time_j - exact_flow_rate(FLOAT_TYPE(j) * dt, s_m, delta_0)) > eps) {
            std::cerr << "Flow rate differs with exact one!";
            throw std::runtime_error("Flow rate differs with exact one!");
        }
#endif
    }

    // 'i' is for x measure, 'j' is for time measure
    for (j = 0; j < M - 1; j++) {
        const auto flow_rate_j = flow_rate_at_t[j];
        const auto concentration_front = concentration_front_at_t[j];
        const auto concentration_front_next = concentration_front_at_t[j + 1];
        for (i = 0; i < N - 1; i++) {
            const auto x_i = FLOAT_TYPE(i) * dx;
            const auto x_i_next = x_i + dx;
            const auto c_i_j = c.coeff(i, j);
            const auto s_i_j = s.coeff(i, j);
            const auto ds_dt = Delta(s_i_j) * c_i_j;
            s.coeffRef(i, j + 1) = s_i_j + dt * ds_dt;
            if (((concentration_front_next - x_i_next) > eps) && ((concentration_front - x_i) > eps)) {
                const auto c_i_j_next = c.coeff(i, j + 1);
                const auto dc_dt = (c_i_j_next - c_i_j) / dt;
                const auto dc_dx = (-ds_dt - dc_dt) / flow_rate_j;

                const auto c_i_next_j_next = c_i_j_next + dx * dc_dx;
                c.coeffRef(i + 1, j + 1) = c_i_next_j_next;

                if (fabs(((c_i_j_next - c_i_j) + (flow_rate_j * (c_i_next_j_next - c_i_j_next)) +
                          (Delta(s_i_j) * c_i_j * dt))) >= std::max(dx, dt)) {
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

                    std::cerr << "x = " << x_i << " Concentration front G " << concentration_front << '\n';
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
                    std::cout << c_i_next_j << " - " << c_i_j << '\n';
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
                   //<< i << ';' << j << ';'
                   << c_i_j << ';'
                   << s_i_j << ';'
                   << flow_rate_at_t[j] << ';'
                   << concentration_front_at_t[j] << ';'
//                   << dc_dx << ';'
//                   << dc_dt << ';'
//                   << ds_dt << ';'
                   << '\n';
        }
    }

    std::cout << "Saved!\n";
    if (output.is_open()) { output.close(); }
#endif
    return 0;
}
