#include <iostream>
#include <json/json.h>
#include <fstream>
#include <filesystem>
#include <cmath>
#include <Eigen/Sparse>

#define out(mes, val) std::cout << "\t" << mes << ": " << val << '\n';

#define TESTING

#define INDEX_TYPE long long

#define IS_FLOAT

#ifndef IS_FLOAT
#define TYPE double
#else
#define TYPE float
#endif


const TYPE eps = TYPE(1e-5);


TYPE v(TYPE t, const std::function<TYPE(TYPE, TYPE, TYPE)> &f_s_0, TYPE s_m, TYPE delta_0) {
    return TYPE(1) - (f_s_0(t, delta_0, s_m) / s_m);
}

TYPE s_0(TYPE t, TYPE delta_0, TYPE s_m) {
    return s_m * (TYPE(1) - exp((TYPE(-1) * t * delta_0) / s_m));
}

TYPE delta(TYPE S, TYPE delta_0, TYPE S_m) {
    return delta_0 * (TYPE(1) - TYPE(S / S_m));
}

TYPE
concentration_front(const std::function<TYPE(TYPE)> &f_flow_rate, TYPE h, INDEX_TYPE n, INDEX_TYPE multiplier = 100,
                    TYPE result_0 = TYPE(0), const TYPE *concentration_front_till_previous = nullptr) {
    if (concentration_front_till_previous != nullptr) {
        TYPE almost_result = *concentration_front_till_previous;
        return almost_result + h * f_flow_rate(TYPE(n - 1) * h);
    }
    TYPE result = result_0;
    n *= multiplier;
    h /= TYPE(multiplier);
    for (INDEX_TYPE i = 0; i < n; i++) {
        result += h * f_flow_rate(h * TYPE(i));
    }
    return result;
}


TYPE exact_flow_rate(TYPE t, TYPE s_m, TYPE delta_0) {
    return exp((TYPE(-1) * delta_0 * t) / s_m);
}


int main(int argc, char **argv) {
    if (argc != 2) {
        std::cout << "Usage: porous_media [configuration_file.json]\n";
        return -1;
    }
    if (!std::filesystem::is_directory("results") || !std::filesystem::exists("results")) {
        std::filesystem::create_directory("results");
    }

    INDEX_TYPE i, j;
    INDEX_TYPE M, N;
    TYPE h;
    TYPE delta_0 = TYPE(1);
    TYPE s_m = TYPE(0.5);
    TYPE c_x_0 = TYPE(0);
    TYPE c_0_t = TYPE(1);
    TYPE s_x_0 = TYPE(0);
    std::cout << "Configuration: " << argv[1] << '\n';
    std::ifstream configuration_json(argv[1], std::ifstream::binary);
    Json::Reader reader;
    Json::Value root;
    reader.parse(configuration_json, root);
    M = root.get("M", TYPE(0)).asInt64();
    N = root.get("N", TYPE(0)).asInt64();
    std::cout << std::setprecision(5);

    {
        h = root.get("h", TYPE(1)).as<TYPE>();
        delta_0 = root.get("delta_0", TYPE(1)).as<TYPE>();
        s_m = root.get("s_m", TYPE(1)).as<TYPE>();
//        a = root.get("a", TYPE(0.5)).as<TYPE>();
//        b = root.get("b", TYPE(1)).as<TYPE>();
        c_x_0 = root.get("C_x_0", TYPE(0)).as<TYPE>();
        c_0_t = root.get("C_0_t", TYPE(1)).as<TYPE>();
        s_x_0 = root.get("S_x_0", TYPE(0)).as<TYPE>();
    }

    std::cout << "Configuration:\n";
    out("M", M)
    out("N", N)
    out("h", h)
    out("c_x_0", c_x_0)
    out("c_0_t", c_0_t)
    out("s_x_0", s_x_0)
    std::cout << '\n';
    const std::function<TYPE(TYPE)> S_0 = [delta_0, s_m](TYPE t) { return s_0(t, delta_0, s_m); };
    const std::function<TYPE(TYPE)> flow_rate = [delta_0, s_m](TYPE t) { return v(t, s_0, s_m, delta_0); };
    const std::function<TYPE(TYPE)> Delta = [delta_0, s_m](TYPE S) { return delta(S, delta_0, s_m); };


    Eigen::SparseMatrix<TYPE> c(N, M);
    Eigen::SparseMatrix<TYPE> s(N, M);

    for (i = 0; i < N; i++) {
        c.coeffRef(i, 0) = c_x_0;
    }

    for (i = 0; i < N; i++) {
        s.coeffRef(i, 0) = s_x_0;
    }

    for (j = 0; j < M; j++) {
        c.coeffRef(0, j) = c_0_t;
    }

    std::vector<TYPE> concentration_front_at_t;
    std::vector<TYPE> flow_rate_at_t;

    for (j = 0; j < M; j++) {
        if (j == 0) { concentration_front_at_t.push_back(concentration_front(flow_rate, h, j)); }
        else {
            concentration_front_at_t.push_back(
                    concentration_front(flow_rate, h, j, 100, 0, &concentration_front_at_t[j - 1]));
        }
//        std::cout << "t = " << TYPE(j) * h << " x(t) = " << concentration_fron_at_t[j] << '\n';
        const auto flow_rate_at_time_j = flow_rate(TYPE(j) * h);
        flow_rate_at_t.push_back(flow_rate_at_time_j);
//        std::cout << "v(" << TYPE(j) * h << ") = " << flow_rate_at_time_j << '\n';
        if (fabs(flow_rate_at_time_j - exact_flow_rate(TYPE(j) * h, s_m, delta_0)) > eps) {
            throw std::runtime_error("Flow rate differs with exact one!");
        }

    }

    // 'i' is for x measure, 'j' is for time measure
    for (j = 0; j < M - 1; j++) {
        const auto flow_rate_j = flow_rate_at_t[j];
        const auto concentration_front = concentration_front_at_t[j];
        for (i = 0; i < N - 1; i++) {
            const auto x_i = TYPE(i) * h;
            if (concentration_front > x_i) {
                const auto c_i_j = c.coeff(i, j);
                const auto c_i_j_next = c.coeff(i, j + 1);
                const auto s_i_j = s.coeff(i, j);

                const auto c_i_j_div_flow_rate = c_i_j / flow_rate_j;
                const auto c_i_j_next_div_flow_rate = c_i_j_next / flow_rate_j;
                c.coeffRef(i + 1, j + 1) = c_i_j_next + c_i_j_div_flow_rate - c_i_j_next_div_flow_rate - c_i_j_div_flow_rate * Delta(s_i_j) * h;
                s.coeffRef(i, j + 1) = s_i_j + h * Delta(s_i_j) * c_i_j;

#ifndef TESTING
                if (c.coeff(i + 1, j) > c_0_t) {
                    std::cerr << "C(" + std::to_string(x_i) + ", " + std::to_string(TYPE(j) * h) +
                                 ") became more than c_0_t!\n";
                    std::cerr << "C(i + 1, j + 1) = " << c.coeffRef(i + 1, j) << '\n';
                    std::cerr << "C(i, j) = " << C_i_j << '\n';
                    std::cerr << "C(i, j + 1) = " << C_i_j_next << '\n';
                    std::cerr << "S(i, j) = " << S_i_j << '\n';
                    std::cerr << "/\\(S) = " << delta_s << '\n';
                    std::cerr << "v(t) = " << flow_rate_j << '\n';

                    std::cerr << "x = " << x_i << " Concentration front G " << concentration_front << '\n';
                    std::cerr << "i = " << i << " j = " << j << '\n';
                    throw std::runtime_error("Calculation error!");
                }
                if (c.coeff(i + 1, j + 1) < 0) {
                    std::cerr << "C(" + std::to_string(x_i) + ", " + std::to_string(TYPE(j + 1) * h) +
                                 ") became less than 0!\n";
                    std::cerr << "C(i + 1, j + 1) = " << c.coeffRef(i + 1, j + 1) << '\n';
                    std::cerr << "C(i, j) = " << C_i_j << '\n';
                    std::cerr << "C(i, j + 1) = " << C_i_j_next << '\n';
                    std::cerr << "S(i, j) = " << S_i_j << '\n';
                    std::cerr << "/\\(S) = " << delta_s << '\n';
                    std::cerr << "v(t) = " << flow_rate_j << '\n';

                    std::cerr << x_i << ' ' << concentration_front << '\n';
                    std::cerr << i << ' ' << j << '\n';
                    throw std::runtime_error("Calculation error!");
                }
                if (s.coeff(i, j + 1) >= s_m) {
                    std::cerr << "S(" + std::to_string(x_i) + ", " + std::to_string(TYPE(j + 1) * h) +
                                 ") became more than s_m!\n";
                    std::cerr << "C(i + 1, j + 1) = " << c.coeffRef(i + 1, j + 1) << '\n';
                    std::cerr << "C(i, j) = " << C_i_j << '\n';
                    std::cerr << "C(i, j + 1) = " << C_i_j_next << '\n';
                    std::cerr << "S(i, j) = " << S_i_j << '\n';
                    std::cerr << "/\\(S) = " << delta_s << '\n';
                    std::cerr << "v(t) = " << flow_rate_j << '\n';

                    std::cerr << x_i << ' ' << concentration_front << '\n';
                    std::cerr << i << ' ' << j << '\n';
                    throw std::runtime_error("Calculation error!");
                }
                if (isnan(c.coeffRef(i + 1, j + 1))) {
                    std::cerr << "C(" + std::to_string(x_i + h) + ", " + std::to_string(TYPE(j + 1) * h) +
                                 ") became nan!\n";
                    std::cerr << "C(i + 1, j + 1) = " << c.coeffRef(i + 1, j + 1) << '\n';
                    std::cerr << "C(i, j) = " << C_i_j << '\n';
                    std::cerr << "C(i, j + 1) = " << C_i_j_next << '\n';
                    std::cerr << "S(i, j) = " << S_i_j << '\n';
                    std::cerr << "/\\(S) = " << delta_s << '\n';
                    std::cerr << "v(t) = " << flow_rate_j << '\n';

                    std::cerr << x_i << ' ' << concentration_front << '\n';
                    std::cerr << i << ' ' << j << '\n';
                    throw std::runtime_error("Calculation error!");
                }
#endif
            }
        }
    }

    std::cout << "Ready!" << std::endl;
    std::ofstream output("results/result.csv");
    output << "h;\n" << h << "\n\n\n";
    output << "x;t;c(i, j);s;v;G;\n";
    std::cout << "Starting writing results!\n";

    for (i = 0; i < N; i++) {
        const auto x = h * TYPE(i);
        for (j = 0; j < M; j++) {
            const auto t = h * TYPE(j);
            auto c_i_j = TYPE(0);
            auto s_i_j = TYPE(0);
            if (concentration_front_at_t[j] >= x) {
                c_i_j = c.coeff(i, j);
                s_i_j = s.coeff(i, j);
            }
            output << x << ';' << t << ';'
                   //<< i << ';' << j << ';'
                   << c_i_j << ';'
                   << s_i_j << ';'
                   << flow_rate_at_t[j] << ';'
                   << concentration_front_at_t[j] << ';' << '\n';
        }
    }

    std::cout << "Saved!\n";
    if (output.is_open()) { output.close(); }
    return 0;
}
