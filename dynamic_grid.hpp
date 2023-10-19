//
// Created by a.pesterev on 12.03.2023.
//

#ifndef POROUS_MEDIA_DYNAMIC_GRID_HPP
#define POROUS_MEDIA_DYNAMIC_GRID_HPP

#include <map>

template<INTEGER_TYPE c>
class DynamicGrid {
public:
    explicit DynamicGrid(const std::function<void(const std::string&, std::ofstream&)>& mapKey, const std::string& filename="", const std::string& path="./") {
        this->xField = std::vector<std::vector<FLOAT_TYPE>>(c, std::vector<FLOAT_TYPE>(c, static_cast<FLOAT_TYPE>(0.)));
        this->yField = std::vector<std::vector<FLOAT_TYPE>>(c, std::vector<FLOAT_TYPE>(c, static_cast<FLOAT_TYPE>(0.)));
        if (!std::filesystem::exists(path)) {
            std::cout << "Target file path does not exists!";
            exit(1);
        }
        if (!std::filesystem::is_directory(path)) {
            std::cout << "Target file path is not directory!";
            exit(1);
        }
        std::filesystem::path path_o = std::filesystem::path(path);
        path_o.append(filename);
        this->output_file = std::ofstream(path_o.string());
        this->mapKey = mapKey;
    }

    FLOAT_TYPE put(INTEGER_TYPE _x, INTEGER_TYPE _y, FLOAT_TYPE _val) {
        this->xField[_x][_y] = _val;
        this->yField[_y][_x] = _val;
        this->mapKey(std::to_string(_x) + ";" + std::to_string(_y), this->output_file);
        this->output_file << _val << ";\n";
        return _val;
    }

    [[nodiscard]] FLOAT_TYPE get(INTEGER_TYPE _x, INTEGER_TYPE _y) const {
        const auto v1 = xField[_x][_y];
        const auto v2 = yField[_y][_x];
        if (std::abs(v1 - v2) > std::numeric_limits<FLOAT_TYPE>::epsilon()) {
            exit(12);
        }
        return v1;
    }

    ~DynamicGrid() {
        if (this->output_file.is_open()) {
            this->output_file.flush();
            this->output_file.close();
        }
    }

    void flush() {
        this->output_file.flush();
    }
private:
    DynamicGrid() = default;
private:
    std::vector<std::vector<FLOAT_TYPE>> xField;
    std::vector<std::vector<FLOAT_TYPE>> yField;
    std::ofstream output_file;
    std::function<void(const std::string&, std::ofstream&)> mapKey;
};

#endif //POROUS_MEDIA_DYNAMIC_GRID_HPP
