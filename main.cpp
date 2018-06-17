#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <cstring>
#include <iomanip>
#include "model.hpp"
#include "fft.hpp"
#include "utils.hpp"
#include "interpolation.hpp"
#include <boost/program_options.hpp>

using namespace std::complex_literals;
namespace po = boost::program_options;

template<typename T, typename C, typename D>
auto generate_interpolation(const C& coords, const D& data) {
    std::vector<T> result; result.reserve(coords.size());
    for (size_t i = 0; i < coords.size(); ++i)
        result.emplace_back(coords[i].begin(), coords[i].end(), data[i].begin(), data[i].end());
    return result;
};

template<bool UseCAP>
void output_result(const isope::model_<UseCAP>& model, const std::string& out_file, const bool binary,
                    const std::vector<std::vector<isope::types::complex_t>>& result, const size_t skip) {
    const uint64_t s1 = result.size(), s2 = result[0].size() / (UseCAP ? 2 : 1);
    if (binary) {
        std::ofstream out(out_file, std::ios_base::binary);
        const auto dx = model.dx(), dz = model.dz();
        out.write(reinterpret_cast<const char*>(&s2), sizeof(s1));
        out.write(reinterpret_cast<const char*>(&s1), sizeof(s2));
        out.write(reinterpret_cast<const char*>(&dx), sizeof(dx));
        out.write(reinterpret_cast<const char*>(&dz), sizeof(dz));
        out.write(reinterpret_cast<const char*>(model.x_coords().data() + (UseCAP ? s2 / 2 : 0)), s2 * sizeof(double));
        for (size_t i = 0; i < model.z_coords().size(); i += skip)
            out.write(reinterpret_cast<const char*>(&model.z_coords()[i]), sizeof(double));
        out.write(reinterpret_cast<const char*>(&model.z_coords().back()), sizeof(double));
        for (const auto& it : result)
            out.write(reinterpret_cast<const char*>(it.data() + (UseCAP ? s2 / 2 : 0)), s2 * sizeof(isope::types::complex_t));
    } else {
        std::ofstream out(out_file);
        out << s2 << '\n';
        out << s1 << '\n';
        out << model.dx() << '\n';
        out << model.dz() << '\n';
        for (size_t i = 0; i < s2; ++i)
            out << model.x_coords()[i + (UseCAP ? s2 / 2 : 0)] << ' ';
        out << '\n';
        for (size_t i = 0; i < model.z_coords().size(); i += skip)
            out << model.z_coords()[i] << ' ';
        out << model.z_coords().back() << '\n';
        for (const auto& it : result) {
            for (size_t i = 0; i < s2; ++i)
                out << std::setprecision(6) << it[i + (UseCAP ? s2 / 2 : 0)] << ' ';
            out << '\n';
        }
    }
}

int main(int argc, char* argv[]) {
    std::string in_file, out_file;
    double k, eps, delta, x0, x1, z0, z1;
    size_t nx, nz, skip, ne, verb;
    po::options_description desc("Allowed options");
    desc.add_options()
            ("help,h", "Print this message")
            ("wavenumber,k", po::value(&k), "Wave number")
            ("eps,e", po::value(&eps), "Kerr constant")
            ("delta,d", po::value(&delta),"Delta parameter for CAP")
            ("initcond,i", po::value(&in_file), "Initial condition data")
            ("outfile,o", po::value(&out_file), "Output file")
            ("nx", po::value(&nx), "Number of points for x")
            ("x0", po::value(&x0), "Left x coordinate")
            ("x1", po::value(&x1),"Right x coordinate")
            ("nz", po::value(&nz),"Number of points for z")
            ("z0", po::value(&z0),"Bottom z coordinate")
            ("z1", po::value(&z1),"Top z coordinate")
            ("n,n", po::value(&ne), "Number of equations to use")
            ("skip,s", po::value(&skip)->default_value(100), "Output each \'skip\' row")
            ("verbose,v", po::value(&verb)->default_value(0), "Show progress")
            ("binaryinput", "Binary file output")
            ("binaryoutput", "Binary file input")
            ("linear,l", "Use linear interpolation")
            ("cap,c", "Use CAP");
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    if (vm.count("help")) {
        std::cout << desc << std::endl;
        return 0;
    }
    std::vector<std::vector<double>> coords;
    std::vector<std::vector<isope::types::complex_t>> data;
    if (vm.count("binaryinput")) {
        std::ifstream in(in_file, std::ios_base::binary);
        uint16_t n = 0;
        in.read(reinterpret_cast<char*>(&n), 16);
        coords.resize(n); data.resize(n);
        for (uint16_t i = 0; i < n; ++i) {
            uint64_t count = 0;
            in.read(reinterpret_cast<char*>(&count), 64);
            coords[i].resize(count); data[i].resize(count);
            in.read(reinterpret_cast<char*>(coords[i].data()), count * sizeof(double));
            in.read(reinterpret_cast<char*>(data[i].data()), count * sizeof(isope::types::complex_t));
        }
    }
    else {
        std::ifstream in(in_file);
        uint16_t n = 0; in >> n;
        coords.resize(n); data.resize(n);
        for (uint16_t i = 0; i < n; ++i) {
            uint64_t count = 0; in >> count;
            coords[i].resize(count); data[i].resize(count);
            for (auto& it : coords[i]) in >> it;
            for (auto& it : data[i]) in >> it;
        }
    }
    if (vm.count("cap")) {
        isope::cap_model model(k, eps, nx, x0, x1, nz, z0, z1, delta);
        isope::types::vector2d_t<isope::types::complex_t> result;
        if (vm.count("linear"))
            result = model.solve(generate_interpolation<isope::interpolation::linear<double, isope::types::complex_t>>(coords, data), ne, skip, verb);
        else
            result = model.solve(generate_interpolation<isope::interpolation::polynomial<double, isope::types::complex_t>>(coords, data), ne, skip, verb);
        output_result(model, out_file, vm.count("binaryoutput") > 0, result, skip);
    }
    else {
        isope::model model(k, eps, nx, x0, x1, nz, z0, z1, delta);
        isope::types::vector2d_t<isope::types::complex_t> result;
        if (vm.count("linear"))
            result = model.solve(generate_interpolation<isope::interpolation::linear<double, isope::types::complex_t>>(coords, data), ne, skip, verb);
        else
            result = model.solve(generate_interpolation<isope::interpolation::polynomial<double, isope::types::complex_t>>(coords, data), ne, skip, verb);
        output_result(model, out_file, vm.count("binaryoutput") > 0, result, skip);
    }
}
