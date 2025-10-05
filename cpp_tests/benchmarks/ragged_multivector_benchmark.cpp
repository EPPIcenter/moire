#include "../src/multivector.h"
#include <iostream>
#include <chrono>
#include <vector>

void benchmark_ragged_construction(const std::array<size_t, 2>& dimensions, const std::vector<size_t>& ragged_dimensions) {
    auto start = std::chrono::high_resolution_clock::now();
    RaggedMultiVector<float, 3> mv(dimensions, ragged_dimensions);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Ragged construction time: " << duration.count() << " seconds\n";
}

void benchmark_ragged_element_assignment(const std::array<size_t, 2>& dimensions, const std::vector<size_t>& ragged_dimensions) {
    RaggedMultiVector<float, 3> mv(dimensions, ragged_dimensions);
    auto start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < dimensions[0]; ++i) {
        for (size_t j = 0; j < dimensions[1]; ++j) {
            for (size_t k = 0; k < ragged_dimensions.size(); ++k) {
                for (size_t l = 0; l < ragged_dimensions[k]; ++l) {
                    mv.at({i, j, k, l}) = i;
                }
            }
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Ragged element assignment time: " << duration.count() << " seconds\n";
}

void benchmark_ragged_iteration(const std::array<size_t, 2>& dimensions, const std::vector<size_t>& ragged_dimensions) {
    RaggedMultiVector<int, 3> mv(dimensions, ragged_dimensions, 1);
    long long target = std::accumulate(ragged_dimensions.begin(), ragged_dimensions.end(), 0) * dimensions[0] * dimensions[1];
    volatile long long tmp = 0;
    auto start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < dimensions[0]; ++i) {
        for (size_t j = 0; j < dimensions[1]; ++j) {
            for (size_t k = 0; k < ragged_dimensions.size(); ++k) {
                for (size_t l = 0; l < ragged_dimensions[k]; ++l) {
                    tmp += mv.at({i, j, k, l});
                }
            }
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Ragged iteration time: " << duration.count() << " seconds\n";
    std::cout << "Ragged iteration result: " << tmp << std::endl;
    std::cout << "Total size: " << mv.total_size() << std::endl;
    std::cout << "Target: " << target << std::endl;
}

void benchmark_ragged_iterator_iteration(const std::array<size_t, 2>& dimensions, const std::vector<size_t>& ragged_dimensions) {
    RaggedMultiVector<int, 3> mv(dimensions, ragged_dimensions, 1);
    volatile long long tmp = 0;
    auto start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < dimensions[0]; ++i) {
        for (size_t j = 0; j < dimensions[1]; ++j) {
            for (size_t k = 0; k < ragged_dimensions.size(); ++k) {
                auto [begin, end] = mv.inner_iterators({i, j, k});
                for (auto it = begin; it != end; ++it) {
                    tmp += *it;
                }
            }
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Ragged iterator iteration time: " << duration.count() << " seconds\n";
    std::cout << "Ragged iterator iteration result: " << tmp << std::endl;
}

void naive_benchmark_ragged_construction(const std::array<size_t, 2>& dimensions, const std::vector<size_t>& ragged_dimensions) {
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<std::vector<std::vector<float>>> nv(dimensions[0], std::vector<std::vector<float>>(dimensions[1]));
    for (size_t i = 0; i < dimensions[0]; ++i) {
        for (size_t j = 0; j < dimensions[1]; ++j) {
            for (size_t d = 0; d < ragged_dimensions.size(); ++d) {
                nv[i][j].resize(ragged_dimensions[d], 0);
            }
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Naive ragged construction time: " << duration.count() << " seconds\n";
}

void naive_benchmark_ragged_element_assignment(const std::array<size_t, 2>& dimensions, const std::vector<size_t>& ragged_dimensions) {
    std::vector<std::vector<std::vector<std::vector<float>>>> nv;
    for (size_t i = 0; i < dimensions[0]; ++i) {
        nv.push_back(std::vector<std::vector<std::vector<float>>>());
        for (size_t j = 0; j < dimensions[1]; ++j) {
            nv[i].push_back(std::vector<std::vector<float>>());
            for (size_t k = 0; k < ragged_dimensions.size(); ++k) {
                nv[i][j].push_back(std::vector<float>());
                nv[i][j][k].resize(ragged_dimensions[k], 1);
            }
        }
    }

    auto start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < dimensions[0]; ++i) {
        for (size_t j = 0; j < dimensions[1]; ++j) {
            for (size_t k = 0; k < ragged_dimensions.size(); ++k) {
                for (size_t l = 0; l < ragged_dimensions[k]; ++l) {
                    nv[i][j][k][l] = i;
                }
            }
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Naive ragged element assignment time: " << duration.count() << " seconds\n";
}

void naive_benchmark_ragged_iteration(const std::array<size_t, 2>& dimensions, const std::vector<size_t>& ragged_dimensions) {
    std::vector<std::vector<std::vector<std::vector<int>>>> nv;
    for (size_t i = 0; i < dimensions[0]; ++i) {
        nv.push_back(std::vector<std::vector<std::vector<int>>>());
        for (size_t j = 0; j < dimensions[1]; ++j) {
            nv[i].push_back(std::vector<std::vector<int>>());
            for (size_t k = 0; k < ragged_dimensions.size(); ++k) {
                nv[i][j].push_back(std::vector<int>());
                nv[i][j][k].resize(ragged_dimensions[k], 1);
            }
        }
    }
    volatile long long tmp = 0;
    auto start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < dimensions[0]; ++i) {
        for (size_t j = 0; j < dimensions[1]; ++j) {
            for (size_t k = 0; k < ragged_dimensions.size(); ++k) {
                for (size_t l = 0; l < ragged_dimensions[k]; ++l) {
                    tmp += nv[i][j][k][l];
                }
            }
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Naive ragged iteration time: " << duration.count() << " seconds\n";
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <outer_dim1> <outer_dim2> <ragged_dim1,ragged_dim2,...>\n";
        return 1;
    }

    std::array<size_t, 2> dimensions = {static_cast<size_t>(std::stoi(argv[1])),
                                        static_cast<size_t>(std::stoi(argv[2]))};

    std::vector<size_t> ragged_dimensions;
    std::string arg(argv[3]);
    size_t pos = 0;
    while ((pos = arg.find(',')) != std::string::npos) {
        ragged_dimensions.push_back(static_cast<size_t>(std::stoi(arg.substr(0, pos))));
        arg.erase(0, pos + 1);
    }
    ragged_dimensions.push_back(static_cast<size_t>(std::stoi(arg)));

    benchmark_ragged_construction(dimensions, ragged_dimensions);
    benchmark_ragged_element_assignment(dimensions, ragged_dimensions);
    benchmark_ragged_iteration(dimensions, ragged_dimensions);
    benchmark_ragged_iterator_iteration(dimensions, ragged_dimensions);

    naive_benchmark_ragged_construction(dimensions, ragged_dimensions);
    naive_benchmark_ragged_element_assignment(dimensions, ragged_dimensions);
    naive_benchmark_ragged_iteration(dimensions, ragged_dimensions);

    return 0;
} 