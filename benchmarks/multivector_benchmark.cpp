#include "../src/multivector.h"
#include <iostream>
#include <chrono>
#include <vector>

void benchmark_construction(const std::array<size_t, 3>& dimensions) {
    auto start = std::chrono::high_resolution_clock::now();
    MultiVector<float, 3> mv(dimensions);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Construction time: " << duration.count() << " seconds\n";
}

void benchmark_element_assignment(const std::array<size_t, 3>& dimensions) {
    MultiVector<float, 3> mv(dimensions);
    auto start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < dimensions[0]; ++i) {
        for (size_t j = 0; j < dimensions[1]; ++j) {
            for (size_t k = 0; k < dimensions[2]; ++k) {
                mv.at({i, j, k}) = i;
            }
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Element assignment time: " << duration.count() << " seconds\n";
}

void benchmark_read_for_loop_iteration(const std::array<size_t, 3>& dimensions) {
    MultiVector<float, 3> mv(dimensions);
    volatile float tmp;
    auto start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < dimensions[0]; ++i) {
        for (size_t j = 0; j < dimensions[1]; ++j) {
            for (size_t k = 0; k < dimensions[2]; ++k) {
                tmp = mv.at({i, j, k});
            }
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Read for loop iteration time: " << duration.count() << " seconds\n";
}

void benchmark_read_for_loop_swapped_iteration(const std::array<size_t, 3>& dimensions) {
    MultiVector<float, 3> mv(dimensions);
    volatile float tmp;
    auto start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < dimensions[0]; ++i) {
        for (size_t j = 0; j < dimensions[2]; ++j) {
            for (size_t k = 0; k < dimensions[1]; ++k) {
                tmp = mv.at({i, k, j});
            }
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Read for loop swapped iteration time: " << duration.count() << " seconds\n";
}

void benchmark_read_for_loop_inverted_iteration(const std::array<size_t, 3>& dimensions) {
    MultiVector<float, 3> mv(dimensions);
    volatile float tmp;
    auto start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < dimensions[2]; ++i) {
        for (size_t j = 0; j < dimensions[1]; ++j) {
            for (size_t k = 0; k < dimensions[0]; ++k) {
                tmp = mv.at({k, j, i});
            }
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Read for loop inverted iteration time: " << duration.count() << " seconds\n";
}


void benchmark_iterator_iteration(const std::array<size_t, 3>& dimensions) {
    MultiVector<float, 3> mv(dimensions);
    auto start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < dimensions[0]; ++i) {
        for (size_t j = 0; j < dimensions[1]; ++j) {
            auto it = mv.inner_begin({i, j});
            auto end = mv.inner_end({i, j});
            for (; it != end; ++it) {
                *it = i;
            }
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Iterator iteration time: " << duration.count() << " seconds\n";
}


void naive_benchmark_construction(const std::array<size_t, 3>& dimensions) {
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<std::vector<std::vector<float>>> nv(dimensions[0], std::vector<std::vector<float>>(dimensions[1], std::vector<float>(dimensions[2], 0)));
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Naive construction time: " << duration.count() << " seconds\n";
}

void naive_benchmark_element_assignment(const std::array<size_t, 3>& dimensions) {
    std::vector<std::vector<std::vector<float>>> nv(dimensions[0], std::vector<std::vector<float>>(dimensions[1], std::vector<float>(dimensions[2])));
    auto start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < dimensions[0]; ++i) {
        for (size_t j = 0; j < dimensions[1]; ++j) {
            for (size_t k = 0; k < dimensions[2]; ++k) {
                nv[i][j][k] = i;
            }
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Naive element assignment time: " << duration.count() << " seconds\n";
}

void naive_benchmark_read_for_loop_iteration(const std::array<size_t, 3>& dimensions) {
    std::vector<std::vector<std::vector<float>>> nv(dimensions[0], std::vector<std::vector<float>>(dimensions[1], std::vector<float>(dimensions[2])));
    volatile float tmp;
    auto start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < dimensions[0]; ++i) {
        for (size_t j = 0; j < dimensions[1]; ++j) {
            for (size_t k = 0; k < dimensions[2]; ++k) {
                tmp = nv[i][j][k];
            }
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Naive read for loop iteration time: " << duration.count() << " seconds\n";
}

void naive_benchmark_read_for_loop_swapped_iteration(const std::array<size_t, 3>& dimensions) {
    std::vector<std::vector<std::vector<float>>> nv(dimensions[0], std::vector<std::vector<float>>(dimensions[1], std::vector<float>(dimensions[2])));
    volatile float tmp;
    auto start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < dimensions[0]; ++i) {
        for (size_t j = 0; j < dimensions[2]; ++j) {
            for (size_t k = 0; k < dimensions[1]; ++k) {
                tmp = nv[i][k][j];
            }
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Naive read for loop swapped iteration time: " << duration.count() << " seconds\n";
}

void naive_benchmark_read_for_loop_inverted_iteration(const std::array<size_t, 3>& dimensions) {
    std::vector<std::vector<std::vector<float>>> nv(dimensions[0], std::vector<std::vector<float>>(dimensions[1], std::vector<float>(dimensions[2])));
    volatile float tmp;
    auto start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < dimensions[2]; ++i) {
        for (size_t j = 0; j < dimensions[1]; ++j) {
            for (size_t k = 0; k < dimensions[0]; ++k) {
                tmp = nv[k][j][i];
            }
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Naive read for loop inverted iteration time: " << duration.count() << " seconds\n";
}

void naive_benchmark_iterator_iteration(const std::array<size_t, 3>& dimensions) {
    std::vector<std::vector<std::vector<float>>> nv(dimensions[0], std::vector<std::vector<float>>(dimensions[1], std::vector<float>(dimensions[2])));
    auto start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < dimensions[0]; ++i) {
        for (size_t j = 0; j < dimensions[1]; ++j) {
            auto it = nv[i][j].begin();
            auto end = nv[i][j].end();
            for (; it != end; ++it) {
                *it = i;
            }
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Naive iterator iteration time: " << duration.count() << " seconds\n";
}

int main(int argc, char* argv[]) {
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <dim1> <dim2> <dim3>\n";
        return 1;
    }

    std::array<size_t, 3> dimensions = {static_cast<size_t>(std::stoi(argv[1])),
                                      static_cast<size_t>(std::stoi(argv[2])),
                                      static_cast<size_t>(std::stoi(argv[3]))};

    benchmark_construction(dimensions);
    benchmark_element_assignment(dimensions);
    benchmark_read_for_loop_iteration(dimensions);
    benchmark_read_for_loop_swapped_iteration(dimensions);
    benchmark_read_for_loop_inverted_iteration(dimensions);
    benchmark_iterator_iteration(dimensions);
    naive_benchmark_construction(dimensions);
    naive_benchmark_element_assignment(dimensions);
    naive_benchmark_read_for_loop_iteration(dimensions);
    naive_benchmark_read_for_loop_swapped_iteration(dimensions);
    naive_benchmark_read_for_loop_inverted_iteration(dimensions);
    naive_benchmark_iterator_iteration(dimensions);
    return 0;
} 