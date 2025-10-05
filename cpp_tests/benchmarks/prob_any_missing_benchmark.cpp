#include "../src/prob_any_missing.h"
#include <iostream>
#include <chrono>
#include <vector>
#include <sstream>

void benchmarkProbAnyMissing(const std::vector<float>& eventProbs, int numEvents) {
    probAnyMissingFunctor functor;
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<float> result;
    for (int k = 0; k < 1000000; ++k) {
        result = functor.vectorized(eventProbs, numEvents);
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    std::cout << "Result = ";
    for (size_t l = 0; l < result.size(); ++l) {
        std::cout << result[l] << " ";
    }
    std::cout << ", Time = " << elapsed.count() << " seconds" << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <eventProb1> <eventProb2> ... <numEvents>" << std::endl;
        return 1;
    }

    std::vector<float> eventProbs;
    for (int i = 1; i < argc - 1; ++i) {
        float prob;
        std::istringstream(argv[i]) >> prob;
        eventProbs.push_back(prob);
    }

    int numEvents;
    std::istringstream(argv[argc - 1]) >> numEvents;

    std::cout << "Benchmarking probAnyMissingFunctor..." << std::endl;
    benchmarkProbAnyMissing(eventProbs, numEvents);
    return 0;
} 