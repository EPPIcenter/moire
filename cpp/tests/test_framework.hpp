#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <chrono>
#include <functional>
#include <exception>
#include <memory>

namespace moire_test {

// Test result structure
struct TestResult {
    std::string name;
    bool passed;
    std::string error_message;
    double duration_ms;
    
    TestResult(const std::string& n, bool p, const std::string& msg = "", double dur = 0.0)
        : name(n), passed(p), error_message(msg), duration_ms(dur) {}
};

// Test suite class
class TestSuite {
private:
    std::string suite_name;
    std::vector<TestResult> results;
    std::vector<std::function<void()>> tests;
    
public:
    TestSuite(const std::string& name) : suite_name(name) {}
    
    // Add a test to the suite
    void add_test(const std::string& test_name, std::function<void()> test_func) {
        tests.push_back([this, test_name, test_func]() {
            auto start = std::chrono::high_resolution_clock::now();
            try {
                test_func();
                auto end = std::chrono::high_resolution_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0;
                results.emplace_back(test_name, true, "", duration);
                std::cout << "✅ " << test_name << " passed (" << duration << " ms)" << std::endl;
            } catch (const std::exception& e) {
                auto end = std::chrono::high_resolution_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0;
                results.emplace_back(test_name, false, e.what(), duration);
                std::cout << "❌ " << test_name << " failed: " << e.what() << " (" << duration << " ms)" << std::endl;
            } catch (...) {
                auto end = std::chrono::high_resolution_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0;
                results.emplace_back(test_name, false, "Unknown error", duration);
                std::cout << "❌ " << test_name << " failed: Unknown error (" << duration << " ms)" << std::endl;
            }
        });
    }
    
    // Run all tests in the suite
    void run_tests() {
        std::cout << "\n🧪 Running test suite: " << suite_name << std::endl;
        std::cout << "=" << std::string(suite_name.length() + 20, '=') << std::endl;
        
        for (auto& test : tests) {
            test();
        }
        
        print_summary();
    }
    
    // Print test summary
    void print_summary() {
        int passed = 0;
        int failed = 0;
        double total_time = 0.0;
        
        for (const auto& result : results) {
            if (result.passed) {
                passed++;
            } else {
                failed++;
            }
            total_time += result.duration_ms;
        }
        
        std::cout << "\n📊 Test Summary for " << suite_name << ":" << std::endl;
        std::cout << "   Passed: " << passed << std::endl;
        std::cout << "   Failed: " << failed << std::endl;
        std::cout << "   Total:  " << (passed + failed) << std::endl;
        std::cout << "   Time:   " << total_time << " ms" << std::endl;
        
        if (failed > 0) {
            std::cout << "\n❌ Failed tests:" << std::endl;
            for (const auto& result : results) {
                if (!result.passed) {
                    std::cout << "   - " << result.name << ": " << result.error_message << std::endl;
                }
            }
        }
    }
    
    // Get test results
    const std::vector<TestResult>& get_results() const { return results; }
    
    // Check if all tests passed
    bool all_passed() const {
        return std::all_of(results.begin(), results.end(), 
                          [](const TestResult& r) { return r.passed; });
    }
};

// Global test registry
class TestRegistry {
private:
    static std::vector<std::unique_ptr<TestSuite>> suites;
    
public:
    static void register_suite(std::unique_ptr<TestSuite> suite) {
        suites.push_back(std::move(suite));
    }
    
    static void run_all_tests() {
        std::cout << "🚀 Starting C++ Test Suite" << std::endl;
        std::cout << "=========================" << std::endl;
        
        int total_passed = 0;
        int total_failed = 0;
        
        for (auto& suite : suites) {
            suite->run_tests();
            
            for (const auto& result : suite->get_results()) {
                if (result.passed) {
                    total_passed++;
                } else {
                    total_failed++;
                }
            }
        }
        
        std::cout << "\n🎯 Overall Test Summary:" << std::endl;
        std::cout << "   Passed: " << total_passed << std::endl;
        std::cout << "   Failed: " << total_failed << std::endl;
        std::cout << "   Total:  " << (total_passed + total_failed) << std::endl;
        
        if (total_failed > 0) {
            std::cout << "\n❌ Some tests failed!" << std::endl;
            exit(1);
        } else {
            std::cout << "\n✅ All tests passed!" << std::endl;
        }
    }
};

// Test assertion macros
#define ASSERT_TRUE(condition) \
    do { \
        if (!(condition)) { \
            throw std::runtime_error("Assertion failed: " #condition " is not true"); \
        } \
    } while(0)

#define ASSERT_FALSE(condition) \
    do { \
        if (condition) { \
            throw std::runtime_error("Assertion failed: " #condition " is not false"); \
        } \
    } while(0)

#define ASSERT_EQ(expected, actual) \
    do { \
        if ((expected) != (actual)) { \
            throw std::runtime_error("Assertion failed: " #expected " == " #actual); \
        } \
    } while(0)

#define ASSERT_NE(expected, actual) \
    do { \
        if ((expected) == (actual)) { \
            throw std::runtime_error("Assertion failed: " #expected " == " #actual); \
        } \
    } while(0)

#define ASSERT_LT(left, right) \
    do { \
        if (!((left) < (right))) { \
            throw std::runtime_error("Assertion failed: " #left " >= " #right); \
        } \
    } while(0)

#define ASSERT_LE(left, right) \
    do { \
        if (!((left) <= (right))) { \
            throw std::runtime_error("Assertion failed: " #left " > " #right); \
        } \
    } while(0)

#define ASSERT_GT(left, right) \
    do { \
        if (!((left) > (right))) { \
            throw std::runtime_error("Assertion failed: " #left " <= " #right); \
        } \
    } while(0)

#define ASSERT_GE(left, right) \
    do { \
        if (!((left) >= (right))) { \
            throw std::runtime_error("Assertion failed: " #left " < " #right); \
        } \
    } while(0)

#define ASSERT_THROWS(statement) \
    do { \
        bool threw = false; \
        try { \
            statement; \
        } catch (...) { \
            threw = true; \
        } \
        if (!threw) { \
            throw std::runtime_error("Assertion failed: " #statement " did not throw"); \
        } \
    } while(0)

#define ASSERT_NO_THROW(statement) \
    do { \
        try { \
            statement; \
        } catch (const std::exception& e) { \
            throw std::runtime_error("Assertion failed: " #statement " threw: " + std::string(e.what())); \
        } catch (...) { \
            throw std::runtime_error("Assertion failed: " #statement " threw unknown exception"); \
        } \
    } while(0)

// Test suite registration macro
#define REGISTER_TEST_SUITE(name) \
    static auto test_suite_##name = std::make_unique<moire_test::TestSuite>(#name); \
    moire_test::TestRegistry::register_suite(std::move(test_suite_##name))

// Test registration macro
#define ADD_TEST(suite, name, func) \
    do { \
        static bool registered = false; \
        if (!registered) { \
            suite->add_test(name, func); \
            registered = true; \
        } \
    } while(0)

} // namespace moire_test

// Define static member
std::vector<std::unique_ptr<moire_test::TestSuite>> moire_test::TestRegistry::suites;

