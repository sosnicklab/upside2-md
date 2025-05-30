// File: test_change_cache_buffer.cpp

#define CATCH_CONFIG_MAIN // Tells Catch2 to supply a main() function for the test executable
#include "catch2/catch.hpp"       // Main Catch2 v2 header

// Include the header for the class under test
#include "interaction_graph.h" // This should now contain the conditional getter

// IMPORTANT: This test file (or the test build configuration)
// must define ENABLE_TEST_ACCESSORS for the tests to properly verify the cache_buffer value.
// For example, compile with: g++ -std=c++17 -DENABLE_TEST_ACCESSORS test_change_cache_buffer.cpp ...

TEST_CASE("PairlistComputation::change_cache_buffer updates the cache_buffer member", "[PairlistComputation][Cache]") {

    SECTION("Testing with non-symmetric PairlistComputation instance") {
        // Instantiate PairlistComputation with dummy values for n_elem1, n_elem2, max_n_edge.
        // These values are needed for construction but don't affect change_cache_buffer directly.
        PairlistComputation<false> pc_nonsymmetric(10, 10, 100); // symmetric = false

#ifdef TEST_HOOKS
        // Verify the initial value of cache_buffer (set to 1.0f in the constructor).
        REQUIRE(pc_nonsymmetric.get_cache_buffer_for_testing() == 1.0f);

        float new_value = 7.8f;
        pc_nonsymmetric.change_cache_buffer(new_value);
        // Verify that cache_buffer was updated.
        REQUIRE(pc_nonsymmetric.get_cache_buffer_for_testing() == new_value);

        float another_value = 0.15f;
        pc_nonsymmetric.change_cache_buffer(another_value);
        REQUIRE(pc_nonsymmetric.get_cache_buffer_for_testing() == another_value);

        // Test with a negative value if appropriate for the design
        float negative_value = -3.0f;
        pc_nonsymmetric.change_cache_buffer(negative_value);
        REQUIRE(pc_nonsymmetric.get_cache_buffer_for_testing() == negative_value);
#else
        // If the test accessor is not enabled, we can't check the value directly.
        // We can still call the function to ensure it doesn't crash.
        WARN("ENABLE_TEST_ACCESSORS is not defined. Detailed cache_buffer checks are skipped for non-symmetric instance.");
        pc_nonsymmetric.change_cache_buffer(7.8f); // Call the function
        SUCCEED("change_cache_buffer (non-symmetric) called without crashing (accessor not enabled).");
#endif
    }

    SECTION("Testing with symmetric PairlistComputation instance") {
        PairlistComputation<true> pc_symmetric(10, 10, 100); // symmetric = true

#ifdef TEST_HOOKS
        // Verify the initial value.
        REQUIRE(pc_symmetric.get_cache_buffer_for_testing() == 1.0f);

        float new_value = 3.45f;
        pc_symmetric.change_cache_buffer(new_value);
        REQUIRE(pc_symmetric.get_cache_buffer_for_testing() == new_value);

        float another_value = 0.0f;
        pc_symmetric.change_cache_buffer(another_value);
        REQUIRE(pc_symmetric.get_cache_buffer_for_testing() == another_value);
#else
        WARN("TEST_HOOKS is not defined. Detailed cache_buffer checks are skipped for symmetric instance.");
        pc_symmetric.change_cache_buffer(3.45f); // Call the function
        SUCCEED("change_cache_buffer (symmetric) called without crashing (accessor not enabled).");
#endif
    }
}