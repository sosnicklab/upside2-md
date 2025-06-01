#include <iostream>
#include <string>
#include <vector>
#include <stdexcept> // For std::stof, std::invalid_argument, std::out_of_range
#include <iomanip>   // For std::fixed, std::setprecision

#include "interaction_graph.h"

// Dummy function for acceptable_id_pair_t signature
// This function always returns a mask indicating all elements are acceptable for testing purposes.
static Int4 dummy_acceptable_id_pair(const Int4& id1_vec, const Int4& id2_vec) {
    return Int4(-1); // All bits set to 1, effectively true for all lanes
}

template <bool symmetric_flag>
class PairlistComputationTestBridge {
public:
    static float get_cache_buffer(const PairlistComputation<symmetric_flag>& pc) {return pc.cache_buffer;}

    // Static method to call the protected ensure_cache_valid
    // Note: ensure_cache_valid takes the function pointer as a non-type template parameter.
    static void call_ensure_cache_valid(
        PairlistComputation<symmetric_flag>& pc,
        float cutoff,
        const float* aligned_pos1, const int pos1_stride, int* id1,
        const float* aligned_pos2, const int pos2_stride, int* id2)
    {
        // Directly use the dummy_acceptable_id_pair function as the template argument
        pc.template ensure_cache_valid<dummy_acceptable_id_pair>(
            cutoff,
            aligned_pos1, pos1_stride, id1,
            aligned_pos2, pos2_stride, id2);
    }

    // Accessors for protected members
    static bool get_cache_valid(const PairlistComputation<symmetric_flag>& pc) { return pc.cache_valid; }
    static float get_cache_cutoff(const PairlistComputation<symmetric_flag>& pc) { return pc.cache_cutoff; }
    static int get_cache_n_edge(const PairlistComputation<symmetric_flag>& pc) { return pc.cache_n_edge; }
};

float stringToFloat(const char* str) {
    try {
        return std::stof(str);
    } catch (const std::invalid_argument& ia) {
        std::cerr << "Error: Invalid argument for float conversion: '" << str << "'" << std::endl;
        exit(EXIT_FAILURE);
    } catch (const std::out_of_range& oor) {
        std::cerr << "Error: Float value out of range: '" << str << "'" << std::endl;
        exit(EXIT_FAILURE);
    }
}

int stringToInt(const char* str) {
    try {
        return std::stoi(str);
    } catch (const std::invalid_argument& ia) {
        std::cerr << "Error: Invalid argument for integer conversion: '" << str << "'" << std::endl;
        exit(EXIT_FAILURE);
    } catch (const std::out_of_range& oor) {
        std::cerr << "Error: Integer value out of range: '" << str << "'" << std::endl;
        exit(EXIT_FAILURE);
    }
}

int main(int argc, char* argv[]) {
    // Disable synchronization with C stdio for potentially faster I/O, if needed.
    // std::ios_base::sync_with_stdio(false);
    // std::cin.tie(NULL);

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <method_name> [args...]" << std::endl;
        std::cerr << "Supported methods:" << std::endl;
        std::cerr << "  change_cache_buffer <new_buffer_value>" << std::endl;
        std::cerr << "  test_ensure_cache_valid <n_elem> <max_n_edge> <pos_stride> <cutoff_value>" << std::endl; // Updated usage
        return 1;
    }

    std::string methodName = argv[1];

    if (methodName == "change_cache_buffer") {
        if (argc != 3) {
            std::cerr << "Usage: " << argv[0] << " change_cache_buffer <new_buffer_value>" << std::endl;
            return 1;
        }

        float new_buffer_value = stringToFloat(argv[2]);

        // Test both symmetric and non-symmetric versions
        PairlistComputation<true> pc_symmetric(10, 10, 45);
        pc_symmetric.change_cache_buffer(new_buffer_value);
        float symmetric_buffer = PairlistComputationTestBridge<true>::get_cache_buffer(pc_symmetric);

        PairlistComputation<false> pc_non_symmetric(10, 10, 45);
        pc_non_symmetric.change_cache_buffer(new_buffer_value);
        float non_symmetric_buffer = PairlistComputationTestBridge<false>::get_cache_buffer(pc_non_symmetric);
        
        // Output both results
        std::cout << "symmetric: " << std::fixed << std::setprecision(8) << symmetric_buffer << std::endl;
        std::cout << "non_symmetric: " << std::fixed << std::setprecision(8) << non_symmetric_buffer << std::endl;

    }
    else if (methodName == "test_ensure_cache_valid") {
        if (argc != 6) {
            std::cerr << "Usage: " << argv[0] << " test_ensure_cache_valid <n_elem> <max_n_edge> <pos_stride> <cutoff_value>" << std::endl;
            return 1;
        }

        int n_elem = stringToInt(argv[2]);
        int max_n_edge = stringToInt(argv[3]);
        int pos_stride = stringToInt(argv[4]);
        float cutoff_value = stringToFloat(argv[5]);

        int num_pos_floats = n_elem * pos_stride;

        std::vector<float> pos1_data(num_pos_floats);
        std::vector<float> pos2_data(num_pos_floats);
        std::vector<int> id1_data(n_elem);
        std::vector<int> id2_data(n_elem);

        // Initialize dummy position and ID data
        for (int i = 0; i < n_elem; ++i) {
            for (int d = 0; d < pos_stride; ++d) {
                pos1_data[i * pos_stride + d] = (float)i * 1.0f + (float)d * 0.1f;
                pos2_data[i * pos_stride + d] = (float)i * 1.0f + (float)d * 0.1f + 0.05f;
            }
            id1_data[i] = i;
            id2_data[i] = i;
        }

        std::cout << "--- Testing Symmetric PairlistComputation ---" << std::endl;
        PairlistComputation<true> pc_symmetric(n_elem, n_elem, max_n_edge);
        // Initially, cache_valid should be false
        std::cout << "symmetric initial cache_valid: " << std::boolalpha << PairlistComputationTestBridge<true>::get_cache_valid(pc_symmetric) << std::endl;
        
        PairlistComputationTestBridge<true>::call_ensure_cache_valid(
            pc_symmetric,
            cutoff_value,
            pos1_data.data(), pos_stride, id1_data.data(),
            pos1_data.data(), pos_stride, id1_data.data() // For symmetric, pos2/id2 are pos1/id1
        );
        std::cout << "symmetric after ensure_cache_valid:" << std::endl;
        std::cout << "  cache_valid: " << std::boolalpha << PairlistComputationTestBridge<true>::get_cache_valid(pc_symmetric) << std::endl;
        std::cout << "  cache_cutoff: " << std::fixed << std::setprecision(8) << PairlistComputationTestBridge<true>::get_cache_cutoff(pc_symmetric) << std::endl;
        std::cout << "  cache_n_edge: " << PairlistComputationTestBridge<true>::get_cache_n_edge(pc_symmetric) << std::endl;
        std::cout << std::endl;


        std::cout << "--- Testing Non-Symmetric PairlistComputation ---" << std::endl;
        PairlistComputation<false> pc_non_symmetric(n_elem, n_elem, max_n_edge);
        // Initially, cache_valid should be false
        std::cout << "non_symmetric initial cache_valid: " << std::boolalpha << PairlistComputationTestBridge<false>::get_cache_valid(pc_non_symmetric) << std::endl;

        PairlistComputationTestBridge<false>::call_ensure_cache_valid(
            pc_non_symmetric,
            cutoff_value,
            pos1_data.data(), pos_stride, id1_data.data(),
            pos2_data.data(), pos_stride, id2_data.data()
        );
        std::cout << "non_symmetric after ensure_cache_valid:" << std::endl;
        std::cout << "  cache_valid: " << std::boolalpha << PairlistComputationTestBridge<false>::get_cache_valid(pc_non_symmetric) << std::endl;
        std::cout << "  cache_cutoff: " << std::fixed << std::setprecision(8) << PairlistComputationTestBridge<false>::get_cache_cutoff(pc_non_symmetric) << std::endl;
        std::cout << "  cache_n_edge: " << PairlistComputationTestBridge<false>::get_cache_n_edge(pc_non_symmetric) << std::endl;

    }
    else {
        std::cerr << "Error: Unknown method name '" << methodName << "'" << std::endl;
        std::cerr << "Supported methods:" << std::endl;
        std::cerr << "  change_cache_buffer <new_buffer_value>" << std::endl;
        std::cerr << "  test_ensure_cache_valid <n_elem> <max_n_edge> <pos_stride> <cutoff_value>" << std::endl; // Updated usage
        return 1;
    }

    return 0;
}