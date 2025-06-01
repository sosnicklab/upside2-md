#include <iostream>
#include <string>
#include <vector>
#include <stdexcept> // For std::stof, std::invalid_argument, std::out_of_range
#include <iomanip>   // For std::fixed, std::setprecision

#include "interaction_graph.h" 

template <bool symmetric_flag>
class PairlistComputationTestBridge {
public:
    static float get_cache_buffer(const PairlistComputation<symmetric_flag>& pc) {return pc.cache_buffer;}
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
        // Add other methods here as they are implemented
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
    // Add else-if blocks here for other methods:
    // else if (methodName == "some_other_method") {
    //    // ... logic for testing some_other_method ...
    // } 
    else {
        std::cerr << "Error: Unknown method name '" << methodName << "'" << std::endl;
        std::cerr << "Supported methods:" << std::endl;
        std::cerr << "  change_cache_buffer <new_buffer_value>" << std::endl;
        return 1;
    }

    return 0;
}