#include <iostream>
#include <string>
#include <vector>
#include <stdexcept> // For std::stof, std::invalid_argument, std::out_of_range
#include <iomanip>   // For std::fixed, std::setprecision
#include <sstream>   // For std::stringstream
#include <algorithm> // For std::transform, std::max, std::fill
#include <numeric>   // For std::iota (for generating id arrays)
#include <limits>    // For std::numeric_limits
#include <initializer_list> // For std::max({})

#include "interaction_graph.h" // Your actual header

// Dummy function for acceptable_id_pair_t signature
// This function always returns a mask indicating all elements are acceptable for testing purposes.
static Int4 dummy_acceptable_id_pair(const Int4& id1_vec, const Int4& id2_vec) {
    return Int4(-1); // All bits set to 1, effectively true for all lanes
}

template <bool symmetric_flag>
class PairlistComputationTestBridge {
public:
    // Accessors for protected members (getters)
    static float get_cache_buffer(const PairlistComputation<symmetric_flag>& pc) { return pc.cache_buffer; }
    static bool get_cache_valid(const PairlistComputation<symmetric_flag>& pc) { return pc.cache_valid; }
    static float get_cache_cutoff(const PairlistComputation<symmetric_flag>& pc) { return pc.cache_cutoff; }
    static int get_cache_n_edge(const PairlistComputation<symmetric_flag>& pc) { return pc.cache_n_edge; }

    // Accessors for public members (for consistency in testing)
    static int get_n_edge(const PairlistComputation<symmetric_flag>& pc) { return pc.n_edge; }
    static const int32_t* get_edge_indices1_ptr(const PairlistComputation<symmetric_flag>& pc) { return pc.edge_indices1.get(); }
    static const int32_t* get_edge_indices2_ptr(const PairlistComputation<symmetric_flag>& pc) { return pc.edge_indices2.get(); }
    static const int32_t* get_edge_id1_ptr(const PairlistComputation<symmetric_flag>& pc) { return pc.edge_id1.get(); }
    static const int32_t* get_edge_id2_ptr(const PairlistComputation<symmetric_flag>& pc) { return pc.edge_id2.get(); }


    // Accessors for protected members (setters)
    static void set_cache_valid(PairlistComputation<symmetric_flag>& pc, bool valid) { pc.cache_valid = valid; }
    static void set_cache_cutoff(PairlistComputation<symmetric_flag>& pc, float cutoff) { pc.cache_cutoff = cutoff; }
    static void set_cache_n_edge(PairlistComputation<symmetric_flag>& pc, int n_edge) { pc.cache_n_edge = n_edge; }
    
    // Static method to call the protected ensure_cache_valid (still useful if we wanted to test it in isolation)
    static void call_ensure_cache_valid(
        PairlistComputation<symmetric_flag>& pc,
        float cutoff,
        const float* aligned_pos1, const int pos1_stride, int* id1,
        const float* aligned_pos2, const int pos2_stride, int* id2)
    {
        pc.template ensure_cache_valid<dummy_acceptable_id_pair>(
            cutoff,
            aligned_pos1, pos1_stride, id1,
            aligned_pos2, pos2_stride, id2);
    }
};

// Helper to convert string to boolean
bool stringToBool(const std::string& str_in) {
    std::string str = str_in;
    std::transform(str.begin(), str.end(), str.begin(), 
                   [](unsigned char c){ return std::tolower(c); });
    size_t at_pos = str.find('@');
    if (at_pos != std::string::npos) {
        str = str.substr(0, at_pos);
    }

    if (str == "true" || str == "1") {
        return true;
    } else if (str == "false" || str == "0") {
        return false;
    }
    std::cerr << "Error: Invalid argument for boolean conversion: '" << str_in << "'" << std::endl;
    exit(EXIT_FAILURE);
}

float stringToFloat(const char* str_cstr) {
    std::string str(str_cstr);
    size_t at_pos = str.find('@');
    if (at_pos != std::string::npos) {
        str = str.substr(0, at_pos);
    }
    try {
        if (str.empty()) {
            std::cerr << "Error: Empty string provided for float conversion." << std::endl;
            exit(EXIT_FAILURE);
        }
        std::string lower_str = str;
        std::transform(lower_str.begin(), lower_str.end(), lower_str.begin(), ::tolower);
        if (lower_str == "inf" || lower_str == "+inf" || lower_str == "infinity" || lower_str == "+infinity") {
            return std::numeric_limits<float>::infinity();
        } else if (lower_str == "-inf" || lower_str == "-infinity") {
            return -std::numeric_limits<float>::infinity();
        } else if (lower_str == "nan") {
            return std::numeric_limits<float>::quiet_NaN();
        }
        return std::stof(str);
    } catch (const std::invalid_argument& ia) {
        std::cerr << "Error: Invalid argument for float conversion: '" << str_cstr << "' (parsed as '" << str << "')" << std::endl;
        exit(EXIT_FAILURE);
    } catch (const std::out_of_range& oor) {
        std::cerr << "Error: Float value out of range: '" << str_cstr << "' (parsed as '" << str << "')" << std::endl;
        exit(EXIT_FAILURE);
    }
}

int stringToInt(const char* str_cstr) {
    std::string str(str_cstr);
    size_t at_pos = str.find('@');
    if (at_pos != std::string::npos) {
        str = str.substr(0, at_pos);
    }
    try {
        if (str.empty()) { 
            std::cerr << "Error: Empty string provided for int conversion." << std::endl;
            exit(EXIT_FAILURE);
        }
        return std::stoi(str);
    } catch (const std::invalid_argument& ia) {
        std::cerr << "Error: Invalid argument for integer conversion: '" << str_cstr << "' (parsed as '" << str << "')" << std::endl;
        exit(EXIT_FAILURE);
    } catch (const std::out_of_range& oor) {
        std::cerr << "Error: Integer value out of range: '" << str_cstr << "' (parsed as '" << str << "')" << std::endl;
        exit(EXIT_FAILURE);
    }
}

std::vector<float> parseFloatArrayString(std::string s) {
    std::vector<float> result;
    if (!s.empty() && s.front() == '{') s.erase(0, 1);
    if (!s.empty() && s.back() == '}') s.pop_back();

    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, ',')) {
        item.erase(0, item.find_first_not_of(" \t\n\r\f\v"));
        item.erase(item.find_last_not_of(" \t\n\r\f\v") + 1);
        if (!item.empty()) {
            result.push_back(stringToFloat(item.c_str()));
        }
    }
    return result;
}

int round_up_util(int val, int multiple) {
    if (val < 0) val = 0; 
    if (multiple == 0) return val; 
    return ((val + multiple - 1) / multiple) * multiple;
}

// Helper to print an integer array for test output
void print_int_array(const std::string& name, const int32_t* arr, int size) {
    std::cout << name << ": [";
    for (int k = 0; k < size; ++k) {
        std::cout << arr[k] << (k == size - 1 ? "" : ", ");
    }
    std::cout << "]" << std::endl;
}


template<bool symmetric_flag>
void run_ensure_cache_valid_test(
    int n_elem1, int n_elem2, int max_n_edge_arg,
    bool initial_cache_valid, float initial_cache_buffer, float initial_cache_cutoff, int initial_cache_n_edge,
    float call_cutoff, int pos1_stride, int pos2_stride,
    std::vector<float>& aligned_pos1_data_arg, std::vector<float>& aligned_pos2_data_arg, 
    std::vector<int>& id1_data_arg, std::vector<int>& id2_data_arg) 
{
    PairlistComputation<symmetric_flag> pc(n_elem1, n_elem2, max_n_edge_arg);

    PairlistComputationTestBridge<symmetric_flag>::set_cache_valid(pc, initial_cache_valid);
    pc.change_cache_buffer(initial_cache_buffer); 
    PairlistComputationTestBridge<symmetric_flag>::set_cache_cutoff(pc, initial_cache_cutoff);
    PairlistComputationTestBridge<symmetric_flag>::set_cache_n_edge(pc, initial_cache_n_edge);
    
    PairlistComputationTestBridge<symmetric_flag>::call_ensure_cache_valid(
        pc,
        call_cutoff,
        aligned_pos1_data_arg.empty() ? nullptr : aligned_pos1_data_arg.data(), pos1_stride, 
        id1_data_arg.empty() ? nullptr : id1_data_arg.data(), 
        aligned_pos2_data_arg.empty() ? nullptr : aligned_pos2_data_arg.data(), pos2_stride, 
        id2_data_arg.empty() ? nullptr : id2_data_arg.data()   
    );

    std::cout << "final_cache_valid: " << std::boolalpha << PairlistComputationTestBridge<symmetric_flag>::get_cache_valid(pc) << std::endl;
    std::cout << "final_cache_cutoff: " << PairlistComputationTestBridge<symmetric_flag>::get_cache_cutoff(pc) << std::endl;
    std::cout << "final_cache_n_edge: " << PairlistComputationTestBridge<symmetric_flag>::get_cache_n_edge(pc) << std::endl; // This is pc.cache_n_edge
    std::cout << "final_cache_buffer: " << PairlistComputationTestBridge<symmetric_flag>::get_cache_buffer(pc) << std::endl;
}

template<bool symmetric_flag>
void run_find_edges_test(
    int n_elem1, int n_elem2, int max_n_edge_arg,
    bool initial_cache_valid, float initial_cache_buffer, float initial_cache_cutoff, int initial_cache_n_edge,
    float call_cutoff, int pos1_stride, int pos2_stride,
    std::vector<float>& aligned_pos1_data_arg, std::vector<float>& aligned_pos2_data_arg, 
    std::vector<int>& id1_data_arg, std::vector<int>& id2_data_arg)
{
    PairlistComputation<symmetric_flag> pc(n_elem1, n_elem2, max_n_edge_arg);

    // Set initial state for the PairlistComputation object (as if ensure_cache_valid was called or state was set)
    PairlistComputationTestBridge<symmetric_flag>::set_cache_valid(pc, initial_cache_valid);
    pc.change_cache_buffer(initial_cache_buffer); 
    PairlistComputationTestBridge<symmetric_flag>::set_cache_cutoff(pc, initial_cache_cutoff);
    PairlistComputationTestBridge<symmetric_flag>::set_cache_n_edge(pc, initial_cache_n_edge);

    // Call the public find_edges method
    pc.template find_edges<dummy_acceptable_id_pair>(
        call_cutoff,
        aligned_pos1_data_arg.empty() ? nullptr : aligned_pos1_data_arg.data(), pos1_stride,
        id1_data_arg.empty() ? nullptr : id1_data_arg.data(),
        aligned_pos2_data_arg.empty() ? nullptr : aligned_pos2_data_arg.data(), pos2_stride,
        id2_data_arg.empty() ? nullptr : id2_data_arg.data()
    );

    // Output results
    std::cout << "final_n_edge: " << PairlistComputationTestBridge<symmetric_flag>::get_n_edge(pc) << std::endl;
    
    // Output internal cache state (affected by ensure_cache_valid call within find_edges)
    std::cout << "final_cache_valid: " << std::boolalpha << PairlistComputationTestBridge<symmetric_flag>::get_cache_valid(pc) << std::endl;
    std::cout << "final_cache_cutoff: " << PairlistComputationTestBridge<symmetric_flag>::get_cache_cutoff(pc) << std::endl;
    std::cout << "final_internal_cache_n_edge: " << PairlistComputationTestBridge<symmetric_flag>::get_cache_n_edge(pc) << std::endl;
    std::cout << "final_cache_buffer: " << PairlistComputationTestBridge<symmetric_flag>::get_cache_buffer(pc) << std::endl;

    // Output edge arrays
    print_int_array("edge_indices1", PairlistComputationTestBridge<symmetric_flag>::get_edge_indices1_ptr(pc), pc.n_edge);
    print_int_array("edge_indices2", PairlistComputationTestBridge<symmetric_flag>::get_edge_indices2_ptr(pc), pc.n_edge);
    print_int_array("edge_id1", PairlistComputationTestBridge<symmetric_flag>::get_edge_id1_ptr(pc), pc.n_edge);
    print_int_array("edge_id2", PairlistComputationTestBridge<symmetric_flag>::get_edge_id2_ptr(pc), pc.n_edge);
}


int main(int argc, char* argv[]) {
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(NULL);
    std::cout << std::fixed << std::setprecision(8);

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <method_name> [args...]" << std::endl;
        std::cerr << "Supported methods:" << std::endl;
        std::cerr << "  change_cache_buffer <new_buffer_value>" << std::endl;
        std::cerr << "  test_ensure_cache_valid <symmetric_type_str> <n_elem1_str> <n_elem2_str> \\" << std::endl;
        std::cerr << "    <initial_cache_valid_str> <initial_cache_buffer_str> <initial_cache_cutoff_str> <initial_cache_n_edge_str> \\" << std::endl;
        std::cerr << "    <call_cutoff_str> <pos1_stride_str> <pos2_stride_str> \\" << std::endl;
        std::cerr << "    <pos_array_1_input_str> <pos_array_2_input_str> \\" << std::endl;
        std::cerr << "    <initial_id1_val_str> <initial_id2_val_str>" << std::endl;
        std::cerr << "  test_find_edges <symmetric_type_str> <n_elem1_str> <n_elem2_str> \\" << std::endl;
        std::cerr << "    <initial_cache_valid_str> <initial_cache_buffer_str> <initial_cache_cutoff_str> <initial_cache_n_edge_str> \\" << std::endl;
        std::cerr << "    <call_cutoff_str> <pos1_stride_str> <pos2_stride_str> \\" << std::endl;
        std::cerr << "    <pos_array_1_input_str> <pos_array_2_input_str> \\" << std::endl;
        std::cerr << "    <initial_id1_val_str> <initial_id2_val_str>" << std::endl;
        return 1;
    }

    std::string methodName = argv[1];

    if (methodName == "change_cache_buffer") {
         if (argc != 3) {
            std::cerr << "Usage: " << argv[0] << " change_cache_buffer <new_buffer_value>" << std::endl;
            return 1;
        }
        float new_buffer_value = stringToFloat(argv[2]);
        PairlistComputation<true> pc_symmetric(10, 10, 100); 
        pc_symmetric.change_cache_buffer(new_buffer_value);
        float symmetric_buffer = PairlistComputationTestBridge<true>::get_cache_buffer(pc_symmetric);
        
        PairlistComputation<false> pc_non_symmetric(10, 10, 100); 
        pc_non_symmetric.change_cache_buffer(new_buffer_value);
        float non_symmetric_buffer = PairlistComputationTestBridge<false>::get_cache_buffer(pc_non_symmetric);
        
        std::cout << "symmetric_cache_buffer: " << symmetric_buffer << std::endl;
        std::cout << "non_symmetric_cache_buffer: " << non_symmetric_buffer << std::endl;

    } else if (methodName == "test_ensure_cache_valid" || methodName == "test_find_edges") {
        if (argc != 16) { 
            std::cerr << "Usage for " << methodName << ": <symmetric_type_str> <n_elem1_str> <n_elem2_str> \\" << std::endl;
            std::cerr << "    <initial_cache_valid_str> <initial_cache_buffer_str> <initial_cache_cutoff_str> <initial_cache_n_edge_str> \\" << std::endl;
            std::cerr << "    <call_cutoff_str> <pos1_stride_str> <pos2_stride_str> \\" << std::endl;
            std::cerr << "    <pos_array_1_input_str> <pos_array_2_input_str> \\" << std::endl;
            std::cerr << "    <initial_id1_val_str> <initial_id2_val_str>" << std::endl;
            std::cerr << "Received " << argc -1 << " arguments instead of expected 14 data arguments." << std::endl;
            for(int k=0; k<argc; ++k) std::cerr << "argv[" << k << "]: " << argv[k] << std::endl;
            return 1;
        }

        const char* symmetric_type_str_arg    = argv[2];
        const char* n_elem1_str_arg             = argv[3];
        const char* n_elem2_str_arg             = argv[4];
        const char* initial_cache_valid_str_arg = argv[5];
        const char* initial_cache_buffer_str_arg= argv[6];
        const char* initial_cache_cutoff_str_arg= argv[7];
        const char* initial_cache_n_edge_str_arg= argv[8];
        const char* call_cutoff_str_arg         = argv[9];
        const char* pos1_stride_str_arg         = argv[10];
        const char* pos2_stride_str_arg         = argv[11];
        std::string pos_array_1_input_str_arg = argv[12];
        std::string pos_array_2_input_str_arg = argv[13];
        const char* initial_id1_val_str_arg     = argv[14];
        const char* initial_id2_val_str_arg     = argv[15];

        bool symmetric = (std::string(symmetric_type_str_arg).find("PairlistComputation<true>") != std::string::npos);
        int n_elem1 = stringToInt(n_elem1_str_arg);
        int n_elem2 = stringToInt(n_elem2_str_arg);
        bool initial_cache_valid = stringToBool(initial_cache_valid_str_arg);
        float initial_cache_buffer = stringToFloat(initial_cache_buffer_str_arg);
        float initial_cache_cutoff = stringToFloat(initial_cache_cutoff_str_arg);
        int initial_cache_n_edge = stringToInt(initial_cache_n_edge_str_arg);
        float call_cutoff = stringToFloat(call_cutoff_str_arg);
        int pos1_stride = stringToInt(pos1_stride_str_arg);
        int pos2_stride = stringToInt(pos2_stride_str_arg);
        
        std::vector<float> aligned_pos1_data = parseFloatArrayString(pos_array_1_input_str_arg);
        if (n_elem1 > 0 && pos1_stride > 0) {
            size_t expected_min_size1 = static_cast<size_t>(n_elem1) * pos1_stride;
            if(aligned_pos1_data.size() < expected_min_size1) {
                aligned_pos1_data.resize(expected_min_size1, 1e20f); 
            }
            size_t target_padded_size1 = static_cast<size_t>(round_up_util(n_elem1, 4)) * pos1_stride;
             if (target_padded_size1 < aligned_pos1_data.size()) target_padded_size1 = aligned_pos1_data.size(); // ensure we don't shrink if already larger
            if (aligned_pos1_data.size() < target_padded_size1) {
                aligned_pos1_data.resize(target_padded_size1, 1e20f); 
            }
        } else if (n_elem1 == 0 && pos1_stride > 0 ) { 
             size_t target_padded_size1 = static_cast<size_t>(round_up_util(n_elem1, 4)) * pos1_stride;
             if (target_padded_size1 > 0) aligned_pos1_data.assign(target_padded_size1, 1e20f);
             else aligned_pos1_data.clear();
        }


        std::vector<float> aligned_pos2_data;
        if (!symmetric) {
            aligned_pos2_data = parseFloatArrayString(pos_array_2_input_str_arg);
            if (n_elem2 > 0 && pos2_stride > 0) {
                 size_t expected_min_size2 = static_cast<size_t>(n_elem2) * pos2_stride;
                 if(aligned_pos2_data.size() < expected_min_size2) {
                    aligned_pos2_data.resize(expected_min_size2, 1e20f);
                 }
                size_t target_padded_size2 = static_cast<size_t>(round_up_util(n_elem2, 4)) * pos2_stride;
                if (target_padded_size2 < aligned_pos2_data.size()) target_padded_size2 = aligned_pos2_data.size();
                if (aligned_pos2_data.size() < target_padded_size2) {
                    aligned_pos2_data.resize(target_padded_size2, 1e20f);
                }
            } else if (n_elem2 == 0 && pos2_stride > 0) {
                 size_t target_padded_size2 = static_cast<size_t>(round_up_util(n_elem2, 4)) * pos2_stride;
                 if (target_padded_size2 > 0) aligned_pos2_data.assign(target_padded_size2, 1e20f);
                 else aligned_pos2_data.clear();
            }
        } else {
            aligned_pos2_data = aligned_pos1_data; 
        }

        int initial_id1_val = stringToInt(initial_id1_val_str_arg);
        std::vector<int> id1_data;
        int padded_n_elem1 = round_up_util(n_elem1, 4); 
        if (n_elem1 > 0) {
            id1_data.resize(padded_n_elem1);
            std::iota(id1_data.begin(), id1_data.begin() + n_elem1, initial_id1_val);
            if (n_elem1 < padded_n_elem1) {
                std::fill(id1_data.begin() + n_elem1, id1_data.end(), (n_elem1 > 0 ? id1_data[n_elem1-1] : 0) );
            }
        } else if (padded_n_elem1 > 0) { 
            id1_data.assign(padded_n_elem1, 0);
        }

        int initial_id2_val = stringToInt(initial_id2_val_str_arg);
        std::vector<int> id2_data;
        if (!symmetric) {
            int padded_n_elem2 = round_up_util(n_elem2, 4);
            if (n_elem2 > 0) {
                id2_data.resize(padded_n_elem2);
                std::iota(id2_data.begin(), id2_data.begin() + n_elem2, initial_id2_val);
                 if (n_elem2 < padded_n_elem2) {
                    std::fill(id2_data.begin() + n_elem2, id2_data.end(), (n_elem2 > 0 ? id2_data[n_elem2-1] : 0) );
                }
            } else if (padded_n_elem2 > 0) {
                id2_data.assign(padded_n_elem2, 0);
            }
        } else {
            id2_data = id1_data; 
        }
        
        int true_max_potential_edges = 0;
        if (symmetric) {
            if (n_elem1 > 1) { 
                true_max_potential_edges = n_elem1 * (n_elem1 - 1) / 2;
            }
        } else {
            if (n_elem1 > 0 && n_elem2 > 0) {
                true_max_potential_edges = n_elem1 * n_elem2;
            }
        }
        int max_n_edge_arg = std::max({100, true_max_potential_edges, initial_cache_n_edge});
        max_n_edge_arg = round_up_util(max_n_edge_arg, 16); 
        if (max_n_edge_arg == 0) max_n_edge_arg = 16; 

        if (methodName == "test_ensure_cache_valid") {
            if (symmetric) {
                run_ensure_cache_valid_test<true>(
                    n_elem1, n_elem2, max_n_edge_arg,
                    initial_cache_valid, initial_cache_buffer, initial_cache_cutoff, initial_cache_n_edge,
                    call_cutoff, pos1_stride, pos2_stride,
                    aligned_pos1_data, aligned_pos1_data, 
                    id1_data, id1_data 
                );
            } else {
                run_ensure_cache_valid_test<false>(
                    n_elem1, n_elem2, max_n_edge_arg,
                    initial_cache_valid, initial_cache_buffer, initial_cache_cutoff, initial_cache_n_edge,
                    call_cutoff, pos1_stride, pos2_stride,
                    aligned_pos1_data, aligned_pos2_data,
                    id1_data, id2_data
                );
            }
        } else { // test_find_edges
             if (symmetric) {
                run_find_edges_test<true>(
                    n_elem1, n_elem2, max_n_edge_arg,
                    initial_cache_valid, initial_cache_buffer, initial_cache_cutoff, initial_cache_n_edge,
                    call_cutoff, pos1_stride, pos2_stride,
                    aligned_pos1_data, aligned_pos1_data, 
                    id1_data, id1_data 
                );
            } else {
                run_find_edges_test<false>(
                    n_elem1, n_elem2, max_n_edge_arg,
                    initial_cache_valid, initial_cache_buffer, initial_cache_cutoff, initial_cache_n_edge,
                    call_cutoff, pos1_stride, pos2_stride,
                    aligned_pos1_data, aligned_pos2_data,
                    id1_data, id2_data
                );
            }
        }

    } else {
        std::cerr << "Error: Unknown method name '" << methodName << "'" << std::endl;
        return 1;
    }

    return 0;
}