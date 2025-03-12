#include "thread_safe_faidx.hpp"
#include <iostream>
#include <vector>
#include <thread>
#include <random>
#include <chrono>
#include <atomic>
#include <string>
#include <cstring>
#include <fstream>
#include <mutex>

void print_usage(const char* program) {
    std::cerr << "Usage: " << program << " [options] <fasta_file>" << std::endl;
    std::cerr << "Options:" << std::endl;
    std::cerr << "  -t <threads>   Number of threads to use (default: 8)" << std::endl;
    std::cerr << "  -n <count>     Number of sequences to read (default: 1000)" << std::endl;
    std::cerr << "  -o <file>      Output file to write sequences to" << std::endl;
    std::cerr << "  -l <length>    Fixed sequence length to use for all queries (default: random 10-100)" << std::endl;
    std::cerr << "  -h             Show this help message" << std::endl;
}

// Worker function that reads sequences in a thread
void worker_thread(const ts_faidx::FastaReader& reader, 
                   const std::vector<std::string>& regions,
                   size_t start_idx, size_t end_idx,
                   std::atomic<size_t>& completed,
                   int thread_id,
                   std::ofstream* output_file,
                   std::mutex* file_mutex) {
    for (size_t i = start_idx; i < end_idx; ++i) {
        try {
            auto sequence = reader.fetch_sequence(regions[i]);
            completed++;
            
            // Write to output file if specified
            if (output_file && file_mutex) {
                std::lock_guard<std::mutex> lock(*file_mutex);
                *output_file << ">" << regions[i] << "\n" 
                           << sequence << "\n";
            }
            
            // Occasionally print progress
            if (completed % 100 == 0) {
                std::cout << "Progress: " << completed << "/" << regions.size() 
                          << " sequences read" << std::endl;
            }
        } catch (const std::exception& e) {
            std::cerr << "Thread " << thread_id << " error: " << e.what() 
                      << " for region " << regions[i] << std::endl;
        }
    }
}

int main(int argc, char* argv[]) {
    // Parse command line arguments
    int num_threads = 8;
    int num_sequences = 1000;
    int fixed_length = 0;  // 0 means random length
    const char* fasta_file = nullptr;
    const char* output_file = nullptr;
    
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "-t") == 0) {
            if (i + 1 >= argc) {
                std::cerr << "Error: -t option requires a value" << std::endl;
                print_usage(argv[0]);
                return 1;
            }
            num_threads = std::stoi(argv[++i]);
        } else if (strcmp(argv[i], "-n") == 0) {
            if (i + 1 >= argc) {
                std::cerr << "Error: -n option requires a value" << std::endl;
                print_usage(argv[0]);
                return 1;
            }
            num_sequences = std::stoi(argv[++i]);
        } else if (strcmp(argv[i], "-o") == 0) {
            if (i + 1 >= argc) {
                std::cerr << "Error: -o option requires a value" << std::endl;
                print_usage(argv[0]);
                return 1;
            }
            output_file = argv[++i];
        } else if (strcmp(argv[i], "-l") == 0) {
            if (i + 1 >= argc) {
                std::cerr << "Error: -l option requires a value" << std::endl;
                print_usage(argv[0]);
                return 1;
            }
            fixed_length = std::stoi(argv[++i]);
            if (fixed_length <= 0) {
                std::cerr << "Error: sequence length must be positive" << std::endl;
                return 1;
            }
        } else if (strcmp(argv[i], "-h") == 0) {
            print_usage(argv[0]);
            return 0;
        } else if (argv[i][0] == '-') {
            std::cerr << "Unknown option: " << argv[i] << std::endl;
            print_usage(argv[0]);
            return 1;
        } else if (fasta_file == nullptr) {
            fasta_file = argv[i];
        } else {
            std::cerr << "Unexpected argument: " << argv[i] << std::endl;
            print_usage(argv[0]);
            return 1;
        }
    }
    
    if (fasta_file == nullptr) {
        std::cerr << "No FASTA file specified" << std::endl;
        print_usage(argv[0]);
        return 1;
    }
    
    try {
        // Verify file exists first
        FILE* test_file = fopen(fasta_file, "r");
        if (!test_file) {
            std::cerr << "Error: Cannot open file " << fasta_file << ": " 
                      << strerror(errno) << std::endl;
            return 1;
        }
        fclose(test_file);
        
        // Create the FASTA reader
        auto start_time = std::chrono::high_resolution_clock::now();
        std::cout << "Loading FASTA file: " << fasta_file << std::endl;
        ts_faidx::FastaReader reader(fasta_file, true);
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        std::cout << "FASTA index loaded in " << duration.count() << " ms" << std::endl;
        
        // Get all sequence names
        auto sequences = reader.get_sequence_names();
        std::cout << "Found " << sequences.size() << " sequences in the file" << std::endl;
        
        if (sequences.empty()) {
            std::cerr << "No sequences found in the FASTA file" << std::endl;
            return 1;
        }
        
        // Generate random regions to query
        std::cout << "Generating " << num_sequences << " random queries..." << std::endl;
        std::vector<std::string> regions;
        regions.reserve(num_sequences);
        
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> seq_dist(0, sequences.size() - 1);
        
        for (int i = 0; i < num_sequences; ++i) {
            // Pick a random sequence
            std::string seq_name = sequences[seq_dist(gen)];
            int64_t seq_len = reader.get_sequence_length(seq_name);
            
            if (fixed_length > 0) {
                // Use fixed length if specified
                if (seq_len > fixed_length) {
                    // Generate a random start position within the sequence
                    std::uniform_int_distribution<> start_dist(0, seq_len - fixed_length);
                    int64_t start = start_dist(gen);
                    int64_t end = start + fixed_length;
                    
                    regions.push_back(seq_name + ":" + std::to_string(start) + "-" + std::to_string(end));
                } else {
                    // Sequence is shorter than requested fixed length, use whole sequence
                    regions.push_back(seq_name);
                }
            } else if (seq_len > 100) {
                // Generate a random region with random length
                std::uniform_int_distribution<> start_dist(0, seq_len - 100);
                std::uniform_int_distribution<> len_dist(10, 100);
                
                int64_t start = start_dist(gen);
                int64_t len = len_dist(gen);
                int64_t end = std::min(start + len, seq_len);
                
                regions.push_back(seq_name + ":" + std::to_string(start) + "-" + std::to_string(end));
            } else {
                // Small sequence, just take the whole thing
                regions.push_back(seq_name);
            }
        }
        
        // Adjust thread count if needed
        num_threads = std::min(num_threads, num_sequences);
        std::cout << "Using " << num_threads << " threads to read " 
                  << num_sequences << " sequences" << std::endl;
        
        // Open output file if specified
        std::ofstream output_file_stream;
        std::mutex file_mutex;
        
        if (output_file) {
            output_file_stream.open(output_file);
            if (!output_file_stream.is_open()) {
                std::cerr << "Error: Could not open output file " << output_file << std::endl;
                return 1;
            }
            std::cout << "Writing sequences to " << output_file << std::endl;
        }
        
        // Assign regions to threads
        std::vector<std::thread> threads;
        std::atomic<size_t> completed(0);
        
        size_t regions_per_thread = num_sequences / num_threads;
        size_t remainder = num_sequences % num_threads;
        
        start_time = std::chrono::high_resolution_clock::now();
        
        size_t start_idx = 0;
        for (int t = 0; t < num_threads; ++t) {
            size_t count = regions_per_thread + (t < remainder ? 1 : 0);
            size_t end_idx = start_idx + count;
            
            threads.emplace_back(worker_thread, std::ref(reader), std::ref(regions), 
                                start_idx, end_idx, std::ref(completed), t,
                                output_file ? &output_file_stream : nullptr,
                                output_file ? &file_mutex : nullptr);
            
            start_idx = end_idx;
        }
        
        // Wait for all threads to finish
        for (auto& thread : threads) {
            thread.join();
        }
        
        end_time = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        
        // Print results
        std::cout << "Completed " << completed << " sequence reads in " 
                  << duration.count() << " ms" << std::endl;
        double reads_per_second = (completed * 1000.0) / duration.count();
        std::cout << "Read rate: " << reads_per_second << " sequences/second" << std::endl;
        
        // Close output file if opened
        if (output_file_stream.is_open()) {
            output_file_stream.close();
            std::cout << "Sequences written to " << output_file << std::endl;
        }
        
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}
