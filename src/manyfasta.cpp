#include "thread_safe_faidx.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <thread>
#include <mutex>
#include <atomic>
#include <future>
#include <algorithm>
#include <memory>
#include <chrono>

// Structure to hold BED entries
struct BedEntry {
    std::string chrom;
    int64_t start;
    int64_t end;
    std::string name;  // Optional name field
    
    BedEntry(const std::string& c, int64_t s, int64_t e, const std::string& n = "")
        : chrom(c), start(s), end(e), name(n) {}
};

// Structure to represent a sequence and its source
struct SequenceSource {
    std::string filename;
    size_t file_index;
    int64_t length;  // Sequence length
};

void print_usage(const char* program) {
    std::cerr << "Usage: " << program << " [options] -b <bed_file> <fasta_file1> [fasta_file2 ...]\n"
              << "Options:\n"
              << "  -b <file>       BED file with regions to extract\n"
              << "  -t <threads>    Number of threads to use (default: number of CPUs)\n"
              << "  -p <prefix>     Prefix for FASTA headers (default: none)\n"
              << "  -n              Add numeric suffix to FASTA headers\n"
              << "  -v              Verbose output\n"
              << "  -h              Show this help message\n";
}

// Read BED file into vector of BedEntry structures
std::vector<BedEntry> read_bed_file(const std::string& filename) {
    std::vector<BedEntry> entries;
    std::ifstream bed_file(filename);
    
    if (!bed_file.is_open()) {
        throw std::runtime_error("Could not open BED file: " + filename);
    }
    
    std::string line;
    while (std::getline(bed_file, line)) {
        // Skip comments and empty lines
        if (line.empty() || line[0] == '#') continue;
        
        std::istringstream iss(line);
        std::string chrom;
        int64_t start, end;
        std::string name;
        
        if (!(iss >> chrom >> start >> end)) {
            continue; // Skip malformed lines
        }
        
        // Optional name field
        if (iss >> name) {
            entries.emplace_back(chrom, start, end, name);
        } else {
            // Create a name in the format "chrom:start-end" if not provided
            name = chrom + ":" + std::to_string(start) + "-" + std::to_string(end);
            entries.emplace_back(chrom, start, end, name);
        }
    }
    
    return entries;
}

// Build an index of sequence names to their source files
std::map<std::string, SequenceSource> build_sequence_index(
    const std::vector<std::string>& fasta_files,
    bool verbose) {
    
    std::map<std::string, SequenceSource> sequence_index;
    std::mutex index_mutex;
    std::vector<std::future<void>> futures;
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    for (size_t file_idx = 0; file_idx < fasta_files.size(); ++file_idx) {
        const std::string& fasta_file = fasta_files[file_idx];
        
        futures.push_back(std::async(std::launch::async, [&, file_idx, fasta_file]() {
            try {
                if (verbose) {
                    std::cerr << "Loading index for " << fasta_file << std::endl;
                }
                
                ts_faidx::FastaReader reader(fasta_file, true);
                auto sequence_names = reader.get_sequence_names();
                
                std::lock_guard<std::mutex> lock(index_mutex);
                for (const auto& seq_name : sequence_names) {
                    // Only add if not already present or if this file has precedence
                    if (sequence_index.find(seq_name) == sequence_index.end()) {
                        SequenceSource source;
                        source.filename = fasta_file;
                        source.file_index = file_idx;
                        source.length = reader.get_sequence_length(seq_name);
                        
                        sequence_index[seq_name] = source;
                    }
                }
            } catch (const std::exception& e) {
                std::cerr << "Error processing " << fasta_file << ": " << e.what() << std::endl;
            }
        }));
    }
    
    // Wait for all threads to complete
    for (auto& f : futures) {
        f.wait();
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    
    if (verbose) {
        std::cerr << "Indexed " << sequence_index.size() << " sequences from " 
                  << fasta_files.size() << " files in " << duration << "ms" << std::endl;
    }
    
    return sequence_index;
}

// Process BED entries and output FASTA sequences
void process_bed_entries(
    const std::vector<BedEntry>& bed_entries,
    const std::map<std::string, SequenceSource>& sequence_index,
    const std::vector<std::string>& fasta_files,
    int num_threads,
    const std::string& header_prefix,
    bool add_numeric_suffix,
    bool verbose) {
    
    // Create a vector of FastaReader objects, one per file
    std::vector<std::unique_ptr<ts_faidx::FastaReader>> readers;
    for (const auto& file : fasta_files) {
        try {
            readers.push_back(std::make_unique<ts_faidx::FastaReader>(file));
        } catch (const std::exception& e) {
            std::cerr << "Error opening " << file << ": " << e.what() << std::endl;
            throw;
        }
    }
    
    // Split work among threads
    std::vector<std::vector<size_t>> thread_work(num_threads);
    for (size_t i = 0; i < bed_entries.size(); ++i) {
        thread_work[i % num_threads].push_back(i);
    }
    
    std::mutex output_mutex;
    std::atomic<size_t> processed_count(0);
    std::vector<std::thread> threads;
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // Process each batch of BED entries in parallel
    for (int t = 0; t < num_threads; ++t) {
        threads.emplace_back([&, t]() {
            for (const auto& idx : thread_work[t]) {
                const auto& entry = bed_entries[idx];
                
                try {
                    // Find the source file for this sequence
                    auto it = sequence_index.find(entry.chrom);
                    if (it == sequence_index.end()) {
                        std::lock_guard<std::mutex> lock(output_mutex);
                        std::cerr << "Warning: Sequence " << entry.chrom << " not found in any FASTA file" << std::endl;
                        continue;
                    }
                    
                    const auto& source = it->second;
                    auto& reader = readers[source.file_index];
                    
                    // Validate coordinates
                    if (entry.start < 0 || entry.end > source.length || entry.start >= entry.end) {
                        std::lock_guard<std::mutex> lock(output_mutex);
                        std::cerr << "Warning: Invalid coordinates for " << entry.chrom 
                                  << ":" << entry.start << "-" << entry.end 
                                  << " (sequence length: " << source.length << ")" << std::endl;
                        continue;
                    }
                    
                    // Fetch the sequence
                    std::string sequence = reader->fetch_sequence(entry.chrom, entry.start, entry.end);
                    
                    // Construct the header
                    std::string header;
                    if (!header_prefix.empty()) {
                        header = header_prefix;
                    }
                    header += entry.name;
                    
                    if (add_numeric_suffix) {
                        header += "_" + std::to_string(idx + 1);
                    }
                    
                    // Output FASTA entry (thread-safe)
                    {
                        std::lock_guard<std::mutex> lock(output_mutex);
                        std::cout << ">" << header << "\n";
                        
                        // Output sequence with line wrapping at 60 characters
                        for (size_t i = 0; i < sequence.length(); i += 60) {
                            std::cout << sequence.substr(i, 60) << "\n";
                        }
                    }
                    
                    // Update progress counter
                    size_t current = ++processed_count;
                    if (verbose && (current % 1000 == 0 || current == bed_entries.size())) {
                        std::lock_guard<std::mutex> lock(output_mutex);
                        std::cerr << "Processed " << current << " / " << bed_entries.size() 
                                  << " regions (" << (current * 100 / bed_entries.size()) << "%)" << std::endl;
                    }
                } catch (const std::exception& e) {
                    std::lock_guard<std::mutex> lock(output_mutex);
                    std::cerr << "Error processing " << entry.chrom << ":" << entry.start << "-" << entry.end 
                              << ": " << e.what() << std::endl;
                }
            }
        });
    }
    
    // Wait for all threads to complete
    for (auto& thread : threads) {
        thread.join();
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    
    if (verbose) {
        std::cerr << "Processed " << bed_entries.size() << " regions in " 
                  << duration << "ms (" << (bed_entries.size() * 1000.0 / duration) << " regions/s)" << std::endl;
    }
}

int main(int argc, char* argv[]) {
    // Parse command line arguments
    std::string bed_file;
    int num_threads = std::thread::hardware_concurrency();
    std::string header_prefix;
    bool add_numeric_suffix = false;
    bool verbose = false;
    std::vector<std::string> fasta_files;
    
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-b" && i + 1 < argc) {
            bed_file = argv[++i];
        } else if (arg == "-t" && i + 1 < argc) {
            num_threads = std::stoi(argv[++i]);
        } else if (arg == "-p" && i + 1 < argc) {
            header_prefix = argv[++i];
        } else if (arg == "-n") {
            add_numeric_suffix = true;
        } else if (arg == "-v") {
            verbose = true;
        } else if (arg == "-h") {
            print_usage(argv[0]);
            return 0;
        } else if (arg[0] != '-') {
            fasta_files.push_back(arg);
        } else {
            std::cerr << "Unknown option: " << arg << std::endl;
            print_usage(argv[0]);
            return 1;
        }
    }
    
    // Check required arguments
    if (bed_file.empty()) {
        std::cerr << "Error: BED file (-b) is required" << std::endl;
        print_usage(argv[0]);
        return 1;
    }
    
    if (fasta_files.empty()) {
        std::cerr << "Error: At least one FASTA file is required" << std::endl;
        print_usage(argv[0]);
        return 1;
    }
    
    try {
        if (verbose) {
            std::cerr << "Loading BED file: " << bed_file << std::endl;
        }
        
        // Read BED file
        std::vector<BedEntry> bed_entries = read_bed_file(bed_file);
        
        if (verbose) {
            std::cerr << "Read " << bed_entries.size() << " regions from BED file" << std::endl;
            std::cerr << "Building sequence index from " << fasta_files.size() << " FASTA files" << std::endl;
        }
        
        // Build index of sequence names to source files
        auto sequence_index = build_sequence_index(fasta_files, verbose);
        
        if (verbose) {
            std::cerr << "Processing regions using " << num_threads << " threads" << std::endl;
        }
        
        // Process BED entries and output sequences
        process_bed_entries(
            bed_entries, 
            sequence_index, 
            fasta_files, 
            num_threads, 
            header_prefix, 
            add_numeric_suffix,
            verbose
        );
        
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}
