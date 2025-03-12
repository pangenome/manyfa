# Thread-Safe FASTA/FASTQ Reader

A header-only C++ library providing thread-safe random access to bgzip compressed and indexed FASTA and FASTQ files.

## Why?

[htslib's faidx_t reader is not thread-safe.](https://github.com/samtools/htslib/issues/663)
To work around this, we reimplement the index reader in a thread-safe manner in C++.

## Features

- Thread-safe access to FASTA/FASTQ files
- Support for both compressed (bgzip) and uncompressed files
- Clean, modern C++ API
- High-performance random access using HTSlib
- Auto-discovery and loading of index files
- Automatic index creation when missing

## Requirements

- C++11 compatible compiler
- HTSlib (for BGZF compression support)
- zlib (for compression)
- CMake (for building)

## Building

```bash
cmake -H. -Bbuild -DCMAKE_BUILD_TYPE=Release && cmake --build build --
```

## Example Usage

```cpp
#include "thread_safe_faidx.hpp"
#include <iostream>
#include <string>

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <fasta_file>" << std::endl;
        return 1;
    }
    
    try {
        // Create a FASTA reader (will build index if needed)
        ts_faidx::FastaReader reader(argv[1], true);
        
        // Get all sequence names
        auto sequences = reader.get_sequence_names();
        std::cout << "Found " << sequences.size() << " sequences" << std::endl;
        
        // Fetch a region
        if (!sequences.empty()) {
            std::string first_seq = sequences[0];
            std::cout << "First sequence: " << first_seq << std::endl;
            
            int64_t length = reader.get_sequence_length(first_seq);
            std::cout << "Length: " << length << std::endl;
            
            // Fetch the first 100 bases (or whole sequence if shorter)
            int fetch_len = std::min(length, (int64_t)100);
            std::string region = first_seq + ":0-" + std::to_string(fetch_len);
            std::string sequence = reader.fetch_sequence(region);
            
            std::cout << "Sequence: " << sequence << std::endl;
        }
        
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}
```

## Multi-threaded Example

The included demo application (`fasta_reader`) demonstrates how to use the library
in a multi-threaded environment. Run with:

```bash
./fasta_reader -t 8 -n 1000 <fasta_file>
```

## API Reference

### FastaReader Class

```cpp
// Constructor
FastaReader(const std::string& fasta_path, bool build_index = false);

// Retrieve sequence by region string (e.g., "chr1:1000-2000")
std::string fetch_sequence(const std::string& region) const;

// Retrieve sequence by coordinates
std::string fetch_sequence(const std::string& contig, int64_t start, int64_t end) const;

// Get all sequence names
std::vector<std::string> get_sequence_names() const;

// Get sequence length
int64_t get_sequence_length(const std::string& contig) const;
```

## License

This library is distributed under the MIT License.
