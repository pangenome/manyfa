# ManyFA

A high-performance tool for extracting sequences from multiple FASTA files in parallel using BED coordinates.

## Overview

ManyFA efficiently extracts genomic sequences from multiple reference genomes or FASTA files using coordinates 
specified in a BED file. It builds a smart index-of-indexes to map sequence names to their source files, allowing 
fast lookup and retrieval even when working with large datasets or multiple genome builds.

## Key Features

- **Multi-file support**: Query multiple FASTA files simultaneously
- **Parallel processing**: Uses multi-threading for both indexing and sequence extraction
- **BED-driven extraction**: Extract thousands of regions using standard BED format
- **Index-of-indexes**: Automatically maps sequence names to source files
- **Smart caching**: Efficiently handles large reference files
- **Flexible output**: Customizable FASTA headers and formatting
- **Memory efficient**: Streams results to avoid excessive memory usage

## Installation

### Requirements

- C++14 compatible compiler
- HTSlib (for FASTA indexing and BGZF compression support)
- zlib
- CMake (3.10+)

### Building from Source

```bash
git clone https://github.com/pangenome/manyfa.git
cd manyfa
cmake -H. -Bbuild -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

For a static build (with all dependencies linked statically):

```bash
cmake -H. -Bbuild -DCMAKE_BUILD_TYPE=Release -DBUILD_STATIC=ON
cmake --build build
```

## Usage

```bash
./manyfa [options] -b <bed_file> <fasta_file1> [fasta_file2 ...]
```

### Required Arguments

- `-b <file>`: BED file with regions to extract
- At least one FASTA file (multiple files can be specified directly or via a list file)

### Options

- `-f <file>`: File containing a list of FASTA files (one path per line)
- `-t <threads>`: Number of threads to use (default: number of CPU cores)
- `-v`: Verbose output with progress information
- `-d`: Debug output (shows detailed file loading information)
- `-h`: Show help message

### Examples

Extract regions from a single genome:
```bash
./manyfa -v -b regions.bed reference.fa > extracted_sequences.fa
```

Extract regions from multiple genome builds (prioritizing the first match):
```bash
./manyfa -v -t 16 -b regions.bed hg19.fa hg38.fa chm13.fa > multi_build_sequences.fa
```

Extract regions using a list of FASTA files:
```bash
# Create a file with FASTA paths
echo "/path/to/hg19.fa" > fasta_list.txt
echo "/path/to/hg38.fa" >> fasta_list.txt
echo "/path/to/chm13.fa" >> fasta_list.txt

# Run with the list file
./manyfa -v -t 16 -b regions.bed -f fasta_list.txt > multi_build_sequences.fa
```

## How It Works

1. ManyFA loads all FASTA indexes in parallel
2. It builds a sequence-to-file mapping (index-of-indexes)
3. It parses the BED file for extraction coordinates
4. It distributes extraction work across multiple threads
5. It retrieves and outputs sequences in FASTA format

## Underlying Library

ManyFA is built on a thread-safe FASTA/FASTQ accessor library that solves concurrency issues 
in htslib's faidx implementation. The library can be used independently in your own projects.

```cpp
// Usage example
ts_faidx::FastaReader reader("reference.fa");
std::string sequence = reader.fetch_sequence("chr1:1000-2000");
```

## API Reference

The underlying `FastaReader` class provides:

```cpp
// Constructor
FastaReader(const std::string& fasta_path, bool build_index = false);

// Retrieve sequence by region string (e.g., "chr1:1000-2000")
std::string fetch_sequence(const std::string& region) const;

// Retrieve sequence by coordinates
std::string fetch_sequence(const std::string& contig, int64_t start, int64_t end) const;

// Get all sequence names and information
std::vector<std::string> get_sequence_names() const;
int64_t get_sequence_length(const std::string& contig) const;
```

## License

This tool is distributed under the MIT License.
