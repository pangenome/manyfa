#!/bin/bash

# Make sure the test data exists
if [ ! -f "test_data/test1.fa" ]; then
  echo "Test data not found. Running generator..."
  bash test_data/generate_test_data.sh
fi

echo "============ Testing with command line arguments ============"
./build/bin/manyfa -v -b test_data/regions.bed test_data/test1.fa test_data/test2.fa

echo
echo "============ Testing with list file ============"
./build/bin/manyfa -v -b test_data/regions.bed -f test_data/fasta_list.txt

echo
echo "============ Testing with different thread counts ============"
./build/bin/manyfa -v -t 4 -b test_data/regions.bed test_data/test1.fa test_data/test2.fa

echo
echo "============ Testing without verbose output ============"
./build/bin/manyfa -b test_data/regions.bed test_data/test1.fa test_data/test2.fa
