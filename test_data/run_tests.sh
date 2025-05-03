#!/bin/bash

# Make sure the test data exists
if [ ! -f "test_data/test1.fa" ]; then
  echo "Test data not found. Running generator..."
  bash test_data/generate_test_data.sh
fi

echo "============ Testing with command line arguments ============"
./build/bin/manyfasta -v -b test_data/regions.bed test_data/test1.fa test_data/test2.fa

echo
echo "============ Testing with list file ============"
./build/bin/manyfasta -v -b test_data/regions.bed -f test_data/fasta_list.txt

echo
echo "============ Testing with prefix ============"
./build/bin/manyfasta -v -b test_data/regions.bed -p "TEST_" test_data/test1.fa test_data/test2.fa

echo
echo "============ Testing with numeric suffixes ============"
./build/bin/manyfasta -v -b test_data/regions.bed -n test_data/test1.fa test_data/test2.fa
