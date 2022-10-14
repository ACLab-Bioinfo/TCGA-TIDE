#!/bin/sh

cat Genelist_config.txt | while read gene
do

Input_file1="./input_TIDE/${gene}/order_data_TUMOR.first_last_15_normalized.${gene}.txt"
Input_file2="./input_TIDE/${gene}/order_data_TUMOR.first_last_25_normalized.${gene}.txt"
Input_file3="./input_TIDE/${gene}/order_data_TUMOR.first_last_30_normalized.${gene}.txt"

Output_file="./output_TIDE/${gene}/TIDE.output.first_last_15_normalized.${gene}.txt"
Output_file="./output_TIDE/${gene}/TIDE.output.first_last_25_normalized.${gene}.txt"
Output_file="./output_TIDE/${gene}/TIDE.output.first_last_30_normalized.${gene}.txt"

tidepy ${Input_file1} -o ${Output_file1} -c Other --ignore_norm &
tidepy ${Input_file2} -o ${Output_file2} -c Other --ignore_norm &
tidepy ${Input_file3} -o ${Output_file3} -c Other --ignore_norm &

done
