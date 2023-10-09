
rm -r output-* 

# find . -executable -type f
./build/src/HPOphenoSemSimComb_run

wait

# combine user input ("Input/gene_variants_table.csv") with output into a separate directory
dd=$(date | sed 's/ /-/g')
mkdir output-$dd
paste Input/gene_variants_table.csv output_col_combined_score_ERIC_ACMG.dat  | sed 's/\t/;/' > output-$dd/out.csv

rm output_col_combined_score_ERIC_ACMG.dat

