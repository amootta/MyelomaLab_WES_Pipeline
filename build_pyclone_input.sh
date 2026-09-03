module unload gcc-libs
module unload mpi
module load r
Rscript build_pyclone_input.R manifest.tsv ./ vi
