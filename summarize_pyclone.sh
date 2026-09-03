module unload gcc-libs
module unload mpi
module load r
Rscript summarise_pyclone_vi.R <results.tsv> <outdir>

