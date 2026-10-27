conda activate pyclone-vi
pyclone-vi fit -i ./COS_003_009.pyclone_vi.tsv -o COS_003_009.h5 -c 40 -d beta-binomial -r 10 --num-threads 2
pyclone-vi write-results-file -i COS_003_009.h5 -o COS_003_009.results.tsv
