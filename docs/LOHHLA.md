# Stage 11: LOHHLA (`./wes.sh lohhla`)

[LOHHLA](https://bitbucket.org/mcgranahanlab/lohhla) estimates allele-specific
copy number of HLA-A, -B and -C in each tumour, and calls loss of
heterozygosity (LOH) of an HLA allele. It runs once per patient, with all of
that patient's tumours against its normal.

The code is in `stages/11_lohhla/` (`run.sh`, `lohhla.job.sh`, `lohhla_env.sh`)
and `bin/{build,summarize}_lohhla.py`. Like every other stage, it reads the
batch config and sample sheet, and takes these earlier outputs:

| From | Used for |
|---|---|
| stage 01 recal bams (normal + tumours, indexed) | reads |
| stage 01 ASCAT `<T>_vs_<N>.purityploidy.txt` | `--CopyNumLoc` (purity/ploidy). Without it, LOHHLA reports allelic imbalance only |
| stage 04 `wes.sh hla collect` → `HLA/calls/<caller>/` | germline class I alleles (OptiType by default) |

## One-time setup (login node, once per `REF_DIR`)

```bash
conda env create -f envs/lohhla.yml
./wes.sh lohhla setup
```

**novoalign is not on conda.** Setup and the job load Myriad's module,
`novocraft/3.04.06` (`LOHHLA_NOVO_MODULE`). Elsewhere, set
`LOHHLA_NOVO_MODULE=""` and point `LOHHLA_NOVO_DIR` at a folder containing
`novoalign` and `novoindex`, downloaded from <https://www.novocraft.com/support/download/>.

`setup` writes to `$REF_DIR/lohhla` (shared across batches):

- `src/`: clone of the original bitbucket LOHHLA.
- `LOHHLAscript.hg38.R`: the script with two changes, applied by `sed`. The run fails if either change did not apply.
  1. The three HLA-A/B/C read windows are hg19 in the original. They are moved to GRCh38 (+32,223 bp).
  2. Calls to `SamToFastq.jar`, `SortSam.jar` and `FilterSamReads.jar` become `picard.jar <Tool>`.
- `imgt/`: IMGT/HLA `hla_gen.fasta` and `hla.dat` from the same release, with the release number in `VERSION`. Also `hla_abc_gen.fasta`, which holds only the class I genomic sequences, renamed to LOHHLA's style (`hla_a_02_01_01_01`).

## Per batch

```bash
./wes.sh lohhla inputs              # writes LOHHLA/<P>/input/
# check LOHHLA/<P>/input/hla_alleles.txt and solutions.txt
DRY_RUN=1 ./wes.sh lohhla run       # see the qsub line
./wes.sh lohhla run                 # all patients
./wes.sh lohhla run COS_001_063     # or just some
```

`inputs` writes, for each patient:

- `input/bams/<S>.bam` (+ `.bam.bai`): links to the recal bams. LOHHLA treats every `*.bam` in this folder as a sample and names it after the file, so the names must be `<sample>.bam` with no other dots.
- `input/hla_alleles.txt`: each 2-field call (e.g. `A*02:01`) mapped to the lowest-numbered IMGT genomic allele in that group (`hla_a_02_01_01_01`). Homozygous genes are listed once, and LOHHLA skips them because LOH cannot be detected there. Calls missing from the fasta are reported in the output.
- `input/solutions.txt`: ASCAT purity and ploidy for each tumour.

To test one patient interactively, without SGE:

```bash
TASK=1 PIPELINE_DIR=$PWD bash stages/11_lohhla/lohhla.job.sh $HOME/Scratch/AMATTA/<batch>/samplesheets/lohhla_patients.txt
```

## Output (`$BATCH_DIR/LOHHLA/<P>/`)

- `<P>.<minCov>.DNA.HLAlossPrediction_CI.xls`: tab-separated, despite the extension. It has one row per tumour × gene. Main columns:
  - `HLA_type1copyNum_withBAFBin`, `HLA_type2copyNum_withBAFBin` (+ `_lower`/`_upper`): copy number of each allele.
  - `PVal_unique`: paired test of allelic imbalance.
  - `LossAllele` / `KeptAllele`.
  - `numMisMatchSitesCov`.
- `<P>.<minCov>.DNA.IntegerCPN_CI.xls`: the per-site values behind the calls.
- `Figures/`: plots for each tumour.
- `running.hla.loh.exome@<date>_log.txt`: look here first if the run fails.

A common LOH call (McGranahan et al. 2017) is copy number < 0.5 for one allele with `PVal_unique` < 0.01.

### Summary across patients

```bash
./wes.sh lohhla summary                  # defaults below
./wes.sh lohhla summary --min-sites 10   # also require >= 10 covered mismatch sites
```

`summarize_lohhla.py` reads every finished patient's prediction table and drops LOHHLA's duplicate rows. It then gives each tumour × gene one of these calls:

| Call | Rule (defaults: `--loss-cn 0.5`, `--pval 0.01`) |
|---|---|
| `LOH` | `PVal_unique` < 0.01 and one allele's `copyNum_withBAFBin` < 0.5. The lost allele is named |
| `imbalance` | `PVal_unique` < 0.01, but no allele below 0.5, or no copy number because the tumour has no ASCAT result |
| `no_change` | `PVal_unique` ≥ 0.01 |
| `insufficient` | no p-value, or fewer mismatch sites than `--min-sites` |

It writes two files to `$BATCH_DIR/LOHHLA/`:
- `LOHHLA_summary.tsv`: one row per tumour × gene, with the numbers behind each call.
- `LOHHLA_summary_wide.tsv`: one row per tumour, with columns HLA-A, HLA-B and HLA-C. Genes that LOHHLA skipped as homozygous show `homozygous/not run`.

## Settings

All of these are `LOHHLA_*` lines in `config/pipeline.config`, so change them there, per batch. An older batch config without these lines still works: `stages/11_lohhla/lohhla_env.sh` fills in the defaults below.

| Variable | Default | |
|---|---|---|
| `LOHHLA_HLA_SOURCE` | `$HLA_SOURCE_FOR_PVACSEQ` (optitype) | optitype, hlahd or t1k |
| `LOHHLA_MIN_COV` | 10 | `--minCoverageFilter`. LOHHLA's own default is 30 |
| `LOHHLA_PLOTS` | TRUE | set FALSE if plotting fails |
| `LOHHLA_CLEANUP` | FALSE | TRUE deletes the intermediate SAM/BAM files |
| `LOHHLA_HLA_FASTA` | IMGT A/B/C genomic | Polysolver's `abc_complete.fasta` can also be used |
| `LOHHLA_DIR` | `$BATCH_DIR/LOHHLA` | per-batch inputs and output |
| `LOHHLA_HOME` | `$REF_DIR/lohhla` | shared install written by `setup` |
| `LOHHLA_NOVO_MODULE` / `LOHHLA_NOVO_DIR` | `novocraft/3.04.06` / `$REF_DIR/novocraft` | where novoalign comes from |

`./wes.sh status` has an `lohhla` column, which shows `+` for a tumour once its patient's prediction table exists. `./wes.sh check` checks the install once `setup` has run.

## Untested on the first run

- **Leftover hg19 contig names:** the script also extracts reads from the hg19 alt-haplotype contigs (`chr6_apd_hap1` …). These do not exist in the no-alt GRCh38 reference, so samtools prints "invalid region" for them. That is harmless, because most reads come from k-mer "fishing" across the whole bam.
- **Tool versions:** the script was written against samtools 1.3, novoalign 3.07, Picard 1.x and R 3.3. If a command fails, the LOHHLA log shows which one.
- **novoalign licence:** I could not confirm whether the current download runs without a licence file. See the setup section.
- **hla.dat:** the plotting step looks up exon boundaries in `hla.dat` by allele name. If it fails on the current IMGT release, use `LOHHLA_PLOTS=FALSE`.
- **Run time:** fishing greps every read in each bam, so allow several hours per patient. The job requests 24 h.
- **Licence:** LOHHLA is free for non-commercial academic research only.
