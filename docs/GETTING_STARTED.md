# Running a batch from scratch

From FASTQs to neoantigens, clonality and signatures, on UCL Myriad. Commands
run on a Myriad login node unless they say otherwise.

- Part A (once per person) gets the software and references in place.
- Part B (once per batch) sets up a batch.
- Part C runs it.

---

## Part A — one-time setup

### A1. Access

- A Myriad account: https://www.rc.ucl.ac.uk/docs/Clusters/Myriad/
- Access to the lab RDSS project (ask Amatta):
  - `/rdss/rd02/ritd-ag-project-rd028n-amira43/reference_files/WES_pipeline`: references
  - `/rdss/rd02/ritd-ag-project-rd028n-amira43/Shared_Scripts`: tool installs

Everything big goes under `~/Scratch` (your home quota is small).

### A2. Get the pipeline

```bash
cd ~/Scratch
git clone https://github.com/amootta/MyelomaLab_WES_Pipeline.git
cd MyelomaLab_WES_Pipeline
```

### A3. Nextflow

```bash
module load java/temurin-17/17.0.2_8
mkdir -p ~/bin && cd ~/bin
curl -s https://get.nextflow.io | bash
chmod +x nextflow
```

Make sure `~/bin` is on your `PATH` in `~/.bashrc`, so SGE jobs find it too.
You only need one install. mgp1000 is pinned to Nextflow 22.10.8 through
`NXF_VER`, and Nextflow downloads that version itself the first time.

### A4. Conda + the four environments

```bash
cd ~/Scratch
curl -LO https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh -b -p ~/Scratch/miniforge3
source ~/Scratch/miniforge3/etc/profile.d/conda.sh

cd ~/Scratch/MyelomaLab_WES_Pipeline
conda env create -f envs/t1k.yml          # stages 04 + 05
conda env create -f envs/hla-la.yml       # stage 04
conda env create -f envs/pyclone-vi.yml   # stage 08
conda env create -f envs/qc.yml           # stage 10: FastQC, Qualimap, MultiQC
```

The jobs load conda themselves (`CONDA_SH`), but your login shell doesn't. If
`conda` is "command not found", run the `source` line above, or add it to
`~/.bashrc` so every session has it.

HLA-LA also needs its reference graph. Copy `PRG_MHC_GRCh38_withIMGT` from
the RDSS into the env:

```bash
cp -r <RDSS>/.../PRG_MHC_GRCh38_withIMGT \
  ~/Scratch/miniforge3/envs/hla-la/opt/hla-la/graphs/
```

### A5. Tools and references from the RDSS

Copy every item below into one folder, `~/Scratch/AMATTA/reference_files`
(`REF_DIR`). It holds programs and references only and is shared by every
batch: each stage is launched from the batch folder, so Nextflow `work/`,
intermediate and final files all go under `$BATCH_DIR`. If you put anything
elsewhere, point its variable in `config/pipeline.config` at it.

| What | Config variable |
|---|---|
| GRCh38 no-alt FASTA + `.fai` | `REF_FASTA`, `REF_FAI` |
| dbSNP, Mills indels, gnomAD AF, panel of normals | `DBSNP`, `KNOWN_INDELS`, `GNOMAD_AF`, `PON` |
| Funcotator data sources | `FUNCOTATOR_DS` |
| ASCAT loci/alleles/GC/RT, mappability, Control-FREEC chr files | `ASCAT_*`, `MAPPABILITY`, `CONTROLFREEC_CHR_DIR` |
| T1K folder (has `hlaidx/`, `kiridx/`) | `T1K_DIR` |
| KIRCLE folder (has `KIRCLE/`, `ref_files/`) | `KIRCLE_DIR` |
| patched nf-core/hlatyping checkout + `hlahd.1.7.1.tar.gz` | `HLATYPING_DIR`, `HLAHD_TARBALL` |
| Myriad-adapted mgp1000 folder | `MGP1000_DIR` |
| IgCaller + `IgCaller_reference_files/` + samtools wrapper | `IGCALLER_HOME`, `IGCALLER_SAMTOOLS` |
| nextflow-pvacseq v1.0.0, IEDB, VEP plugins | `PVACSEQ_PIPELINE`, `PVACSEQ_IEDB`, `PVACSEQ_VEP_PLUGINS` |
| exome capture BED (Agilent V6 padded, hg38), for Qualimap; optional | `QC_TARGET_BED` |

Qualimap needs a 6-column BED (chrom, start, end, name, score, strand). If
yours has 3 or 4 columns, pad it once:

```bash
awk 'BEGIN { OFS = "\t" } !/^(#|track|browser)/ { print $1, $2, $3, "t" NR, 0, "+" }' \
    Exome-Agilent_V6_Padded_hg38.bed > Exome-Agilent_V6_Padded_hg38.6col.bed
```

For R (PyClone summaries, MMSig), `module load r`, then install `data.table`
and `ggplot2` once from R.

---

## Part B — set up a batch

### B1. A config for the batch

```bash
cd ~/Scratch/MyelomaLab_WES_Pipeline
config/pipeline.config
cp config/samples.tsv     [BATCH WORKING DIRECTORY]/samples.tsv
```

In `config/pipeline.config`, change:

```bash
BATCH_NAME=BATCH WORKING DIRECTORY NAME
SAMPLES=[BATCH WORKINGD IRECTORY]/samples.tsv
```

Every output then goes under `~/Scratch/[BATCH WORKING IRECTORY]`. Then, for this shell:

```bash
export PIPELINE_CONFIG=~/Scratch/MyelomaLab_WES_Pipeline/config/pipeline.config
```

### B2. The sample sheet

Fill samples.tsv` with **one row per FASTQ pair (lane)**.
It's tab-separated and keeps the header row.

| patient | sample | status | sex | lane | fastq_1 | fastq_2 | cram | bam |
|---|---|---|---|---|---|---|---|---|
| COS_005_001 | COS_005_001_PBMC | 0 | XX | 3 | /path/..._1.fq.gz | /path/..._2.fq.gz | | |
| COS_005_001 | COS_005_001_BL | 1 | XX | 3 | ... | ... | | |
| COS_005_001 | COS_005_001_BL | 1 | XX | 4 | ... | ... | | |

- `status`: 0 is the normal (PBMC/germline) and 1 is a tumour. Each patient needs one normal.
- Two lanes of one sample share the `sample` name with different `lane`
  values. Sarek merges them, so there's no need to merge them by hand.
- Use short sample names (`<patient>_<timepoint>`). Every output is named after them.
- Leave `cram`/`bam` empty when starting from FASTQs.

### B3. Check before submitting

```bash
./wes.sh check
```

This lists every config path as `ok` or `MISSING`, plus any badly formed
sample-sheet rows. Fix everything it reports before going on.

---

## Part C — run it

Each command below submits SGE jobs and returns straight away. Before it
submits, it checks that the previous step's files exist, so running one too
early just stops with a "missing" message.

| Command | Use |
|---|---|
| `DRY_RUN=1 ./wes.sh <stage> <step>` | shows the `qsub` command without submitting |
| `qstat` | your queued and running jobs |
| `./wes.sh status` | which outputs exist so far, per sample |
| `$BATCH_DIR/logs/<job>.o<id>.<task>` | each job's log |

The Nextflow steps (sarek, hlatyping, mgp1000, pVACseq) use `-resume`. If one
fails, fix the cause and rerun the same command; it carries on from where it stopped.

### C1. Sarek: alignment and variant calling (about 1–2 days)

```bash
./wes.sh sarek fastq
```

Wait for the `sarek_fastq` job to finish, then:

```bash
./wes.sh sarek bams             # recal cram -> bam (HLA/KIR/IgCaller/mgp1000 read bams)
./wes.sh sarek filter_mutect2   # lab Mutect2 filters + Funcotator MAF per pair
```

### C2. These run in parallel once the bams exist

```bash
./wes.sh mgp1000 somatic        # Mutect2 + Strelka + VarScan (+ Manta, Delly, IgCaller)
./wes.sh mgp1000 germline       # DeepVariant germline calls + fastNGSadmix ancestry, per normal
./wes.sh hla t1k
./wes.sh hla hlatyping          # OptiType + HLA-HD
./wes.sh hla hlala
./wes.sh kir t1k
./wes.sh kir kircle
./wes.sh igcaller run
./wes.sh qc fastqc              # FastQC on FASTQs + bams
./wes.sh qc qualimap            # Qualimap bamqc per bam
```

Once the two QC jobs finish:

```bash
./wes.sh qc multiqc             # -> QC/multiqc/multiqc_report.html
```

Or queue it straight away, held until both QC jobs are done:

```bash
HOLD=fastqc,qualimap ./wes.sh qc multiqc
```

### C3. Gather the HLA calls (after the C2 HLA jobs)

```bash
./wes.sh hla collect            # one folder per caller under HLA/calls/
./wes.sh hla concordance        # optional: needs HLA_TRUTH set in the config
./wes.sh kir concordance        # optional
```

### C4. Consensus SNVs (after mgp1000)

```bash
./wes.sh consensus vcf          # 2-of-3 consensus VCF per tumour (pVACseq, PyClone)
./wes.sh consensus isec         # 2-of-3 variant tables (MMSig)
```

### C5. Downstream analyses

```bash
./wes.sh pvacseq hla            # writes pvacseq/hla/<patient>.csv from OptiType
#   -> open the csv files and check the alleles look right, then:
./wes.sh pvacseq run

./wes.sh pyclone run            # needs C4 vcf + sarek ASCAT

./wes.sh signatures input       # -> consensus/isec/mmsig_input.tsv
```

For signatures, run MMSig in R on `mmsig_input.tsv`
(https://github.com/evenrus/mmsig). These calls are **hg38** and mmsig
defaults to hg19, so set the genome when you run it.

### Order at a glance

```
C1 sarek fastq ─► sarek bams ─┬─► mgp1000 ─► consensus vcf ─┬─► pvacseq hla ─► pvacseq run
                              │                │            └─► pyclone run
                              │                └─► consensus isec ─► signatures
                              ├─► hla (t1k, hlatyping, hlala) ─► hla collect ─┘(feeds pvacseq hla)
                              ├─► kir (t1k, kircle)
                              ├─► igcaller
                              └─► qc (fastqc, qualimap) ─► qc multiqc
               sarek filter_mutect2  (only needs sarek; for the MAFs)
```

---

## Starting from existing crams/bams instead of FASTQs

Put the paths in the `cram` (or `bam`) column and leave the FASTQ columns
empty. Then either:

- run `./wes.sh sarek cram`, which repeats variant calling only from those crams, or
- skip sarek entirely if `variant_calling/` already exists, and start at
  `./wes.sh sarek bams`.

## Things to know

- **One mgp1000 install, one output folder per batch.** The install
  (`MGP1000_DIR`) is shared; each batch is launched from `$BATCH_DIR/mgp1000`
  (`MGP1000_RUN_DIR`) and writes to its `output/somatic/`. Both `consensus` steps take only the pairs
  whose tumour is in the current sample sheet, and log how many others they
  skipped. Tumour names with library suffixes from older runs
  (`COS_003_009_PRG_EKDN…_L7`) still match their short sample-sheet name. If a
  tumour has pairs from two mgp1000 runs, `consensus vcf` warns, so you can
  remove the one you don't want.
- mgp1000 can also preprocess from FASTQs itself instead of using sarek's
  bams. That command is kept, commented, in `stages/02_mgp1000/mgp1000.job.sh`,
  but `wes.sh` uses sarek's bams.
