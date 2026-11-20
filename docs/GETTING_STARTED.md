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

### A4. Conda + the environments

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
conda env create -f envs/lohhla.yml       # stage 11: LOHHLA
```

The jobs load conda themselves (`CONDA_SH`), but your login shell doesn't. If
`conda` is "command not found", run the `source` line above, or add it to
`~/.bashrc` so every session has it.

`envs/t1k.yml` pins only T1K on purpose: pinning every library stops conda
from finding a samtools 0.1.19 build that works with T1K.

HLA-LA also needs its reference graph. Copy `PRG_MHC_GRCh38_withIMGT` from
the RDSS into the env:

```bash
cp -r <RDSS>/.../PRG_MHC_GRCh38_withIMGT \
  ~/Scratch/miniforge3/envs/hla-la/opt/hla-la/graphs/
```

### A5. Tools and references from the RDSS

Everything the pipeline needs is on the lab RDSS, in two folders:

- `/rdss/rd02/ritd-ag-project-rd028n-amira43/reference_files/WES_pipeline`: genome, variant resources, ASCAT, Control-FREEC, capture BED, container images
- `/rdss/rd02/ritd-ag-project-rd028n-amira43/Shared_Scripts`: tools (mgp1000, T1K, KIRCLE, hlatyping, HLA-HD, IgCaller, IEDB, VEP plugins)

**Jobs don't read the RDSS.** Copy both into one folder on Scratch,
`~/Scratch/AMATTA/reference_files` (`REF_DIR`), and unpack the zips **before
running anything**. `REF_DIR` holds programs and references only, and every
batch shares it. Each stage is launched from the batch folder, so Nextflow
`work/` and all intermediate and final files go under `$BATCH_DIR`.

```bash
RDSS=/rdss/rd02/ritd-ag-project-rd028n-amira43
REF=~/Scratch/AMATTA/reference_files
mkdir -p $REF && cd $REF

# 1. references, as they are (the Funcotator folder is already unzipped there)
rsync -a --info=progress2 --exclude 'funcotator_dataSources*.zip' $RDSS/reference_files/WES_pipeline/ $REF/

# 2. tools: folders as they are, zips unpacked (annovar and pyclone are not needed)
rsync -a --info=progress2 --exclude '*.zip' --exclude annovar --exclude pyclone \
      $RDSS/Shared_Scripts/ $REF/
for z in mgp1000 T1K KIRCLE IgCaller samtools_wrapper VEP_plugins; do
    unzip -qn $RDSS/Shared_Scripts/$z.zip -d $REF
done
chmod +x $REF/samtools_wrapper

# 3. reference archives the config reads as folders
unzip -qn ASCAT.zip                  # -> ASCAT/ with the 4 battenberg/GC/RT zips (leave those zipped)
unzip -qn chrfiles_controlfreec.zip  # -> chrfiles_controlfreec/
```

`hlahd.1.7.1.tar.gz` stays a tarball. The ASCAT files inside `ASCAT/` stay
zipped too, because sarek reads them that way.

Then **check that everything is where the config expects**. Fix each
`MISSING` line, either by renaming the folder (for example, `T1K.zip` might
unpack as `T1K/` instead of `T1K-1.0.10/`) or by changing that line in
`config/pipeline.config`:

```bash
cd ~/Scratch/MyelomaLab_WES_Pipeline
./wes.sh check
```

| What | Config variable | From |
|---|---|---|
| GRCh38 no-alt FASTA + `.fai` + `.dict` | `REF_FASTA`, `REF_FAI` | `reference_files` |
| dbSNP, Mills indels, gnomAD AF, panel of normals | `DBSNP`, `KNOWN_INDELS`, `GNOMAD_AF`, `PON` | `reference_files` |
| Funcotator data sources (folder) | `FUNCOTATOR_DS` | `reference_files` |
| `ASCAT/` (battenberg loci/alleles, GC, RT zips) | `ASCAT_*` | `reference_files/ASCAT.zip` |
| mappability `out100m2_hg38.gem`, `chrfiles_controlfreec/` | `MAPPABILITY`, `CONTROLFREEC_CHR_DIR` | `reference_files` (unzip the chr files) |
| Agilent V6 capture BED, 6 columns, for Qualimap | `QC_TARGET_BED` | `reference_files` |
| pre-pulled apptainer images (`*.img`; Nextflow pulls anything missing) | `CONTAINER_LIBRARY` | `reference_files` |
| `T1K-1.0.10/` (has `hlaidx/`, `kiridx/`) | `T1K_DIR` | `Shared_Scripts/T1K.zip` |
| `KIRCLE/` (has `KIRCLE/`, `ref_files/`) | `KIRCLE_DIR` | `Shared_Scripts/KIRCLE.zip` |
| `hlatyping-patched/` + `hlahd.1.7.1.tar.gz` | `HLATYPING_DIR`, `HLAHD_TARBALL` | `Shared_Scripts` |
| `mgp1000/` (Myriad-adapted; has `references/`, `containers/`) | `MGP1000_DIR` | `Shared_Scripts/mgp1000.zip` |
| `IgCaller/` (has `IgCaller_reference_files/`) + `samtools_wrapper` | `IGCALLER_HOME`, `IGCALLER_SAMTOOLS` | `Shared_Scripts` zips |
| `iedb/`, `VEP_plugins/` | `PVACSEQ_IEDB`, `PVACSEQ_VEP_PLUGINS` | `Shared_Scripts` |
| `nextflow-pvacseq-v.1.0.0/` (has `main.nf`) | `PVACSEQ_PIPELINE` | RDSS |
| `vep_cache/` (VEP 102, has `homo_sapiens/102_GRCh38/`). Without it, every patient's pVACseq run downloads ~20 GB | `PVACSEQ_VEP_CACHE` | RDSS |
| `lohhla/` | `LOHHLA_HOME` | made by `./wes.sh lohhla setup` (A6) |

The pVACseq pipeline and the VEP cache are large. If either is a `.zip` or
`.tar.gz` on the RDSS, the commands above don't unpack it. Unpack it into
`$REF` the same way, so `$REF/nextflow-pvacseq-v.1.0.0/main.nf` and
`$REF/vep_cache/homo_sapiens/` exist.

Qualimap needs a 6-column BED (chrom, start, end, name, score, strand). If
yours has 3 or 4 columns, pad it once:

```bash
awk 'BEGIN { OFS = "\t" } !/^(#|track|browser)/ { print $1, $2, $3, "t" NR, 0, "+" }' \
    Exome-Agilent_V6_Padded_hg38.bed > Exome-Agilent_V6_Padded_hg38.6col.bed
```

### A6. LOHHLA install (once per `REF_DIR`)

```bash
./wes.sh lohhla setup   # clones LOHHLA, patches it for hg38, downloads IMGT HLA files into $REF_DIR/lohhla
```

### A7. R packages

For R (PyClone summaries, MMSig), `module load r`, then install `data.table`
and `ggplot2` once from R.

---

## Part B — set up a batch

### B1. A config for the batch

Pick a name for the batch (for example `WES_tumor_batch2`). That name is the
batch's working directory: `~/Scratch/AMATTA/<batch name>/`.

In `config/pipeline.config`, set:

```bash
BATCH_NAME=WES_tumor_batch2
```

`SAMPLES` already points at `$HOME/Scratch/AMATTA/$BATCH_NAME/samples.tsv`,
so make that folder and start the sample sheet from the template:

```bash
cd ~/Scratch/MyelomaLab_WES_Pipeline
mkdir -p ~/Scratch/AMATTA/WES_tumor_batch2
cp config/samples.tsv ~/Scratch/AMATTA/WES_tumor_batch2/samples.tsv
```

Every output then goes under `~/Scratch/AMATTA/WES_tumor_batch2/`. Then, for this shell:

```bash
export PIPELINE_CONFIG=~/Scratch/MyelomaLab_WES_Pipeline/config/pipeline.config
```

### B2. The sample sheet

Fill `samples.tsv` with **one row per FASTQ pair (lane)**.
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

### C6. HLA loss of heterozygosity (LOHHLA, stage 11)

Needs sarek's bams and ASCAT (C1) and `hla collect` (C3). It can run alongside
C5. Once per `REF_DIR` (not per batch), install it with `./wes.sh lohhla setup`.
Then, per batch:

```bash
./wes.sh lohhla inputs     # LOHHLA/<P>/input/: bams, hla_alleles.txt, solutions.txt
#   -> check hla_alleles.txt and solutions.txt, then:
./wes.sh lohhla run        # one SGE task per patient (~hours each)
./wes.sh lohhla summary    # LOHHLA/LOHHLA_summary{,_wide}.tsv
```

Details, settings and how to read the output: [LOHHLA.md](LOHHLA.md).

### Order at a glance

```
C1 sarek fastq ─► sarek bams ─┬─► mgp1000 ─► consensus vcf ─┬─► pvacseq hla ─► pvacseq run
                              │                │            └─► pyclone run
                              │                └─► consensus isec ─► signatures
                              ├─► hla (t1k, hlatyping, hlala) ─► hla collect ─┘(feeds pvacseq hla)
                              │                                     └─► lohhla (+ ASCAT)
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
