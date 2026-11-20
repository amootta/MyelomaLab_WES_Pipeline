# MyelomaLab_WES_Pipeline

Tools and pipelines for paired tumour–normal whole-exome sequencing in the
Myeloma Lab, set up for UCL Myriad (SGE + apptainer).

Germline calls (HLA, KIR) come from each patient's normal. Somatic calls and
everything downstream (consensus SNVs, IgCaller, neoantigens, clonality,
signatures) are per tumour/normal pair. HLA loss of heterozygosity (LOHHLA)
is per patient, across all of its tumours.

```
FASTQ ──► 01 sarek ──► recal crams/bams, Mutect2, Strelka, ASCAT, CNVkit, Manta
                │
                ├──► 02 mgp1000 ──► 03 consensus ──► 07 pVACseq ◄── 04 HLA
                │                       └──────────► 09 signatures (MMSig)
                ├──► 04 HLA   (T1K, OptiType/HLA-HD, HLA-LA)  on normals
                ├──► 05 KIR   (T1K, KIRCLE)                   on normals
                ├──► 06 IgCaller                              per pair
                ├──► 08 PyClone-VI  (03 consensus + ASCAT)    per patient
                ├──► 10 QC      (FastQC, Qualimap, MultiQC)   every sample
                └──► 11 LOHHLA  (04 HLA calls + ASCAT)        per patient
```

**New here? [docs/GETTING_STARTED.md](docs/GETTING_STARTED.md)** walks through
installing everything, setting up a batch and running it from FASTQs.

**[docs/STAGES.md](docs/STAGES.md)** lists, for every step, what it reads,
what it writes and where. [docs/LOHHLA.md](docs/LOHHLA.md) covers stage 11
in more depth.

## Layout

| Path | What's there |
|---|---|
| `wes.sh` | the one command you run |
| `config/pipeline.config` | **every path** (references, tools, output folders) |
| `config/samples.tsv` | **every sample**: one row per sequencing lane |
| `config/nextflow/` | Myriad fixes for sarek, hlatyping and nextflow-pvacseq |
| `stages/NN_name/` | `run.sh` (checks inputs, submits) + `*.job.sh` (the SGE jobs) |
| `bin/` | Python/R helpers (consensus filtering, VCF conversion, concordance, PyClone) |
| `lib/common.sh` | shared helpers: sample lookups, file locations, `qsub` wrapper |
| `envs/` | conda environments (T1K, HLA-LA, PyClone-VI, QC, LOHHLA) |

## Quick start (on a Myriad login node)

```bash
git clone https://github.com/amootta/MyelomaLab_WES_Pipeline.git
cd MyelomaLab_WES_Pipeline
# 1. edit config/pipeline.config (paths) and config/samples.tsv (samples)
./wes.sh check                    # every path and the sample sheet look right?
DRY_RUN=1 ./wes.sh sarek fastq    # see the qsub line without submitting
./wes.sh sarek fastq              # submit
./wes.sh status                   # what exists so far, per sample
```

Then work down [docs/STAGES.md](docs/STAGES.md). Each step checks the
previous step's output exists before it queues anything, so a missing file
fails in seconds on the login node instead of hours into a job. Job logs go
to `$BATCH_DIR/logs/`.

A new batch: copy `config/pipeline.config` and `config/samples.tsv`, change
`BATCH_NAME` / `SAMPLES`, and run with
`PIPELINE_CONFIG=/path/to/new.config ./wes.sh ...`.

### The sample sheet

`config/samples.tsv` is tab-separated, one row per lane:

| column | |
|---|---|
| `patient` | e.g. `COS_001_037`. HLA/KIR results are keyed by this |
| `sample` | unique per sample, e.g. `COS_001_037_BL`; lanes of one sample share it |
| `status` | `0` normal, `1` tumour |
| `sex` | `XX` / `XY`, used by sarek and PyClone. Defaults to `XX`; edit it per patient where known |
| `lane`, `fastq_1`, `fastq_2` | for `wes.sh sarek fastq`. A sample sequenced on several lanes gets one row per lane (same `sample`); sarek aligns each lane and merges them into one bam per sample. A repeated lane number within a sample is renamed `3_2` etc. automatically |
| `cram` | optional: use this cram instead of sarek's (e.g. a manually merged one) |
| `bam` | optional: use this bam instead of sarek's (e.g. germline bams from elsewhere) |

Every tumour is paired with its patient's normal automatically. The
sarek, hlatyping, mgp1000 and PyClone sheets are generated from this file.

## Installing the tools

Shared installs and references live on the RDSS (ask Amatta for access):
`/rdss/rd02/ritd-ag-project-rd028n-amira43/Shared_Scripts` (packages) and
`/rdss/rd02/ritd-ag-project-rd028n-amira43/reference_files` (references).
Point the matching variable in `config/pipeline.config` at each one.

| Tool | Install | Config variable |
|---|---|---|
| Nextflow | https://docs.seqera.io/nextflow/install (mgp1000 needs 22.10.8; handled via `NXF_VER`) | on `PATH` |
| conda | miniforge, see https://docs.conda.io | `CONDA_SH` |
| T1K | `conda env create -f envs/t1k.yml` + T1K folder from RDSS (has `hlaidx/`, `kiridx/`) | `T1K_DIR` |
| KIRCLE | KIRCLE zip from RDSS | `KIRCLE_DIR` |
| HLA-LA | `conda env create -f envs/hla-la.yml`, then copy `PRG_MHC_GRCh38_withIMGT` from RDSS into `<env>/opt/hla-la/graphs/` | — |
| nf-core/hlatyping | the patched checkout (in `$REF_DIR/hlatyping-patched`), launched from `$HLA_DIR` | `HLATYPING_DIR` |
| HLA-HD | `hlahd.1.7.1.tar.gz` from RDSS (passed to hlatyping as `--hlahd_path`) | `HLAHD_TARBALL` |
| MGP1000 | the Myriad-adapted mgp1000 folder from RDSS | `MGP1000_DIR` |
| IgCaller | IgCaller + `IgCaller_reference_files/` + samtools wrapper (`$REF_DIR/samtools_wrapper`) | `IGCALLER_HOME`, `IGCALLER_SAMTOOLS` |
| nextflow-pvacseq | v1.0.0 + IEDB + VEP plugins + VEP cache (v102; without it every patient downloads ~20 GB) | `PVACSEQ_*` |
| PyClone-VI | `conda env create -f envs/pyclone-vi.yml` | — |
| FastQC, Qualimap, MultiQC | `conda env create -f envs/qc.yml` (stage 10); Agilent V6 capture BED for Qualimap | `QC_TARGET_BED` |
| LOHHLA | `conda env create -f envs/lohhla.yml`, then `./wes.sh lohhla setup` (once per `REF_DIR`); novoalign from Myriad's `novocraft` module | `LOHHLA_*` |
| Pre-pulled container images | `.sif`/`.img` files Nextflow uses instead of pulling | `CONTAINER_LIBRARY` |
| Funcotator data sources | https://console.cloud.google.com/storage/browser/gcp-public-data--broad-references/funcotator (also on RDSS) | `FUNCOTATOR_DS` |
| MMSig | in R: https://github.com/evenrus/mmsig (calls are **hg38**; mmsig defaults to hg19) | — |

## Analysis choices

- PyClone, MMSig and pVACseq all use the mgp1000 2-of-3 consensus SNVs
  (`CONSENSUS_MIN_CALLERS`). Mutect2-only filtered calls + Funcotator MAFs are
  still produced by `wes.sh sarek filter_mutect2`.
- pVACseq takes OptiType's class I typing (`HLA_SOURCE_FOR_PVACSEQ`).
- LOHHLA (stage 11) takes the same typing by default
  (`LOHHLA_HLA_SOURCE`) and ASCAT purity/ploidy from sarek.

## Myriad workarounds the scripts apply for you

- **Nextflow config parser:** newer Nextflow rejects nf-core's `ucl_myriad`
  config, so every Nextflow job sets `NXF_SYNTAX_PARSER=v1` (`lib/common.sh`).
- **mgp1000 reference:** each batch's `mgp1000/references/` links to the
  shared install, except the FASTA/`.fai`/`.dict`, which point at `REF_FASTA`
  (the reference sarek aligned to). The GATK hg38 FASTA has contigs the bams
  lack, and Manta, Strelka and GATK refuse that.
- **mgp1000 sample names:** sarek writes read-group `SM` as
  `<patient>_<sample>`, but mgp1000 expects the bare sample name. The staged
  bams get a new header with `SM:<sample>`. These are new files, so sarek's
  own bams are not changed.
- **mgp1000 somatic + germline together:** they share one launch folder, so
  each one resumes its own last Nextflow session instead of a bare `-resume`.

## Still to check (marked `TODO` in the code)

- A few mgp1000 output folder/file names in stage 03 are assumed, not seen.
- `grep -rn TODO .` lists the rest.
