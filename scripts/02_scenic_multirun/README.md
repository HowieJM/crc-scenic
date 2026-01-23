# SCENIC multirun (VSN pipelines) — runbook


**Purpose:** run pySCENIC **multirun** via the VSN Nextflow pipeline, using the fork: https://github.com/HowieJM/vsn-pipelines

.

## Prerequisites

This runbook assumes you have produced a TAS-only loom file using the scripts in "01_loom_prep":

1) Build the TAS-only SCENIC input loom (script **01**):
   `scripts/01_loom_prep/01_tas_seurat_filtered_to_scenic_loom.R`
   - Loom should be at: `resources/pyscenic_resources/Joanito_StromaOnly_Filtered_Seurat_Object_All_Cohorts_FullDataSet_SCopeLoomR.loom`

2) (Optional) Quick loom check:
   - Python checker: `scripts/01_loom_prep/02_check_loom_all.py`
   - Jupyter tunnel note: `scripts/01_loom_prep/03_jupyter_tunnel.md`


## Set up environment (on the Virtual Machine, Work Station or equivalent)

1) Create/activate a conda env (you can reuse any sensible env)
```bash
conda create -n vsn_pipelines_env -c conda-forge -c bioconda -y
conda activate vsn_pipelines_env
```

2) Install Singularity (or Apptainer)
```bash
conda install -y conda-forge::singularity
```

3) Install Nextflow (from binary, and move to conda bin -> exact version)
```bash
wget https://github.com/nextflow-io/nextflow/releases/download/v21.04.3/nextflow-21.04.3-all
mv nextflow-21.04.3-all nextflow && chmod +x nextflow
mkdir -p "$CONDA_PREFIX/bin" && mv nextflow "$CONDA_PREFIX/bin/"
```

4) Sanity checks
```bash
singularity --version     # expect 3.8.6 (or site-equivalent)
nextflow -version         # expect 21.04.3 -> this is **crucial**
```

5) Session manager (optional, but can be wise)
```bash
tmux new-session -s VSN_pySCENIC   # later: tmux attach -t VSN_pySCENIC
```

6) Locale (some clusters require this)
```bash
export LANG=C
export LC_ALL=C
```

---

7) Pull the VSN pipeline (**use HowieJM fork or equivalent by default**) 

```bash
nextflow pull HowieJM/vsn-pipelines -r master
ls -l ~/.nextflow/assets/HowieJM/vsn-pipelines
```
**VSN-pipelines upstream is not maintained, so this fork allows:**
- up-to-date **cisTarget v10** motif databases
- multirun **parallel** mode works
- but, skips visual reports by default (`skipReports = true`) to avoid broken reports 

**Alternative (not recommended)**
  - upstream: `vib-singlecell-nf/vsn-pipelines` (last maintained before these fixes)
  - works with **older** motif DBs only
  - useful if using v9 motif DBs and wanting visual reports -> may run with manual patching




### IN PROGRESS - BELOW NOT FIXED CODE


## Acquire motif/track resources (if not already present)

With the the environment and pipeline prepared and the loom filed added to the resources, we can now prepare the further data resources needed to run SCENIC. To do so:


First, choose a resources directory (default in this repo layout):

```bash
export RESOURCES_DIR="${PWD}/resources/pyscenic_resources"
mkdir -p "$RESOURCES_DIR" && cd "$RESOURCES_DIR"
```

Then obtain the following files:

1) TF list (hg38) -> a list of transcription factors for GRN inference via GRNBoost2
```bash
wget https://resources.aertslab.org/cistarget/tf_lists/allTFs_hg38.txt
```

2) Motif rankings database (v10, hg38, ±10 kb) -> per-gene motif enrichment rankings for cisTarget pruning 
```bash
FEATHER_DB_URL='https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather'
wget "${FEATHER_DB_URL}"
```

3) Motif→TF "motif2TF" mapping file (cisTarget v10, hg38) → mapping table linking motifs to TFs (must match v10)

The TF list is used by GRNBoost2 to infer GRNs based on TF-gene coexpression across cells. cisTarget then uses the matched v10 motif ranking statistics database and motif2TF annotation files to prune putative TF network modules, retaining only those whose target genes are enriched for the TF’s cognate motifs within ±10 kb of the TSS. In this repo, we use **cisTarget v10 rankings database with ±10 kb around TSS** (full transcript + matching v10 motif→TF table). Note that Aerts Lab host *alternatives, including:* promoter-centric (**±500 bp**), not used here. 


#   Optional checksum
wget https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather.sha1sum.txt
feather_database='hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather'
awk -v feather_database="${feather_database}" '$2==feather_database' hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather.sha1sum.txt | sha1sum -c -

# 3) Motif→TF mapping (v10)
wget https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl

# 4) (Optional) ENCODE track rankings & mapping (only if using track pruning)
TRACK_DB_URL='https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/tc_v1/gene_based/encode_20190621__ChIP_seq_transcription_factor.hg38__refseq-r80__10kb_up_and_down_tss.max.genes_vs_tracks.rankings.feather'
wget "${TRACK_DB_URL}"

wget https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/tc_v1/gene_based/encode_20190621__ChIP_seq_transcription_factor.hg38__refseq-r80__10kb_up_and_down_tss.max.genes_vs_tracks.rankings.feather.sha1sum.txt
track_db="$(basename "${TRACK_DB_URL}")"
awk -v feather_database="${track_db}" '$2==feather_database' encode_20190621__ChIP_seq_transcription_factor.hg38__refseq-r80__10kb_up_and_down_tss.max.genes_vs_tracks.rankings.feather.sha1sum.txt | sha1sum -c -

wget https://resources.aertslab.org/cistarget/track2tf/encode_project_20190621__ChIP-seq_transcription_factor.homo_sapiens.hg38.bigwig_signal_pvalue.track_to_tf_in_motif_to_tf_format.tsv

# Return to your run directory
cd -





--------------------------------------------------------------------------------
Get the pipeline (fork vs upstream) and apply tiny patch if needed
--------------------------------------------------------------------------------
# Nextflow caches pipelines under ~/.nextflow/assets. Pull the pipeline:
nextflow pull vib-singlecell-nf/vsn-pipelines -r master
# If you used/need a fork (e.g., restored multirun parallelism), pull that instead:
# nextflow pull <your-fork-owner>/vsn-pipelines -r master

# (Historical patch) if a ModuleNotFoundError on pyscenic.genesig occurs:
#   replace pyscenic.genesig with ctxcore.genesig in the pipeline helper.
#   Only needed for certain snapshots.
PIPE_DIR="${HOME}/.nextflow/assets/vib-singlecell-nf/vsn-pipelines/src/scenic/bin"
# or: PIPE_DIR="${HOME}/.nextflow/assets/<your-fork-owner>/vsn-pipelines/src/scenic/bin"
if [ -d "$PIPE_DIR" ]; then
  sed -i 's/from pyscenic.genesig import GeneSignature/from ctxcore.genesig import GeneSignature/' "$PIPE_DIR/utils.py" || true
fi

--------------------------------------------------------------------------------
Config (use the committed one, or generate then edit)
--------------------------------------------------------------------------------
# Option A (recommended): use the repo config you committed/edited
#   config/nf_CPUopt_scenic_multirun_25runs_minGene5.config
#   - numRuns = 25
#   - min_genes_regulon = 5
#   - min_regulon_gene_occurrence = 5
#   - rank_threshold = 5000; nes_threshold = 3.0; auc_threshold = 0.05; max_similarity_fdr = 0.001; min_genes = 20
#   - container = 'aertslab/pyscenic_scanpy:0.12.0_1.9.1'
#   - skipReports = true
#   - (optional tracks) comment tracksDb/tracksAnnotation if not used
#
# Option B (advanced): draft a config and then edit
# nextflow config vib-singlecell-nf/vsn-pipelines \
#   -profile scenic,scenic_multiruns,scenic_use_cistarget_motifs,scenic_use_cistarget_tracks,hg38,singularity \
#   > nf_draft.config
# # then edit fields as above (container, numRuns, resources paths, skips, etc.)

--------------------------------------------------------------------------------
Run pySCENIC multirun
--------------------------------------------------------------------------------
# Choose a clean run directory (work/ and out/ will be created here)
RUN_DIR="${PWD}/runs/02_scenic_multirun"
mkdir -p "$RUN_DIR" && cd "$RUN_DIR"

# Locale (if needed for reproducibility on your cluster)
export LANG=C
export LC_ALL=C

# Launch
nextflow -C "${PWD}/../../config/nf_CPUopt_scenic_multirun_25runs_minGene5.config" \
  run vib-singlecell-nf/vsn-pipelines -profile scenic_multiruns

# - You can add `-resume` to continue from a previous partial run:
# nextflow -C ... run vib-singlecell-nf/vsn-pipelines -profile scenic_multiruns -resume

# Notes:
# - Ensure Singularity/Apptainer cache and /tmp have enough space.
# - Consider tmux for long sessions: tmux attach -t VSN_pySCENIC
# - work/ and out/ live under RUN_DIR unless you override NXF_* env vars.

--------------------------------------------------------------------------------
Outputs (where to find the aggregate loom)
--------------------------------------------------------------------------------
# The Nextflow pipeline writes an aggregate SCENIC loom under out/scenic/…/data/
# Example (path shape; names/timestamps vary):
# out/scenic/<project>__25Runs__YYYY_MM_DD/data/<project>__25Runs__YYYY_MM_DD.SCENIC.loom
#
# That loom is the input for the downstream import step:
# scripts/03_downstream/01_import_SCENIC_to_Seurat.R

--------------------------------------------------------------------------------
Real run settings (as in the paper)
--------------------------------------------------------------------------------
# numRuns = 25
# cpus = 32; num_workers = 32
# min_genes_regulon = 5
# min_regulon_gene_occurrence = 5
# rank_threshold = 5000; nes_threshold = 3.0; auc_threshold = 0.05
# skipReports = true
# Container: aertslab/pyscenic_scanpy:0.12.0_1.9.1

# CRUCIAL:
# - Keep tracks* commented if you’re not using track pruning
# - Ensure /tmp and cache dirs have space (failures can be opaque)
# - If you patched utils.py as above, verify it remains patched in the cache

--------------------------------------------------------------------------------
Next step
--------------------------------------------------------------------------------
# Import SCENIC results into Seurat (AUC assays & embeddings):
# scripts/03_downstream/01_import_SCENIC_to_Seurat.R
