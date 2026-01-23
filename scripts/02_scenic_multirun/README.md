# SCENIC multirun (VSN pipelines) — runbook


**Purpose:** run pySCENIC **multirun** via the VSN Nextflow pipeline, using the fork: 

https://github.com/HowieJM/vsn-pipelines

.

## Prerequisites

This runbook assumes you have produced a TAS-only loom file using the scripts in "01_loom_prep":

1) Build the TAS-only SCENIC input loom (script **01**):
   `scripts/01_loom_prep/01_tas_seurat_filtered_to_scenic_loom.R`
   - Loom should be at: `resources/pyscenic_resources/Joanito_StromaOnly_Filtered_Seurat_Object_All_Cohorts_FullDataSet_SCopeLoomR.loom`

2) (Optional) Quick loom check:
   - Python checker: `scripts/01_loom_prep/02_check_loom_all.py`
   - Jupyter tunnel note: `scripts/01_loom_prep/03_jupyter_tunnel.md`


## Set up environment (on the Virtual Machine, Workstation or equivalent)

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


## Acquire motif/track resources (if not already present)

With the environment and pipeline prepared and the loom file added to the resources, we can now prepare the further data resources needed to run SCENIC. To do so:


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
```bash
wget https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl
```
For clarity, the TF list is used by **GRNBoost2** to infer GRNs based on TF–gene co-expression across cells. **cisTarget** then uses the matched v10 motif-ranking database and motif→TF annotation files to prune putative TF network modules, to retain only those where a TF’s target motifs are enriched within ±10 kb of the target gene’s TSS. In this repo, we use **cisTarget v10 files**. Note that Aerts Lab host alternatives (not used here), including promoter-centric (±500 bp). 

Activity scores per regulon are calculated later via **AUCell**.

**Optional tracks.** In addition to motif files, you can choose to include the ENCODE **track rankings** database and **track→TF** mapping so cisTarget can prune modules by **track enrichment** (experimental TF-binding signal from ChIP/ATAC) rather than—or in addition to—sequence PWMs. Tracks are complementary to motifs; if you use them, download the matching track rankings and track-to-tf, and set `tracksDb` / `tracksAnnotation` paths in the config.

4) Track rankings database
```bash
TRACK_DB_URL='https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/tc_v1/gene_based/encode_20190621__ChIP_seq_transcription_factor.hg38__refseq-r80__10kb_up_and_down_tss.max.genes_vs_tracks.rankings.feather'
wget "${TRACK_DB_URL}"
```

5) Track→TF file
```bash
wget https://resources.aertslab.org/cistarget/track2tf/encode_project_20190621__ChIP-seq_transcription_factor.homo_sapiens.hg38.bigwig_signal_pvalue.track_to_tf_in_motif_to_tf_format.tsv
```

.

To check completeness of the motif rankings database, and, if used, the track database, run these optional checksums:

Motif rankings (cisTarget v10) checksum
```bash
wget https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather.sha1sum.txt
feather_database='hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather'
awk -v feather_database="${feather_database}" '$2==feather_database' \
  hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather.sha1sum.txt | sha1sum -c -
```

Track rankings (ENCODE) checksum
```bash
wget https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/tc_v1/gene_based/encode_20190621__ChIP_seq_transcription_factor.hg38__refseq-r80__10kb_up_and_down_tss.max.genes_vs_tracks.rankings.feather.sha1sum.txt
feather_database="$(basename "${TRACK_DB_URL}")"
awk -v feather_database="${feather_database}" '$2==feather_database' \
  encode_20190621__ChIP_seq_transcription_factor.hg38__refseq-r80__10kb_up_and_down_tss.max.genes_vs_tracks.rankings.feather.sha1sum.txt | sha1sum -c -
```

**With the resources in place, return to your run directory**
```bash
cd -
```

## Prepare the config file

**Option A (recommended):** use the committed config 

```bash
`config/nf_CPUopt_scenic_multirun_25runs_minGene5.config`
```

Required fields already set:
```bash
- `numRuns = 25`
- `min_genes_regulon = 5`, `min_regulon_gene_occurrence = 5`
- `rank_threshold = 5000`, `nes_threshold = 3.0`, `auc_threshold = 0.05`, `max_similarity_fdr = 0.001`, `min_genes = 20`
- `container = 'aertslab/pyscenic_scanpy:0.12.0_1.9.1'`
- `skipReports = true`
- *(tracks optional)* keep `tracksDb` / `tracksAnnotation` **commented** unless using tracks
```

**Option B (advanced):** draft then edit
```bash
nextflow config HowieJM/vsn-pipelines \
  -profile scenic,scenic_multiruns,scenic_use_cistarget_motifs,scenic_use_cistarget_tracks,hg38,singularity \
  > nf_draft.config
```
If using option B, set the same fields as Option A (container, numRuns, resource paths, skips, etc.)


## Run pySCENIC multirun

Now, we can run pySCENIC multirun via VSN-pipelines, using the fork for v10 files and parallel processing

This will run 25 times, and aggregate results 

> For detailed run instructions see the fork README: https://github.com/HowieJM/vsn-pipelines

To run:
```bash
# 0) (Recommended) run inside tmux
# attach if created earlier, or create a new session
tmux attach -t VSN_pySCENIC || tmux new -s VSN_pySCENIC

# 1) Choose a clean run directory (Nextflow will create work/ and out/ here)
RUN_DIR="${PWD}/outputs/02_scenic_multirun"
mkdir -p "$RUN_DIR" && cd "$RUN_DIR"

# 2) Locale (already exported above; re-export here only if needed)
# export LANG=C
# export LC_ALL=C

# 3) Launch (uses committed config + VSN fork)
nextflow -C "${PWD}/../../config/nf_CPUopt_scenic_multirun_25runs_minGene5.config" \
  run HowieJM/vsn-pipelines -entry scenic -r master

# Resume a partial run if needed:
# nextflow -C ... run HowieJM/vsn-pipelines -entry scenic -r master -resume
```

### Outputs
The aggregated SCENIC loom will be written under:
```
out/scenic/<project>__25Runs__YYYY_MM_DD/data/<project>__25Runs__YYYY_MM_DD.SCENIC.loom
```

.
## End

Use the output loom in step 03 (import to Seurat):  
```
scripts/03_downstream/01_import_SCENIC_to_Seurat.R
```
.
