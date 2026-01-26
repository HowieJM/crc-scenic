# Here is a list of required resources - files - to run this analysis (not bundled)

> This repo **starts with the TAS-derived loom** file (as defined below). 
> It then runs pySCENIC multirun + downstream analyses.

> **Third-party** resources are **not** bundled, due to GitHub space & licensing.
> But, scripts to create the core loom file from the scRNA Seurat data are included.

##

> To run, place the following files under `${projectDir}/pyscenic_resources/`
> (default in this repo: `resources/pyscenic_resources/`), or run with `--resources_dir /ABS/PATH`.

# Files to include:

| Path under `pyscenic_resources/` | File (what it is) | Source / link | Notes |
|---|---|---|---|
| `Joanito_StromaOnly_Filtered_Seurat_Object_All_Cohorts_FullDataSet_SCopeLoomR.loom` | **SCENIC input (TAS-only) loom.** Derived from **Joanito et al.** CRC scRNA-seq; **TAS re-annotation/selection** by M. Frank; **SCENIC pre-filtering + loom export** by J. Howie. | Raw scRNA-seq: Synapse **syn26844071** → https://www.synapse.org/#!Synapse:syn26844071/ • For access to the TAS subset Seurat or filtered loom, contact **Dietmar Herndler-Brandstetter** <dietmar.herndler-brandstetter@meduniwien.ac.at>. • For **loom-generation**, see scripts or contact **James Howie**. |
| `allTFs_hg38.txt` | HG38 TF list (pySCENIC) | Aerts lab / pySCENIC docs | 
| `hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather` | cisTarget motif rankings (hg38 v10) | Aerts lab resources | 
| `motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl` | Motif→TF mapping | Aerts lab resources | 
| *(optional)* `encode_20190621__ChIP_seq_transcription_factor.hg38__refseq-r80__10kb_up_and_down_tss.max.genes_vs_tracks.rankings.feather` | ENCODE TF tracks | ENCODE / pySCENIC | 
| *(optional)* `encode_project_20190621__ChIP-seq_transcription_factor.homo_sapiens.hg38.bigwig_signal_pvalue.track_to_tf_in_motif_to_tf_format.tsv` | ENCODE TF track→TF mapping | ENCODE / pySCENIC | 

> Note: we do **not** here redistribute the TAS Seurat object or derived loom (size and licensing). This repository begins at the TAS-only loom file, see below. However, the scripts to regenerate this loom from the TAS Seurat object are provided in `scripts/prefiltering/`.

## Main scRNA CRC-TAS Loom File Provenance:
- **Raw scRNA-seq**: Joanito et al. (Nat Genet 2022), Synapse **syn26844071** (raw scRNA CRC) 
- **TAS subset**: stromal re-annotation & selection by M. Frank (8 TAS classes, TAS-only CRC Seurat file)  
- **SCENIC pre-filtering + loom**: performed by J. Howie; scripts in `scripts/prefiltering/` (Seurat → loom file).
