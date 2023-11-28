## Reference Construction Module
Scripts to construct both diploid and graph-based references.

### Scripts
- `gen_diploid_ref.sh`: Constructs a personalized diploid reference genome from VCF files.
- `gen_bwa_index.sh`: Creates BWA indexes for the Haplotype_1 and Haplotype_2 fasta files.
- `gen_vg_autoindex.sh`: Utilizes vg autoindex for indexing genomic data, essential for constructing personalized reference graphs.

### Personalized Graph Reference
The `gen_vg_autoindex.sh` script is crucial for generating personalized reference graphs. This process involves using vg autoindex with the giraffe workflow, as detailed in the study for HG001. The script uses the GRCh38 reference FASTA and personalized diploid variant calls, embedding unique variants into the graph. 

### Installation Prerequisites
- **VG Autoindex and VG Giraffe (v1.46.0)**:
  - For installation, visit the [VG GitHub Repository](https://github.com/vgteam/vg/releases).
- **BWA (0.7.17)**:
  - Install via Bioconda: `conda install -c bioconda bwa`
  - More details are available at the [BWA GitHub Repository](https://github.com/lh3/bwa).

### Usage
- Follow the installation guidelines for vg and BWA.
- Execute `gen_diploid_ref.sh` to create diploid references from VCF files.
- Use `gen_bwa_index.sh` and `gen_vg_autoindex.sh` for indexing the references. 

Refer to the comments within each script for detailed usage instructions. The indexes generated are prerequisites for the processes in the alignment and variant calling phases.

