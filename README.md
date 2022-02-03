# RADAR: RNA sensing using ADAR

Code, analyses, and other resources accompanying the RADAR manuscript.

`radar.py` contains functions to generate new sensor designs, and to analyze whole genomes based
on Biomart transcript exports.

You can access Biomart exports for mouse and human transcriptomes at https://rnasensing.bio/biomart-files

You additionally need blast-indexed, Fasta-formatted genome files.
E.g., for mouse, download `GCF_000001635.27_GRCm39_genomic.fna.gz` from https://www.ncbi.nlm.nih.gov/assembly/GCF_000001635.27/, uncompress and index it for blast:
```bash
zcat GCF_000001635.27_GRCm39_genomic.fna.gz > GRCm39.fa
makeblastdb -in GRCm39.fa -dbtype nucl
```

`RADAR human transcripts.ipynb` and `RADAR mouse transcripts.ipynb` contain analysis of transcriptomes.
The generated files are available at https://rnasensing.bio/


# Example: generating sensors for murine Bdnf

```bash
python radar.py design_based_on_mart_export \
    --gene Bdnf --variant CCA \
    --mart-export musmusculus_3UTR_mart_export.txt.gz --fa-genome GRCm39.fa \
    --verbose 5 \
    > Bdnf_designs
less Bdnf_designs
```

The `--verbose` flag can be omitted for less output.
