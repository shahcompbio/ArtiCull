# Extract Features

## Usage:
```
python extract_features.py variant_maf mappability_bedgraph output bam_dir1 bam_dir2 ...
```

Extract features will create a tsv file, containing the following columns

- `chrm`
- `pos`
- `f_var_length`. Average length over all variant reads at position
- `f_aligned_length`. Average aligned length over all variant reads at position
- `f_nm`. Average number of mismatches over all variant reads at position
- `f_softclip`. Average number of softclip bases over all variant reads at position
- `f_mapq`. Average mapping quality over all variant reads at position
- `f_tlen`. Average template length over all variant reads at position
- `f_tot`. Total number of reads at position
- `f_map`. Mappability in 300bp window around position
