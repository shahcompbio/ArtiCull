# Features

## Genome Features

### `f_map`  Mappability

The mean mappability score over a 300bp window surrounding the variant. The mappability scores are based on the [ENCODE 50mer alignability track](https://www.genome.ucsc.edu/cgi-bin/hgTrackUi?g=wgEncodeMapability&hgsid=1597917733_IMCnJiLIKMGrBE61XIGRFdlP4nCo). 

## Read Features

Let $R$ be the set of reads containing the variant allele at a locus.   

### `f_var_length` Mean Variant Read Length

$$ \frac{1}{|R|} \sum_{r \in R} L(r) $$

