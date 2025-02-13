# Reported Features

## Genome Features

### Mappability

`f_map` is the mean mappability score over a 300bp window surrounding the variant. The mappability scores are based on the [ENCODE 50mer alignability track](https://www.genome.ucsc.edu/cgi-bin/hgTrackUi?g=wgEncodeMapability&hgsid=1597917733_IMCnJiLIKMGrBE61XIGRFdlP4nCo). 

## Read Features

Let $R$ be the set of reads containing the variant allele at a locus.   

### Mean Length

 `f_mean_length` = $ \frac{1}{|R|} \sum_{r \in R} L(r) $
where $L(r)$ is the length of read $r$. 

### Mean Template Length

`f_mean_tlen` = $ \frac{1}{|R|} \sum_{r \in R} T(r) $
where $T(r)$ is the template length for the read pair containing $r$. 

### Mean Num. of Mismatches

`f_mean_nm` = $ \frac{1}{|R|} \sum_{r \in R} M(r) $
where $M(r)$ is the number of mismatches on read $r$. 

###  Mean Mapping Quality

`f_mean_mapq` = $ \frac{1}{|R|} \sum_{r \in R} Q(r) $
where $Q(r)$ is the mapping quality for read $r$. 

### Mean Variant Distance to Read End

`f_mean_distreadend` = $ \frac{1}{|R|} \sum_{r \in R} \min(v(r), L(r) - v(r))$
where $v(r)$ is the position of the variant in read $r$. 

### Proportion of Reads with Softclipping

`f_p_softclip`=$ \frac{1}{|R|} \sum_{r \in R} S(r)$
where $S(r)$ is an indicator if read $r$ contains softclipped bases. 

###  Proportion of Reads with Insertions
`f_p_ins` = $ \frac{1}{|R|} \sum_{r \in R} I(r)$
where $I(r)$ is an indicator if read $r$ contains inserted bases. 

### Proportion of Reads with Deletions
 `f_p_del`=$ \frac{1}{|R|} \sum_{r \in R} D(r)$
where $D(r)$ is an indicator if read $r$ contains deleted bases. 

### Proportion of Reads with Properly Mapped Mates
 `f_p_matemapped` = $ \frac{1}{|R|} \sum_{r \in R} P(r)$
where $P(r)$ is an indicator if read $r$ has a properly mapped mate. 

### Directionality
`f_directionality` = $ \left|0.5 - \frac{1}{|R|} \sum_{r \in R} F(r)\right| $
where $F(r)$ is an indicator if read $r$ is mapped on the forward strand.  

### Standard Deviation of Start Positions

`f_std_start` is the standard deviation of the start positions of reads $R$. 

### Standard Deviation of End Positions

`f_std_end` is the standard deviation of the end positions of reads $R$. 




