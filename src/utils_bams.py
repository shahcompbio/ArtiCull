import functools
import os
import pysam

def match_variants_to_filenames(df, data_dirs):
    # Could be shared with extract features?
    """
    Bam files are assumed to be merged and in region format. Thus, for each variant
    we match it to the appropriate file that contains that locus
    """

    files = [(f, f.strip().split('.')[0].split('-')) for f in os.listdir(data_dirs[0]) if f.endswith('.bam')]
    files = [(filename, (region[0], int(region[1]), int(region[2]))) for filename, region in files]

    def get_file(chrm, pos):
        for filename, region in files:
            if region[0] == chrm and region[1] <= pos and region[2] > pos:
                return filename
        else:
            raise Error('Variant does not fit in any region')

    df['filename'] = df.apply(lambda x: get_file(x['chrm'], x['pos']), axis=1)
    return df

@functools.lru_cache # memoizes
def get_sam(data_dir, filename):
    # Could be shared with extract features?
    return pysam.AlignmentFile(os.path.join(data_dir, filename), "rb")

def generate_reads(samfile, x):
    # Could be shared with extract features?
    chrm, pos, ref, alt = x['chrm'], x['pos'], x['ref_allele'], x['alt_allele']

    for pileupcolumn in samfile.pileup(chrm, pos-1, pos+1, min_base_quality=20, min_mapping_quality=0):
        if pileupcolumn.pos != pos-1: continue
        for i, pileupread in enumerate(pileupcolumn.pileups):
            if pileupread.is_del or pileupread.is_refskip: continue
            allele = pileupread.alignment.query_sequence[pileupread.query_position]
            if allele != alt and allele != ref: continue

            yield pileupread, allele == alt
