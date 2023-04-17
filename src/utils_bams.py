import functools
import os
import pysam


def match_variants_to_filenames(df, data_dirs):
    # Could be shared with extract features?
    """
    Bam files are assumed to be merged and in region format. Thus, for each variant
    we match it to the appropriate file that contains that locus
    """
    def get_region(f):
        v = f.strip().split('.')[0].split('_')[-1].split('-')
        return (v[0], int(v[1]), int(v[2]))

    files = [(f, get_region(f)) for f in os.listdir(data_dirs[0]) if f.endswith('.bam')]

    #files = [(filename, (region[0], int(region[1]), int(region[2]))) for filename, region in files]

    def get_file(chrm, pos):
        for filename, region in files:
            if region[0] == chrm and region[1] <= pos and region[2] > pos:
                return "{}-{}-{}.bam".format(region[0], str(region[1]), str(region[2]))
        raise RuntimeError('Variant does not fit in any region: {} {}'.format(chrm, pos))

    try:
        temp = df.apply(lambda x: get_file(x['chrm'], x['pos']), axis=1)
        df['filename'] = temp
    except:
        print(temp)
    return df

@functools.lru_cache # memoizes
def get_sam(data_dir, filename):
    pysam.set_verbosity(0)
    ## TODO This is not a great way to handle things. Should probably switch to taking as input a metadata file which
    ## maps regions to filenames and all that logic can be handled externally in a (potentially) sample specific way
    try:
        return pysam.AlignmentFile(os.path.join(data_dir, filename), "rb")
    except:
        real_filename = [f for f in os.listdir(data_dir) if f.endswith("_"+filename)][0]
        return pysam.AlignmentFile(os.path.join(data_dir, real_filename), "rb")

def generate_reads(samfile, x):
    chrm, pos, ref, alt = x['chrm'], x['pos'], x['ref_allele'], x['alt_allele']

    for pileupcolumn in samfile.pileup(chrm, pos-1, pos+1, min_base_quality=20, min_mapping_quality=0, ignore_orphans=False):
        if pileupcolumn.pos != pos-1: continue
        for i, pileupread in enumerate(pileupcolumn.pileups):
            if pileupread.is_del or pileupread.is_refskip: continue
            allele = pileupread.alignment.query_sequence[pileupread.query_position]
            if allele != alt and allele != ref: continue

            yield pileupread, allele == alt

def generate_reads_w_normal(samfile, x, labels):
    chrm, pos, ref, alt = x['chrm'], x['pos'], x['ref_allele'], x['alt_allele']
    for pileupcolumn in samfile.pileup(chrm, pos-1, pos+1, min_base_quality=20, min_mapping_quality=0, ignore_orphans=False):
        if pileupcolumn.pos != pos-1: continue
        for i, pileupread in enumerate(pileupcolumn.pileups):
            if pileupread.is_del or pileupread.is_refskip: continue
            allele = pileupread.alignment.query_sequence[pileupread.query_position]
            if allele != alt and allele != ref: continue

            RG = pileupread.alignment.get_tag('RG')
            cell_id = '_'.join(RG.split('_')[1:-2])

            try:
                is_hq_normal = labels.loc[cell_id].is_high_quality_normal_cell
            except:
                is_hq_normal = False

            yield pileupread, allele == alt, is_hq_normal
