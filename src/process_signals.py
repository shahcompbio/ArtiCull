

from os import path
import pandas as pd
import subprocess

signals_dir="/data1/shahs3/junobackup/isabl_data_lake/analyses/75/59/37559/results/"
signals_result = path.join(signals_dir, 'signals.Rdata')
signals_cns = path.join(signals_dir, 'hscn.csv.gz')
output_dir="temp"

###
# From signals.Rdata extracts the cell to clone mapping
# There's not a great way to read Rdata in python, so I call out to an R script to do this
###
c2c_file="{}/cell2clone.tsv".format(output_dir)
# Gets the directory of the articull src folder
directory = path.dirname(path.realpath(path.expanduser(__file__)))
command = "Rscript {}/extract_cell2clone.R {} {}".format(directory, signals_result, c2c_file)
subprocess.check_call(command.split(' '))
ids_map = pd.read_table(c2c_file)['clone_id']

###
# From hscn, gets per cell copy number, then uses clonemap to map cells to clones
###
fields = ['chr','start', 'end', 'clone', 'A', 'B']
try:
    # Depending on the version of signals, either 'A'/'B' and 'Maj'/'Min' may be the columns of interest
    fields = ['chr','start', 'end', 'clone', 'A', 'B']
    cn_df = pd.read_table(signals_cns, sep=',', fields=fields)
    cn_df['Maj'] = cn_df['A']
    cn_df['Min'] = cn_df['B']
except:
    fields = ['chr','start', 'end', 'clone', 'Maj', 'Min']
    cn_df = pd.read_table(signals_cns, sep=',', fields=fields)

# map cells to clones
cn_df['clone']= cn_df['cell_id'].map(ids_map)

# get clone copy number
clone_cns['CN'] = clone_cns['Maj'] + clone_cns['Min']
clone_cns = cn_df.groupby(['chr', 'start', 'end', 'clone'])[['CN']].median() \
    .reset_index().pivot(index=['chrm', 'start', 'end'], columns = 'clone', values = 'CN')

clone_cns.to_tsv("{}/clone_cns.tsv".format(output_dir))
