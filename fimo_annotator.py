"""
Associate filtered FIMO pickle with gene annotations

Removes matches that are not part of gene or promoter region.
"""
import argparse
import re
from functools import partial

import pandas as pd
from tqdm import tqdm

from utils import read_fimo


def get_upstream_loc(group, length=1000):
    strand, chrom = group.name
    if strand == "+":
        return (group[3] - length).clip(lower=1)
    elif strand == "-":
        return group[4] + length


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("gene_annotation", type=str, help="gff3 file of gene annotations")
    parser.add_argument("matches", type=str, help="FIMO matches")
    parser.add_argument("-p", "--pseudo", action="store_true", help="include pseudogenes")
    parser.add_argument("-o", "--output", type=argparse.FileType('wb'), help="output file", required=True)
    parser.add_argument("-u", "--upstream-length", type=int, default=1000,
                        help="define promoter region (default: 1000bp upstream)")

    args = parser.parse_args()

    matches = read_fimo(args.matches)
    match_group = matches.groupby(['strand', 'sequence name'])

    gff = pd.read_csv(
        args.gene_annotation,
        sep="\t", header=None, engine='c', comment='#', dtype={0: str})

    if args.pseudo:
        genes = gff[gff[2].str.match(r'^(?:pseudo)?gene$', flags=re.I)].copy()
    else:
        genes = gff[gff[2] == 'gene'].copy()

    genes[8] = genes[8].str.extract(r"ID=([^;]+)", expand=False).str.strip().str.upper()
    genes[9] = genes.groupby([6, 0], as_index=False).apply(
        partial(get_upstream_loc, length=args.upstream_length)).reset_index(level=0, drop=True).sort_index()
    genes = genes.reset_index(drop=True).drop([1, 2, 5, 7], axis=1)


    def annotate_filtered():
        with tqdm(total=genes.shape[0], desc="Genes Filtered") as pbar:
            for name, group in genes.groupby([6, 0]):
                mg = match_group.get_group(name)
                for row in group.itertuples(index=False, name=None):
                    if row[3] == "+":
                        temp_mg = mg[(mg["start"] >= row[5]) & (mg["stop"] <= row[2])].copy(False)
                        temp_mg['dist'] = (temp_mg['start'] - row[1]).where(lambda x: x < 0, lambda x: x + 1)
                    elif row[3] == "-":
                        temp_mg = mg[(mg["start"] >= row[1]) & (mg["stop"] <= row[5])].copy(False)
                        temp_mg['dist'] = (row[2] - temp_mg['stop']).where(lambda x: x < 0, lambda x: x + 1)
                    else:
                        raise Exception("undefined strand {}".format(row[3]))

                    temp_mg['gene'] = row[4]
                    pbar.update(1)
                    yield temp_mg


    annotated: pd.DataFrame = (pd.concat(annotate_filtered())
                               .reset_index()
                               .rename(columns={'index': 'match_id'})
                               .sort_values(['gene', 'match_id'])
                               .set_index('gene'))

    annotated.to_pickle(args.output, compression='gzip')
