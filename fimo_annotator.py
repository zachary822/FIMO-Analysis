"""
Associate filtered FIMO pickle with gene annotations

Removes matches that are not part of gene or promoter region.
"""
import re

import pandas as pd


def get_upstream_loc(group, length=1000):
    strand, chrom = group.name
    if strand == "+":
        g = group[3] - length
        return g.clip_lower(1)
    elif strand == "-":
        return group[4] + length


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("gene_annotation", type=str, help="gff3 file of gene annotations")
    parser.add_argument("filtered_matches", type=str, help="FIMO matches that are filtered for duplicated matches")
    parser.add_argument("--output", type=str, default="annotated.pickle.gz")

    args = parser.parse_args()

    matches = pd.read_pickle(args.filtered_matches)
    match_group = matches.groupby(['strand', 'sequence name'])

    gff = pd.read_csv(
        args.gene_annotation,
        sep="\t", header=None, engine='c', comment='#')

    genes = gff[gff[2] == "gene"].copy(False)

    genes[8] = genes[8].str.extract(r"id=(.+?);", re.I, expand=True)
    genes[9] = genes.groupby([6, 0], as_index=False).apply(
        get_upstream_loc).reset_index(level=0, drop=True).sort_index()
    genes.reset_index(drop=True, inplace=True)
    genes.drop([1, 2, 5, 7], axis=1, inplace=True)


    def annotate_filtered():
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
                yield temp_mg


    annotated: pd.DataFrame = pd.concat(annotate_filtered())
    annotated.reset_index(inplace=True)
    annotated.rename(columns={'index': 'match_id'}, inplace=True)
    annotated.sort_values(['gene', 'match_id'], inplace=True)
    annotated['gene'] = annotated['gene'].astype('category')
    annotated.set_index('gene', inplace=True)

    annotated.to_pickle(args.output, compression='gzip')
