import sys
from itertools import zip_longest
from typing import List

import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import fdrcorrection


def promoter_size(s):
    value = int(s)
    if value > 1000:
        raise argparse.ArgumentTypeError("Promoter size has to be less than or equal to 1000.")

    return value


def probability(s):
    value = float(s)
    if value < 0 or value > 1:
        raise argparse.ArgumentTypeError("Probability cannot be less than 0 or greater than 1.")

    return value


def get_lists(path: str) -> List[pd.Series]:
    data = pd.read_csv(path, header=None, squeeze=True)
    data = data.str.strip()

    names = data.index[data.str.startswith(">")]

    return [data[start + 1:end].rename(data[start].strip(">")) for start, end in zip_longest(names, names[1:])]


if __name__ == "__main__":
    import argparse
    import shutil
    import os

    width, height = shutil.get_terminal_size()

    parser = argparse.ArgumentParser()
    parser.add_argument("genelist", type=str, help="gene list in FASTA format")
    parser.add_argument("-A", "--annotated", type=str, help="annotated matches", default="annotated.pickle.gz")
    parser.add_argument("-P", "--promoter", type=promoter_size, help="limit promoter size")
    parser.add_argument("-p", "--p-value", type=probability, help="p-value cutoff")
    parser.add_argument("-a", "--alpha", type=probability, help="alpha for enrichment", default=0.05)
    parser.add_argument("--hide-rejected", action="store_false", help="only display significant results")
    parser.add_argument("-o", "--output", type=str, help="output folder name")

    group = parser.add_argument_group(title="limit search")
    mutex = group.add_mutually_exclusive_group()
    mutex.add_argument("--promoter-only", action="store_true", help="limit search to promoter")
    mutex.add_argument("--promoter-overlap", action="store_true", help="search matches that overlap promoter")
    mutex.add_argument("--gene-only", action="store_true", help="limit search to gene body")
    mutex.add_argument("--gene-overlap", action="store_true", help="search matches that overlap gene body")

    args = parser.parse_args()

    annotated = pd.read_pickle(args.annotated)

    if args.p_value:
        annotated = annotated[annotated['p-value'] < args.p_value]

    if args.promoter:
        annotated = annotated[annotated['dist'] >= -args.promoter]

    if args.promoter_only:
        annotated = annotated[(annotated['stop'] - annotated['start'] + annotated['dist']) < 0]

    if args.promoter_overlap:
        annotated = annotated[annotated['dist'] < 0]

    if args.gene_only:
        annotated = annotated[annotated['dist'] > 0]

    if args.gene_overlap:
        annotated = annotated[(annotated['stop'] - annotated['start'] + annotated['dist']) > 0]

    ann_dedup = annotated.drop_duplicates('match_id')
    cluster_size = ann_dedup.groupby('#pattern name').size()

    print("total matches: {0[0]:,d}\n".format(ann_dedup.shape), file=sys.stderr)


    def get_list_enrichment(gene_list, alpha: float = 0.05, show_rejected: bool = True) -> pd.DataFrame:
        list_cluster_dedup = annotated[annotated.index.isin(gene_list)].drop_duplicates('match_id')
        list_cluster_size = list_cluster_dedup.groupby('#pattern name').size()

        def cluster_fisher(row):
            return fisher_exact(
                [[row[0], row[1] - row[0]],
                 [list_cluster_dedup.shape[0] - row[0],
                  ann_dedup.shape[0] - list_cluster_dedup.shape[0] - row[1] + row[0]]],
                alternative='greater')[1]

        p_values = pd.concat([list_cluster_size, cluster_size],
                             axis=1).fillna(0).apply(cluster_fisher, axis=1).sort_values()
        reject, adj_p = fdrcorrection(p_values, alpha=alpha, is_sorted=True)

        if not show_rejected:
            p_values = p_values[reject]
            adj_p = adj_p[reject]

        adj_p = pd.Series(adj_p, index=p_values.index)
        return pd.concat([p_values, adj_p], axis=1).rename(columns={0: 'p', 1: 'adj_p'})


    if args.output:
        os.mkdir(args.output)

    for test_list in get_lists(args.genelist):
        result = get_list_enrichment(test_list, alpha=args.alpha, show_rejected=args.hide_rejected)
        if args.output:
            result.to_csv(os.path.join(args.output, "{}.csv".format(test_list.name)))
        else:
            print(test_list.name)
            result.to_csv(sys.stdout, sep="\t")
            print("-" * (width if width <= 80 else 80))
