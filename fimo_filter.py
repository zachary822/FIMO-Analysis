import csv
import multiprocessing as mp
from typing import Generator, Tuple

import numpy as np
import pandas as pd

from utils import read_fimo


def remove_duplicate(name: Tuple[str], group: pd.DataFrame) -> pd.DataFrame:
    group = group.sort_values(["p-value", "start", "stop"])

    def get_min_unique() -> Generator[pd.Series, None, None]:
        nonlocal group
        while not group.empty:
            p_min = group.iloc[0]
            yield p_min
            group = group.iloc[1:]
            selector: pd.Series = (group.start <= p_min.stop) & (group.stop >= p_min.start)
            if selector.any():
                group = group[~selector]

    total = pd.concat(get_min_unique(), axis=1).T
    total[['#pattern name', 'sequence name', 'strand']] = total[['#pattern name', 'sequence name', 'strand']].astype(
        str)
    total[['start', 'stop']] = total[['start', 'stop']].astype(np.int64)
    total[['score', 'p-value']] = total[['score', 'p-value']].astype(np.float64)

    return total


if __name__ == "__main__":
    import argparse
    from pathlib import Path

    parser = argparse.ArgumentParser(description="filter FIMO output for overlapping matches")
    parser.add_argument('-p', '--pickle', help='output pickle file', action='store_true')
    parser.add_argument('fimo', type=Path, help='FIMO file or data frame.')
    parser.add_argument('output', type=Path, help='output file')

    args = parser.parse_args()

    fimo = read_fimo(args.fimo)

    with mp.Pool() as pool:
        filtered: pd.DataFrame = pd.concat(
            pool.starmap(remove_duplicate, fimo.groupby(["strand", "#pattern name", "sequence name"]), chunksize=100))

    filtered.sort_index(inplace=True)

    if args.pickle:
        filtered.to_pickle(args.output)
    else:
        filtered.to_csv(args.output, sep='\t', index=False, quoting=csv.QUOTE_NONE)
