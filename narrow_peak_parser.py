from concurrent.futures import ProcessPoolExecutor
from functools import partial
from pathlib import Path
from typing import Union

import pandas as pd

from cluster_parser import get_cluster_dap_files
from fimo_adjust import fimo_coord_translate


def read_narrow_peak(narrow_peak_file: Union[str, Path], top: Union[int, None] = 600) -> pd.DataFrame:
    peaks = pd.read_csv(
        narrow_peak_file,
        sep="\t", header=None,
        names=['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'signalValue', 'pValue', 'qValue', 'peak'],
        engine='c'
    )

    if top is not None:
        peaks = peaks.sort_values('score').head(top)

    return peaks


def p_values_from_peaks(narrow_peak_file: str, grouped_fimo: pd.DataFrame, cluster_name: str,
                        top: Union[int, None] = 600) -> pd.Series:
    peaks = read_narrow_peak(narrow_peak_file, top)

    hit_pvalues = []

    # which peaks are actually in the promoter
    # boxplot or violin plot for each cluster, using all files in cluster, not just random.
    # look for DAP peaks in promoters
    for name, group in peaks.groupby(['chrom']):
        gc = grouped_fimo.get_group((cluster_name, name))
        for idx, row in group.iterrows():
            hit = gc[(row['chromStart'] <= gc.start_x) & (row['chromEnd'] >= gc.stop)]
            if not hit.empty:
                hit_pvalues.append(hit['p-value'])

    return pd.concat(hit_pvalues)
