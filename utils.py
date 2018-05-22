from pathlib import Path
from typing import IO, Union

import pandas as pd


def read_fimo(fimo_path: Union[str, IO, Path]) -> pd.DataFrame:
    try:
        return pd.read_pickle(fimo_path)
    except KeyError:
        return pd.read_csv(fimo_path, sep="\t", header=0,
                           usecols=['#pattern name', 'sequence name', 'start', 'stop', 'strand', 'score',
                                    'p-value'], engine='c')
