import argparse
import gzip
import pickle
from pathlib import Path
from typing import Union

import pandas as pd

from extract_coordinates import cut_promoter, non_negative_num, promoter_file
from utils import read_fimo


def gzip_pickle(obj, name: Union[str, Path]):
    with open(name, 'wb') as f:
        with gzip.open(f, 'wb') as g:
            pickle.dump(obj, g, protocol=pickle.HIGHEST_PROTOCOL)


def fimo_coord_translate(fimo_path: Union[str, Path], promoter_path: Union[str, Path],
                         cut: Union[int, None] = None) -> pd.DataFrame:
    """
    Translate relative coordinates to genome coordinates
    :param fimo_path:
    :param promoter_path:
    :param cut:
    :return:
    """
    fimo = read_fimo(fimo_path)

    promoter = promoter_file(promoter_path)

    if cut is not None:
        promoter = cut_promoter(promoter, cut)

    merged = fimo.merge(promoter, left_on="sequence name", right_on="id", how="left")
    merged["start_x"] = merged["start_y"] + merged["start_x"] - 1
    merged["stop"] = merged["stop"] + merged["start_x"] - 1
    merged.drop(["id", "start_y", "end", "direction", "length"], inplace=True, axis=1)

    return merged


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("fimo", type=Path, help="FIMO result file")
    parser.add_argument("promoter", type=Path, help="file with promoters")
    parser.add_argument("-c", "--cut", type=non_negative_num, help="cut promoter to length")

    args = parser.parse_args()

    print(fimo_coord_translate(args.fimo, args.promoter, args.cut))
