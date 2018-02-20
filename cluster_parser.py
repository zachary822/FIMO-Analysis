import glob
import os.path as path
import re
from typing import Generator, List, Tuple

DAP_REGEX = re.compile(r"^Collection_\d_m\d+_(\w+)\.(\w+?_col(?:amp)?(?:_a)?)", re.I)


def parse_cluster_file(cluster_path: str):
    with open(cluster_path) as f:
        for line in f:
            cluster, motifs = line.split("\t")
            motifs = list(map(str.strip, motifs.split(",")))
            yield cluster, motifs


def get_dap_motifs(motifs: List) -> List:
    return [m.groups() for m in filter(None, map(DAP_REGEX.match, motifs))]


def get_dap_paths(motifs: List, dap_base_path):
    for m, n in get_dap_motifs(motifs):
        yield from glob.glob(path.join(dap_base_path, m, n, "**/*.narrowPeak"), recursive=True)


def get_cluster_dap_files(tab: str, dap_base_path: str) -> Generator[Tuple[str, List[str]], None, None]:
    """
    Get DAP peak files from cluster-to-gene-family file
    :param tab: file with cluster-to-gene-family relationships
    :param dap_base_path: folder with DAP peaks
    :return:
    """
    for cluster, motifs in parse_cluster_file(tab):
        motif_files = list(get_dap_paths(motifs, dap_base_path))
        if motif_files:
            yield cluster, motif_files


if __name__ == '__main__':
    import argparse
    from pathlib import Path

    parser = argparse.ArgumentParser()
    parser.add_argument('clusters', type=Path, help='file with cluster members')
    parser.add_argument('peaks', type=Path, help='narrowPeak folder')
    args = parser.parse_args()

    for cluster, motif_files in get_cluster_dap_files(args.clusters,
                                                      args.peaks):
        print(cluster, motif_files)
