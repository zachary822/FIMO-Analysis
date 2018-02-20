import argparse
import re
from typing import Generator, Tuple
import sys

import pandas as pd

COMMENT_LINE_REGEX = re.compile(
    r'^>(?P<id>AT[1-5MC]G\d{5}) \| (?P<chr>chr[1-5MC]):(?P<start>\d+)-(?P<end>\d+) (?P<direction>REVERSE|FORWARD) '
    r'LENGTH=(?P<length>\d+)$',
    re.I)


def non_negative_num(s):
    value = int(s)
    if value < 0:
        raise argparse.ArgumentTypeError("{} is not a non-negative integer".format(s))
    return value


def parse_promoter_file(path: str) -> Generator[Tuple[str], None, None]:
    with open(path) as f:
        for l in f:
            if l.startswith(">"):
                yield COMMENT_LINE_REGEX.match(l).groups()


def promoter_file(path: str) -> pd.DataFrame:
    data = pd.DataFrame(parse_promoter_file(path),
                        columns=['id', 'chr', 'start', 'end', 'direction', 'length'])
    data[['start', 'end', 'length']] = data[['start', 'end', 'length']].astype(int)

    return data


def cut_promoter(data: pd.DataFrame, cut_size: int = 500) -> pd.DataFrame:
    """
    Adjust start, end coordinates and length to desired promoter size
    :param data:
    :param cut_size:
    :return:
    """
    data.loc[((data.end - data.start) >= cut_size), 'length'] = cut_size
    data.loc[((data.end - data.start) >= cut_size) & (data.direction == 'FORWARD'), 'end'] = \
        data.loc[((data.end - data.start) >= cut_size) & (data.direction == 'FORWARD'), 'start'] + cut_size - 1
    data.loc[((data.end - data.start) >= cut_size) & (data.direction == 'REVERSE'), 'start'] = \
        data.loc[((data.end - data.start) >= cut_size) & (data.direction == 'REVERSE'), 'end'] - cut_size + 1
    return data


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("promoter", type=str, help="file with promoters")
    parser.add_argument("-c", "--cut", type=non_negative_num, help="cut promoter to length")

    args = parser.parse_args()

    data = promoter_file(args.promoter)

    if args.cut is not None:
        data = cut_promoter(data, args.cut)

    data.to_csv(sys.stdout)
