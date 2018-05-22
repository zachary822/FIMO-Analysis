# FIMO-Tools

## Introduction

Collection of scripts to filter [FIMO](http://meme-suite.org/tools/fimo) output.

## Installation

```bash
pip3 install -r requirements.txt
```

## Usage

### FIMO Filter

Removes overlapping matches, only keeping the one with the lowest p\-value.

```bash
python3 fimo_filter.py input.fimo output.fimo
```

## LICENSE

[MIT](LICENSE)
