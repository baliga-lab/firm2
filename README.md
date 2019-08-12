# firm2 - integrated FIRM

## Description

This is a version of FIRM that was created with the goal of being
easier to use for the general user.
The FIRM pipeline is made available through a single command line
command.


## Installation

Install the firm from pypi.org here:

```
  $ pip install firm
```

FIRM relies on a customized version of Weeder 1.4.2 that is included in
the project source code.

Obtain the weeder-1.4.3.tar.gz from

https://github.com/baliga-lab/firm2

which is distributed as a autoconf project, so you can use the
common

```
configure / make / make install
```

sequence to install the Weeder suite.


## Usage

### Main firm tool

```
  $ firm -h
```

```
usage: firm [-h] [-ue] [-t TMPDIR] expdir outdir

firm - Run FIRM pipeline

positional arguments:
  expdir                expression input directory
  outdir                output directory

optional arguments:
  -h, --help            show this help message and exit
  -ue, --use_entrez     input file uses entrez IDs instead of RefSeq
  -t TMPDIR, --tmpdir TMPDIR
                        temporary directory

```


Typically the firm pipeline will expect an input directory ("expdir")
and an output directory.
The input directory should contain one or more tab-separated files with
the suffix ".sgn". Their format should be

```
Gene<Tab>Group
<Gene identifier><Tab><Group number>
<Gene identifier><Tab><Group number>
...
```

**Note**

The -ue / --use_entrez switch is used if the input gene identifiers
are in Entrez format, otherwise Refseq is used by default.

### firm-convertminer tool

This is a tool to conver MINER regulon files in JSON format to
FIRMs expected input format.

```
  $ firm-convertminer -h
```

```
usage: firm-convertminer [-h] regulons mappings outdir

firm-convertminer - convert MINER input files to a FIRM input directory

positional arguments:
  regulons    regulons file (JSON format)
  mappings    mappings file
  outdir      output directory

optional arguments:
  -h, --help  show this help message and exit
```
