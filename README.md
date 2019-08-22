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
the suffix ".sgn"

```
inputdir
  |
  file1.sgn
  file2.sgn
  ...
```


Their format should be

```
Gene<Tab>Group
<Gene identifier><Tab><Group number>
<Gene identifier><Tab><Group number>
...
```

**Note**

The -ue / --use_entrez switch is used if the input gene identifiers
are in Entrez format, otherwise Refseq is used by default.

**Pipeline Stages**

![Pipeline overview](images/pipeline_stages.png)

The FIRM pipeline consists of 3 stages:

  1. Motif discovery: Based on the genes contained in the input clusters, finds
     common motifs using Weeder
  2. miRNA discovery: Tries to find miRNAs from mirbase.org using miRvestigator
     an HMM (Hidden-Markov-Model) based approach
  3. Matching with prediction databases: This step tries to find the miRNAs obtained
     in step 2 within target prediction databases, and provides a score. By default this will be
     PITA and TargetScan.


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


## An Example Run

For this example, we will use the data in this projects "exp" directory.
It looks like this

```
exp
  |
  AD_Lung_Beer.sgn
```

and the contents of the "AD_Lung_Beer.sgn" file is:

```
Gene    Group
NM_000014       32
NM_000015       23
NM_000016       6
NM_000017       10
NM_000018       56
...
```

The first column of the input files are gene identifiers and the second column contain the
regulon/cluster numbers.

Assuming FIRM has been successfully installed on the user's system the pipeline can
be run by

```
  $ firm exp example-out
```

This command runs all stages of the pipeline, which can take a few hours.
The final and intermediate results will be found in the "example-out" directory
after successfully running the pipeline.

```
  $ ls example-out

combinedResults.csv  filtered_TargetScan.csv  mergedResults_PITA.csv        miRNA       seqs.txt
filtered_PITA.csv    m2m_standalone.pkl       mergedResults_TargetScan.csv  pssms.json
```

The combinedResults.csv file is the most important file of the pipeline output.

For our example data set it would look something like this:

```
Dataset,signature,miRvestigator.miRNA,miRvestigator.model,miRvestigator.mature_seq_ids,PITA.miRNA,PITA.percent_targets,PITA.P_Value,PITA.mature_seq_ids,TargetScan.miRNA,TargetScan.percent_targets,TargetScan.P_Value,TargetScan.mature_seq_ids
AD_Lung_Beer,19,hsa-miR-487b-3p,7mer_a1,MIMAT0003180,hsa-mir-371-3,0.02,0.0006879498,,hsa-mir-199b,0.09859154929577464,0.0001291377,
AD_Lung_Beer,31,hsa-miR-29a-3p_hsa-miR-29b-3p_hsa-miR-29c-3p,8mer,MIMAT0000086 MIMAT0000100 MIMAT0000681,hsa-mir-29 hsa-mir-29 hsa-mir-29b,0.2222222222222222,5.801169e-06,,hsa-mir-29b,0.12658227848101267,1.594762e-06,
AD_Lung_Beer,29,hsa-miR-8071,8mer,MIMAT0030998,hsa-mir-661,0.125,0.0002108273,MIMAT0003324,hsa-mir-100,0.03333333333333333,0.0001221282,
...
```
