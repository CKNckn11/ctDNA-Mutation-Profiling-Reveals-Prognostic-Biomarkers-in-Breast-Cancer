# ctDNA-Mutation-Profiling-Reveals-Prognostic-Biomarkers-in-Breast-Cancer
ctDNA mutation profiling identifies plasma-derived prognostic biomarkers in breast cancer and supports a seven-gene survival risk model for liquid biopsy-based risk stratification.
## Workflow

### cfDNA Workflow

![cfDNA study workflow](Figure1/workflow_cfDNA.png)

This workflow summarizes the overall cfDNA study design, including sample collection, mutation profiling, clinical integration, biomarker screening, and survival risk model development.

### Mutect2 Analysis Workflow

![Mutect2 analysis workflow](Figure1/workflow_Mutect2.png)

This workflow summarizes the Mutect2-based somatic mutation calling process used for cfDNA BAM files, including sample metadata preparation, BAM indexing, Mutect2 calling, contamination estimation, filtering, ANNOVAR annotation, and final mutation table aggregation.

## cfDNA Mutect2 Pipeline Generator

The cfDNA Mutect2 pipeline generator used in this project is maintained as a separate repository:

[cfDNA_Mutect2_generator](https://github.com/CKNckn11/cfDNA_Mutect2_generator)

It provides a lightweight shell-script generator for tumor-only cfDNA somatic mutation calling, including Mutect2 calling, contamination estimation, filtering, ANNOVAR annotation, and final mutation table aggregation from a YAML configuration file.

This tool does not run all analysis jobs automatically inside Python. Instead, it generates an independent shell-based project that can be submitted and monitored manually on a local server or computing cluster.

```bash
$ python generate_pipeline.py -h
usage: cfdna-mutect2-generate [-h] [-v] --config CONFIG [--outdir OUTDIR]

cfDNA Mutect2 Pipeline Generator (Version = 0.1.0): Generate tumor-only cfDNA Mutect2 shell scripts from a YAML config.

optional arguments:
  -h, --help       show this help message and exit
  -v, --version    show the version of cfdna-mutect2-generate and exit.
  --config CONFIG  YAML config file.
  --outdir OUTDIR  Override project.outdir from the YAML config.
```
