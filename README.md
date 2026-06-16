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
