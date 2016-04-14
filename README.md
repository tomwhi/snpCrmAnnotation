# Summary
This is a repository of the source code that was used to perform annotation of transcription factor binding effects by prostate cancer risk SNPs, published in the research article "Gene regulatory mechanisms underpinning prostate cancer susceptibility", Whitington et al., Nature Genetics 2016 (doi:10.1038/ng.3523).

This software is provided with the goal of make the source code available to researchers wishing to replicate or extend upon this work. However, please note that the annotation script itself will not perform annotation without the relevant ChIP-seq resources and the PWM resource, which includes a MySQL database storing metadata for the PWMs.

# Required resources
To function, the SNP annotation tool requires a ChIP-seq resource in a structured set of folders, and a PWM resource, with metadata for the PWMs stored in a MySQL database.

Neither of these resources are currently included in this repository due to the size limitations. To produce a ChIP-seq resource with the required structure, one would need to place the mapped bam files for the ChIP-seq data into a nested folder heirarchy similar structured according to <StudyName>/<CellType>/<AntibodyTarget>/<ExperimentID>. E.g.:
Huang2014_ERP004190/VCaP/FOXA1/ERX332517/merged_ERR359745.bam

One would then need to use the tool "runPeakCalling.py" to generate called peaks and tdf files as required for downstream analysis as was performed in the manuscript.

The schema of the MySQL database (without the large amount of data included) is included as a database dump, in "PWM_db_schema.sql".

# Dependencies
The resource requires a hacked version of FIMO from the MEME package, in which an additional command line option "--allow-overlap" has been added. A git repository containing this hack is available at https://github.com/tomwhi/meme_WithFimoHack

# Usage
"snpCrmAnnotation.py" is used to run annotation of a set of SNPs of interest
"runPeakCalling.py" is used to process a given ChIP-seq dataset such that it can then be used by the pipeline

# Citation
If you publish work using these algorithms, please cite the following paper:
"Gene regulatory mechanisms underpinning prostate cancer susceptibility", Whitington et al., Nature Genetics 2016 (doi:10.1038/ng.3523)

# Future plans
This project is not currently being developed by the original author, but feel free to fork the project.
