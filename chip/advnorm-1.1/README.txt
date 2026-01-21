VERSION
-------
v1.1


INTRODUCTION
------------
Advnorm is a command-line tool that performs plate normalization (a.k.a advanced normalization) 
of probesets with outlier plates and genotyping calling analysis using the normalized probeset signals.

This tool does 3 things:
    1. identifies probesets with outlier plates
    2. adjusts the signals of these problematic probesets using regression
    3. performs genotype calling analysis using the normalized signals

Please note that this algorithm assumes the signal distributions of a probeset in different plates are the same.
Therefore samples must be randomized across plates. Otherwise, the genotyping result may become worse.
Advnorm requires a genotyping batch with at least 3 plates to work.

The directory structure of this package is as follows:

    advnorm/
    ├── README.txt
    ├── advnorm.sh
    ├── bin
    │   ├── buildPlatemap.py
    │   ├── exclude_comment_header_mapper.py
    │   ├── find_outliers.R
    │   ├── PlateNorm_reduce.R
    │   ├── plot_per_plate.R
    │   ├── subsetNormalized.py
    │   └── transpose.pl
    └── demo
        ├── demo.sh
        ├── demo_input
        ├── demo_output
        └── library_files


The advnorm.sh script is the entry point run while the bin folder contains the other auxiliary scripts.
The demo folder contains the demo data for testing purpose. The demo.sh script contains the example command,
The demo_input and demo_output folder contain the input data and the expected output.
library_files contains the required analysis files to perform genotype calling analysis of the demo data.


REQUIREMENTS
------------
This tool has the following dependencies
* Perl 5.20 or later
* Python 2.7.0 or later
* R 3.1.0 or later
* SNPolisher 2.0.3 or later  (https://www.thermofisher.com/tw/zt/home/life-science/microarray-analysis/microarray-analysis-partners-programs/affymetrix-developers-network/affymetrix-devnet-tools.html)
* APT 2.11.3 or later (https://www.thermofisher.com/tw/zt/home/life-science/microarray-analysis/microarray-analysis-partners-programs/affymetrix-developers-network/affymetrix-power-tools.html)


USAGE
-----
Run advnorm.sh without options will print the version number and usage as follows:

VERSION:1.0
USAGE:
        bash advnorm.sh --summary-file                  SUMMARY_FILE \
                        --calls-file                    CALLS_FILE \
                        --report-file                   REPORT_FILE \
                        --trustcheck-file               TRUSTCHECK_FILE \
                        --analysis-files-path           LIBRARY_DIR \
                        --snp-priors-file               SNP_SPECIFIC_PRIORS \
                        --snp-specific-param-file       SNP_SPECIFIC_PARAMETERS \
                        --special-snps-file             SPECIAL_SNPS \
                        --ps2snp-file                   PS2SNP \
                        --output-dir                    OUTPUT_DIR
                        --plot

REQUIRED PARAMETER

    --summary-file      
            The summary file in text format

    --calls-file        
            The calls.txt file

    --report-file       
            The report.txt file

    --trustcheck-file               
            The trustcheck file produced by artifact-reduction

    --analysis-files-path           
            Default directory to search for analysis library files
            
    --snp-priors-file               
            The snp priors to use when genotyping

    --snp-specific-param-file       
            The SNP-specific parameters overriding standard parameters in genotyping

    --special-snps-file             
            The file specifying non-autosomal chromosome information

    --ps2snp-file                   
            The file containing the mapping of probeset IDs to SNP IDs

    --output-dir                    
            The output directory for results files


OPTIONAL PARAMETER

    --plot      Whether to draw the contrast-size plots of the probesets before and after normalization


OUTPUT 
------

The contents of the output folder are as follows:

demo_output/
├── apt-summary-genotype-axiom.errors
├── apt-summary-genotype-axiom.log
├── AxiomGT1.calls.txt
├── AxiomGT1.confidences.txt
├── AxiomGT1.report.txt
├── AxiomGT1.snp-posteriors.txt
├── AxiomGT1.summary.norm.txt
├── AxiomGT1.summary.txt
├── AxiomGT1.trustcheck.txt
├── Axiom_PMRA.r1.generic_prior.txt
├── Axiom_PMRA.r1.snp_specific_parameters.txt
├── NoMinorHom.ps
├── Other.ps
├── plot
│   ├── AX-59056814.after.png
│   ├── AX-59056814.before.png
│   ├── AX-98050496.after.png
│   ├── AX-98050496.before.png
│   ├── AX-98081668.after.png
│   ├── AX-98081668.before.png
│   ├── AX-98295632.after.png
│   ├── AX-98295632.before.png
│   └── pidFile.ps
├── PolyHighResolution.ps
├── ps-classification.log
├── ps-metrics.log
├── ps-metrics.txt
├── Ps.performance.txt
├── Recommended.ps
└── tmp
    ├── gender.tsv
    ├── plate_map.tsv
    ├── probeset.ps
    ├── ps_status.tsv
    └── sample_order.tsv


AxiomGT1.summary.txt is a subset of the original summary file containing the signals of probesets to be normalized.
AxiomGT1.summary.norm.txt is the normalized version of AxiomGT1.summary.txt.
AxiomGT1.calls.txt, AxiomGT1.confidences.txt, AxiomGT1.report.txt, AxiomGT1.snp-posteriors.txt, etc. are the genotype calling result of these normalized probesets.
The probeset.ps in the tmp folder is the list of probeset IDs with outlier plates.
The plot folder contains the cluster plot before and after normalization.
