#!/bin/bash
VERSION='1.1'

main(){

    ##
    ## parse arguments
    ##
    parse_args "$@"


    ##
    ## Preparation
    ##
    init_env
    gen_plate_map 
    gen_sample_order
    gen_gender

    ##
    ## Find Outliers
    ##
    find_outlier

    ##
    ## Normalization
    ##
    subset_data
    normalize_outlier

    ##
    ## Genotyping
    ##
    genotype
    ps_metrics
    ps_classification

    ##
    ## Cluster Plot
    ##
    if [[ $PLOT = 'TRUE' ]]; then
        plot_clusters
    fi
}


plot_clusters() {
    BEFORE_DIR=$(dirname $SUMMARY_FILE)
    AFTER_DIR=$OUTPUT_DIR
    PLOT_DIR="$OUTPUT_DIR/plot"

    mkdir -p $PLOT_DIR

    echo 'probeset_id' > $PLOT_DIR/pidFile.ps;
    cat $PROBESETS >> $PLOT_DIR/pidFile.ps;

    CMD="
        Rscript $BIN/plot_per_plate.R
            $PLATE_MAP
            $BEFORE_DIR/AxiomGT1.calls.txt
            $BEFORE_DIR/AxiomGT1.summary.txt
            $PLOT_DIR/pidFile.ps
            $PLOT_DIR
            'before'
    "

    execute "$CMD"

    CMD="
        Rscript $BIN/plot_per_plate.R
            $PLATE_MAP
            $AFTER_DIR/AxiomGT1.calls.txt
            $AFTER_DIR/AxiomGT1.summary.norm.txt
            $PLOT_DIR/pidFile.ps
            $PLOT_DIR
            'after'
    "

    execute "$CMD"
    
}


parse_args() {
    if [[ $# -eq 0 ]]; then
        usage
    fi

    ARGUMENTS=$(
        getopt -a                            \
               -n advnorm.sh                 \
               -o ''                         \
               -l summary-file:              \
               -l calls-file:                \
               -l report-file:               \
               -l trustcheck-file:           \
               -l analysis-files-path:       \
               -l snp-priors-file:           \
               -l snp-specific-param-file:   \
               -l special-snps-file:         \
               -l ps2snp-file:               \
               -l output-dir:                \
               -l plot                       \
               -- "$@"                                 
    )

    EXIT_STATUS=$?

    if [[ $EXIT_STATUS -ne 0 ]]; then
        usage
    fi

    eval set -- "$ARGUMENTS"
    unset ARGUMENTS

    while true; do
        case "$1" in 
            '--summary-file') SUMMARY_FILE=$2
                shift 2
                ;;
            '--calls-file') CALLS_FILE=$2
                shift 2
                ;;
            '--report-file') REPORT_FILE=$2
                shift 2
                ;;
            '--trustcheck-file') TRUSTCHECK_FILE=$2
                shift 2
                ;;
            '--analysis-files-path') LIBRARY_DIR=$2
                shift 2
                ;;
            '--snp-priors-file') SNP_SPECIFIC_PRIORS=$2
                shift 2
                ;;
            '--snp-specific-param-file') SNP_SPECIFIC_PARAMETERS=$2
                shift 2
                ;;
            '--special-snps-file') SPECIAL_SNPS=$2
                shift 2
                ;;
            '--ps2snp-file') PS2SNP=$2
                shift 2
                ;;
            '--output-dir') OUTPUT_DIR=$2
                shift 2
                ;;
            '--plot') PLOT='TRUE'
                shift 1
                ;;
            '--') 
                shift; break 
                ;;
            *) echo "Unexpected option: $1"
                usage
                ;;
        esac
    done

    [[ -z $SUMMARY_FILE ]] && usage
    [[ -z $CALLS_FILE ]] && usage
    [[ -z $REPORT_FILE ]] && usage
    [[ -z $TRUSTCHECK_FILE ]] && usage

    [[ -z $LIBRARY_DIR ]] && usage
    [[ -z $SNP_SPECIFIC_PRIORS ]] && usage
    [[ -z $SNP_SPECIFIC_PARAMETERS ]] && usage
    [[ -z $SPECIAL_SNPS ]] && usage
    [[ -z $PS2SNP ]] && usage

    [[ -z $OUTPUT_DIR ]] && usage
    [[ -z $PLOT ]] && PLOT='FALSE'
    TMP_DIR="$OUTPUT_DIR/tmp"

    BIN="$(dirname $0)/bin"

    PLATE_MAP="$TMP_DIR/plate_map.tsv"
    SAMPLE_ORDER="$TMP_DIR/sample_order.tsv"
    PS_STATUS="$TMP_DIR/ps_status.tsv"
    PROBESETS="$TMP_DIR/probeset.ps"
    GENDER="$TMP_DIR/gender.tsv"
}

usage() {
    printf "\nUSAGE:\n\tbash advnorm.sh --summary-file\t\t\tSUMMARY_FILE \\
                        --calls-file\t\t\tCALLS_FILE \\
                        --report-file\t\t\tREPORT_FILE \\
                        --trustcheck-file\t\tTRUSTCHECK_FILE \\
                        --analysis-files-path\t\tLIBRARY_DIR \\
                        --snp-priors-file\t\tSNP_SPECIFIC_PRIORS \\
                        --snp-specific-param-file\tSNP_SPECIFIC_PARAMETERS \\
                        --special-snps-file\t\tSPECIAL_SNPS \\
                        --ps2snp-file\t\t\tPS2SNP \\
                        --output-dir\t\t\tOUTPUT_DIR\n\n"
  exit 2
}

init_env() {
    rm -rf $OUTPUT_DIR
    mkdir -p $OUTPUT_DIR
    mkdir -p $TMP_DIR
}

ps_metrics() {
    CMD="
        ps-metrics
            --call-file $OUTPUT_DIR/AxiomGT1.calls.txt
            --summary-file $OUTPUT_DIR/AxiomGT1.summary.norm.txt
            --report-file $OUTPUT_DIR/AxiomGT1.report.txt
            --output-dir $OUTPUT_DIR
            --log-file $OUTPUT_DIR/ps-metrics.log
            --posterior-file $OUTPUT_DIR/AxiomGT1.snp-posteriors.txt
            --metrics-file $OUTPUT_DIR/ps-metrics.txt
            --special-snps $SPECIAL_SNPS
    "

    execute "$CMD"
}

ps_classification() {
    CMD="
        ps-classification
            --metrics-file $OUTPUT_DIR/ps-metrics.txt
            --output-dir $OUTPUT_DIR
            --log-file $OUTPUT_DIR/ps-classification.log
            --species-type human
            --ps2snp-file $PS2SNP
    "

    execute "$CMD"
}

genotype() {
    CMD="
        apt-summary-genotype-axiom
            --analysis-files-path $LIBRARY_DIR
            --genotyping-node:brlmmp-clustertype 2
            --genotyping-node:brlmmp-MS 0.15
            --genotyping-node:brlmmp-SB 0.75
            --genotyping-node:brlmmp-CSepPen 0.1
            --genotyping-node:brlmmp-CSepThr 4
            --genotyping-node:brlmmp-copytype -1
            --genotyping-node:brlmmp-ocean 0.00001
            --out-dir $OUTPUT_DIR
            --log-file $OUTPUT_DIR/apt-summary-genotype-axiom.log
            --snp-posteriors-output true
            --genotyping-node:snp-posteriors-output-file $OUTPUT_DIR/AxiomGT1.snp-posteriors.txt
            --summary-input-file $OUTPUT_DIR/AxiomGT1.summary.norm.txt 
            --artifact-reduction-trustcheck-file $OUTPUT_DIR/AxiomGT1.trustcheck.txt
            --read-genders $GENDER
            --genotyping-node:snp-priors-input-file $OUTPUT_DIR/$(basename $SNP_SPECIFIC_PRIORS)
            --snp-specific-param-file $OUTPUT_DIR/$(basename $SNP_SPECIFIC_PARAMETERS)
            --special-snps $SPECIAL_SNPS
    "

    execute "$CMD"


}

gen_gender() {

    CMD="
        echo $'cel_files\tgender' > $GENDER;
        grep -v '^#' $REPORT_FILE 
            | sed 1d
            | cut -f1,2 
            >> $GENDER
    "

    execute "$CMD"
}

normalize_outlier() {
    CMD="
        grep -m 1 ^probeset_id
            $SUMMARY_FILE
        > $OUTPUT_DIR/AxiomGT1.summary.norm.txt
    "

    execute "$CMD"
 
    CMD="
        cat $OUTPUT_DIR/AxiomGT1.summary.txt
            | Rscript $BIN/PlateNorm_reduce.R
                $PROBESETS
                $PLATE_MAP
                $SAMPLE_ORDER
            >> $OUTPUT_DIR/AxiomGT1.summary.norm.txt
    "

    execute "$CMD"
}

subset_data() {

    # summary
    CMD="
        python2 $BIN/subsetNormalized.py
            $PROBESETS
            $SUMMARY_FILE
            > $OUTPUT_DIR/$(basename $SUMMARY_FILE)
    "

    execute "$CMD"

    # trustcheck
    CMD="
        python2 $BIN/subsetNormalized.py
            $PROBESETS
            $TRUSTCHECK_FILE
            > $OUTPUT_DIR/$(basename $TRUSTCHECK_FILE)
    "

    execute "$CMD"

    # snp specific priors
    CMD="
        python2 $BIN/subsetNormalized.py
            $PROBESETS
            $SNP_SPECIFIC_PRIORS
            > $OUTPUT_DIR/$(basename $SNP_SPECIFIC_PRIORS)
    "

    execute "$CMD"

    # snp specific parameters
    CMD="
        python2 $BIN/subsetNormalized.py
            $PROBESETS
            $SNP_SPECIFIC_PARAMETERS
            > $OUTPUT_DIR/$(basename $SNP_SPECIFIC_PARAMETERS)
    "

    execute "$CMD"

}


find_outlier() {
    CMD="
        python2 $BIN/exclude_comment_header_mapper.py < $CALLS_FILE
            | Rscript $BIN/find_outliers.R
                $PLATE_MAP
                $SAMPLE_ORDER
            > $PS_STATUS
    "

    execute "$CMD"
 
    CMD="
        egrep 'OUTLIER PLATE DETECTED' $PS_STATUS
            | grep -v 'NO OUTLIER PLATE DETECTED'
            | cut -f 1
            > $PROBESETS
    "

    execute "$CMD"
}



gen_sample_order() {
    CMD="
        grep '^probeset_id' -m1 $CALLS_FILE
            | $BIN/transpose.pl 2>/dev/null
            | tail -n +2 
         > $SAMPLE_ORDER
    "

    execute "$CMD"
}

gen_plate_map() {
    CMD="
        python2 $BIN/buildPlatemap.py
            $REPORT_FILE
        > $PLATE_MAP
    "

    execute "$CMD"
}

execute() {
    echo "$1"
    eval $1

    EXIT_STATUS=$?

    if [[ $EXIT_STATUS -ne 0 ]]; then
        exit $EXIT_STATUS
    fi

}


main "$@"
