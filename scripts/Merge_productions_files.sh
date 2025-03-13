#!/bin/bash

PRODUCTION_START=$1
PRODUCTION_STOP=$2

FOLDER_PRODUCTIONS="/storage/gpfs_data/neutrino/users/gi/SAND-DRIFT-STUDY/geometry/production_antinumucc_01"
echo "FOLDER_PRODUCTIONS : ${FOLDER_PRODUCTIONS}"

timestamp=$(date +"%Y%m%d_%H%M%S")

OUTPUT_DIR="/storage/gpfs_data/neutrino/users/gi/sand-physics/production_antinumucc_01/data_preunfold/"
MERGE_RECAP="${OUTPUT_DIR}merge_recap.txt"

LIST_OF_PREUNFOLD_FILES="${OUTPUT_DIR}list_of_files_preufold.${PRODUCTION_START}.to.${PRODUCTION_STOP}.txt"
LIST_OF_OUTPUT_ANALYSIS_FILES="${OUTPUT_DIR}list_of_files_output_analysis.${PRODUCTION_START}.to.${PRODUCTION_STOP}.txt"

echo -e "****** TIME STAMP : $(date +"%A, %d %B %Y, %H:%M:%S")" >> ${MERGE_RECAP}
echo -e "FOLDER_PRODUCTIONS=${FOLDER_PRODUCTIONS}" >> ${MERGE_RECAP}
echo -e "LIST_OF_PREUNFOLD_FILES=${LIST_OF_PREUNFOLD_FILES}" >> ${MERGE_RECAP}
echo -e "LIST_OF_OUTPUT_ANALYSIS_FILES=${LIST_OF_OUTPUT_ANALYSIS_FILES}" >> ${MERGE_RECAP}

total_productions=0
total_files_report_root=0
total_files_report_preunfold=0
total_processed_events=0

for (( PRODUCTION=PRODUCTION_START; PRODUCTION<=PRODUCTION_STOP; PRODUCTION++ )); do
    PRODUCTION_FORMATTED=$(printf "%04d" $PRODUCTION)
    PRODUCTION_PATH="${FOLDER_PRODUCTIONS}/production_${PRODUCTION_FORMATTED}"
    FOLDER_EVENT_PREUNFOLD="${PRODUCTION_PATH}/preunfold"
    
    count_files_preunfold=$(find "${FOLDER_EVENT_PREUNFOLD}/" -type f -name "*preunfold.root"| wc -l)

    total_productions=$((total_productions + 1))
    total_files_report_root=$((total_files_report_root + count_files_report_root))
    total_files_report_preunfold=$((total_files_report_preunfold +count_files_preunfold))

    ls ${FOLDER_EVENT_PREUNFOLD}/*preunfold.root >> ${LIST_OF_PREUNFOLD_FILES}
    ls ${FOLDER_EVENT_PREUNFOLD}/*output_analysis.root >> ${LIST_OF_OUTPUT_ANALYSIS_FILES}

    echo "production_${PRODUCTION_FORMATTED} | report txt files : ${count_files_report_root} |files preunfold ${count_files_preunfold}"
done

# success_processed=$((total_files_report_preunfold * 1000))
echo "__________________________________________________________________________________"
echo "TOTAL | productions ${total_productions} | report.root files : ${total_files_report_root} | preunfold.root files ${total_files_report_preunfold} "

# for preunfold_file in "${FOLDER_EVENT_PREUNFOLD}"/*preunfold.root; do
#     echo "Checking file: $preunfold_file"
#     if [[ -f "$preunfold_file" ]]; then
#         entries=$(get_entries "$preunfold_file")
#         total_processed_events=$((total_entries + entries))
#     else
#         echo "$preunfold_file file not found"
#     fi
# done

# for run in {0..999}; do
#     run_formatted=$((PRODUCTION * 1000 + run))
#     echo "run_formatted ${run_formatted}"
#     THIS_FOLDER_RUN="${PRODUCTION_PATH}/run_${run_formatted}"
#     # file_unfold="${FOLDER_EVENT_PREUNFOLD}"
# done