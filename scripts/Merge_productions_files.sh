#!/bin/bash

NOF_PRODUCTIONS=59

PRODUCTION_START=$1
PRODUCTION_STOP=$2

FOLDER_PRODUCTIONS="/storage/gpfs_data/neutrino/users/gi/SAND-DRIFT-STUDY/geometry/production_antinumucc"
echo "FOLDER_PRODUCTIONS : ${FOLDER_PRODUCTIONS}"

timestamp=$(date +"%Y%m%d_%H%M%S")
LIST_OF_PREUNFOLD_FILES="/storage/gpfs_data/neutrino/users/gi/sand-physics/production_antinumucc/antinumu_CCQE_on_H_like/preunfold/list_of_files_production.${PRODUCTION_START}.to.${PRODUCTION_STOP}.txt"
LIST_OF_PREUNFOLD_SELECTED_EVENTS="/storage/gpfs_data/neutrino/users/gi/sand-physics/production_antinumucc/antinumu_CCQE_on_H_like/events_selection/list_of_files_production.${PRODUCTION_START}.to.${PRODUCTION_STOP}.txt"

total_productions=0
total_files_report_root=0
total_files_report_preunfold=0
total_processed_events=0

for (( PRODUCTION=PRODUCTION_START; PRODUCTION<=PRODUCTION_STOP; PRODUCTION++ )); do
    PRODUCTION_FORMATTED=$(printf "%04d" $PRODUCTION)
    PRODUCTION_PATH="${FOLDER_PRODUCTIONS}/production_${PRODUCTION_FORMATTED}"
    FOLDER_EVENT_SELECTION="${PRODUCTION_PATH}/event_selection"
    FOLDER_EVENT_PREUNFOLD="${PRODUCTION_PATH}/preunfold"
    
    count_files_report_root=$(find "${FOLDER_EVENT_SELECTION}/" -type f -name "*.root" | wc -l)
    count_files_preunfold=$(find "${FOLDER_EVENT_PREUNFOLD}/" -type f -name "*.root"| wc -l)

    total_productions=$((total_productions + 1))
    total_files_report_root=$((total_files_report_root + count_files_report_root))
    total_files_report_preunfold=$((total_files_report_preunfold +count_files_preunfold))

    ls ${FOLDER_EVENT_PREUNFOLD}/*root >> ${LIST_OF_PREUNFOLD_FILES}
    ls ${FOLDER_EVENT_SELECTION}/*root >> ${LIST_OF_PREUNFOLD_SELECTED_EVENTS}

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