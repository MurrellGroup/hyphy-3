if(NOT DEFINED DIFFUBAR OR NOT DEFINED DATA_DIR OR NOT DEFINED WORK_DIR)
    message(FATAL_ERROR "DIFFUBAR, DATA_DIR, WORK_DIR are required")
endif()
file(MAKE_DIRECTORY "${WORK_DIR}")
execute_process(COMMAND "${DIFFUBAR}" --help RESULT_VARIABLE status OUTPUT_VARIABLE help)
if(NOT status EQUAL 0 OR NOT help MATCHES "difFUBAR")
    message(FATAL_ERROR "difFUBAR help failed")
endif()
execute_process(COMMAND "${DIFFUBAR}"
    --alignment "${DATA_DIR}/demo.fasta" --tree "${DATA_DIR}/demo.nwk"
    --fitted-model "${DATA_DIR}/fixed-model.json"
    --foreground-grid 2 --background-grid 1 --iterations 80 --seed 12 --threads 1
    --output "${WORK_DIR}/result.json" --grid-likelihoods "${WORK_DIR}/grid.tsv"
    --save-fitted-model "${WORK_DIR}/model.json" --quiet
    RESULT_VARIABLE status ERROR_VARIABLE error)
if(NOT status EQUAL 0)
    message(FATAL_ERROR "difFUBAR CLI failed: ${error}")
endif()
file(READ "${WORK_DIR}/result.json" result)
string(JSON name GET "${result}" analysis name)
string(JSON count LENGTH "${result}" sites)
string(JSON grid_count GET "${result}" grid category_count)
if(NOT name STREQUAL "difFUBAR" OR NOT count EQUAL 12 OR NOT grid_count EQUAL 128)
    message(FATAL_ERROR "Incorrect difFUBAR result shape")
endif()
file(STRINGS "${WORK_DIR}/grid.tsv" grid_lines)
list(LENGTH grid_lines grid_line_count)
if(NOT grid_line_count EQUAL 128)
    message(FATAL_ERROR "Incorrect category-major likelihood export")
endif()
foreach(bad_option IN ITEMS --iterations --foreground-grid --concentration)
    execute_process(COMMAND "${DIFFUBAR}" --alignment "${DATA_DIR}/demo.fasta"
        --tree "${DATA_DIR}/demo.nwk" "${bad_option}" -1
        RESULT_VARIABLE status OUTPUT_QUIET ERROR_QUIET)
    if(status EQUAL 0)
        message(FATAL_ERROR "Accepted invalid ${bad_option}")
    endif()
endforeach()
# The previous result JSON is itself a valid fixed-model input.
execute_process(COMMAND "${DIFFUBAR}"
    --alignment "${DATA_DIR}/demo.fasta" --tree "${DATA_DIR}/demo.nwk"
    --fitted-model "${WORK_DIR}/result.json" --foreground-grid 2 --background-grid 1
    --iterations 80 --seed 12 --threads 1 --output "${WORK_DIR}/replayed.json" --quiet
    RESULT_VARIABLE status ERROR_VARIABLE error)
if(NOT status EQUAL 0)
    message(FATAL_ERROR "Result-model replay failed: ${error}")
endif()
file(READ "${WORK_DIR}/replayed.json" replayed)
string(JSON first_sites GET "${result}" sites)
string(JSON second_sites GET "${replayed}" sites)
if(NOT first_sites STREQUAL second_sites)
    message(FATAL_ERROR "Seeded replay changed site posteriors")
endif()
message(STATUS "difFUBAR CLI tests passed")
