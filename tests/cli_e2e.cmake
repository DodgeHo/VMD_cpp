if(NOT DEFINED CLI OR NOT DEFINED OUTPUT_DIR)
    message(FATAL_ERROR "CLI and OUTPUT_DIR are required")
endif()

file(REMOVE_RECURSE "${OUTPUT_DIR}")
execute_process(
    COMMAND "${CLI}" --demo --sample-rate 1200 --modes 8 --alpha 50
            --output-dir "${OUTPUT_DIR}" --report
    RESULT_VARIABLE result
    OUTPUT_VARIABLE stdout
    ERROR_VARIABLE stderr
)
if(NOT result EQUAL 0)
    message(FATAL_ERROR "vmd_cli failed with ${result}\nstdout:\n${stdout}\nstderr:\n${stderr}")
endif()

foreach(name modes.csv spectra.csv omega.csv summary.json report.html)
    set(path "${OUTPUT_DIR}/${name}")
    if(NOT EXISTS "${path}")
        message(FATAL_ERROR "missing CLI output: ${path}")
    endif()
    file(SIZE "${path}" size)
    if(size EQUAL 0)
        message(FATAL_ERROR "empty CLI output: ${path}")
    endif()
endforeach()

file(READ "${OUTPUT_DIR}/modes.csv" modes LIMIT 256)
string(FIND "${modes}" "index,time_s,signal,imf_1" modes_header)
if(modes_header LESS 0)
    message(FATAL_ERROR "modes.csv header is invalid")
endif()

file(READ "${OUTPUT_DIR}/spectra.csv" spectra LIMIT 256)
string(FIND "${spectra}" "frequency_hz,imf_1_real,imf_1_imag,imf_1_magnitude" spectra_header)
if(spectra_header LESS 0)
    message(FATAL_ERROR "spectra.csv header is invalid")
endif()

file(READ "${OUTPUT_DIR}/omega.csv" omega LIMIT 256)
string(FIND "${omega}" "iteration,imf_1_hz" omega_header)
if(omega_header LESS 0)
    message(FATAL_ERROR "omega.csv header is invalid")
endif()

file(READ "${OUTPUT_DIR}/summary.json" summary)
foreach(fragment "\"schema_version\"" "\"converged\": true" "\"reconstruction_error\"" "\"report\": \"report.html\"")
    string(FIND "${summary}" "${fragment}" found)
    if(found LESS 0)
        message(FATAL_ERROR "summary.json is missing ${fragment}")
    endif()
endforeach()


set(index_input "${OUTPUT_DIR}/index-input.csv")
file(WRITE "${index_input}" "label,value\na,1\nb,2\n")
execute_process(
    COMMAND "${CLI}" --input "${index_input}" --column 1 --sample-rate 2 --modes 1
            --alpha 50 --dc --output-dir "${OUTPUT_DIR}/index-output"
    RESULT_VARIABLE index_result
    ERROR_VARIABLE index_stderr
)
if(NOT index_result EQUAL 0)
    message(FATAL_ERROR "index-selected mixed text/numeric input failed: ${index_stderr}")
endif()
file(STRINGS "${OUTPUT_DIR}/index-output/modes.csv" index_rows)
list(LENGTH index_rows index_row_count)
if(NOT index_row_count EQUAL 3)
    message(FATAL_ERROR "index-selected input lost a data row; modes.csv has ${index_row_count} rows")
endif()

set(invalid_input "${OUTPUT_DIR}/invalid.csv")
file(WRITE "${invalid_input}" "value\n1\nNaN\n")
execute_process(
    COMMAND "${CLI}" --input "${invalid_input}" --column value
            --output-dir "${OUTPUT_DIR}/invalid-output"
    RESULT_VARIABLE invalid_result
    ERROR_VARIABLE invalid_stderr
)
if(NOT invalid_result EQUAL 3)
    message(FATAL_ERROR "invalid numeric input returned ${invalid_result}, expected exit code 3")
endif()
string(FIND "${invalid_stderr}" "invalid.csv:3:1" invalid_location)
if(invalid_location LESS 0)
    message(FATAL_ERROR "invalid input error did not include line and column: ${invalid_stderr}")
endif()

execute_process(
    COMMAND "${CLI}" --demo --seed -1 --output-dir "${OUTPUT_DIR}/negative-seed"
    RESULT_VARIABLE negative_seed_result
    ERROR_VARIABLE negative_seed_stderr
)
if(NOT negative_seed_result EQUAL 2)
    message(FATAL_ERROR "negative seed returned ${negative_seed_result}, expected exit code 2")
endif()
string(FIND "${negative_seed_stderr}" "--seed expects a non-negative integer" negative_seed_message)
if(negative_seed_message LESS 0)
    message(FATAL_ERROR "negative seed error was unclear: ${negative_seed_stderr}")
endif()

execute_process(
    COMMAND "${CLI}" --demo --max-iterations 1 --output-dir "${OUTPUT_DIR}/nonconverged"
    RESULT_VARIABLE nonconverged_result
)
if(NOT nonconverged_result EQUAL 4)
    message(FATAL_ERROR "nonconverged run returned ${nonconverged_result}, expected exit code 4")
endif()
if(NOT EXISTS "${OUTPUT_DIR}/nonconverged/summary.json")
    message(FATAL_ERROR "nonconverged run did not preserve output files")
endif()

file(READ "${OUTPUT_DIR}/report.html" report)
foreach(fragment "<!doctype html>" "createElementNS" "VMD decomposition report" "XMLSerializer" "mode-toggle")
    string(FIND "${report}" "${fragment}" found)
    if(found LESS 0)
        message(FATAL_ERROR "report.html is missing ${fragment}")
    endif()
endforeach()
