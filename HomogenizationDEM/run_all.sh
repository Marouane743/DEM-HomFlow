#!/usr/bin/env bash

# Ensure the script exits immediately if any command fails
set -e

# Function to calculate elapsed time
calculate_elapsed_time() {
    local start=$1
    local end=$2
    # Force the C locale so that printf accepts the dot as the decimal separator
    LC_NUMERIC=C printf "%.3f" "$(echo "$end - $start" | bc)"
}


# 1. Record the overall start time (optional, if you want total script time)
# overall_start_time=$(date +%s.%N)

# 2) User supplies a new RC:
RC_VALUE="$1"
if [ -z "$RC_VALUE" ]; then
   echo "Usage: $0 <RC_VALUE>"
   exit 1
fi

# 3) Compare with last_rc.txt
LAST_RC_FILE="last_rc.txt"
LAST_RC="none"
if [ -f "$LAST_RC_FILE" ]; then
   LAST_RC=$(cat "$LAST_RC_FILE")
fi

# Initialize elapsed time variables
elapsed_computing_X_v3="Not Executed"
elapsed_homogenization_nv1="Not Executed"

# 4) If different => sed in shared_config.h, remove integrals_csv/*.csv, re-run computing_X
if [ "$RC_VALUE" != "$LAST_RC" ]; then
   echo "RC changed from $LAST_RC to $RC_VALUE. Re-generating integrals."

   # Update shared_config.h with the new RC value
   sed -i "s/^#define RC .*/#define RC $RC_VALUE/" shared_config.h

   # Remove old CSV files
   rm -f integrals_csv/*.csv

   # Compile computing_X_v3.c
   echo "Compiling computing_X_v3.c..."
   gcc computing_X_v3.c -o computing_X_v3 -fopenmp -lm
   echo "Compilation of computing_X_v3.c succeeded."

   # Record start time for computing_X_v3
   start_computing_X_v3=$(date +%s.%N)

   # Run computing_X_v3 with the new RC value
   echo "Running computing_X_v3 with RC=$RC_VALUE..."
   ./computing_X_v3 "$RC_VALUE"
   echo "Execution of computing_X_v3 completed."

   # Record end time for computing_X_v3
   end_computing_X_v3=$(date +%s.%N)

   # Calculate elapsed time for computing_X_v3
   elapsed_computing_X_v3=$(calculate_elapsed_time "$start_computing_X_v3" "$end_computing_X_v3")

   # Update last_rc.txt with the new RC value
   echo "$RC_VALUE" > "$LAST_RC_FILE"
else
   echo "RC has not changed ($RC_VALUE). Skipping integrals recompute."
fi

# 5) Compile and run homogenization_nv1.c
echo "Compiling homogenization_nv1.c..."
gcc homogenization_nv1.c -o homogenization_nv1 -lm
echo "Compilation of homogenization_nv1.c succeeded."

echo "Running homogenization_nv1 with RC=$RC_VALUE..."

# Record start time for homogenization_nv1
start_homogenization_nv1=$(date +%s.%N)

# Run homogenization_nv1
./homogenization_nv1 "$RC_VALUE"
echo "Execution of homogenization_nv1 completed."

# Record end time for homogenization_nv1
end_homogenization_nv1=$(date +%s.%N)

# Calculate elapsed time for homogenization_nv1
elapsed_homogenization_nv1=$(calculate_elapsed_time "$start_homogenization_nv1" "$end_homogenization_nv1")

# 6. Display the elapsed times
echo "----------------------------------------"
echo "Execution Time Summary:"
if [ "$elapsed_computing_X_v3" != "Not Executed" ]; then
    echo "computing_X_v3: $elapsed_computing_X_v3 seconds"
else
    echo "computing_X_v3: Not Executed (RC unchanged)"
fi
echo "homogenization_nv1: $elapsed_homogenization_nv1 seconds"
# overall_end_time=$(date +%s.%N)
# elapsed_overall=$(calculate_elapsed_time "$overall_start_time" "$overall_end_time")
# echo "Total script execution time: $elapsed_overall seconds"
echo "----------------------------------------"
