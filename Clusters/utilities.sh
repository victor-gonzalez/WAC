#!/bin/bash

# Function to extract a JSON array, multiply values by 10, convert to integers, and store in a variable
extract_json_scaled_num_array() {
    local json_file="$1"
    local key_name="$2"
    local var_name="$3"  # Name of the variable to store the result
    local output

    # Check if python3 is installed
    if ! command -v python3 >/dev/null 2>&1; then
        echo "Error: python3 is not installed." >&2
        return 1
    fi

    # Check if the JSON file exists
    if [ ! -f "$json_file" ]; then
        echo "Error: JSON file '$json_file' not found." >&2
        return 1
    fi

    # Check if var_name is provided
    if [ -z "$var_name" ]; then
        echo "Error: Variable name not provided." >&2
        return 1
    fi

    # Use Python to extract the array, multiply by 10, convert to integers, and output as space-separated list
    output=$(python3 -c "
import json
import sys
try:
    with open('$json_file', 'r') as f:
        data = json.load(f)
    if '$key_name' not in data:
        print(f\"Error: Key '$key_name' not found in $json_file.\", file=sys.stderr)
        sys.exit(1)
    if not isinstance(data['$key_name'], list):
        print(f\"Error: Key '$key_name' is not an array in $json_file.\", file=sys.stderr)
        sys.exit(1)
    values = [int(x * 10) for x in data['$key_name'] if isinstance(x, (int, float))]
    if not values:
        print(f\"Error: No valid numeric values found in '$key_name' array.\", file=sys.stderr)
        sys.exit(1)
    print(' '.join(map(str, values)))
except json.JSONDecodeError:
    print(f\"Error: Invalid JSON format in $json_file.\", file=sys.stderr)
    sys.exit(1)
" 2>&1)

    # Check if Python returned an error
    if [ $? -ne 0 ]; then
        echo "$output" >&2
        return 1
    fi

    # Store the output in the specified variable
    printf -v "$var_name" "%s" "$output"
}

extract_json_strings_array() {
    local json_file="$1"
    local key_name="$2"
    local var_name="$3"  # Name of the variable to store the result
    local output

    # Check if python3 is installed
    if ! command -v python3 >/dev/null 2>&1; then
        echo "Error: python3 is not installed." >&2
        return 1
    fi

    # Check if the JSON file exists
    if [ ! -f "$json_file" ]; then
        echo "Error: JSON file '$json_file' not found." >&2
        return 1
    fi

    # Check if var_name is provided
    if [ -z "$var_name" ]; then
        echo "Error: Variable name not provided." >&2
        return 1
    fi

    # Use Python to extract the array
    output=$(python3 -c "
import json
import sys
try:
    with open('$json_file', 'r') as f:
        data = json.load(f)
    if '$key_name' not in data:
        print(f\"Error: Key '$key_name' not found in $json_file.\", file=sys.stderr)
        sys.exit(1)
    if not isinstance(data['$key_name'], list):
        print(f\"Error: Key '$key_name' is not an array in $json_file.\", file=sys.stderr)
        sys.exit(1)
    values = data['$key_name']
    for val in values :
        val.strip()
    print(' '.join(map(str, values)))
except json.JSONDecodeError:
    print(f\"Error: Invalid JSON format in $json_file.\", file=sys.stderr)
    sys.exit(1)
" 2>&1)

    # Check if Python returned an error
    if [ $? -ne 0 ]; then
        echo "$output" >&2
        return 1
    fi

    # Store the output in the specified variable
    printf -v "$var_name" "%s" "$output"
}
