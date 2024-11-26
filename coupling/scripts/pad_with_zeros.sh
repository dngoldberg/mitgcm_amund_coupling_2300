#!/bin/sh

pad_with_zeros() {
    number="$1"
    length="$2"

    # Calculate the length of the number
    num_length=$(echo -n "$number" | wc -c)

    # Calculate the number of zeros to pad
    zeros_to_pad=$((length - num_length))

    # Construct the padded number
    padded_number=""
    for i in $(seq 1 $zeros_to_pad); do
        padded_number="${padded_number}0"
    done
    padded_number="${padded_number}${number}"

    echo "$padded_number"
}

