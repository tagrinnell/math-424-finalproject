#!/bin/bash

for filename in *; do
    if [ -f "$filename" ]; then
        echo "on $filename"
        cat $filename | grep -v export > final_data/$filename
    fi
done