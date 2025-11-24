#!/bin/bash

for filename in *; do
    if [ -f "$filename" ]; then
        echo "on $filename"
    fi
done
