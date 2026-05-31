#!/bin/bash
export ASAN_OPTIONS=detect_leaks=0
rm "g1v.txt"
touch "g1v.txt"
for x in {500..8500..40};
do
    y=$((x-40))
    sed -i "90s/${y}/${x}"/ input.txt
    echo "line ${x}"
    echo "=============================="
    ./main.x
    done

sed -i "90s/8500/500"/ input.txt
