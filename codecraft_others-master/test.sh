#!/bin/bash

dir1=SDK-gcc
dir2=case_example
datetime=`date +%m%d%H%M`
logfile="score"${datetime}

cd $dir1
./build.sh

cd ..
cp "$dir1/bin/cdn" "$dir2"

cd "$dir2"

g++ test.cpp -o test
date +%m%d%H%M >> "$logfile"
for file in `ls *.txt`; do
	file2=$(ls "$file" | sed -e 's/txt/out/')
    echo "$file" ":" >> "$logfile"
    ./cdn "$file" "$file2"
    ./test "$file" "$file2" >> "$logfile"
done
#cat "$logfile"
