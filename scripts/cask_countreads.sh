#!/usr/bin/bash

filein="$(readlink -f ${1})"
ext="${filein##*.}"

if [ $ext == "gz" ]
then
  x=$(zcat "${filein}" | wc -l | cut -d " " -f1-1)
else
  x=$(wc -l "${filein}" | cut -d " " -f1-1)
fi
y=$((x/4))
echo "${filein} : ${y}"

echo "${y}" >"${filein}.counts.txt"

