#!/usr/bin/bash

FF="$(readlink -f ${1})"
ext="${FF##*.}"

if [ $ext == "gz" ]
then
  x=$(zcat "${FF}" | wc -l | cut -d " " -f1-1)
else
  x=$(wc -l "${FF}" | cut -d " " -f1-1)
fi
y=$((x/4))
echo "${FF} : ${y}"

echo "${y}" >"${FF}.counts.txt"

