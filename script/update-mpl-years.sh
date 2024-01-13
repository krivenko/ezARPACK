#!/bin/env bash

#
# Update copyright years in files
#

FIND=$(which find)
SED=$(which sed)

# Collect paths of files to be updated
FILES="README.md doc/Doxyfile.in doc/conf.py.in"
for PATTERN in "*.cpp" "*.hpp" "CMakeLists.txt" "*.cmake" "*.cmake.in"
do
    FILES="$FILES $($FIND . -name $PATTERN)"
done

# Write the new year string
YEARS_REGEX="2016-[0-9][0-9][0-9][0-9]"
YEARS_UPDATED="2016-$(date +%Y)"

COPYRIGHT_REGEX="Copyright (C) ${YEARS_REGEX} Igor Krivenko"
COPYRIGHT_UPDATED="Copyright (C) ${YEARS_UPDATED} Igor Krivenko"
for FILE in $FILES
do
    $SED -i "s/${COPYRIGHT_REGEX}/${COPYRIGHT_UPDATED}/g" $FILE
done

# doc/conf.py.in requires special treatment
sed -i "s/${YEARS_REGEX}, Igor Krivenko/${YEARS_UPDATED}, Igor Krivenko/g" \
       "doc/conf.py.in"
