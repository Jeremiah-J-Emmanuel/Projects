#!/bin/bash
curl https://celestrak.org/NORAD/elements/gp.php?GROUP=stations -o tle.txt

echo "This is the list of the available satellites: "
awk 'NR % 3 == 1' ./tle.txt #This print all the satellites in the tle file.

read -p "Which satellite would you like to llok at thier orbital position: " sat_name

awk -v name="$sat_name" 'index($0, name) > 0 {print; getline; print; getline; print}' tle.txt > selected_tle

if (cat selected_tle | wc -l) > 3; then
    echo "Please, enter the exact name of the satellite"
    rm -rf selected_tle
else;
    :

./fetch_tle.sh