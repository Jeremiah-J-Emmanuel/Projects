#!/bin/bash

clear
echo "________________Welcome to the DesJerIva Orbital Tracking System🛰️_______________"
echo "Loading..."

# Download the TLE file
if curl -s https://celestrak.org/NORAD/elements/gp.php?GROUP=stations -o tle.txt; then
    clear

    echo "================== DesJerIva Orbital Tracking System🛰️ =================="
    echo "Available Satellites:"
    echo "---------------------------------------------------------"
    awk 'NR % 3 == 1' tle.txt  # Print every 3rd line starting from line 1 (satellite names)

    echo
    read -p "Which satellite do you want to track? Enter the exact name: " sat_name

    # Extract satellite TLE lines into selected_tle
    awk -v name="$sat_name" 'index($0, name) > 0 {print; getline; print; getline; print}' tle.txt > selected_tle

    # Check if we successfully got 3 lines
    line_count=$(wc -l < selected_tle)
    if [ "$line_count" -ne 3 ]; then
        echo "❌ Satellite not found or name not entered exactly. Please try again."
        rm -f selected_tle
    else
        echo "✅ Satellite '$sat_name' selected."
        echo "---------------------------------------------------------"
        cat selected_tle
    fi

else
    echo "❌ Download failed! Please check your internet connection and try again."
fi
