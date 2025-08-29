#!/usr/bin/python3
import requests
from predictor import predict_passes
from sgp4.api import Satrec
from sgp4.api import jday
#from date import datetime, timedelta    
#from time import datetime

#You will need to add at the beginning of this file, choose from favourites or not.

"""
This is a sample ouput of the json file from the ipinfo.io/json.com
{
  "ip": "102.22.168.128",
  "city": "Kigali",
  "region": "Kigali",
  "country": "RW",
  "loc": "-1.9500,30.0588",
  "org": "AS36924 GVA Cote d'Ivoire SAS",
  "timezone": "Africa/Kigali",
  "readme": "https://ipinfo.io/missingauth"
}
"""


results = requests.get("https://ipinfo.io/json")
data = results.json()
lat, long = data["loc"].split(",") #Stored as a tuple
(city, region, country) = (data["city"], data["region"], data["country"])


print("These are your location and coordinates.")
#Will add option for turning off autolocater and allow manual entry.
print(f"{city}, {region}, {country}")
print("Latitude:", lat)
print("Longitude:", long)

ans = input("Do you want to save this location as a favourite location? [Y/N]").strip().lower()
try:
    while True:
        if ans == "y":
            #Check how to save somthing from a python to a txt file.
            print("saved succesfully")
            break
        elif ans == "n":
            print("Location not added to favourites")
            break
        else:
            print("Please select Y or N")
            continue
            #You can also add a feauture that allows users to save a location as a favourite location
except KeyboardInterrupt:
    print("Location Bookmark closed")

with open("selected_tle.txt") as file:
    lines = file.readlines()
    name, line1, line2 = lines[0], lines[1], lines[2]

predict_passes(line1, line2, lat, long)