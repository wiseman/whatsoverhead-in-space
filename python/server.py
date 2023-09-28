import csv
import dataclasses
import datetime
import io
import logging
import math
import os
import sys
from typing import Optional

import flask
import requests
import ephem

logger = logging.getLogger(__name__)


def get_url(url: str, path: str, max_age: datetime.timedelta):
    if (
        os.path.exists(path)
        and (
            datetime.datetime.utcnow()
            - datetime.datetime.utcfromtimestamp(os.path.getmtime(path))
        ).total_seconds()
        < max_age.total_seconds()
    ):
        with open(path, "r") as f:
            return f.read()
    else:
        sys.stderr.write(f"Downloading {url}...\n")
        response = requests.get(url)
        with open(path, "w") as f:
            f.write(response.text)
        return response.text


satcat = None


def get_satcat():
    global satcat
    if satcat is None:
        load_satcat()
    return satcat


# Load the catcat.tsv.
def load_satcat():
    global satcat
    logger.info("Loading satcat.tsv...")
    satcat_str = get_url(
        "https://planet4589.org/space/gcat/tsv/cat/satcat.tsv",
        "satcat.tsv",
        datetime.timedelta(days=1),
    )
    satcat_stream = io.StringIO(satcat_str)
    reader = csv.DictReader(satcat_stream, delimiter="\t")
    satcat = {row["Satcat"]: row for row in reader}


ORBIT_DESCS = {
    "ATM": "atmospheric orbit",
    "SO": "suborbital",
    "TA": "trans-atmospheric orbit",
    "LLEO/E": "equatorial lower LEO",
    "LLEO/I": "intermediate lower LEO",
    "LLEO/P": "polar lower LEO",
    "LLEO/S": "sun-synchronous lower LEO",
    "LLEO/R": "retrograde lower LEO",
    "LEO/E": "equatorial upper LEO",
    "LEO/I": "intermediate upper LEO",
    "LEO/P": "polar upper LEO",
    "LEO/S": "sun-synchronous upper LEO",
    "LEO/R": "retrograde upper LEO",
    "MEO": "medium Earth orbit",
    "HEO": "highly elliptical orbit",
    "HEO/M": "Molniya orbit",
    "GTO": "geotransfer orbit",
    "GEO/S": "stationary GEO",
    "GEO/I": "inclined GEO",
    "GEO/T": "synchronous GEO",
    "GEO/D": "drift GEO",
    "GEO/SI": "inclined GEO",
    "GEO/ID": "inclined drift GEO",
    "GEO/NS": "near-sync GEO",
    "VHEO": "very high Earth orbit",
    "DSO": "deep space orbit",
    "CLO": "cislunar/translunar orbit",
    "EEO": "Earth escape orbit",
    "HCO": "heliocentric orbit",
    "PCO": "planetocentric orbit",
    "SSE": "solar system escape orbit",
}


@dataclasses.dataclass
class Orbit:
    category: str

    def description(self):
        return ORBIT_DESCS.get(self.category, "unknown orbit")


def orbit_class(catno: int) -> Optional[Orbit]:
    satcat = get_satcat()
    record = satcat.get(str(catno))
    if record and record["OpOrbit"]:
        return Orbit(record["OpOrbit"])
    else:
        return None


def get_tles():
    tles_str = get_url(
        "https://www.celestrak.com/NORAD/elements/active.txt",
        "active.txt",
        datetime.timedelta(days=1),
    )
    tles = tles_str.split("\n")[:-1]
    return tles


def leo_meo_or_geo(alt_miles):
    if alt_miles < 621:
        return "LEO"
    elif alt_miles > 22000:
        return "GEO"
    else:
        return "MEO"


# Flask request handler that handles requests that look like
# "/?loc=34.0522,-118.2437".
app = flask.Flask(__name__)


def rad2deg(rad: float) -> float:
    return rad * 180.0 / math.pi


@app.route("/")
def index():
    loc = flask.request.args.get("loc")
    if loc is None:
        # Return a 400.
        return "Location not specified", 400
    lat = float(loc.split(",")[0])
    lon = float(loc.split(",")[1])
    observer = ephem.Observer()
    observer.lat = str(lat)
    observer.lon = str(lon)
    observer.elevation = 0
    observer.date = datetime.datetime.utcnow()
    tles = get_tles()
    objs = []
    for i in range(0, len(tles), 3):
        line1 = tles[i]
        line2 = tles[i + 1]
        line3 = tles[i + 2]
        name = line1.strip()
        sat = ephem.readtle(name, line2, line3)
        sat.compute(observer)
        objs.append(sat)
    # Group the objects by orbital class, and find the closest object in each
    # class.
    groups = {}
    for obj in objs:
        orbit_group = leo_meo_or_geo(obj.elevation * 3.28084 / 5280)
        if orbit_group not in groups:
            groups[orbit_group] = []
        groups[orbit_group].append(obj)
    closest_objs = []
    for orbit_group in ["LEO", "MEO", "GEO"]:
        closest_objs.append(
            min(
                groups[orbit_group],
                key=lambda obj: haversine(
                    rad2deg(obj.sublat), rad2deg(obj.sublong), lat, lon
                ),
            )
        )
    response = []
    for closest in closest_objs:
        print(closest.catalog_number)
        print(orbit_class(closest.catalog_number))
        range_miles = haversine(
            rad2deg(closest.sublat), rad2deg(closest.sublong), lat, lon
        )
        velocity_mph = closest.range_velocity * 2.23694
        elevation_miles = closest.elevation * 3.28084 / 5280
        closest_bearing = bearing(
            lat, lon, rad2deg(closest.sublat), rad2deg(closest.sublong)
        )
        cardinal = cardinal_direction(closest_bearing)
        orbit_desc = ""
        orbit = orbit_class(closest.catalog_number)
        if orbit:
            orbit_desc = f" in {orbit.description()}"
        response.append(
            f"{closest.name} ({closest.catalog_number}) is {range_miles:.0f} miles {cardinal}, {elevation_miles:.0f} miles up{orbit_desc}, and moving at {with_sigdigs(velocity_mph, 3)} mph."
        )
    return "\n".join(response)


def haversine(lat1, lon1, lat2, lon2):
    R = 3958.8  # Earth's radius in miles
    dlat = math.radians(lat2 - lat1)
    dlon = math.radians(lon2 - lon1)
    a = (
        math.sin(dlat / 2) ** 2
        + math.cos(math.radians(lat1))
        * math.cos(math.radians(lat2))
        * math.sin(dlon / 2) ** 2
    )
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    distance = R * c
    return distance


def bearing(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    dlat = math.radians(lat2 - lat1)
    dlon = math.radians(lon2 - lon1)
    y = math.sin(dlon) * math.cos(math.radians(lat2))
    x = math.cos(math.radians(lat1)) * math.sin(math.radians(lat2)) - math.sin(
        math.radians(lat1)
    ) * math.cos(math.radians(lat2)) * math.cos(dlon)
    return math.degrees(math.atan2(y, x))


def normalize_bearing(bearing: float) -> float:
    while bearing < 0.0:
        bearing += 360.0
    while bearing >= 360.0:
        bearing -= 360.0
    return bearing


def cardinal_direction(bearing: float) -> str:
    bearing = normalize_bearing(bearing)
    if bearing < 22.5:
        return "north"
    elif bearing < 67.5:
        return "northeast"
    elif bearing < 112.5:
        return "east"
    elif bearing < 157.5:
        return "southeast"
    elif bearing < 202.5:
        return "south"
    elif bearing < 247.5:
        return "southwest"
    elif bearing < 292.5:
        return "west"
    elif bearing < 337.5:
        return "northwest"
    else:
        return "north"


# Returns a number with at most the requested number of significant digits.
# (1801.4, 3) => "1800"
# (1.8, 3) => "1.8"

def with_sigdigs(x: float, n: int) -> str:
    if x == 0:
        return "0"
    order_of_magnitude = math.floor(math.log10(abs(x)))
    decimal_places = n - 1 - order_of_magnitude
    power_of_ten = 10 ** decimal_places
    rounded_number = round(x * power_of_ten) / power_of_ten
    # Format the number to a string with the required decimal places
    return format(rounded_number, f'.{max(decimal_places, 0)}f')


if __name__ == "__main__":
    app.run(host="0.0.0.0", port=8080, debug=True)
