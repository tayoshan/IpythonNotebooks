import pandas as pd
import urllib2 as url
from BeautifulSoup import BeautifulSoup
import geopandas as gp
import numpy as np
from math import radians, cos, sin, asin, sqrt

def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a))

    # 6367 km is the radius of the Earth
    km = 6367 * c
    return km


csvIn = "C:\\Users\\Taylor\\Desktop\\barclayscyclehireusagestats2012\\uniqueStations.csv"

csvOut = "C:\\Users\\Taylor\\Desktop\\barclayscyclehireusagestats2012\\uniqueStationsDistances.csv"

stations = pd.DataFrame(pd.read_csv(csvIn, index_col=False, names=[0,1,2], dtype={0:int, 1:float, 2:float}))



combos = []
for station in stations[0].values:
    for pair in stations[0].values:
        #print station, pair
        lon1 = stations[1][stations[0]== station].values
        lat1 = stations[2][stations[0]== station].values
        lon2 = stations[1][stations[0]== pair].values
        lat2 = stations[2][stations[0]== pair].values
        distance = haversine(lon1, lat1, lon2, lat2)
        #print distance
        combos.append((station,pair,distance))


pd.DataFrame(combos).to_csv(csvOut, header=False, index=False)