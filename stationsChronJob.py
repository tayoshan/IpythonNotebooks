import pandas as pd
import urllib2 as url
from BeautifulSoup import BeautifulSoup
import geopandas as gp
import numpy as np
import os

#Where to store station data
print '1'
csv = "C:\\Users\\Taylor\\Desktop\\barclayscyclehireusagestats2012\\uniqueStations.csv"

#Scrape url for station xml data and process it into a pandas dataframe

stations = []
xml = "http://www.tfl.gov.uk/tfl/syndication/feeds/cycle-hire/livecyclehireupdates.xml"
site = url.urlopen(xml)
soup = BeautifulSoup(site)
for station in soup.findAll('station'):
    stations.append([int(station.id.text), float(station.long.text), float(station.lat.text)])
stats = pd.DataFrame(stations)




#Lastly add stations and associated data to station csv file
if not os.path.isfile(csv):
    stats.to_csv(csv, header=False, index=False)
    print stats


else:
    data = pd.DataFrame(pd.read_csv(csv, index_col=False, names=[0,1,2], dtype={0:int, 1:float, 2:float}))
    stats =  pd.concat([data,stats])
    stats = stats.drop_duplicates(0)
    print stats
    os.remove(csv)
    stats.to_csv(csv, header=False, index=False)

