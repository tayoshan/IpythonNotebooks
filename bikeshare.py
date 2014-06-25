import pandas as pd
import datetime as dt
import time
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import igraph as ig
import pySI as SI

#Turn on interactive plotting
#plt.ion()

#Get latest stations and their coordinates - aggregated over time so may change slightly
stations = pd.DataFrame(pd.read_csv("C:\\Users\\Taylor\\Desktop\\barclayscyclehireusagestats2012\\stations.csv"))

#Get pre-calculated distances for all stations pairs
dists = pd.DataFrame(pd.read_csv("C:\\Users\\Taylor\\Desktop\\barclayscyclehireusagestats2012\\stationDistances.csv"))

#Read in bike share trip data
data = pd.DataFrame(pd.read_csv("C:\\Users\\Taylor\\Desktop\\barclayscyclehireusagestats2012\\All\\allBikeTrips.csv"))


data['month'] = data['Start Date'].map(lambda x: str(x)[3:5])

#change trip start date/time to a datetime object (dayfirst = True for European date style)
data['Start Date'] = pd.to_datetime(data['Start Date'], dayfirst=True)
data

#Sort data by date/time
data = data.sort(['Start Date'])

#Set start date as index for fast slicing
data = data.set_index("Start Date")

test = pd.DataFrame(data['month'].resample('M', how='count'))
months = test.apply(lambda x: pd.to_datetime(x.name).strftime("%B"), axis=1)
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlabel('Month')
ax.set_ylabel('Trip Count')
test.plot(kind='bar',x=months, title='Number of Cycle Trips By Month', ax=ax, legend=False)
plt.show()




data['day'] = data.index.weekday
groups = data.groupby('day').count()


fig = plt.figure()
ax = fig.add_subplot(111)
groups.day.plot(kind='bar', legend=False, ax=ax, title='Cycle Trips By Day of the Week')
days = ["Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday"]
ax.set_xticklabels(days)
ax.set_xlabel('Day of the Week')
ax.set_ylabel('Trip Count')

plt.show()