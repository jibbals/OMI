
# read hdfeos5 module
import h5py

# write averaged time series to csv
import csv

# regular expressions to get date from name of file
import re

# use maths and dates, and glob for files list
import numpy as np
from datetime import datetime
from glob import glob

def read_avg_hcho_subset(fpath, lllat,lllon,urlat,urlon):
    '''
    Read the average HCHO from a file
    molecules/cm2
    '''
    HCHOPath='/HDFEOS/GRIDS/OMI Total Column Amount HCHO/Data Fields/ColumnAmountHCHO'
    QFPath='/HDFEOS/GRIDS/OMI Total Column Amount HCHO/Data Fields/MainDataQualityFlag'
    lonpath='/HDFEOS/GRIDS/OMI Total Column Amount HCHO/Data Fields/Longitude'
    latpath='/HDFEOS/GRIDS/OMI Total Column Amount HCHO/Data Fields/Latitude'
    out_avg=np.NaN
    out_avg_pos=np.NaN
    
    with h5py.File(fpath,'r') as inf:
        # entries x lats x lons
        # 15 x I x J
        
        flags=inf[QFPath].value
        hcho=inf[HCHOPath].value
        hcho[~(flags == 0)]=np.NaN
        lats, lons = inf[latpath].value, inf[lonpath].value
        lonsub = (lons >= lllon) * (lons <= urlon)
        latsub = (lats >= lllat) * (lats <= urlat)
        hcho[~lonsub] = np.NaN
        hcho[~latsub] = np.NaN
        out_avg=np.nanmean(hcho)
        out_avg_pos=np.nanmean(hcho[hcho>0.0])
        #print(lllat,lllon,urlat,urlon)
        #print(lats[latsub])
        #print(lons[lonsub])
    # names like: OMI-Aura_L2G-OMHCHOG_2004m1231_v003-2014m0707t161734.SUB
    # use regex to grab YYYYmMMDD
    m = re.search('OMHCHOG_(.+?)_v003', fpath)
    assert m, "Couldn't figure out the date from "+fpath
    out_date = datetime.strptime(m.group(1), "%Ym%m%d")
    return (out_avg, out_avg_pos, out_date)


files=glob("data/AusSubset/*OMHCHOG*.he5")
files.sort()
#GEOS-Chem grid box: -36, 147.5, -32, 152.5
#Larger NSW region: -38, 145, -30, 153
#Sydney region: -35.5, 150, -33.5, 151.5
subsets=[ [-36, 147.5, -32, 152.5],[-38,145,-30,153],[-35.5, 150, -33.5, 151.5] ]
outnames=['TS_GC.csv','TS_LargerNSW.csv', 'TS_Sydney.csv']
for subset, outname in zip(subsets, outnames):
    hcho=[]
    hcho_pos=[]
    times=[]
    for ff in files:
        h, hp, t = read_avg_hcho_subset(ff, subset[0],subset[1],subset[2],subset[3])
        hcho.append(h)
        hcho_pos.append(hp)
        times.append(t.strftime("%Y%m%d"))
    print('writing %s'%outname)
    print(times[0:5],hcho[0:5])
    with open(outname, 'w') as outf:
        writer=csv.writer(outf,quoting=csv.QUOTE_NONE)
        writer.writerows(zip(times,hcho, hcho_pos))

__testing__=True

if __testing__:
    import matplotlib.pyplot as plt
    colours=['k','b','m']
    for i,outcsv in enumerate(outnames):
        with open(outcsv,'r') as inf:
            reader=csv.reader(inf)
            data=list(reader)
            t = [ d[0] for d in data ]
            h = [ d[1] for d in data ]
            hp = [ d[2] for d in data ]
        
        plt.plot(h,colours[i], label=outcsv)
        
    plt.xlabel('Days since '+t[0])
    plt.ylabel('molecules/cm2')
    plt.title('Average OMI HCHO (daily Gridded V3) VC subset to three regions')
    plt.ylim([1.0e14, 2.0e17])
    plt.legend()
    savename="TS_AllSubsets.png"
    print("saving %s"%savename)
    plt.savefig(savename)
    plt.close()
    
    # zoom in a bit on the newer stuff:
    colours=['k','b','m']
    for i,outcsv in enumerate(outnames[0:3]):
        with open(outcsv,'r') as inf:
            reader=csv.reader(inf)
            data=list(reader)
            t = [ d[0] for d in data ]
            h = [ d[1] for d in data ]
            hp = [ d[2] for d in data ]
        
        plt.plot(h[-1000:],colours[i], label=outcsv)
        
    plt.xlabel('Days since '+t[-1000])
    plt.ylabel('molecules/cm2')
    plt.title('Average OMI HCHO (daily Gridded V3) VC subset to three regions')
    plt.ylim([0, 2.0e17])
    plt.legend()
    savename="TS_Subsetslast3years.png"
    print("saving %s"%savename)
    plt.savefig(savename)
    plt.close()