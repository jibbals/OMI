
# read hdfeos5 module
import h5py

# write averaged time series to csv
import csv

# regular expressions to get date from name of file
import re

# use maths and dates, and glob for files list
import numpy as np
from datetime import datetime, timedelta
from glob import glob

#_swathesfolder="/media/jesse/My Book/jwg366/OMI/OMHCHOSubset/"
_swathesfolder="data/"

def read_omi_swath(path):
    '''
    Read info from a single swath file
    NANify entries with main quality flag not equal to zero
    Returns:
        (hcho, hcho_corrected, lats, lons)
    '''
    # Total column amounts are in molecules/cm2
    datafields='/HDFEOS/SWATHS/OMI Total Column Amount HCHO/Data Fields/'
    geofields='/HDFEOS/SWATHS/OMI Total Column Amount HCHO/Geolocation Fields/'
    
    #useful fields
    field_hcho = datafields+'ColumnAmount'
    field_hcho_rsc = datafields+'ReferenceSectorCorrectedVerticalColumn'
    field_qf    = datafields+'MainDataQualityFlag'
    #field_clouds= datafields+'AMFCloudFraction'
    #field_xqf   = geofields +'XtrackQualityFlags'
    field_lon   = geofields +'Longitude'
    field_lat   = geofields +'Latitude'
    
    ## read in file:
    with h5py.File(path,'r') as in_f:
        ## get data arrays
        lats    = in_f[field_lat].value     #[ 1644, 60 ]
        lons    = in_f[field_lon].value     #
        hcho    = in_f[field_hcho].value    #
        hcho_rsc= in_f[field_hcho_rsc].value#
        qf      = in_f[field_qf].value      #
        
        ## remove missing values and bad flags: 
        # QF: missing<0, suss=1, bad=2
        suss = qf != 0
        hcho[suss]=np.NaN
        lats[suss]=np.NaN
        lons[suss]=np.NaN
    
    #return hcho, lats, lons, amf, amfg, w, apri, plevs
    return {'HCHO':hcho,'lats':lats,'lons':lons,'HCHO_rsc':hcho_rsc,'qualityflag':qf}

def read_day_avg(day, lllat,lllon,urlat,urlon):
    '''
    Read the average HCHO from a days worth of swathes
    molecules/cm2
    Returns both the ColumnAmountHCHO and the refsectorcorrected amount
    '''
    
    hchos=0
    hcount=0
    hcho_rscs=0
    hrsccount=0
    # files look like this:
    #  OMI-Aura_L2-OMHCHO_2009m1230t0156-o29035_v003-2014m0626t164117.SUB.he5
    pattern="*OMHCHO_%st*"%day.strftime("%Ym%m%d")
    swaths=glob(_swathesfolder+pattern)
    #print("reading %d swaths like %s"%(len(swaths),pattern))
    # for each swath, grab the entries within our lat/lon bounds
    for fpath in swaths:
        swath=read_omi_swath(fpath)
        # rows x sensors
        # I x 60
        hcho=swath['HCHO']
        hcho_rsc=swath['HCHO_rsc']
        lats, lons = swath['lats'], swath['lons']
        # subset to region of interest
        lonsub = (lons >= lllon) * (lons <= urlon)
        latsub = (lats >= lllat) * (lats <= urlat)
        hcho[~lonsub] = np.NaN
        hcho[~latsub] = np.NaN
        hcho_rsc[~lonsub] = np.NaN
        hcho_rsc[~latsub] = np.NaN
        hchos=hchos+np.nansum(hcho)
        hcount=hcount+np.sum(np.isnan(hcho))
        hcho_rscs=hcho_rscs+np.nansum(hcho_rsc)
        hrsccount=hrsccount+np.sum(np.isnan(hcho_rsc))
    
    out_avg=hchos/float(hcount)
    out_avg_corr=hcho_rscs/float(hrsccount)
    # use regex to grab YYYYmMMDD
    m = re.search('OMHCHO_(.+?)t', swaths[0])
    assert m, "Couldn't figure out the date from "+swaths[0]
    out_date = datetime.strptime(m.group(1), "%Ym%m%d")
    return (out_avg, out_avg_corr, out_date, hrsccount)

# Dates where we have data:
ndays=(datetime(2015,1,1)-datetime(2005,1,1)).days
dates=[datetime(2005,1,1)+timedelta(days=d) for d in range(ndays)]
#GEOS-Chem grid box: -36, 147.5, -32, 152.5
#Larger NSW region: -38, 145, -30, 153
#Sydney region: -35.5, 150, -33.5, 151.5
subsets=[ [-36, 147.5, -32, 152.5],[-38,145,-30,153],[-35.5, 150, -33.5, 151.5] ]
outnames=['TS_GC.csv','TS_LargerNSW.csv', 'TS_Sydney.csv']
# subset,outname=subsets[0],outnames[0]
for subset, outname in zip(subsets, outnames):
    hcho=[]
    hcho_rsc=[]
    times=[]
    counts=[]
    for day in dates:
        h, hc, t, c = read_day_avg(day, subset[0],subset[1],subset[2],subset[3])
        hcho.append(h)
        hcho_rsc.append(hc)
        times.append(t.strftime("%Y%m%d"))
        counts.append(c)
    print('writing %s'%outname)
    print(times[0:5],hcho[0:5])
    with open(outname, 'w') as outf:
        writer=csv.writer(outf,quoting=csv.QUOTE_NONE)
        writer.writerows(zip(times,hcho, hcho_rsc,counts))

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
            c = [ d[3] for d in data ]
        plt.plot(h,colours[i], label=outcsv)
        
    plt.xlabel('Days since '+t[0])
    plt.ylabel('molecules/cm2')
    plt.title('Average OMI HCHO (daily Gridded V3) VC subset to three regions')
    plt.ylim([1.0e14, 2.0e17])
    ax=plt.Axes
    newax=ax.twinx()
    newax.set_ylabel('satellite readings count')
    newax.plot(c,'cyan',label=outnames[3]+' sat count')
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