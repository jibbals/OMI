# read hdfeos5 module
import h5py

# write averaged time series to csv
import csv

# use maths and dates, and glob for files list
import numpy as np
from datetime import datetime, timedelta
from glob import glob

_swathesfolder="/media/jesse/My Book/jwg366/Satellite/Aura/OMI/OMHCHOSubset"
#_swathesfolder="data/"

def read_omi_swath(path,removerowanomaly=True, cloudy=0.4, screen=[-0.5e16, 1e17], szamax=60):
    '''
    Read info from a single swath file
    NANify entries with main quality flag not equal to zero
    Filtering: remove pixels with following properties
        cloud frac > cloudy
        Col Density outside screen range
        solar zenith angle > szamax
        quality flag not 0
        xtrack flag not 0
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
    field_clouds= datafields+'AMFCloudFraction'
    field_xqf   = geofields +'XtrackQualityFlags'
    field_lon   = geofields +'Longitude'
    field_lat   = geofields +'Latitude'
    field_sza   = geofields +'SolarZenithAngle'
    
    ## read in file:
    with h5py.File(path,'r') as in_f:
        ## get data arrays
        lats    = in_f[field_lat].value      #[ 1644, 60 ]
        lons    = in_f[field_lon].value      #
        hcho    = in_f[field_hcho].value     #
        hcho_rsc= in_f[field_hcho_rsc].value #
        qf      = in_f[field_qf].value       #
        xqf     = in_f[field_xqf].value      #
        cld     = in_f[field_clouds].value   #
        sza     = in_f[field_sza].value      #
        
        ## remove missing values and bad flags: 
        # QF: missing<0, suss=1, bad=2
        suss = (qf != 0)
        
        hcho[suss]=np.NaN
        lats[suss]=np.NaN
        lons[suss]=np.NaN
        hcho_rsc[suss]=np.NaN
        
        # XQF:
        if removerowanomaly: 
            xsuss=(xqf != 0)
            hcho[xsuss]=np.NaN
            lats[xsuss]=np.NaN
            lons[xsuss]=np.NaN
            hcho_rsc[xsuss]=np.NaN
        
        # remove cloudiness
        rmcloud=cld > cloudy
        hcho[rmcloud]=np.NaN
        lats[rmcloud]=np.NaN
        lons[rmcloud]=np.NaN
        hcho_rsc[rmcloud]=np.NaN
        
        # remove range outside of screen values
        if screen is not None:
            rm = (hcho < screen[0]) + (hcho > screen[1])
            rmrsc= (hcho_rsc < screen[0]) + (hcho_rsc > screen[1])
            hcho[rm]=np.NaN
            lats[rm]=np.NaN
            lons[rm]=np.NaN
            hcho_rsc[rmrsc]=np.NaN
        
        if szamax is not None:
            rm = (sza > szamax)
            hcho[rm]=np.NaN
            lats[rm]=np.NaN
            lons[rm]=np.NaN
            hcho_rsc[rm]=np.NaN
        
    #return hcho, lats, lons, amf, amfg, w, apri, plevs
    return {'HCHO':hcho,'lats':lats,'lons':lons,'HCHO_rsc':hcho_rsc,'qualityflag':qf,'xqf':xqf,'sza':sza}

def read_day_avg(day, subsets):
    '''
    Read the average HCHO from a days worth of swathes
    molecules/cm2
    Returns both the ColumnAmountHCHO and the refsectorcorrected amount
    Pass in a list of subsets to pull the daily average out for each subset
    '''
    YYYYmMMDD=day.strftime("%Ym%m%d")
    nsubs=len(subsets)
    hchos=np.zeros(nsubs)
    hcount=np.zeros(nsubs)
    hcho_rscs=np.zeros(nsubs)
    hrsccount=np.zeros(nsubs)
    
    # files look like this:
    #  OMI-Aura_L2-OMHCHO_2009m1230t0156-o29035_v003-2014m0626t164117.SUB.he5
    pattern="*OMHCHO_%st*"%YYYYmMMDD
    swaths=glob(_swathesfolder+pattern)
    # if no swaths!?
    if len(swaths) == 0:
        print("Warning: %s missing"%YYYYmMMDD)
        nans=np.repeat(np.NaN,nsubs)
        return(nans,nans,np.repeat(0,nsubs))
    
    #print("reading %d swaths like %s"%(len(swaths),pattern))
    # for each swath, grab the entries within our lat/lon bounds
    for fpath in swaths:
        swath=read_omi_swath(fpath)
        # rows x sensors
        # I x 60
        hcho=swath['HCHO']
        hcho_rsc=swath['HCHO_rsc']
        lats, lons = swath['lats'], swath['lons']
        # subset to region of interest, ignoring warnings
        for si, subset in enumerate(subsets):
            lllat,lllon,urlat,urlon = subset
            with np.errstate(invalid='ignore'):
                lonsub = (lons >= lllon) * (lons <= urlon)
                latsub = (lats >= lllat) * (lats <= urlat)
                sub=lonsub*latsub
            hchos[si]=hchos[si]+np.nansum(hcho[sub])
            hcount[si]=hcount[si]+np.sum(sub)
            hcho_rscs[si]=hcho_rscs[si]+np.nansum(hcho_rsc[sub])
            hrsccount[si]=hrsccount[si]+np.sum(sub)
    
    out_avg=np.zeros(nsubs)
    out_avg_corr=np.zeros(nsubs)
    out_counts=hrsccount
    for si in range(nsubs):
        if hcount[si]==0:
            assert hrsccount[si]==0, "Both counters should be zero if one is..."
            out_avg[si]=np.NaN
            out_avg_corr[si]=np.NaN
        else:
            out_avg[si]=hchos[si]/float(hcount[si])
            out_avg_corr[si]=hcho_rscs[si]/float(hrsccount[si])
    return (out_avg, out_avg_corr, out_counts)

# Dates where we have data:
enddate=datetime(2016,1,1)
startdate=datetime(2005,1,1)
ndays=(enddate-startdate).days
dates=[startdate+timedelta(days=d) for d in range(ndays)]
#GEOS-Chem grid box: -36, 147.5, -32, 152.5
#Larger NSW region: -38, 145, -30, 153
#Sydney region: -35.5, 150, -33.5, 151.5
# Massive comparison region: -50,110,-10,160

#subsets=[ [-36, 147.5, -32, 152.5],[-38,145,-30,153],[-35.5, 150, -33.5, 151.5] , [-50,110,-10, 160]]
#outnames=['TS_GC.csv','TS_LargerNSW.csv', 'TS_Sydney.csv','TS_Aus.csv']
# subset,outname=subsets[0],outnames[0]

# second set of subsets for Kaitlyn Jan2017
# 1. A box like the South Coast region, but with Sydney excluded (what’s the resolution of the raw product? I can give you lat/lon based on that resolution but it’ll be around -34,150,-33.5,151)
#2. The South Coast region but only ocean cells (Jenny said that this should be easy for you to filter?)
#3. The South Coast region but with only land cells
#4. Like #1 but with only ocean cells
#5. Like #1 but with only land cells


# list of lists
n_subs=len(subsets)
hcho=[[] for i in range(n_subs)]
hcho_rsc=[[] for i in range(n_subs)]
times=[[] for i in range(n_subs)]
counts=[[] for i in range(n_subs)]

# for every day, read the averages and count the entries
st=datetime.now()
for day in dates:
    try:
        h, hc, c = read_day_avg(day, subsets)
    except Exception as e:
        print("WARNING: day %s file is bad?"%day.strftime("%Y%m%d"))
        print("WARNING: Skipping this day, printing error message:")
        print(e.message)
        h, hc, c = (np.repeat(np.NaN,n_subs), np.repeat(np.NaN,n_subs), np.repeat(0,n_subs))
    if day == startdate+timedelta(days=100):
        check=(datetime.now()-st).total_seconds()
        print("~ %3.2f seconds per 100 days"%check)
        print("~ %4.2f minutes left..."%(check/100./60*ndays))
        
    ymd=day.strftime("%Y%m%d")
    # store them in lists for each subset
    for i in range(n_subs):
        hcho[i].append(h[i])
        hcho_rsc[i].append(hc[i])
        times[i].append(ymd)
        counts[i].append(c[i])

for i, outname in enumerate(outnames):
    print('writing %s'%outname)
    print(times[i][0:n_subs],hcho[i][0:n_subs])
    with open(outname, 'w') as outf:
        writer=csv.writer(outf,quoting=csv.QUOTE_NONE)
        writer.writerows(zip(times[i],hcho[i], hcho_rsc[i],counts[i]))


