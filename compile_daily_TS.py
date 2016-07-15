

# read hdfeos5 module
import h5py

# write averaged time series to csv
import csv

# use maths and dates, and glob for files list
import numpy as np
from datetime import datetime, timedelta
from glob import glob

_swathesfolder="/media/jesse/My Book/jwg366/OMI/OMHCHOSubset/"
#_swathesfolder="data/"

def read_omi_swath(path,removerowanomaly=True):
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
    field_xqf   = geofields +'XtrackQualityFlags'
    field_lon   = geofields +'Longitude'
    field_lat   = geofields +'Latitude'
    
    ## read in file:
    with h5py.File(path,'r') as in_f:
        ## get data arrays
        lats    = in_f[field_lat].value      #[ 1644, 60 ]
        lons    = in_f[field_lon].value      #
        hcho    = in_f[field_hcho].value     #
        hcho_rsc= in_f[field_hcho_rsc].value #
        qf      = in_f[field_qf].value       #
        xqf     = in_f[field_xqf].value      #
        ## remove missing values and bad flags: 
        # QF: missing<0, suss=1, bad=2
        suss = qf != 0
        # XQF:
        if removerowanomaly: suss = suss * (xqf != 0)
        
        hcho[suss]=np.NaN
        lats[suss]=np.NaN
        lons[suss]=np.NaN
    
    #return hcho, lats, lons, amf, amfg, w, apri, plevs
    return {'HCHO':hcho,'lats':lats,'lons':lons,'HCHO_rsc':hcho_rsc,'qualityflag':qf,'xqf':xqf}

# checking xqf flag
fpath=glob('data/*2009m1230t0156*')[0]
swath=read_omi_swath(fpath)

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
        nans=np.repeat(np.NaN,3)
        return(nans,nans,np.repeat(0,3))
    
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
#enddate=datetime(2010,1,1)
enddate=datetime(2016,1,1)
startdate=datetime(2005,1,1)
ndays=(enddate-startdate).days
dates=[startdate+timedelta(days=d) for d in range(ndays)]
#GEOS-Chem grid box: -36, 147.5, -32, 152.5
#Larger NSW region: -38, 145, -30, 153
#Sydney region: -35.5, 150, -33.5, 151.5
subsets=[ [-36, 147.5, -32, 152.5],[-38,145,-30,153],[-35.5, 150, -33.5, 151.5] ]
outnames=['TS_GC.csv','TS_LargerNSW.csv', 'TS_Sydney.csv']
# subset,outname=subsets[0],outnames[0]

# list of lists
hcho=[[] for i in range(len(subsets))]
hcho_rsc=[[] for i in range(len(subsets))]
times=[[] for i in range(len(subsets))]
counts=[[] for i in range(len(subsets))]

# for every day, read the averages and count the entries
st=datetime.now()
for day in dates:
    try:
        h, hc, c = read_day_avg(day, subsets)
    except Exception as e:
        print("WARNING: day %s file is bad?"%day.strftime("%Y%m%d"))
        print("WARNING: Skipping this day, printing error message:")
        print(e.message)
        h, hc, c = (np.repeat(np.NaN,3), np.repeat(np.NaN,3), np.repeat(0,3))
    if day == startdate+timedelta(days=100):
        check=(datetime.now()-st).total_seconds()
        print("~ %3.2f seconds per 100 days"%check)
        print("~ %4.2f minutes left..."%(check/100./60*ndays))
        
    ymd=day.strftime("%Y%m%d")
    # store them in lists for each subset
    for i in range(len(subsets)):
        hcho[i].append(h[i])
        hcho_rsc[i].append(hc[i])
        times[i].append(ymd)
        counts[i].append(c[i])

for i, outname in enumerate(outnames):
    print('writing %s'%outname)
    print(times[i][0:3],hcho[i][0:3])
    with open(outname, 'w') as outf:
        writer=csv.writer(outf,quoting=csv.QUOTE_NONE)
        writer.writerows(zip(times[i],hcho[i], hcho_rsc[i],counts[i]))

## Testing stuff
if True:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    colours=['k','b','m']
    f,ax = plt.subplots(2,1,sharex=True,figsize=(16,14))
    for i,outcsv in enumerate(outnames):
        with open(outcsv,'r') as inf:
            reader=csv.reader(inf)
            data=list(reader)
            t = [ d[0] for d in data ] # dates
            h = [ d[1] for d in data ] # old averages
            hcor = [ d[2] for d in data ] # corrected averages
            c = [ d[3] for d in data ] # entries averaged
        ax[0].plot(hcor,'.'+colours[i], label=outcsv)
        ax[1].plot(h,'.'+colours[i], label=outcsv)
        
    ax[1].set_xlabel('Days since '+t[0])
    ax[0].set_ylabel('molecules/cm2')
    ax[1].set_title('non corrected')
    f.suptitle('Average OMI HCHO (daily Gridded V3) VC subset to three regions')
    [ ax[i].set_ylim([-3e16, 5e16]) for i in [0,1] ]
    ax[0].legend(loc=0)
    
    # plot average count for every 30 days
    avgs=[]
    c=np.array(c).astype(float)
    for i in np.arange(0,len(t),30):
        mean30=np.mean(c[i:(i+30)])
        avgs.append(mean30)
    newax=plt.twinx(ax[0])
    newax.set_ylabel('good entries')
    newax.plot(np.arange(0,len(t),30), np.array(avgs),'cyan',label='Sydney good entries(30 day mean)')
    newax.set_ylim([20,70])
    newax.legend(loc=0)
    savename="TS_AllSubsets.png"
    print("saving %s"%savename)
    plt.savefig(savename)
    plt.close()
    
    # zoom in a bit and compare:
    plt.cla()
    colours=['k','b','m']
    for i,outcsv in enumerate(outnames[0:3]):
        with open(outcsv,'r') as inf:
            reader=csv.reader(inf)
            data=list(reader)
            t = [ d[0] for d in data ]
            h = [ d[1] for d in data ]
            hcor = [ d[2] for d in data ]
        plt.plot(hcor[-700:],'.'+colours[i], label=outcsv)
        
    plt.xlabel('Days since '+t[-700])
    plt.ylabel('molecules/cm2')
    plt.title('Average OMI HCHO VC from swaths, subset in three regions')
    plt.ylim([-5e16, 7.0e16])
    plt.legend()
    savename="TS_Subsetslast2years.png"
    print("saving %s"%savename)
    plt.savefig(savename)
    plt.close()
