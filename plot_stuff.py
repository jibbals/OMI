#libraries (I run in an environment set up by anaconda which includes these packages)
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter # for label formatting
from matplotlib.colors import LogNorm
# read hdfeos5 module
import h5py
import csv
import numpy as np
from mpl_toolkits.basemap import Basemap
import glob
from datetime import datetime, timedelta

_swathesfolder="/media/jesse/My Book/jwg366/OMI/OMHCHOSubset/"
_fullswathexample="/home/jesse/Desktop/Repos/OMI/data/omhcho_full_20080129/"
_datafields='/HDFEOS/SWATHS/OMI Total Column Amount HCHO/Data Fields/'
_geofields='/HDFEOS/SWATHS/OMI Total Column Amount HCHO/Geolocation Fields/'
_rscfield=_datafields+'ReferenceSectorCorrectedVerticalColumn'
_vcfield =_datafields+'ColumnAmount'        

class Swath:
    def __init__(self,hcho,lats,lons,flags,xflags,clouds):
        self.hcho=hcho
        self.lats=lats
        self.lons=lons
        self.flags=flags
        self.xflags=xflags
        self.clouds=clouds

def read_swath(fname, rsc=True, cloudy=0.4, cutlons=[80,200]):
    # open file
    fh = h5py.File(fname,"r")
    
    # pull out variables of interest
    lons=fh[_geofields+'Longitude'][:]
    lats=fh[_geofields+'Latitude'][:]
    try:
        hcho = fh[[_vcfield,_rscfield][rsc]][:]
    except KeyError:
        print("Warning: file with missing RSC field, using VC field instead: %s "%fname)
        hcho = fh[_vcfield][:]
    
    # flags and cloud fraction
    qf = fh[_datafields+'MainDataQualityFlag'][:]
    try:
        xqf  = fh[_geofields +'XtrackQualityFlags'][:]
    except KeyError:
        print("Warning: file missing xtrack quality flags: %s"%fname)
        xqf = np.zeros(qf.shape)
    cld = fh[_datafields+'AMFCloudFraction'][:] > cloudy
    
    # remove non optimal flagged data
    hcho[qf != 0]=np.NaN
    hcho[xqf != 0]=np.NaN
    lons[qf != 0]=np.NaN
    lons[xqf != 0]=np.NaN
    lats[qf != 0]=np.NaN
    lats[xqf != 0]=np.NaN
    
    # Make sure lons don't jump to -180
    lons[lons<0]=lons[lons<0]+360.0
    # cut away some chaff
    lons[lons<cutlons[0]] = np.NaN
    lons[lons>cutlons[1]] = np.NaN
    lats[lons<cutlons[0]] = np.NaN
    lats[lons>cutlons[1]] = np.NaN
    hcho[lons<cutlons[0]] = np.NaN
    hcho[lons>cutlons[1]] = np.NaN
    
    # remove cloudy bits
    lons[cld]=np.NaN
    lats[cld]=np.NaN
    hcho[cld]=np.NaN
    
    # close file
    fh.close()
    
    return Swath(hcho,lats,lons,qf,xqf,cld)

class Day:
    def __init__(self, day,rsc=True,cloudy=0.4,verbose=False):
        pattern="%s*%4dm%02d%02d*.he5"%(_swathesfolder,day.year,day.month,day.day)
        filenames = glob.glob(pattern) # grab files matching pattern
        filenames.sort()
        
        hcho,lats,lons,flags,xflags,clouds=[],[],[],[],[],[]
        self.cloudy=cloudy
        self.day=day
        self.filenames=filenames
        
        for fname in filenames:
            swath=(read_swath(fname, rsc=rsc,cloudy=cloudy))
            if verbose:
                print("Reading %s"%fname)
                print("Shape: "+str(np.shape(swath.hcho)))
                print("non-nans: %d"%np.sum(np.isnan(swath.hcho)))
            hcho.append(swath.hcho)
            lats.append(swath.lats)
            lons.append(swath.lons)
            flags.append(swath.flags)
            xflags.append(swath.xflags)
            clouds.append(swath.clouds)
        self.hcho=np.vstack(hcho)
        self.lats=np.vstack(lats)
        self.lons=np.vstack(lons)
        self.flags=np.vstack(flags)
        self.xflags=np.vstack(xflags)
        self.clouds=np.vstack(clouds)

def check_flags(day=datetime(2012,1,14)):
    datac = Day(day)
    f, axes = plt.subplots(2,2)
    
    hcho=datac.hcho
    logbins=np.logspace(12,19,50)
    
    #remove nans
    filtrd=~np.isnan(datac.hcho)
    hcho=hcho[filtrd]    
    # look at plus and minus seperately
    plus=hcho[hcho > 0]
    minus=-1*hcho[hcho < 0]
    
    # First do plus and minus log bins!
    ax=axes[0,0]
    plt.sca(ax)
    plt.hist(minus,bins=logbins)
    plt.xscale("log")
    ax.invert_xaxis()
    plt.title("Negative hcho")
    plt.sca(axes[0,1])
    plt.hist(plus,bins=logbins)
    plt.xscale("log")
    plt.title('Positive hcho')
    
    # Then do flag checks
    flags=datac.flags[filtrd]
    plt.sca(axes[1,0])
    plt.scatter(flags,hcho)
    plt.title('flags')
    
    axes[1,0].xaxis.set_major_formatter(FormatStrFormatter('%3.2f'))
    xflags=datac.xflags[filtrd]
    plt.sca(axes[1,1])
    plt.scatter(xflags,hcho)
    plt.title('xflags')
    axes[1,1].xaxis.set_major_formatter(FormatStrFormatter('%3.2f'))
    ymdstr=str(day).split(' ')[0]
    plt.suptitle(ymdstr, fontsize=25)
    plt.savefig("images/DayFlags_%s"%ymdstr)
    plt.close()
    assert np.sum(np.abs(xflags))+np.sum(np.abs(flags)) == 0.0, "Non zero flag somewhere"

def plot_25_days(start=datetime(2012,1,1),rsc=True,cloudy=0.4):
    # set figure window giving x, y sizes
    fig,axes=plt.subplots(5,5,figsize=(20,20))
    colour_range=[-1e16,1e16]
    datetimes=[datetime(2012,1,1)+timedelta(days=d) for d in range(25)]
    mins,means=[],[]
    for i,day in enumerate(datetimes):
        plt.sca(axes[np.floor(i/5),np.mod(i,5)])
        pattern="%s*%4dm%02d%02d*.he5"%(_swathesfolder,day.year,day.month,day.day)
        filename = glob.glob(pattern) # grab files matching pattern
        print("looking at "+ pattern)
        ymdstr="%4d%02d%02d"%(day.year,day.month,day.day)
        # basemap over area of interest
        m=Basemap(llcrnrlat=-55,  urcrnrlat=-5,
                  llcrnrlon=110, urcrnrlon=175,
                  resolution='c',projection='merc')
        dayhcho=[]
        for fname in filename:
            hcho,lats,lons=read_swath(fname, rsc=rsc,cloudy=cloudy)
            dayhcho.append(hcho)
            # basemap gridpoints
            xi, yi = m(lons,lats)
            
            # draw the total column onto the map
            cs = m.pcolormesh(xi,yi,hcho,vmin=colour_range[0], vmax=colour_range[1])#, norm=LogNorm())
        dayhcho=np.vstack(dayhcho)
        mins.append(np.nanmin(dayhcho))
        means.append(np.nanmean(dayhcho))
        # add coastlines and equator
        m.drawcoastlines()
        
        #add title, colorbar
        if np.floor(i/5)==2 and np.mod(i,5) == 2:
            cb=m.colorbar(cs,"right",size="5%", pad="2%")
            cb.set_label('HCHO')
        
        #ax.set_title('Total Column Ascending '+str(date))
        plt.title(ymdstr)
    print("averages")
    print(means)
    print("minimums")
    print(mins)
    plt.tight_layout()
    plt.savefig("images/ExampleDays.png")
    plt.close()
    

def examine_single_day(day=datetime(2012,1,14),cloudy=0.4,rsc=True):
    
    f,axes=plt.subplots(2,2,figsize=(13,13))
    plt.sca(axes[0,0])
    data=Day(day,rsc=rsc,cloudy=cloudy,verbose=True)
    ymdstr="%4d%02d%02d"%(day.year,day.month,day.day)
    
    # basemap over area of interest
    m=Basemap(llcrnrlat=-50, urcrnrlat=10, 
              llcrnrlon=90, urcrnrlon=180,
              resolution='c',projection='merc')
    
    # look at plus and minus seperately
    hcho=data.hcho
    plus=hcho.copy()
    plus[hcho < 0]=np.NaN
    minus=hcho.copy()
    minus[hcho > 0]=np.NaN
    minus=-1*minus
    
    nlats=data.lats.copy()
    nlons=data.lons.copy()
    plats=data.lats.copy()
    plons=data.lons.copy()
    nlats[hcho>0]=np.NaN
    nlons[hcho>0]=np.NaN
    plats[hcho<0]=np.NaN
    plons[hcho<0]=np.NaN
    nxi,nyi=m(nlons,nlats) # negative lats/lons
    pxi,pyi=m(plons,plats) # positive lats/lons
    vmin,vmax=1e14,1e17
    
    # draw the negative CO total column onto the map
    cs = m.pcolormesh(nxi,nyi,minus,vmin=vmin,vmax=vmax, norm=LogNorm())
    cs.cmap.set_over('fuchsia')
    cs.set_clim(vmin, vmax)
    cb=m.colorbar(cs,"right",size="5%", pad="2%")
    cb.set_label('molecules/cm2')
    m.drawcoastlines()
    plt.title('Negative')
    
    # draw positive HCHO without a colourbar
    plt.sca(axes[0,1])
    cs = m.pcolormesh(pxi,pyi,plus,vmin=vmin,vmax=vmax, norm=LogNorm())
    m.drawcoastlines()
    plt.title('Positive')
    
    # now do plus and minus log bins!
    minus=-1*hcho[hcho < 0]
    plus=hcho[hcho>0]
    logbins=np.logspace(12,20,50)
    plt.sca(axes[1,0])
    plt.hist(minus,bins=logbins)
    plt.xscale("log")
    axes[1,0].invert_xaxis()
    plt.title("Negative hcho")
    plt.sca(axes[1,1])
    plt.hist(plus,bins=logbins)
    plt.xscale("log")
    plt.title('Positive hcho')
    plt.suptitle('OMI HCHO %s \n cloudfrac>%3.1f removed, non zero qflags and xtrackflags removed \n %s'%(['ColumnAmount','ReferenceSectorCorrected'][rsc], cloudy,ymdstr),fontsize=20)
    plt.savefig("images/Swath%s_%s"%(['','rsc'][rsc],ymdstr))
    plt.close()

def negative_swath(day=datetime(2012,1,14)):
    
    f,axes=plt.subplots(2,1,figsize=(13,13))
    
    pattern="%s*%4dm%02d%02d*.he5"%(_swathesfolder,day.year,day.month,day.day)
    filename = glob.glob(pattern) # grab files matching pattern
    
    ymdstr="%4d%02d%02d"%(day.year,day.month,day.day)
    
    for j,ax in enumerate(axes):
        # basemap over area of interest
        m=Basemap(llcrnrlat=-50,  urcrnrlat=0,
              llcrnrlon=120, urcrnrlon=190,
              resolution='c',projection='merc')
        plt.sca(ax)
        for fname in filename:
            # open file
            hcho,lats,lons=read_swath(fname)
            
            cbarstr="molecs/cm2"
            if j==1:
                lats[hcho>0]=np.NaN
                lons[hcho>0]=np.NaN
                hcho[hcho>0]=np.NaN
                hcho[hcho<0]=-hcho[hcho<0]
                cbarstr="-ve molecs/cm2"
            
            
            # lat lon are 1D, basemap uses 2D mesh
            #lon,lat = np.meshgrid(lons,lats)
            xi, yi = m(lons,lats)
            
            # draw the CO total column onto the map
            cs = m.pcolormesh(xi,yi,hcho,vmin=1e14, vmax=1e19, norm=LogNorm())
        cb=m.colorbar(cs,"right",size="5%", pad="2%")
        cb.set_label(cbarstr)
        m.drawcoastlines()
        plt.title(["HCHO","-ve HCHO"][j])
    
    
    plt.suptitle(ymdstr,fontsize=25)
    plt.savefig("images/NegativeSwath_%s"%ymdstr)
    plt.close()

def compare_to_non_subset(rsc=True, cloudy=0.4):
    f,axes=plt.subplots(2,1,figsize=(13,13))
    day=datetime(2008,1,29)
    vlims=[-4e16,5e16]
    # read subset data
    patternsub="%s*%4dm%02d%02d*.he5"%(_swathesfolder,day.year,day.month,day.day)
    filenamesub = glob.glob(patternsub) # grab files matching pattern
    
    patternfull="%s*%4dm%02d%02d*.he5"%(_fullswathexample,day.year,day.month,day.day)
    filenamefull = glob.glob(patternfull)
    ymdstr="%4d%02d%02d"%(day.year,day.month,day.day)
    
    for j,ax in enumerate(axes):
        # basemap over area of interest
        m=Basemap(llcrnrlat=-60,  urcrnrlat=5,
              llcrnrlon=90, urcrnrlon=190,
              resolution='c',projection='merc')
        plt.sca(ax)
        for fname in [filenamesub,filenamefull][j]:
            # open file
            hcho,lats,lons=read_swath(fname)
            
            xi, yi = m(lons,lats)
            cs = m.pcolormesh(xi,yi,hcho,vmin=vlims[0], vmax=vlims[1])
        cb=m.colorbar(cs,"right",size="5%", pad="2%")
        cb.set_label('molecs/cm2')
        m.drawcoastlines()
        plt.title(["Subset","Full"][j])
    
    
    plt.suptitle(ymdstr,fontsize=25)
    plt.savefig("images/SubsetVsFull_%s"%ymdstr)
    plt.close()

def plot_time_series(corrected=False):
    outnames=['TS_Aus.csv','TS_Sydney.csv']
    colours=['k','m']
    f = plt.figure(figsize=(16,14))
    for i,outcsv in enumerate(outnames):
        with open(outcsv,'r') as inf:
            reader=csv.reader(inf)
            data=list(reader)
            t = [ d[0] for d in data ] # dates
            h = [ float(d[1]) for d in data ] # old averages
            hcor = [ float(d[2]) for d in data ] # corrected averages
            c = [ int(float(d[3])) for d in data ] # how many entries averaged
            if corrected:
                plt.plot(hcor,'.'+colours[i], label=outcsv)
            else:
                plt.plot(h,'.'+colours[i], label=outcsv)
                print("Minimum entry:%4.2e"%np.nanmin(h))
        
    plt.xlabel('Days since '+t[0])
    plt.ylabel('molecules/cm2')
    f.suptitle('Average OMI HCHO (daily Gridded V3) VC subset to three regions')
    plt.ylim([-1e16, 6e16])
    plt.legend(loc=2)
    
    # plot average counts for every 30 days
    newax=plt.twinx(plt.gca())
    newax.set_ylabel('good entries')
    avgs=[]
    c=np.array(c).astype(float)
    for i in np.arange(0,len(t),30):
        mean30=np.mean(c[i:(i+30)])
        avgs.append(mean30)
    newax.plot(np.arange(0,len(t),30), np.array(avgs),'cyan',label='Sydney good entries(30 day mean)')
    newax.set_ylim([0,60])
    newax.legend(loc=1)
    savename="images/TS_AllSubsets.png"
    print("saving %s"%savename)
    plt.savefig(savename)
    plt.close()

if __name__=="__main__":
    print("running")
    #examine_single_day(cloudy=0.4, rsc=True)
    #examine_single_day(cloudy=0.4, rsc=False)
    #examine_single_day(day=datetime(2008,1,8),rsc=False,cloudy=0.4)
    #examine_single_day(day=datetime(2009,6,19),rsc=False,cloudy=0.3)
    #negative_swath()
    #plot_25_days(rsc=True,cloudy=0.1)
    #compare_to_non_subset()
    plot_time_series()
    #check_flags()
