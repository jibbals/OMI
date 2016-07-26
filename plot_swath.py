#libraries (I run in an environment set up by anaconda which includes these packages)
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import gridspec
# read hdfeos5 module
import h5py
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
    hcho[qf+xqf != 0]=np.NaN
    lons[qf+xqf != 0]=np.NaN
    lats[qf+xqf != 0]=np.NaN
    
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
    return hcho,lats,lons

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
    plt.savefig("ExampleDays.png")
    plt.close()


def examine_single_day(day=datetime(2012,1,14),cloudy=0.4):
    f=plt.figure(figsize=(13,13))
    
    gs=gridspec.GridSpec(3,2)
    ax1=f.add_subplot(gs[0,:])
    plt.sca(ax1)
    
    pattern="%s*%4dm%02d%02d*.he5"%(_swathesfolder,day.year,day.month,day.day)
    filename = glob.glob(pattern) # grab files matching pattern
    
    ymdstr="%4d%02d%02d"%(day.year,day.month,day.day)
    # basemap over area of interest
    m=Basemap(llcrnrlat=-50,  urcrnrlat=10,
              llcrnrlon=90, urcrnrlon=190,
              resolution='c',projection='merc')
    
    dayhcho=[]
    daylats=[]
    daylons=[]
    
    for fname in filename:
        hcho,lats,lons=read_swath(fname,cloudy=cloudy)
        dayhcho.append(hcho)
        daylats.append(lats)
        daylons.append(lons)
        
        # basemap grid
        xi, yi = m(lons,lats)
        
        # draw the CO total column onto the map
        cs = m.pcolormesh(xi,yi,hcho,vmin=1e14, vmax=1e17, norm=LogNorm())
    cb=m.colorbar(cs,"right",size="5%", pad="2%")
    cb.set_label('HCHO')
    m.drawcoastlines()
    
    # turn into numpy arrays and remove nans
    dayhcho=np.vstack(dayhcho)
    dayhcho=dayhcho[~np.isnan(dayhcho)]
    # look at plus and minus seperately
    dayhchoplus=dayhcho[dayhcho > 0]
    dayhchominus=-1*dayhcho[dayhcho<0]
    daylats=np.vstack(daylats)
    daylats=daylats[~np.isnan(daylats)]
    daylons=np.vstack(daylons)
    daylons=daylons[~np.isnan(daylons)]
    
    bins=np.logspace(12,19,50)
    ax2=f.add_subplot(gs[1,0])
    plt.sca(ax2)
    # reversed log bin histogram of negative hcho values
    plt.hist(dayhchominus,bins=bins)
    plt.xscale("log")
    ax2.invert_xaxis()
    plt.title("Negative entries")
    
    ax3=f.add_subplot(gs[1,1])
    plt.sca(ax3)
    plt.hist(dayhchoplus,bins=bins)
    plt.xscale("log")
    plt.title('Positive entries')
    
    ax4=f.add_subplot(gs[2,0])
    plt.sca(ax4)
    plt.hist(daylats,bins=30)
    plt.title('latitude hist')
    
    ax5=f.add_subplot(gs[2,1])
    plt.sca(ax5)
    plt.hist(daylons,bins=30)
    plt.title('longitude hist')
    plt.suptitle(ymdstr,fontsize=25)
    plt.savefig("Swath_%s"%ymdstr)
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
    plt.savefig("NegativeSwath_%s"%ymdstr)
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
    plt.savefig("SubsetVsFull_%s"%ymdstr)
    plt.close()

if __name__=="__main__":
    print("running")
    #examine_single_day(cloudy=0.1)
    #negative_swath()
    #plot_25_days(rsc=True,cloudy=0.1)
    compare_to_non_subset()