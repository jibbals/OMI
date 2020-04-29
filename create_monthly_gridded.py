# read hdfeos5 module
import h5py
import h5netcdf # Write netcdf output


# write averaged time series to csv
import csv

# use maths and dates, and glob for files list
import numpy as np
from datetime import datetime, timedelta
from mpl_toolkits.basemap import maskoceans
from glob import glob
from compile_daily_TS import get_ocean_mask, read_omi_swath

_swathesfolder="/media/jesse/My Book/jwg366/Satellite/Aura/OMI/OMHCHOSubset/"
#_swathesfolder="data/"

__LATRES__ = 0.25
__LONRES__ = 0.3125

def GMAO_lats(dy=__LATRES__):
    '''
        GMAO structure for latitudes (half spaced top and bottom)
        returns lats_mids, lats_edges
        mids are [-90+dy/4, -90+dy/4+dy, -90+dy/4+2dy, ...]
        edges are [-90, -90+dy/2, -90+dy/2 + dy, -90+dy/2+2dy, ...]
    '''
    ny=int(180/dy) # how many latitudes
    late=np.ndarray(ny+2)
    #edges of latitude bounds
    late[0]=-90
    late[-1]=90
    # the rest is regular from half spaced polar boxes
    late[1:-1] = np.linspace(-90+dy/2.0, 90-dy/2.0, ny)

    # Latitude midpoints
    latm=(late[0:-1]+late[1:]) / 2.0
    return latm,late

def GMAO_lons(dx=__LONRES__):
    '''
        GMAO structured longitudes (this one is regularly spaced)
        returns lons_mids, lons_edges
        mids are [-180, -180+dx, ...]
    '''
    nx=int(360/dx) # how many lons
    lone=np.linspace(-180-dx/2.0, 180-dx/2.0, nx+1)
    lonm=(lone[0:-1]+lone[1:]) / 2.0

    return lonm, lone

def list_days(day0,dayn=None):
    '''
        return list of days from day0 to dayn, or just day0
        if month is True, return [day0,...,end_of_month]
    '''
    if dayn is None: return [day0,]
    day0 = datetime(day0.year,day0.month,day0.day)
    dayn = datetime(dayn.year,dayn.month,dayn.day)
    numdays = (dayn-day0).days + 1 # timedelta
    return [day0 + timedelta(days=x) for x in range(0, numdays)]


def list_months(day0,dayn):
    '''
        Return list of months (day=1) included between day0 and dayN
    '''
    # first get list of days
    days=list_days(day0,dayn)
    # Just pull out entries with day==1
    months=[d for d in days if d.day==1]
    return months

def save_to_hdf5(outfilename, arraydict, fillvalue=np.NaN,
                 attrdicts={}, fattrs={},
                 verbose=False):
    '''
        Takes a bunch of arrays, named in the arraydict parameter, and saves
        to outfilename as hdf5 using h5py (with fillvalue specified), and gzip compression

        INPUTS:
            outfilename: name of file to save
            arraydict: named arrays of data to save using given fillvalue and attributes
            attrdicts is an optional dictionary of dictionaries,
                keys should match arraydict, values should be dicts of attributes
            fattrs: extra file level attributes
    '''
    print("saving "+outfilename)
    with h5py.File(outfilename,"w") as f:
        if verbose:
            print("Inside fio.save_to_hdf5()")
            print(arraydict.keys())

        # attribute creation
        # give the HDF5 root some more attributes
        f.attrs['Filename']        = outfilename.split('/')[-1]
        f.attrs['creator']          = 'Jesse Greenslade'
        f.attrs['HDF5_Version']     = h5py.version.hdf5_version
        f.attrs['h5py_version']     = h5py.version.version
        f.attrs['Fill_Value']       = fillvalue
        # optional extra file attributes from argument
        for key in fattrs.keys():
            if verbose:
                print(key,fattrs[key], type(fattrs[key]))
            f.attrs[key] = fattrs[key]


        for name in arraydict.keys():
            # create dataset, using arraydict values
            darr=np.array(arraydict[name])
            if verbose:
                print((name, darr.shape, darr.dtype))

            # handle datetime type and boolean types
            # turn boolean into int8
            if darr.dtype == np.dtype('bool'):
                darr=np.int8(darr)
                attrdicts={}
                if not name in attrdicts:
                    attrdicts[name]={}
                attrdicts[name]['conversion']={'from':'boolean','to':'int8','by':'fio.save_to_hdf5()'}
            # harder to test for datetime type...

            # Fill array using darr,
            #
            dset=f.create_dataset(name,fillvalue=fillvalue,
                                  data=darr, compression_opts=9,
                                  chunks=True, compression="gzip")
            # for VC items and RSC, note the units in the file.
            if name in attrdicts:
                for attrk, attrv in attrdicts[name].items():
                    dset.attrs[attrk]=attrv

        # force h5py to flush buffers to disk
        f.flush()
    print("Saved "+outfilename)

def make_monthly_average_gridded(month):
    '''
    Read the average HCHO from a days worth of swathes
    molecules/cm2
    save monthly averaged HCHO, HCHO_Corrected, pixel counts with dims: [lats, lons]
    '''
    YYYYmMM=month.strftime("%Ym%m")
    outfilename="data/gridded_monthly_%s.he5"%YYYYmMM
    # files look like this:
    #  OMI-Aura_L2-OMHCHO_2009m1230t0156-o29035_v003-2014m0626t164117.SUB.he5
    pattern="*OMHCHO_%s*"%YYYYmMM
    swaths=glob(_swathesfolder+pattern)
    # if no swaths!?
    assert len(swaths) > 0, "ERROR: %s missing"%YYYYmMM
    
    print("reading %d swaths like %s"%(len(swaths),pattern))
    # for each swath, grab the entries within our lat/lon bounds
    allhcho=None
    allhcho_rsc=None
    alllats=None
    alllons=None
    
    for fpath in swaths:
        try:
            swath=read_omi_swath(fpath,mask_ocean=False, mask_land=False)
        except KeyError as KE:
            print(KE)
            print("Continuing without this swath")
            continue
        # rows x sensors
        # I x 60
        hcho=swath['HCHO']
        hcho_rsc=swath['HCHO_rsc']
        lats, lons = swath['lats'], swath['lons']
        
        # read in whole month to one array
        if allhcho is None:
            allhcho=hcho
            allhcho_rsc=hcho_rsc
            alllats=lats
            alllons=lons
        else:
            allhcho=np.append(allhcho, hcho, axis=0)
            allhcho_rsc=np.append(allhcho_rsc, hcho_rsc, axis=0)
            alllats=np.append(alllats, lats, axis=0)
            alllons=np.append(alllons, lons, axis=0)
        del swath
    
    
    ## GRID INTO LATS/LONS BINS:
    # GMAO FOR 0.25 x 0.3125:
    lons_m,lons_e = GMAO_lons()
    lats_m,lats_e = GMAO_lats()
    
    limits=[-50,100,-1,165] # lat,lon, lat,lon  : lower left cornter to upper right
    
    # how many lats, lons
    ny,nx = len(lats_m), len(lons_m)


    # Skip all this if there is no omhcho data on this day
    VC = np.ndarray([ny, nx])+np.NaN # Vertical columns 
    VCC = np.ndarray([ny, nx])+np.NaN # corrected
    counts = np.zeros([ny, nx]) # pixel counts
    
    ## TODO: FASTER WAY USING
    ## from scipy.stats import binned_statistic_2d
    ## no2_gridded = binned_statistic_2d(alllons, alllats, allno2, bins=[outlon_e,outlat_e],statistic=‘mean’).statistic
    # Just looking at subset, can skip most grid squares
    for i in range(ny):
        if (lats_m[i] < limits[0]) or (lats_m[i]>limits[2]) :continue
        
        for j in range(nx):
            if (lons_m[j]< limits[1]) or (lons_m[j]>limits[3]): continue
            
            # how many pixels within this grid box
            matches=(alllats >= lats_e[i]) & (alllats < lats_e[i+1]) & (alllons >= lons_e[j]) & (alllons < lons_e[j+1])
            
            counts[i,j]= np.sum(matches)
            
            # Save the means of each good grid pixel
            if counts[i,j] > 0:
                VC[i,j]      = np.nanmean(allhcho[matches])
                VCC[i,j]     = np.nanmean(allhcho_rsc[matches])
    
    # subset to save space
    loni = (lons_m > limits[1]) & (lons_m < limits[3])
    lati = (lats_m > limits[0]) & (lats_m < limits[2])
    counts = counts[lati,:]
    counts = counts[:,loni]
    VC     = VC[lati,:]
    VC     = VC[:,loni]
    VCC    = VCC[lati,:]
    VCC    = VCC[:,loni]
    lats   = lats_m[lati]
    lons   = lons_m[loni]
    
    
    datadict={"VC":VC, "VCC":VCC, "counts":counts, "lats":lats, "lons":lons}
    attrdicts = {"VC":{"desc":"Vertical column amount from OMHCHO good pixels","units":"molec/cm2"},
                 "VCC":{"desc":"reference sector corrected vertical columns from OMHCHO good pixels","units":"molec/cm2"},
                 "counts":{"desc":"count of good pixels averaged into latlon grid square"}}
    save_to_hdf5(outfilename, arraydict=datadict, attrdicts=attrdicts)
    
def gregorian_from_dates(dates):
    ''' gregorian array from datetime list
        gregorian is hours since 1985,1,1,0,0

    '''
    d0=datetime(1985,1,1,0,0,0)
    return np.array([(date-d0).total_seconds()/3600.0 for date in dates ])

def read_hdf5(filename):
    '''
        Should be able to read hdf5 files created by my method above...
        Returns data dictionary and attributes dictionary
    '''
    retstruct={}
    retattrs={}
    with h5py.File(filename,'r') as in_f:
        #print('reading from file '+filename)

        # READ DATA AND ATTRIBUTES:
        for key in in_f.keys():
        #    if __VERBOSE__: print(key)
            retstruct[key]=in_f[key].value
            attrs=in_f[key].attrs
            retattrs[key]={}
            # print the attributes
            for akey,val in attrs.items():
        #        if __VERBOSE__: print("reading %s(attr)   %s:%s"%(key,akey,val))
                retattrs[key][akey]=val

        # ADD FILE ATTRIBUTES TO ATTRS
        retattrs['file']={}
        for fkey in in_f.attrs.keys():
            retattrs['file'][fkey] = in_f.attrs[fkey]

    return retstruct, retattrs
    
def convert_to_netcdf():
    '''
        Read all the he5 data, output as [lon, lat, time(REC)] netcdf file.
    '''
    months = list_months(datetime(2005,1,1), datetime(2015,12,1))
    data,attrs = read_hdf5('data/gridded_monthly_2005m01.he5')
    lats=data['lats']
    lons=data['lons']
    z       = gregorian_from_dates(months)
    VC      = np.zeros([len(lons),len(lats),len(z)])
    VCC     = np.zeros([len(lons),len(lats),len(z)])
    counts  = np.zeros([len(lons),len(lats),len(z)])
    
    
    
    # Save whole gridded product
    VC[:,:,0]       =np.transpose(data['VC'])
    VCC[:,:,0]      =np.transpose(data['VCC'])
    counts[:,:,0]   =np.transpose(data['counts'])
    for i,month in enumerate(months[1:]):
        data,attrs = read_hdf5(month.strftime('data/gridded_monthly_%Ym%m.he5'))
        VC[:,:,i]=np.transpose(data['VC'])
        VCC[:,:,i]=np.transpose(data['VCC'])
        counts[:,:,i]=np.transpose(data['counts'])
        
        
    
    # WRITE TO NETCDF FILE:
    with h5netcdf.File('GriddedColumns.nc', 'w') as f:
        # set dimensions with a dictionary
        f.dimensions = {'lon':len(lons),'lat':len(lats), 'time':len(z)} # month should be unlimited, not sure how to do that
        # and update them with a dict-like interface
        # f.dimensions['x'] = 5
        # f.dimensions.update({'x': 5})
    
        var1 = f.create_variable('VC', ('lon','lat','time',), data=VC)
        var2 = f.create_variable('VCC', ('lon','lat','time',), data=VCC)
        var3 = f.create_variable('counts', ('lon','lat','time',), data=counts)
        dim1 = f.create_variable('lons', ('lon',), data=lons)
        dim2 = f.create_variable('lats', ('lat',), data=lats)
        dim3 = f.create_variable('times', ('time',), data=z)
        
        # access and modify attributes with a dict-like interface
        for k,v in attrs['VC'].items():
            var1.attrs[k] = v
        for k,v in attrs['VCC'].items():
            var2.attrs[k] = v
        for k,v in attrs['counts'].items():
            var3.attrs[k] = v
        dim3.attrs['desc'] = 'Gregorian dates, hours since 1985 01 01 00:00:00'
        
if __name__=='__main__':
    print("HELLO")
    months=list_months(datetime(2015,1,1),datetime(2015,12,1))
    #for month in months:
    #    start=datetime.now()
    #    make_monthly_average_gridded(month)
    #    check=(datetime.now()-start).total_seconds()
    #    print("~ %3.2f minutes for month"%(check/60.0))
    
    convert_to_netcdf()
    
        
        
