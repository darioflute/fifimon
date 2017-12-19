# Several modules used for reading/reducing FIFI-LS data
import multiprocessing as mp
from lmfit.models import LinearModel
import numpy as np


def exploreDirectory(path):
    ''' Explore FITS files in a directory '''
    # Should be smarter and open only files not yet considered
    
    from astropy.io import fits
    import os, glob
    import numpy as np
    from datetime import datetime
    
    files = sorted(glob.glob(path+'*w.fits'), key=os.path.getctime)
    files = np.array(files)

    # Wavelenghts
    start=[]
    fgid=[]
    ch=[]
    
    for infile in files:
        hlf = fits.open(infile)
        try:
            h = hlf[0].header
            proc = h['PROCSTAT']
            #obstype = h['OBSTYPE']
            if proc == 'LEVEL_1':
                #if obstype == 'OBJECT':
                start.append(h['DATE-OBS'])
                fgid.append(h['FILEGPID'])
                ch.append(h['DETCHAN'])
        except:
            pass
        hlf.close()

    start=np.array(start)
    fgid =np.array(fgid)
    ch   =np.array(ch)

    # Ideally I should order by date and time, now only by time
    a = [datetime.strptime(s, '%Y-%m-%dT%H:%M:%S').time() for s in start]
    a = np.array(a)
    s = np.argsort(a)

    return files[s], start[s], fgid[s], ch[s]



def readData(fitsfile):
    from astropy.io import fits
    import numpy as np

    hdulist = fits.open(fitsfile)
    scidata = hdulist[1].data
    header = hdulist[0].header
    data = scidata.DATA
    hdulist.close()
    procstat = header['PROCSTAT']
    obstype = header['OBSTYPE']
    if obstype != 'OBJECT':
        print (" This program works only with type OBJECT")
        print (fitsfile, "is a ", obstype, " file")
        exit
    elif procstat != 'LEVEL_1':
        print ("This program works only with raw FIFI-LS files (Level 1)")
        exit
    else:
        detchan = header['DETCHAN']
        obsdate = header['DATE-OBS']
        dichroic = header['DICHROIC']
        if detchan == 'RED':
            #wvc = header['G_WAVE_R']
            ncycles = header['C_CYC_R']
            start = header['G_STRT_R']
            step = header['G_SZUP_R']
            ngrat = header['G_PSUP_R']
            order = 1
        else:
            #wvc = header['G_WAVE_B']
            ncycles = header['C_CYC_B']
            start = header['G_STRT_B']
            step = header['G_SZUP_B']
            ngrat = header['G_PSUP_B']
            order = header['G_ORD_B']

        filegpid=header['FILEGPID']    
        nodbeam = header['NODBEAM']
        # Position
        xmap = float(header['DLAM_MAP'])
        ymap = float(header['DBET_MAP'])
        xoff = float(header['DLAM_OFF'])
        yoff = float(header['DBET_OFF'])
        ra   = float(header['OBSLAM'])
        dec  = float(header['OBSBET'])
        dx = (xmap+xoff)/3600.
        dy = (ymap+yoff)/3600.
        # House keeping
        alti_sta = header['ALTI_STA']
        alti_end = header['ALTI_END']
        za_sta = header['ZA_START']
        za_end = header['ZA_END']
        wv_sta = header['WVZ_STA']
        wv_end = header['WVZ_END']
        angle = header['DET_ANGL']
        filename = header['FILENAME']
        filenum = int(filename[:5])            
        data = np.float32(data)+2**15  # signed integer to float
        data *= 3.63/65536.            # ADU to V
        nramps = np.size(data[:,0,0])
        if nramps < (ncycles*4*ngrat*32):
            print ("WARNING: Number of ramps does not agree with header for ",fitsfile)
        else:
            data = data[:ncycles*4*ngrat*32,1:17,:25]
            flux = data.reshape(ngrat,ncycles*4*32,16,25)
            gratpos = start+step*np.arange(ngrat)
            aor = (detchan, order, dichroic, ncycles, nodbeam, filegpid, filenum)
            hk  = (obsdate, (ra,dec), (dx,dy), angle, (za_sta,za_end), (alti_sta,alti_end), (wv_sta,wv_end))
            return aor, hk, gratpos, flux


            

def computeSlope(i,data):
    ''' Module called by  multiSlope to compute ramp slopes '''   
    import numpy as np
    from lmfit.models import LinearModel
    linmodel = LinearModel()

    ds = np.shape(data)
    nramps = ds[1] // 32
    satlim = 2.7
    dtime = 1./250.  ## 250 Hz
    x = dtime * np.arange(32)
    slopes = []
    for ngramps in data:
        ngslopes = []
        for i2 in range(16):
            # Reorganize data in ramps of 32 readouts
            ramps = ngramps[:,i2].reshape(nramps,32)
            # select ramps
            ramp1 = ramps[1::4,:]
            ramp3 = ramps[3::4,:]
            # mask saturated values
            ramp1[ramp1 > satlim] = np.nan
            ramp3[ramp3 > satlim] = np.nan
            # Differential ramps
            ramp = ramp1-ramp3
            # Condense ramp
            m = np.isfinite(ramp)
            if np.sum(m) > 0:
                y = np.nanmedian(ramp,axis=0)
                y[0]=np.nan  # mask 1st value
                # Consider only finite values
                mask = np.isfinite(y)
                # Linear part
                if np.sum(mask) > 10:  # At least 10 points to fit a line
                    pars = linmodel.guess(y[mask],x=x[mask])
                    out = linmodel.fit(y[mask],pars,x=x[mask])
                    ngslopes.append(out.params['slope'].value)
                else:
                    ngslopes.append(np.nan)
            else:
                ngslopes.append(np.nan)
        slopes.append(ngslopes)

    return i,slopes
    
        
def collectResults(results):
    results.extend(results)    

def multiSlopes(data):
    ''' Compute slopes for each pixel and grating position using multiprocessing '''    
    #import multiprocessing as mp
    #import numpy as np
    # If I don't import LinearModel here, then the multiprocessing stops after one file
    #from lmfit.models import LinearModel

        
    with mp.Pool(processes=mp.cpu_count()) as pool:
        res = [pool.apply_async(computeSlope, (i,data[:,:,:,i])) for i in range(25)]
        results = [r.get() for r in res]
    #results.sort()   not needed ...
        
    ds = np.shape(data)
    ng = ds[0]
    spectra = np.zeros((ng,16,25))
    for i,slopes in results:
        for ig,ngslopes in zip(range(ng),slopes):
            spectra[ig,:,i] = ngslopes

    return spectra

def waveCal(gratpos,dichroic,obsdate,array,order):
    import numpy as np
    import pandas as pd
    import os
    
    '''
    Usage:
    l,lw = waveCal( gratpos=1496600, order=1, array='RED',dichroic=105,obsdate='2015-03-12T04:41:33')
    '''

    if array == 'RED':
        channel = 'R'
    else:
        if order == '1':
            channel = 'B1'
        else:
            channel = 'B2'

    # Extract month and year from date
    year = obsdate.split('-')[0] 
    month = obsdate.split('-')[1]   
    odate = int(year[2:]+month)
        
    path0,file0 = os.path.split(__file__)
    wvdf = pd.read_csv(path0+'/CalibrationResults.csv', header=[0, 1])
    ndates = (len(wvdf.columns)-2)//5
    dates = np.zeros(ndates)
    for i in range(ndates):
        dates[i]= wvdf.columns[2+i*5][0]

    # Select correct date
    i = 0
    for date in dates:
        if date < odate:
            i+=1
        else:
            pass
    cols = range(2+5*(i-1),2+5*(i-1)+5)
    w1=wvdf[wvdf.columns[cols]].copy()
    if channel == 'R':
        if dichroic == 105:
            co = w1.columns[0]
        else:
            co = w1.columns[1]            
    elif channel == 'B1':
        co = w1.columns[2]
    else:
        if dichroic == 105:
            co = w1.columns[3]
        else:
            co = w1.columns[4]
    g0 = w1.ix[0][co]
    NP = w1.ix[1][co]
    a =  w1.ix[2][co]
    ISF = w1.ix[3][co]
    gamma = w1.ix[4][co]
    PS = w1.ix[5][co]
    QOFF = w1.ix[6][co]
    QS = w1.ix[7][co]
    ISOFF = w1.ix[8:][co].values

    pix = np.arange(16)+1.
    result = np.zeros((25,16))
    result_dwdp = np.zeros((25,16))
    for module in range(25):
        phi = 2.*np.pi*ISF*(gratpos+ISOFF[module])/2.0**24
        sign = np.sign(pix - QOFF) 
        delta = (pix - 8.5)*PS+sign*(pix-QOFF)**2 * QS
        slitPos = 25 - 6 * (module//5) + module%5
        g = g0*np.cos(np.arctan2(slitPos-NP,a))    # Careful here with arctan ...
        lambd  = 1000. * (g/order) * (np.sin(phi+gamma+delta)+np.sin(phi-gamma))
        dwdp   = 1000. * (g/order) * (PS + 2.*sign*QS*(pix-QOFF))*np.cos(phi+gamma+delta)
        result[module,:]=lambd
        result_dwdp[module,:]=dwdp
    
    return result, result_dwdp

class Obs(object):
    """ Single observation """
    def __init__(self,spectrum,coords,offset,angle,altitude,zenithAngle,waterVapor,nodbeam,fileGroupID,filenum,gratingPosition,channel,order,obsdate,dichroic):
        self.spec = spectrum
        self.ra = coords[0]
        self.dec = coords[1]
        self.x = offset[0]
        self.y = offset[1]
        self.angle = angle
        self.alt = altitude
        self.za = zenithAngle
        self.gp = gratingPosition
        self.wv = waterVapor
        self.fgid = fileGroupID
        self.n = filenum
        self.nod = nodbeam
        self.ch=channel
        self.order=order
        self.obsdate=obsdate
        self.dichroic=dichroic
