# Several modules used for reading/reducing FIFI-LS data
import multiprocessing as mp
#from lmfit.models import LinearModel
import numpy as np

def exploreDirectory(path):
    ''' Explore FITS files in a directory '''
    # Should be smarter and open only files not yet considered
    
    from astropy.io import fits
    import os, glob
    import numpy as np
    from datetime import datetime
    
    files = sorted(glob.glob(os.path.join(path,'*w.fits')), key=os.path.getctime)
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


def readFlats(channel, order, dichroic, obsdate, silent=False):
    ''' Read flats '''
    wflat, specflat, especflat = readSpecFlats(channel, order, dichroic, silent=silent)
    spatflat = readSpatFlats(channel, obsdate, silent=silent)
    return wflat, specflat, especflat, spatflat

def readSpecFlats(channel, order, dichroic, silent=False):
    ''' Read flats '''
    import os
    from astropy.io import fits
    path0, file0 = os.path.split(__file__)
    if channel == 'RED':
        infile = os.path.join(path0,'data','spectralFlatsR1D'+str(dichroic)+'.fits.gz')
    else:
        infile = os.path.join(path0, 'data', 'spectralFlatsB'+str(order)+'D'+str(dichroic)+'.fits.gz')
    hdl = fits.open(infile)
    if silent == False:
        hdl.info()
    wflat = hdl['WAVE'].data
    specflat = hdl['SPECFLAT'].data
    especflat = hdl['ESPECFLAT'].data
    hdl.close()
    return wflat, specflat, especflat

def readSpatFlats(channel, obsdate, silent=False):
    ''' Read spatial flats.'''
    import os, re
    import numpy as np
    path0, file0 = os.path.split(__file__)
    if channel == 'RED':
        infile = os.path.join(path0, 'data', 'spatialFlatR.txt')
    else:
        infile = os.path.join(path0, 'data', 'spatialFlatB.txt')
    data = np.genfromtxt(infile,dtype='str',skip_header=1)
    dates = data[:,0].astype(int)
    spatflats = data[:, 1:].astype(float)
    # Extract month, year, and day from date
    parts = re.split('-|T|:', obsdate)
    odate = int(parts[0]+parts[1]+parts[2])
    # Select correct date
    for date, spatflat in zip(dates, spatflats):
        if date < odate:
            pass
        else:
            return spatflat

def applyFlats(waves, fluxes, channel, order, dichroic, obsdate):
    ''' Apply flats to fluxes '''
    
    wflat, specflat, especflat, spatflat= readFlats(channel, order, dichroic, obsdate, silent=True)
    for i in range(16):
        for j in range(25):
            sf = np.interp(waves[:,i,j], wflat, specflat[:,j,i])
            fluxes[:,i,j] /= sf
    for j in range(25):
        fluxes[:,:,j] /= spatflat[j]
    # Apply bad pixel mask
    import os
    path0, file0 = os.path.split(__file__)
    if channel == 'RED':
        bads = np.loadtxt(os.path.join(path0,'data','badpixels_2019_r.txt'))
    else:
        bads = np.loadtxt(os.path.join(path0,'data','badpixels_2019_b.txt'))
    for bad in bads:
        j,i = bad
        fluxes[:,np.int(i)-1,np.int(j)-1] = np.nan
    return fluxes



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
            #print("file ", fitsfile, " read")
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
    onslopes = []
    offslopes = []
    for ngramps in data:
        onngslopes = []
        offngslopes = []
        for i2 in range(16):
            # Reorganize data in ramps of 32 readouts
            ramps = ngramps[:,i2].reshape(nramps,32)
            # select ramps
            ramp0 = ramps[0::4,:]
            ramp1 = ramps[1::4,:]
            ramp2 = ramps[2::4,:]
            ramp3 = ramps[3::4,:]
            rampOn = np.concatenate((ramp0, ramp1), axis=0)
            rampOff = np.concatenate((ramp2, ramp3), axis=0)
            # mask saturated values
            rampOn[rampOn > satlim] = np.nan
            rampOff[rampOff > satlim] = np.nan
            # Differential ramps
            # ramp = ramp1-ramp3
            # Condense ramp
            m = np.isfinite(rampOn)
            if np.sum(m) > 0:
                # Get rid of first ramp
                y1 = np.nanmedian(rampOn[1:,:],axis=0)
                y2 = np.nanmedian(rampOff[1:,:],axis=0)
                # mask 1st value
                y1[0]=np.nan
                y1[1]=np.nan
                y1[-1]=np.nan
                y2[0]=np.nan
                y2[1]=np.nan
                y2[-1]=np.nan
                dtime = 1./250.  ## 250 Hz
                x = dtime * np.arange(32)
                # mask ramp
                mask = np.isfinite(y1)
                # Linear part
                if np.sum(mask) > 10:  # At least 10 points to fit a line
                    pars = linmodel.guess(y1[mask],x=x[mask])
                    out = linmodel.fit(y1[mask],pars,x=x[mask])
                    slope1 = out.params['slope'].value
                else:
                    slope1 = np.nan
                mask = np.isfinite(y2)
                # Linear part
                if np.sum(mask) > 10:  # At least 10 points to fit a line
                    pars = linmodel.guess(y2[mask],x=x[mask])
                    out = linmodel.fit(y2[mask],pars,x=x[mask])
                    slope2 = out.params['slope'].value
                else:
                    slope2 = np.nan
                #slope = slope1 - slope2 
                #ratio = (slope1 - slope2) / (slope1 + slope2) * 2
                #slope = 1/ratio - 0.5  # Ratio Flux/Sky
                # median sky
                #slope = 0.5 * (slope1 + slope2)
                onngslopes.append(slope1)
                offngslopes.append(slope2)
            else:
                onngslopes.append(np.nan)
                offngslopes.append(np.nan)
        onslopes.append(onngslopes)
        offslopes.append(offngslopes)

    return i, onslopes, offslopes
    
        
def collectResults(results):
    results.extend(results)    

def multiSlopes(data):
    ''' Compute slopes for each pixel and grating position using multiprocessing '''    
    #import multiprocessing as mp
    #import numpy as np
    # If I don't import LinearModel here, then the multiprocessing stops after one file
    #from lmfit.models import LinearModel

    # To avoid forking error in MAC OS-X
    try:
        mp.set_start_method('spawn')
    except RuntimeError:
        pass
        
    with mp.Pool(processes=mp.cpu_count()) as pool:
        res = [pool.apply_async(computeSlope, (i,data[:,:,:,i])) for i in range(25)]
        results = [r.get() for r in res]
    #results.sort()   not needed ...
        
    ds = np.shape(data)
    ng = ds[0]
    onspectra = np.zeros((ng,16,25))
    offspectra = np.zeros((ng,16,25))
    for i, onslopes, offslopes in results:
        for ig, onngslopes, offngslopes in zip(range(ng), onslopes, offslopes):
            onspectra[ig,:,i] = onngslopes
            offspectra[ig,:,i] = offngslopes

    return onspectra, offspectra

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
    obsdate = int(year[2:]+month)
        
    path0,file0 = os.path.split(__file__)
    wvdf = pd.read_csv(os.path.join(path0, 'data', 'CalibrationResults.csv'), header=[0, 1])
    ndates = (len(wvdf.columns)-2)//5
    dates = np.zeros(ndates)
    for i in range(ndates):
        dates[i]= wvdf.columns[2+i*5][0]

    # Select correct date
    # Select correct date
    for i, date in enumerate(dates):
        if date < int(obsdate):
            pass
        else:
            break
    cols = range(2 + 5 * i , 2 + 5 * i + 5)
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
    def __init__(self,wave,spectrum, skyspectrum, coords,offset,angle,altitude,zenithAngle,waterVapor,nodbeam,fileGroupID,filenum,gratingPosition,channel,order,obsdate,dichroic):
        self.wave = wave
        self.spec = spectrum
        self.sky  = skyspectrum
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
