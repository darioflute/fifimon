# Several modules used for reading/reducing FIFI-LS data


def exploreDirectory(path):
    ''' Explore FITS files in a directory '''
    # Should be smarter and open only files not yet considered
    
    from astropy.io import fits
    import os, glob
    import numpy as np

    print "path is ",path    
    files = sorted(glob.glob(path+'*sw.fits'), key=os.path.getctime)
    files = np.array(files)

    # Wavelenghts
    start=[]
    fgid=[]
    
    for infile in files:
        hlf = fits.open(infile)
        try:
            h = hlf[0].header
            proc = h['PROCSTAT']
            if proc == 'LEVEL_1':
                start.append(h['FIFISTRT'])
                fgid.append(h['FILEGPID'])
        except:
            pass
        hlf.close()

    start=np.array(start)
    fgid =np.array(fgid)    
    return files, start, fgid


def readData(fitsfile):
    from astropy.io import fits
    import numpy as np

    hdulist = fits.open(fitsfile)
    scidata = hdulist[1].data
    header = hdulist[0].header
    data = scidata.DATA
    hdulist.close()
    procstat = header['PROCSTAT']
    if procstat != 'LEVEL_1':
        print "This program works only with raw FIFI-LS files (Level 1)"    
        exit
    else:
        detchan = header['DETCHAN']
        obsdate = header['DATE-OBS']
        dichroic = header['DICHROIC']
        if detchan == 'RED':
            wvc = header['G_WAVE_R']
            ncycles = header['C_CYC_R']
            start = header['G_STRT_R']
            step = header['G_SZUP_R']
            ngrat = header['G_PSUP_R']
            order = 1
        else:
            wvc = header['G_WAVE_B']
            ncycles = header['C_CYC_B']
            start = header['G_STRT_B']
            step = header['G_SZUP_B']
            ngrat = header['G_PSUP_B']
            order = header['G_ORD_B']

        filegpid=header['FILEGPID']    
        nodbeam = header['NODBEAM']
        # Position
        xmap = header['DLAM_MAP']
        ymap = header['DBET_MAP']
        xoff = header['DLAM_OFF']
        yoff = header['DBET_OFF']
        ra   = header['OBSLAM']
        dec  = header['OBSBET']
        dx = xmap+xoff
        dy = ymap+yoff
        # House keeping
        altitude = header['ALTI_STA']
        za = header['ZA_START']
        wv = header['WVZ_STA']
        angle = header['TEL_ANGL']
        filename = header['FILENAME']
        filenum = int(filename[:5])            
        data = np.float32(data)+2**15  # signed integer to float
        data *= 3.63/65536.            # ADU to V
        nramps = np.size(data[:,0,0])
        if nramps < (ncycles*4*ngrat*32):
            flux = 0
            print "WARNING: Number of ramps does not agree with header for ",fitsfile
        else:
            data = data[:ncycles*4*ngrat*32,1:17,:25]
            flux = data.reshape(ngrat,ncycles*4*32,16,25)
            gratpos = start+step*np.arange(ngrat)
            aor = (detchan, order, dichroic, ncycles, nodbeam, filegpid, filenum)
            hk  = (obsdate, (ra,dec), (dx,dy), angle, za, altitude, wv)
        return aor, hk, gratpos, flux


def computeSlope(data,i):
    ''' Module called by  multiSlope to compute a ramp slope '''   
    from lmfit.models import LinearModel
    import numpy as np
    ds = np.shape(data)
    ng = ds[0]
    i1=i%ng
    i2=i/25/5
    i3=(i/5)%25
    ramps = data[i1,:,i2,i3]
    # Reorganize them in groups of 32
    nramps = np.size(ramps)/32
    ramps = ramps.reshape(nramps,32)
    # select ramps
    ramp1 = ramps[1::4,:]
    ramp3 = ramps[3::4,:]
    # mask saturated values
    satlim = 2.7
    m1 = ramp1 > satlim
    ramp1[m1] = np.nan
    m3 = ramp3 > satlim
    ramp3[m3] = np.nan
    # Differential ramps
    ramp = ramp1-ramp3
    # Condense ramp
    y = np.nanmedian(ramp,axis=0)
    # mask 1st value
    y[0]=np.nan
    dtime = 1./250.  ## 250 Hz
    x = dtime * np.arange(32)
    # mask ramp
    mask = np.isfinite(y)
    # Linear part
    base = LinearModel()
    if np.sum(mask) > 10:  # At least 10 points to fit a line
        pars = base.guess(y[mask],x=x[mask])
        out = base.fit(y[mask],pars,x=x[mask])
        slope =out.params['slope'].value
    else:
        slope = np.nan
    return i,slope

def multiSlopes(data):
    ''' Compute slopes for each pixel and grating position using multiprocessing '''    
    import multiprocessing as mp
    import numpy as np
    from fifitools import computeSlope

    ds = np.shape(data)
    ng = ds[0]
    ncpu = mp.cpu_count()
    pool = mp.Pool(processes=ncpu)
    results = [pool.apply_async(computeSlope, args=(data,i)) for i in range(ng*16*25)]
    results = [p.get() for p in results]
    results.sort()
    
    spectra = np.zeros((ng,16,25))
    for r in results:
        i = r[0]
        i1=i%5
        i2=i/25/5
        i3=(i/5)%25
        spectra[i1,i2,i3] = r[1]

    return spectra

def waveCal(gratpos,dichroic,obsdate,array,order):
    import numpy as np
    import pandas as pd
    
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
    odate = year[2:]+month
        
    path = '/Users/dfadda/Pipeline/Repository/fifi-ls/data/wave_cal/CalibrationResults.csv'
    wvdf = pd.read_csv(path, header=[0, 1])
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
    cols = range(i-1,i+4)
    
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
    def __init__(self,spectrum,coords,offset,angle,altitude,zenithAngle,waterVapor,nodbeam,fileGroupID,filenum,gratingPosition):
        self.spec = spectrum
        self.ra,self.dec = coords
        self.x,self.y = offset
        self.angle = angle
        self.alt = altitude
        self.za = zenithAngle
        self.gp = gratingPosition
        self.wv = waterVapor
        self.fgid = fileGroupID
        self.n = filenum
        self.nod = nodbeam
