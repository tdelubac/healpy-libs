#
# Set of personal useful functions for healpy
#


def MaskBorders(map,nest=True):
    #
    # Masks pixels in the border of the footprint
    #
    import healpy as hp
    import copy

    nside = hp.npix2nside(len(map))
    masked_map = copy.copy(map)
    for ipix,ival in enumerate(map):
        if ival == 0:
            continue
        else:
            neighb = hp.get_all_neighbours(nside,ipix,nest=nest)
            for ineighb in neighb:
                if ineighb == -1:
                    continue
                if map[ineighb] == 0:
                    masked_map[ipix]=0
                    break
    return masked_map

def PlotMap(map,nest=True,clabel='density',s=6,savefile=None,zrange=None):
    #
    # Plot healpy map 
    #
    import numpy as np
    import matplotlib as mpl
    import healpy as hp
    from matplotlib import pyplot as p
    
    nside = hp.npix2nside(len(map))

    pixels = np.arange(hp.nside2npix(nside))
    pixTh,pixPh = hp.pix2ang(nside,pixels,nest)

    pixRA=[]
    for el in pixPh:
        res=180*el/np.pi
        if res<180:
            pixRA.append(res)
        else:
            pixRA.append(res-360)

    pixRA=np.array(pixRA)
    pixDEC=90.-180*pixTh/np.pi

    '''
    pixRA = pixRA[map>0]
    pixDEC = pixDEC[map>0]
    map = map[map>0]
    '''
    pixRA = pixRA[map!=0]
    pixRA[pixRA<-60]+=360# For Review
    pixDEC = pixDEC[map!=0]
    map = map[map!=0]
    # Set z-axis range 
    if zrange!=None:
        map[map<zrange[0]] = zrange[0] 
        map[map>zrange[1]] = zrange[1] 

    p.figure()
    p.scatter(pixRA,pixDEC,c = map,s=s,edgecolors='none',marker='s')
    p.clim(zrange)
    cbar = p.colorbar()
    cbar.set_label(clabel,fontsize=20)
    p.xlabel('RA',fontsize=20)
    p.ylabel('Dec',fontsize=20)
    if(savefile!=None):
        p.savefig(savefile)
    else:
        p.show()


def PlotHistMap(map,xlabel='density',savefile=None,bins=100,percentile=True,cumulative=False,normed=1):
    #
    # Plot the histogram of a given map
    #
    from matplotlib import pyplot as p

    map = map[map>0]
    
    p.hist(map,bins=bins,cumulative=cumulative,normed=normed,histtype='step')
    if percentile:
        mean = map.mean()
        p.axvline(mean*1.075,linestyle='dashed',color='k')
        p.axvline(mean*0.925,linestyle='dashed',color='k')
    p.xlabel(xlabel)
    if savefile != None:
        p.savefig(savefile)
    else:
        p.show()


def MakeMap(fitsfile,output=None,nside=256,nest=True,norm=True,masked=False):
    #
    # Save map as healpy format from fits input
    #
    import numpy as np
    import pyfits as pf
    import healpy as hp
    
    hdulist = pf.open(fitsfile)
    Cat = hdulist[1].data
    hdulist.close()

    pixarea = hp.nside2pixarea(nside,degrees=True)
    print 'nside = ',nside,' --> Pixel area (deg2) = ',pixarea

    tiles  = hp.ang2pix(nside,-Cat['dec']*np.pi/180.+np.pi/2.,Cat['ra']*np.pi/180.,nest)
    npix = hp.nside2npix(nside)
    n_hit_selec = np.zeros(npix)

    for itile in tiles:
        if norm:
            n_hit_selec[itile]+=1./pixarea
        else:
            n_hit_selec[itile]+=1
    if masked == True:
        n_hit_selec = MaskBorders(n_hit_selec)
    if output == None:
        return n_hit_selec
    else:
        hp.write_map(output,n_hit_selec,nest)   
        return

def MakeMapWCS(fitsfile,output=None,nside=256,nest=True,norm=True,masked=True):
    #
    # Save map as healpy format from WCS fits input
    #
    from astropy.io import fits as pf 
    from astropy.wcs import WCS
    import numpy as np
    import healpy as hp
    
    hdulist = pf.open(fitsfile)
    Cat = hdulist[0].data
    w = WCS(fitsfile)
    hdulist.close()

    pixarea = hp.nside2pixarea(nside,degrees=True)
    print 'nside = ',nside,' --> Pixel area (deg2) = ',pixarea

    x = np.arange(len(Cat[0]))
    y = np.zeros(len(Cat[0]))

    for i in np.arange(len(Cat)):
        if (i==0):
            continue
        x = np.concatenate((x,np.arange(len(Cat[0]))))
        y = np.concatenate((y,np.ones(len(Cat[0]))*i))
        print len(x),len(y)
    radec = w.all_pix2world(x,y,1)
    
    tiles  = hp.ang2pix(nside,-radec[1]*np.pi/180.+np.pi/2.,radec[0]*np.pi/180.,nest)
    npix = hp.nside2npix(nside)
    n_hit_selec = np.zeros(npix)
    val = np.zeros(npix)
    for ix,iy,itile in zip(x,y,tiles):
        n_hit_selec[itile]+=1
        val[itile]+=Cat[iy][ix]
    
    val[n_hit_selec>0]/=n_hit_selec[n_hit_selec>0]
    if masked == True:
        val = MaskBorders(val)
    if output == None:
        return val
    else:
        hp.write_map(output,val,nest)   
        return

def MakeMapWCS2(fitsfile,output=None,nside=256,nest=True,norm=True,masked=True):
    #
    # Save map as healpy format from WCS fits input
    #
    from astropy.io import fits as pf 
    from astropy.wcs import WCS
    import numpy as np
    import healpy as hp
    
    hdulist = pf.open(fitsfile)
    Cat = hdulist[0].data
    w = WCS(fitsfile)
    hdulist.close()

    pixarea = hp.nside2pixarea(nside,degrees=True)
    print 'nside = ',nside,' --> Pixel area (deg2) = ',pixarea

    npix = hp.nside2npix(nside)
    n_hit_selec = np.zeros(npix)
    val = np.zeros(npix)

    x = np.arange(len(Cat[0]))
    for i in np.arange(len(Cat)):
        y = np.ones(len(Cat[0]))*i
        radec = w.all_pix2world(x,y,1)
        tiles  = hp.ang2pix(nside,-radec[1]*np.pi/180.+np.pi/2.,radec[0]*np.pi/180.,nest)

        for ix,iy,itile in zip(x,y,tiles):
            n_hit_selec[itile]+=1
            val[itile]+=Cat[iy][ix]
        print i
    val[n_hit_selec>0]/=n_hit_selec[n_hit_selec>0]
    if masked == True:
        val = MaskBorders(val)
    if output == None:
        return val
    else:
        hp.write_map(output,val,nest)   
        return

def MakeMapColumn(fitsfile,column,output=None,nside=256,nest=True,masked=False):
    #
    # Save map as healpy format from fits input for column in input
    #
    import numpy as np
    import pyfits as pf
    import healpy as hp
    
    hdulist = pf.open(fitsfile)
    Cat = hdulist[1].data
    hdulist.close()

    pixarea = hp.nside2pixarea(nside,degrees=True)
    print 'nside = ',nside,' --> Pixel area (deg2) = ',pixarea

    tiles  = hp.ang2pix(nside,-Cat['dec']*np.pi/180.+np.pi/2.,Cat['ra']*np.pi/180.,nest)
    npix = hp.nside2npix(nside)
    n_hit_selec = np.zeros(npix)
    val_column = np.zeros(npix)

    for i,itile in enumerate(tiles):
            n_hit_selec[itile]+=1
            if isinstance(column, basestring):
                val_column[itile]+=Cat[column][i]
            else:
                val_column[itile]+=column[i]

    for i,val in enumerate(val_column):
        if n_hit_selec[i] > 0:
            val_column[i]/=n_hit_selec[i]
            
    if masked == True:
        val_column = MaskBorders(val_column)
    if output == None:
        return val_column
    else:
        hp.write_map(output,val_column,nest)   
        return


def MakePredictedDensityMap(syst_maps,syst_params,syst_limits,mean_dens,output,nside=256,nest=True):
    #
    # Create a healpy map of the expected density per pixel given the value of the syst in the pixels
    # and the value of the linear parameters for each systmatics. Must also give mean depth. 
    #
    import healpy as hp
    import numpy as np
    import copy
    cop_maps = copy.deepcopy(syst_maps)
    maps = []
    for i,map in enumerate(cop_maps):
        map[map<syst_limits[i][0]] = 0
        map[map>syst_limits[i][1]] = 0
        map/=map[map>0].mean()
        maps.append(map)
    
    predicted_dens = np.zeros(len(maps[0]))
    predicted_dens+= mean_dens
    for i,ipar in enumerate(syst_params):
        predicted_dens+= ipar*(maps[i]-1)

    for imap in maps:
        predicted_dens[imap==0] = 0

    hp.write_map(output,predicted_dens,nest=nest)


def ComputeFullSystCorrection(syst_maps,syst_params,syst_limits=None,output=None,nside=256,nest=True):
    #
    # Create a healpy map of the expected fluctuation density per pixel given the value of the syst in the pixels
    # and the value of the linear parameters for each systmatics.
    #
    import healpy as hp
    import numpy as np
    import copy
    cop_maps = copy.deepcopy(syst_maps)
    maps = []
    for i,map in enumerate(cop_maps):
        if syst_limits is not None:
            map[map<syst_limits[i][0]] = 0
            map[map>syst_limits[i][1]] = 0
        map/=map[map>0].mean()
        maps.append(map)
    
    predicted_fluc = np.zeros(len(maps[0]))
    for i,ipar in enumerate(syst_params):
        predicted_fluc+= ipar*(maps[i]-1)

    for imap in maps:
        predicted_fluc[imap==0] = 0

    if output == None:
        return predicted_fluc
    else:
        hp.write_map(output,predicted_fluc,nest=nest)
        return


def ComputeSystCorrectionM1(syst_maps,syst_params,syst_limits,syst_i,output=None,nside=256,nest=True):
    #
    # Create a healpy map of the expected fluctuation density per pixel but for syst_i, given the value of the syst in the pixels
    # and the value of the linear parameters for each systmatics.
    #
    import healpy as hp
    import numpy as np

    maps = []
    for i,isyst in enumerate(syst_maps):
        map = hp.read_map(isyst,nest=nest)
        map[map<syst_limits[i][0]] = 0
        map[map>syst_limits[i][1]] = 0
        map/=map[map>0].mean()
        maps.append(map)
    
    predicted_fluc = np.zeros(len(maps[0]))
    for i,ipar in enumerate(syst_params):
        if i==syst_i:
            continue
        predicted_fluc+= ipar*maps[i]

    for i,imap in enumerate(maps):
        if i==syst_i:
            continue
        predicted_fluc[imap==0] = 0
    
    if output == None:
        return predicted_fluc
    else:
        hp.write_map(output,predicted_fluc,nest=nest)
        return

def ChangeNside(map,nside):
    #
    # Convert map from nside0 to nside
    # if nside0 < nside call RebinMap
    # if nside0 > nside call ExtrapolMap
    #
    import healpy as hp
    nside0 = hp.npix2nside(len(map))

    if nside0 == nside:
        print 'ChangeNside warning: nside0 == nside, returning map unchanged'
        return map
    elif nside0 < nside:
        nmap = ExtrapolMap(map,nside)
        return nmap
    else:
        nmap = RebinMap(map,nside)
        return nmap

def RebinMap(map,nside,nest=True):
    #
    # Convert map from nside0 to nside for nside0 < nside
    #
    import healpy as hp
    import numpy as np

    nside0 = hp.npix2nside(len(map))
    
    if nside0 <= nside:
        print 'ERROR in RebinMap : nside0 >= nside, should call ChangeNside instead.'
        return
    
    pix = np.arange(len(map))
    ang = hp.pix2ang(nside0,pix,nest=nest)
    nmap = np.zeros(hp.nside2npix(nside))
    nentries = np.zeros(hp.nside2npix(nside))
    old2newpix = hp.ang2pix(nside,ang[0],ang[1],nest=nest)

    for ip,ival in zip(old2newpix,map):
#        if ival!=0:
        nmap[ip] += ival
        nentries[ip] += 1

    nmap[nentries>0]/=nentries[nentries>0]
    return nmap

def ExtrapolMap(map,nside,nest=True):
    #
    # Convert map from nside0 to nside for nside > nside0
    #
    import healpy as hp
    import numpy as np

    nside0 = hp.npix2nside(len(map))
    
    if nside0 >= nside:
        print 'ERROR in ExtrapolMap: nside0 <= nside, should call ChangeNside instead.'
        return
    
    pix = np.arange(hp.nside2npix(nside))
    ang = hp.pix2ang(nside,pix,nest=nest)
    nmap = np.zeros(hp.nside2npix(nside))
    old2newpix = hp.ang2pix(nside0,ang[0],ang[1],nest=nest)

    for ipnew,ipold in zip(pix,old2newpix):
        nmap[ipnew] = map[ipold]
    return nmap
        
def EqualAreaSamples(ra,dec,nside,diff_max=0.05,it_max=100,npix_ini=10,nest=True,verbose=True,optimize=True):
    #
    # Divide a ra,dec list into equal area samples. Elements == 0 are put in the same sample '-1'.
    # Try to optimise the variance in number of entries per pixels normalised by the mean to be below diff_max
    #
    # diff_max : Normalised difference in number of pixels (max(abs(nentries - mean(nentries)))/mean(nentries)) below wich optimization is considered succesful
    # it_max  : Maximum number of iteration in optimisation process
    # npix_ini : Initial number of pix to exchange at each step (must be small to obtain satisfactory samples) 
    #
    # Note : Variance exclude isolated samples
    #
    # Return two lists:
    #     - pixel number (in healpy standard)
    #     - Sample to which the pixel belongs
    # Both lists have len = len(map)
    #
    
    import healpy as hp
    import numpy as np
    from sklearn.neighbors import BallTree
#    import random
    theta,phi = RaDecToThetaPhi(ra,dec,degree=True)
    samples = hp.ang2pix(nside,theta,phi,nest=nest)
    '''If optimize == false return those samples'''
    set_samples = list(set(samples))
    if optimize == False:
        set_samples_sort = np.sort(set_samples)
        for i,isample in enumerate(set_samples_sort):
            samples[samples==isample] = i
        return samples
    
    '''Get nb. entries per sample'''
    nentries = [len(samples[samples==i]) for i in set_samples]
    '''Get neighbours of samples'''
    #neighb = np.transpose(hp.get_neighbours(nside,set_samples,nest=nest)[0])
    neighb = np.transpose(hp.get_all_neighbours(nside,set_samples,nest=nest))
    good_neighb = []
    for i,ineighb in enumerate(neighb):
        good_neighb_i = []
        i_theta,i_phi = hp.pix2ang(nside,set_samples[i],nest=nest)
        for isamp in ineighb:
            if isamp in set_samples:
                if set_samples.index(isamp) == i:
                    continue
                isamp_theta,isamp_phi = hp.pix2ang(nside,isamp,nest=nest)
                if (isamp_theta == i_theta) | (isamp_phi == i_phi):
                    continue
                good_neighb_i.append(set_samples.index(isamp))
        good_neighb.append(good_neighb_i)
    neighb = good_neighb

    '''Remove sample with no neighbours'''
    set_samples_n = [x for x,y in zip(set_samples,neighb) if len(y)>0] 
    nentries_n = [x for x,y in zip(nentries,neighb) if len(y)>0] 
    neighb_n = [x for x,y in zip(neighb,neighb) if len(y)>0] 
    '''Get RA DEC of samples'''
    samp_theta,samp_phi = hp.pix2ang(nside,set_samples_n,nest=nest)
    samp_ra,samp_dec = ThetaPhiToRaDec(samp_theta,samp_phi,degree=False)
    samp_radec = np.transpose([samp_ra,samp_dec])
    '''Get starting mean and diff'''
    mean = np.mean(nentries_n)
    diff = np.max(np.abs(nentries_n-mean))/float(mean)
    if verbose:
        print 'Initial diff :',diff

    '''Start optimization'''
    n_it = 0
    ind = np.arange(len(ra))

    np.random.seed(1)
    inf_loop = 0
    change_samp = True
    while float(diff) > diff_max:
        if diff*mean > 2*npix_ini:
            npix_ex = npix_ini 
        else:
            npix_ex = int(diff*mean/2.)
        if inf_loop >= 100:
            npix_ex = int(npix_ex/2)
            inf_loop = 0
        '''Get random sample'''
        if change_samp:
            chage_samp = False
            iprev = isamp
            while isamp == iprev:
                isamp = np.random.choice(np.arange(len(set_samples_n)))
        '''Look for bigger neighbours'''
        samp_neighbs = [set_samples_n[i] for i in neighb_n[isamp] if nentries_n[i]>nentries_n[isamp]+npix_ex]
        if len(samp_neighbs)==0:
            change_samp = True
            inf_loop+=1
            continue
        '''Get objects in those neighbouring samples'''
        good_ra = np.array([ra[i] for i in ind if samples[i] in samp_neighbs])
        good_dec = np.array([dec[i] for i in ind if samples[i] in samp_neighbs])
        good_ra*=np.pi/180.
        good_dec*=np.pi/180.
        good_ind = [i for i in ind if samples[i] in samp_neighbs]
        '''Construct BallTree with pixels of biggest samples'''
        radec = np.transpose([good_ra,good_dec])
        tree = BallTree(radec, leaf_size = 1, metric='haversine')
        exchange_radec = tree.query(samp_radec[isamp],k=npix_ex,return_distance=False)[0]
        '''Exchange ra,dec from bigger samples to smaller sample'''
        exchange_ind = [good_ind[i] for i in exchange_radec]
        samples[exchange_ind] = set_samples_n[isamp]
        nentries_n = [len(samples[samples==i]) for i in set_samples_n]
        '''Recompute mean and diff'''
        diff = np.max(np.abs(nentries_n-mean))/float(mean)
        if verbose:
            print 'nentries :',nentries_n,' mean :',mean,' n_it :',n_it
        '''Increase iter'''
        n_it+= 1
        if n_it >= it_max:
            if verbose:
                print 'Aborting optimization after',n_it,'iterations'
            break
    '''Correction'''
    

    if verbose:
        print 'Final diff :',diff
    print nentries_n,mean
    set_samples_sort = np.sort(set_samples)
    for i,isample in enumerate(set_samples_sort):
        samples[samples==isample] = i

    return samples

def RaToPhi(ra,degree=True):
    #
    # Convert from RA to healpy standard phi
    #
    import copy
    import numpy as np
    from numbers import Number
    cra = copy.copy(ra)
    if isinstance(ra,Number):
        cra = np.array(cra)
    if degree == True:
        cra*=np.pi/180.
    cra[cra<=0]+=2*np.pi
    phi = cra
    if isinstance(ra,Number):
        phi=float(phi)
    return phi

def PhiToRa(phi,degree=True):
    #
    # Convert from healpy standard phi to Ra
    #
    import copy
    import numpy as np
    from numbers import Number
    ra = copy.copy(phi)
    if isinstance(phi,Number):
        ra = np.array(ra)
    ra[ra>=np.pi]-=2*np.pi
    if degree == True:
        ra*=180./np.pi
    if isinstance(phi,Number):
        ra = float(ra)
    return ra

def DecToTheta(dec,degree=True):
    #
    # Convert from Dec to healpy standard theta
    #
    import copy
    import numpy as np
    cdec = copy.copy(dec)
    if degree == True:
        cdec*=np.pi/180.
    theta = -cdec+np.pi/2.  
    return theta

def ThetaToDec(theta,degree=True):
    #
    # Convert from healpy standard theta to Dec
    #
    import copy
    import numpy as np
    dec = -copy.copy(theta)+np.pi/2.  
    if degree == True:
        dec*=180./np.pi
    return dec

def RaDecToThetaPhi(ra,dec,degree=True):
    #
    # Convert from Ra,Dec to healpy standards theta,phi
    #
    phi = RaToPhi(ra,degree)
    theta = DecToTheta(dec,degree)
    return theta,phi

def ThetaPhiToRaDec(theta,phi,degree=True):
    #
    # Convert from healpy standards theta,phi to Ra,Dec
    #
    ra = PhiToRa(phi,degree)
    dec = ThetaToDec(theta,degree)
    return ra,dec
