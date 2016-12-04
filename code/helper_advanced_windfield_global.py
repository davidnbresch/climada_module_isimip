
##################################################################################

# Copyright Tobias Geiger 2016, geiger@pik-potsdam.de

##################################################################################

#### function to determine advanced winfdield
#### USE with programme advanced_windfield_global.pyplot

from __future__ import division
from pandas import Series, DataFrame
import pandas as pd
import numpy as np
import csv
import os, glob
from datetime import datetime
from time import time


### so far adapted to windfield by holland1980, holland2008, holland2010 cf peduzzi2012
### can be easily expanded to different wind-fields eg chavez2015

def windfield_holland(msize, res, cgps, ngps,penv,pcen,rmax,vmax,tint,prepcen, model):
    # msize=masksize in gridpoints, squared shape assumed
    # res=resolution in degrees, regular grid assumed
    # cgps=coordinate of grid center cgps=(xc,yc)
    # ngps = next gps, consecutive track point
    # penv= pressure around cyclone ie pressure at outermost closed isobar
    # pcen= central pressure in cyclone
    # rmax= radius of maximum winds
    # vmax= max observed wind in knots
    # tint= time interval between grid points
    # prepcen= previous time step central pressure in cyclone
    # static windfield model to choose
    # H80 = holland1980 including gradient wind
    # H08 = holland2008 as H80 but with bs instead of b
    # H10 = holland2010 as H08 so far (x_exp=0.5) but no gradient wind included!
    
    from math import radians, cos, sin, asin, sqrt, asin, acos

    def haversine(lon1, lat1, lon2, lat2):
        """
        Calculate the great circle distance between two points 
        on the earth (specified in decimal degrees)
        """
        # convert decimal degrees to radians 
        lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

        # haversine formula 
        dlon = lon2 - lon1 
        dlat = lat2 - lat1 
        a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
        c = 2 * asin(sqrt(a)) 

        # 6367 km is the radius of the Earth
        km = 6367 * c
        return km 
    
    def final_coord(i,j,xc,yc):
        xf,yf=[(j-temp)*res+xc, (temp-i)*res+yc]
        return (xf,yf)
    
    # holland symmetric and static wind field according to Holland1980 and with b=bs according to Holland2008
    def stat_holland(r_arr,rmax,b,penv,pcen,lat):
        lat=np.abs(lat)
        f=2*0.0000729 * sin(lat*np.pi/180.)
        # units are m/s
        v_arr=np.sqrt(((100*b/rho*(rmax/r_arr)**b*(penv-pcen)*np.exp(-(rmax/r_arr)**b))**0.5)**2+(1000*0.5*r_arr*f)**2) - 0.5*1000*r_arr*f
        v_arr[np.isnan(v_arr)]=0
        # translate to knots
        v_arr=v_arr/0.51444444444
        return v_arr
    
    def stat_holland_alt(r_arr,rmax,vms,b,lat):
        lat=np.abs(lat)
        f=2*0.0000729 * sin(lat*np.pi/180.)
        # units are m/s
        vms=vms*0.51444444444
        v_arr=np.sqrt(vms**2*((rmax/r_arr)**b * exp(1-(rmax/r_arr)**b))+(1000*0.5*r_arr*f)**2) - 0.5*1000*r_arr*f
        v_arr[np.isnan(v_arr)]=0
        # translate to knots
        v_arr=v_arr/0.51444444444
        return v_arr
    
    def stat_holland2010(r_arr,rmax,penv,pcen,bs,x_exp):        
        v_arr= ((100*bs*(penv-pcen)*(rmax/r_arr)**bs) / (rho*exp((rmax/r_arr)**bs)))**x_exp
        # translate to knots
        v_arr=v_arr/0.51444444444
        return v_arr
        
    def x_value(penv,pcen):
        xx=0.6*(1.-(penv-pcen)/215.)
        return xx
    
    #def stat_chavaz():
    
    def v_trans_speed():
        # calculate distance to next grid point and translational cyclone speed
        dist=haversine(lon1=xc,lat1=yc,lon2=xn,lat2=yn)
        dist= dist/1.852 # dist to nautical miles
        # hours = hours between track coordiantes
        hours=tint
        vtrans=dist/hours
        #print vtrans
        if vtrans>30:
            print 'max bound of vtrans (30knts) reached:', vtrans
            vtrans=30
        return vtrans,dist
    
    def b_value(vtrans,vmax,penv,pcen):
        if vmax>0:
            b=rho*np.exp(1)*(0.51444444444*(vmax-vtrans))**2/(100*(penv-pcen))
            if b < 1:
                b=1.0
            elif b > 2.5:
                b=2.5
        else:
            b=1.0
        return b
    
    def bs_value(vtrans,penv,pcen,prepcen,lat,xx,tint):
        vt_ms=vtrans*0.51444444444
        bs=-4.4e-5*(penv-pcen)**2 + 0.01*(penv-pcen) + 0.03*(pcen-prepcen)/tint - 0.014*np.abs(lat) + 0.15*vt_ms**xx + 1.0
        return bs
    
    def extra_v(pcen,lat,lon):
        # determined in knots 10min sustained wind
        #v_from_p = 1131.816+0.064*lat-0.031*lon-1.1000*pcen  # peduzzi
        v_from_p = 1142.000 + 0.056*lat - 0.032*lon - 1.111*pcen  # ibtracs 1980-2013 (r2=0.91)
        return v_from_p

    def extra_p(vmax,lat,lon):
        # determined in hPa
        #p_from_v = 1024.688+0.055*lat-0.028*lon-0.815*vmax # peduzzi
        p_from_v = 1024.388 + 0.047*lat - 0.029*lon - 0.818*vmax # ibtracs 1980 -2013 (r2=0.91)
        return p_from_v
    
    def extra_rmax(pcen, lat, lon):
        #cubic fit , ibtracs 1980-2013 (r2=0.22)
        #rmax_from_p= -45201.105207 + 146.043463*pcen - 0.157263*pcen**2 + 0.000056*pcen**3 + 0.097515*lat + 0.016056*lon
        # rmax extrapolated from historic distribution
        # rmax-thresholds in nm
        rmax1=15; rmax2=25; rmax3=50
        p1=950; p2=980; p3=1020
        if pcen <= p1:
            val=rmax1
        elif  (pcen > p1) & (pcen <= p2):
            val= (pcen-p1)*(rmax2-rmax1)/(p2-p1) + rmax1
        elif (pcen > p2):
            val= (pcen-p2)*(rmax3-rmax2)/(p3-p2) + rmax2
        rmax_from_p=val        
        return rmax_from_p
    
    def angular_array(cgps,msize):
        # define general array that contains angular wind-speed dependence for translational wind speed of cyclone
        # northward movement of cyclone assumed
        angarr=np.empty((msize, msize))
        temp=msize//2
        xc,yc=cgps
        for i in range(msize):
            for j in range(msize):
                xf,yf=[(j-temp), (temp-i)]
                coss = xf/np.sqrt(yf**2+xf**2)
                #northern hemisphere
                if yc >0:            
                    if yf >= 0:
                        angarr[i,j]=np.arccos(coss)*180/np.pi
                    else:
                        angarr[i,j]=360 - np.arccos(coss)*180/np.pi
                #southern hemisphere
                else:
                    if yf >= 0:
                        angarr[i,j]=np.arccos(coss)*180/np.pi
                    else:
                        angarr[i,j]=360 - np.arccos(coss)*180/np.pi
        return angarr
    
    def dynamic_windfield(vwind,vtrans,angarr,dist):
        # calculate angle of direction using vtrans and next track point            
        # determine direction, northern hemisphere deviation in degree from northward movement
        londiff=xn-xc
        latdiff=yn-yc  # if latdiff ... add or subract 90 degrees, test all possible configurations
        londist=haversine(lon1=xc,lat1=yc,lon2=xn,lat2=yc)  # distance along same latitude
        londist= londist/1.852
        latdist=haversine(lon1=xc,lat1=yc,lon2=xc,lat2=yn)  # distance along same longitude
        latdist= latdist/1.852
        # londiff for values abs(londiff) > 300 happen at point where 180 turns to -179 etc 
        if (londiff < 0):
            if londiff < -300:
                pref=1
            else:
                pref=-1                
        else:
            if londiff > 300:
                pref=-1
            else:
                pref=1   
        if dist>0:
            if latdiff < 0:
                direct=pref*(180 - np.arcsin(londist/dist)*180/np.pi)
                #direct=pref*(180 - np.arccos(latdist/dist)*180/np.pi)
            else:
                direct=pref*np.arcsin(londist/dist)*180/np.pi
                #direct=pref*np.arccos(latdist/dist)*180/np.pi
        else:
            print 'no movement of cyclone'
            direct=0
        #print dist, vtrans, direct, londist, londiff, latdiff
    
        angarr=angarr + direct
        angarr=np.cos(angarr*np.pi/180.)
        angarr[np.isnan(angarr)]=0
        if yc < 0:
            angarr=np.flipud(angarr)
            angarr=np.fliplr(angarr)

        # calculate windfield plus translation speed
        vfull=np.sqrt(vwind**2 + vtrans**2 + 2*vwind*vtrans*angarr)        
        return vfull
    
    #### start calcualtion of windfield
    #### static
    # center coordinates from which to measure distance
    xc,yc=cgps
    # next track point
    xn,yn=ngps
    temp=msize//2
    r_arr=np.empty((msize, msize))

    for i in range(msize):
        for j in range(msize):
            xf,yf=final_coord(i=i,j=j,xc=xc,yc=yc)
            r_arr[i,j]=haversine(xc,yc,xf,yf)
            #print xf,yf
            
    # data for windfield calculation
    rho=1.15
    
    vtrans,dist=v_trans_speed()
    
    angarr=angular_array(cgps=cgps,msize=msize)
    
     # extrapolate vmax and pcen from peduzzi analysis
    if (pcen > 850) & (vmax > 5):
        if (np.abs(prepcen -pcen) < 10*tint ):
            pass
        else:
            pcen=extra_p(vmax=vmax,lat=yc,lon=xc)
            #print 'neu pcen', pcen
    elif (pcen < 850) & (vmax > 5):
        pcen=extra_p(vmax=vmax,lat=yc,lon=xc)
        #print 'neu pcen', pcen
    elif (vmax <= 5) & (pcen >= 850):
        vmax=extra_v(pcen=pcen,lat=yc,lon=xc)
        #print 'neu vmax', vmax
    else:
        print 'both v and p undefined!'
    # make sure that penvp always larger than pcen
    if pcen > penv:
        pcen=penv
    #adjust also pressure at previous track point
    if prepcen < 850:
        prepcen=pcen
    # extrapolate rmax from pressure if not given and trnasfer to kilometers for analysis
    if rmax > 0:
        rmax=rmax*1.852  # nm to km
    else:
        rmax=extra_rmax(pcen=pcen, lat=yc, lon=xc)
        rmax=rmax*1.852  # nm to km
    #print rmax
        
    b=b_value(vtrans=vtrans,vmax=vmax,penv=penv,pcen=pcen)
    #print 'b=',b
    xx=x_value(penv=penv,pcen=pcen)
    bs=bs_value(vtrans=vtrans,penv=penv,pcen=pcen,prepcen=prepcen,lat=yc,xx=xx,tint=tint)
    #print 'bs=',bs, 'xx:',xx
    
    if model=='H80':
        vwind=stat_holland(r_arr=r_arr,rmax=rmax,b=b,penv=penv,pcen=pcen,lat=yc)
    elif model=='H08':
        # same as H80 but with bs instead of b, additional info required previous track point central pressure for pressure gradient
        vwind=stat_holland(r_arr=r_arr,rmax=rmax,b=bs,penv=penv,pcen=pcen,lat=yc)
    elif model=='H08_v':
        # same as H08 but uses max-wind speed instead of pressure
        vwind=stat_holland_alt(r_arr=r_arr,rmax=rmax,vms=vmax,b=bs,lat=yc)
    elif model=='H10':
        # in holland2010 x_exp varies with distance, not implemented
        # so far identical to holland2008 but without gradient wind
        ## note: NO gradient wind implemented for holland2010!!!
        x_exp=0.5
        vwind=stat_holland2010(r_arr=r_arr,rmax=rmax,penv=penv,pcen=pcen,bs=bs,x_exp=x_exp)
        
    if np.amin(vwind)< 0:
        print 'negative wind speed encountered, set to zero!'
        vwind[vwind < 0]=0
        
    # calculate vtrans wind field array assuming that effect of vtrans decreases with distance from eye
    r_arr_normed=rmax/r_arr
    r_arr_normed[r_arr_normed>1]=1
    vtrans_arr=vtrans*r_arr_normed

    #### dynamic
    vfull=dynamic_windfield(vwind=vwind,vtrans=vtrans_arr,angarr=angarr,dist=dist)
    
    #print 'penv:', penv, 'pcen:', pcen, 'vmax:', np.amax(vfull)
    
    vfull=np.round(vfull).astype(int)
    
    
    return vfull, pcen, vmax
