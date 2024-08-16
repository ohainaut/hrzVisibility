#!/usr/bin/python3
# Visibility.py
#
# Usage:
#  Visibility.py -f FILEroot -o OBJECT    -s START -e END
#  Visibility.py -f plot    -o "2024 FS5" -s 2024-01-01 -e 2024-12-31
#
# Output: FILEroot.pdf, with the plot
#
#
#.version 2019-01-25T17:19:11 Fri
#.history
#
import matplotlib
matplotlib.use('Agg')  # to avoid Xdisplay issues in remote
import matplotlib.pyplot as plt
import numpy as np
import argparse


from astropy.io import ascii
from astropy.time import Time
from astropy.coordinates import Angle
import astropy.units as u

from astroquery.jplhorizons import Horizons


# Some init:
alphabeta = 0.04 # solar phase effect, mag/deg
ZDmax = Angle(60.*u.deg)  # maximal observable zenithal distance
JD0 =  2454911.50 # an equinox, when RA_sun = 0. (good enough for our purpose)
epsr = Angle(-23.40*u.deg) #  inclination of the Earth orbit
# max zenith distance accepted:



my_dir = './' # working directory where I write files

##start, stop: YYYY-MM-DD
##step: 1d

#------------------------------------------------------------------------------

def computeObs():
    '''compute the object visibility
    Return:
        obs: duration of the observability (obj < ZDmax, sun < -12),
        T: HA of object at ZDmax = 1/2 duration of object < ZDmax
        nightLength: duration of sun < -12
    '''

    DECr = np.radians( ephall['DEC'])

    # Hour Angle when the object reaches ZDmax
    Tcos = (np.cos(ZDmax)-np.sin(LatOBSr)*np.sin( DECr ) ) / (np.cos(LatOBSr)*np.cos( DECr ))
    # Tcos < -1: always visible
    #      > +1: never visible

    # erradicate arccos(>1)
    Tcos = (Tcos <= 1) * Tcos + (Tcos > 1)*1.
    Tcos = (Tcos >= -1)* Tcos - (Tcos <-1)*1.
    T = Angle(np.arccos( Tcos )  )

    # sidereal time during which the object is observable
    # brings in -12, 12h

    RAngle = Angle( ephall['RA'] )

    HAmin = (RAngle -blackRA - T).wrap_at(180.*u.deg)
    #HAmin = (Tcos >=  0.999) *(-90.)*u.deg  + (Tcos <  0.999)*  HAmin
    #HAmin = (Tcos <= -0.999) *  90. *u.deg  + (Tcos > -0.999)*  HAmin

    HAmax = (RAngle -blackRA + T).wrap_at(180.*u.deg)
    #HAmax = (Tcos >=  0.999) *(-90.)*u.deg  + (Tcos <  0.999)*  HAmax
    #HAmax = (Tcos <= -0.999) *  90. *u.deg  + (Tcos > -0.999)*  HAmax



    # when the sun reaches -12deg = antisun reaches +12 = zd 78
    TSuncos =   (np.cos(78.*u.deg) -np.sin(LatOBSr)*np.sin(blackDec)) / (np.cos(LatOBSr)*np.cos(blackDec))

    # erradicate arccos(>1)
    TSuncos = (TSuncos <=  1)*TSuncos + (TSuncos > 1)*1.
    TSuncos = (TSuncos >= -1)*TSuncos - (TSuncos <-1)*1.
    TSun = Angle(np.arccos( TSuncos )  )

    blackHAmin = (-TSun).wrap_at(180.*u.deg)
    blackHAmax =  (TSun).wrap_at(180.*u.deg)


    # length of the night
    nightLength = 2*TSun
    #dbg# print( nightLength.deg/15. ## OK)
    #dbg# plt.plot(ephall['datetime_jd'],nightLength.deg/15.,label="O.Plane", c=mycol) ## OK



    # find the max of the 2 beginning ST = begin obs
    # done trigonometrically to deal with funny angles.

    # begin of night or object raise?
    w = HAmin-blackHAmin
    h1 = (np.sin(w)>0.)*HAmin + (np.sin(w)<=0.)* blackHAmin

    # find the min of the 2 end ST = end obs
    # end of night or object set?
    w = HAmax-blackHAmax
    h2 = (np.sin(w)<0.)*HAmax + (np.sin(w)>=0.)*blackHAmax

    # obs is the observability (h2-h1)
    w = h2 - h1
    obs = (np.sin(w)>0.)*w + 0.

    # erradicate artefacts at solar conjuction
    obs = (Tcos < 0.999) * obs
    obs = (obs < T*2.001) * obs

    return obs, T, nightLength
#------------------------------------------------------------------------------

def densGlx(l,b):
    '''Return a star density flag in the range [0, 1+bulgex]
        with bulgex=2:
        good: [0, 0.5]
        med: [0.5,1]
        bad: [1,2]
        verybag [2,[

    This is done using a very basic model of star density
    purely on geometric parameters.
    '''

    lwidth = 30. # bulge width in longitude
    bthick = 10. # disk thickness in latitude; also used as 2x bulge thickness
    bulgex = 2.  # bulge intensity scaling

    l =  (l+180.)%360 -180
    return (np.exp(-( b/ bthick) **2   ) +
            bulgex*np.exp(-( (l/lwidth)**2 + (b/2./bthick)**2   )))

#------------------------------------------------------------------------------
def doXticks():
    # deal with the X-ticks for the plots
    if tickMonthFlag:
        _ = ax1.set_xticks(tickmonths.jd)
    else:
        _ = ax1.set_xticks(tickyears.jd)
        _ = ax1.set_xticks(tickmonths.jd, minor=True)
    _ = ax1.set_xticklabels([])
    _ = ax1.grid(axis='x',color='k',linewidth=0.1,alpha=0.1,linestyle='-', visible=True)
    _ = ax1.set_xlim(limityears.jd)

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

parser = argparse.ArgumentParser(description='Generate a visibility plot for a solar system object')
parser.add_argument('-f','--outFile',
                        help='Root of the output file (no extension)')
parser.add_argument('-o','--object',
                        help='''Designation of the object;
                        must be resolved by Horizon;
                        in case of doubt use the Unique JPL ID
                        (in the ephem header, 1st line Rec #:''')
parser.add_argument('-s','--start',
                        help='Start time, YYYY-MM-DD')
parser.add_argument('-e','--end',
                        help='End time, a YYYY-MM-DD')
parser.add_argument('-d','--deltaMax', default=0.,
                        help='Max value of Delta for plot; leave to 0 for auto.')

myargs = parser.parse_args()

my_file = myargs.outFile
comet = myargs.object

ystart = 2019.
yend = 2040.5
cystart =  myargs.start
cyend =    myargs.end
Deltamax = myargs.deltaMax



# Epochs required:
Ts = Time(cystart)
Te = Time(cyend)
if (Te-Ts).value > 1000.:
    epochs = {'start': Ts.value,
              'stop' : Te.value,
              'step' : '10d' }
else:
    epochs = {'start': Ts.value,
              'stop' : Te.value,
              'step' : '1d' }


# READ DATA
ephall = Horizons( id=comet, location=500, epochs=epochs ).ephemerides()
print('Ephemerides in')

# absolute mag:
#TBD


# coordinates of the sun and antisun
lsunr = Angle((ephall['datetime_jd']- JD0)/365.25* 360. *u.deg) # RA of the sun, good enough for our purpose
sundec = np.arcsin( -np.sin(lsunr)*np.sin(epsr) )
sunraw = np.arccos(  np.cos(lsunr)/np.cos(sundec) )
sunra =  ((np.sin(lsunr)*np.cos(epsr)/np.cos(sundec) > 0)*2-1)  * sunraw
# coordinates of the antisun
blackRA = Angle(sunra + 180.*u.deg)
blackDec = Angle(-sundec)


# years for the labels; all these variables must be synchronized:
oliyears = np.arange(2000.,2061.)
tickyears = Time(oliyears, format='jyear',scale='utc')
yearlabels = ['2000','01','02','03','04','05','06','07','08','09','2010','11','12','13','14','2015','16','17','18','19','2020','21','22','23','24','2025','26','27','28','29','2030','31','32','33','34','2035','36','37','38','39','2040','41','42','43','44','45','46','47','48','49','2050','51','52','53','54','55','56','57','58','59','2060']
#months for the ticks
tickmonths = tickyears[0] + np.arange(0,500)*30.4375*24.*u.h # that's 1/12 of a yr
monthlabel =  np.arange(0,500) %12 +1

#- magnitude correction
#ephall['magRaw']  = ephall['M1'] + 5*np.log10(ephall['r']*ephall['delta'])
#ephall['magPh']   = ephall['magRaw']  + alphabeta*ephall['alpha']

ephall['magPh'] = ephall['V']
ephall['magRaw'] = ephall['V'] - alphabeta*ephall['alpha']


plt.rc('grid', linestyle="-", color='black')



##=====MAIN PLOTS==========

mycol1 = 'k'
mycol2 = 'b'
nfig = 7  # how many plots

limityears = Time([cystart,cyend], scale='utc')

fig,axAll = plt.subplots(nfig,1, figsize=(8.,11.))
fig.subplots_adjust(hspace=0, wspace=0)
print( limityears.jd)

tickMonthFlag = (limityears[1].jd - limityears[0].jd) < 400.
# true: show the months, only if plot is short enough

#ELT?
ELTflag = ( limityears[0].jd > 2462106.)

ithisplot = 0

#=PLOT 1=============== r, Delta, error, active, not active
ax1 = axAll[ithisplot]
ax1.cla()
ax1.set_title(comet, loc='left', fontsize=15)
ax1.set_title(cystart+' - '+cyend, loc='right')

ax1.set_yticks(np.arange(0.,100.,5.))
ax1.set_yticks(np.arange(0.,100.,1.),  minor=True)



mycol = mycol1
active = (np.cos( np.radians( ephall['true_anom']) )   >0.)*(ephall['r'] < 3.)

ractive = ephall['r']*active
ax1.fill_between(ephall['datetime_jd'],ractive, 0, facecolor='c', alpha=0.25, label="Active")
rinactive = ephall['r']*(1-active)
ax1.fill_between(ephall['datetime_jd'],rinactive, 0, facecolor='r', alpha=0.025, label="Not active")

ax1.legend(loc=1, ncol=3, framealpha=0.2)
ax1.plot(ephall['datetime_jd'],ephall['r'], c=mycol)



ax1.tick_params('y', colors=mycol)
ax1.set_xticks(tickyears.jd)
_ = ax1.set_xticklabels([])


ax1.set_ylabel("r (black)\n$\Delta$ (blue) [AU]",color=mycol)
ax1.tick_params('y', colors=mycol)

resolved = (ephall['delta']< Deltamax)  #(ephall['delta'] <= Deltamax )
Dresolved =ephall['delta']* resolved

ax1.plot(ephall['datetime_jd'],ephall['delta'], c=mycol2)
ax1.fill_between(ephall['datetime_jd'],Dresolved, 0, lw=0, facecolor='r', alpha=0.5)

rmin = min(ephall['delta'])
rmax = max(ephall['delta'])
wmin = rmin -0.05*(rmax - rmin)
wmax = rmax +0.40*(rmax - rmin)
ax1.set_ylim(wmin, wmax)
ax1.grid(axis='both',color='k',linewidth=0.1,alpha=0.1,linestyle='-')

ax2 = ax1.twinx()
ax2.set_yscale('log')

ax2.plot(ephall['datetime_jd'],ephall['SMAA_3sigma'], "-r")
ax2.plot(ephall['datetime_jd'],ephall['SMIA_3sigma'], ":r")
ax2.set_ylabel("Position error [arcsec]",color='r')
ax2.tick_params('y', colors='r')

doXticks()

ithisplot += 1

#=PLOT 2================RA,DEC
ax1 = axAll[ithisplot]
ax1.cla()

mycol = mycol1
ax1.plot(ephall['datetime_jd'],ephall['RA']/15., c=mycol)
ax1.set_ylabel("R.A.",color=mycol)
ax1.tick_params('y', colors=mycol)
_ = ax1.set_yticks(np.arange(0.,24.,6.))
_ = ax1.set_yticks(np.arange(0.,25.,1.),  minor=True)
ax1.set_ylim(-0.5,24.5)

doXticks()



mycol = mycol2
ax2 = ax1.twinx()
_ = ax2.set_yticks(np.arange(-90.,91.,30.))
_ = ax2.set_yticks(np.arange(-90.,91.,10.),  minor=True)
ax2.plot(ephall['datetime_jd'],ephall['DEC'], c=mycol)
ax2.plot(ephall['datetime_jd'],ephall['DEC']*0.,  c=mycol, linewidth=0.2)
ax2.set_ylabel("Dec",color=mycol)
ax2.tick_params('y', colors=mycol)


ax2.set_xlim(limityears.jd)

ithisplot += 1



#=PLOT 3===========Long, Lat, TA
ax1 = axAll[ithisplot]
ax1.cla()

mycol = mycol1
ax1.plot(ephall['datetime_jd'],ephall['EclLon'], ':', c=mycol)
ax1.plot(ephall['datetime_jd'],ephall['true_anom'], c=mycol)
ax1.set_ylabel("Ecl.Long.(:)/TA(-)",color=mycol)

ax1.tick_params('y', colors=mycol)

ax1.set_ylim(-5.,365.)
_ = ax1.set_yticks(np.arange(90.,271.,90.))
_ = ax1.set_yticks(np.arange(0.,361.,30.),  minor=True)

doXticks()

mycol = mycol2
ax2 = ax1.twinx()
ax2.plot(ephall['datetime_jd'],ephall['EclLat'], ':', c=mycol)
ax2.plot(ephall['datetime_jd'],ephall['GlxLat'], c='r')
ax2.plot(ephall['datetime_jd'],ephall['GlxLat']*0., c=mycol, linewidth=0.2)

ax2.tick_params('y', colors=mycol)
_ = ax2.set_yticks(np.arange(-90,91,30))
_ = ax2.set_yticks(np.arange(-90,91,10),  minor=True)

ax2.set_ylabel("Latitude\nEcl.(:) Glx (red)",color=mycol)




ax2.set_xlim(limityears.jd)

ithisplot += 1
#=PLOT 4=============== PlAng, solar phase
ax1 = axAll[ithisplot]
ax1.cla()
ax2 = ax1.twinx()
ax2.cla()

elongNotOk = ephall['elong'] < 85.
elongJwstOk = ( ephall['elong'] >= 85.)*( ephall['elong'] <= 135.)
elongGroundOk =  ephall['elong'] > 135.
elongOk = ephall['elong'] >= 85.

SOTNotOk = ephall['elong']*elongNotOk
SOTJwstOk = ephall['elong']*elongJwstOk
SOTGroundOk = ephall['elong']*elongGroundOk

dustactive = active*  ephall['OrbPlaneAng']
mycol = mycol1

ax1.plot(ephall['datetime_jd'],ephall['OrbPlaneAng'],label="PlAng", c=mycol)
ax1.plot(ephall['datetime_jd'],ephall['OrbPlaneAng']*0., c=mycol, linewidth=0.2)
ax1.fill_between(ephall['datetime_jd'],dustactive, 0, facecolor='c', alpha=0.25)

plmax = max( 11., max(ephall['OrbPlaneAng']))
plmin = max( -11., min(ephall['OrbPlaneAng']))
ax1.set_ylim(plmin, plmax)
_ = ax1.set_yticks(np.arange( int((plmin-1.)/10.)*10., int((plmax+1)/10.)*10., 5), minor=True)

ax1.set_ylabel("Pl. Ang.",color=mycol)
ax1.tick_params('y', colors=mycol)

doXticks()

mycol = mycol2
ax2.plot(ephall['datetime_jd'],ephall['alpha'], c=mycol, alpha=0.5)
ax2.set_ylabel("Solar Phase",color=mycol)
ax2.tick_params('y', colors=mycol)

if True :  #Solar Elong
    ax2.plot(ephall['datetime_jd'],ephall['elong'],  c=mycol)

    ax2.fill_between(ephall['datetime_jd'],SOTNotOk, 0, facecolor='r', alpha=0.25, label='Low Elong.')
    ax2.fill_between(ephall['datetime_jd'],SOTJwstOk, 0, facecolor='limegreen', alpha=0.5, label='JWST+Ground')
    ax2.fill_between(ephall['datetime_jd'],SOTGroundOk, 0, facecolor='g', alpha=0.5, label='Ground')
    ax2.set_ylabel("Solar Elong.",color=mycol)
    _ = ax2.set_yticks(np.arange(0,181,30))
    _ = ax2.set_yticks(np.arange(0,180,10),  minor=True)
    ax2.legend(loc=1, ncol=5, framealpha=0.0)

    ax2.set_ylabel("Solar Elong.\n/Phase",color=mycol)
    ax2.set_ylim(0,215)


doXticks()

ithisplot += 1

#=PLOT 5=============== Mag
ax1 = axAll[ithisplot]

ax1.cla()
mycol = mycol1

#- groundbased limit
maglim = np.array([22.5, 24., 26.])
if ELTflag:
    maglim += 3.5

okspectr = (ephall['magPh'] < maglim[0])
okphot   = (ephall['magPh'] < maglim[1])*(ephall['magPh']>maglim[0])
okastrom = (ephall['magPh'] < maglim[2])*(ephall['magPh']>maglim[1])

magspectr = ephall['magPh']*okspectr + (1-okspectr)*99
magphot   = ephall['magPh']*okphot   + (1-okphot)*99
magastrom = ephall['magPh']*okastrom + (1-okastrom)*99

#- JWST lim
magjwst = 27.
okjwst =  (ephall['magPh'] < magjwst)

mx = max(ephall['magPh'])
mn = min(ephall['magRaw'])


wl = 'Spct (<{:.1f})'.format(maglim[0])
ax1.fill_between(ephall['datetime_jd'],magspectr, 40, facecolor='b', alpha=0.25, label=wl)
wl = 'Phot (<{:.1f})'.format(maglim[1])
ax1.fill_between(ephall['datetime_jd'],magphot,   40, facecolor='g', alpha=0.25, label=wl)
wl = 'Astr (<{:.1f})'.format(maglim[2])
ax1.fill_between(ephall['datetime_jd'],magastrom, 40, facecolor='y', alpha=0.25, label=wl)

wl = 'JWST (<{:.1f})'.format(magjwst)
ax1.fill_between(ephall['datetime_jd'],
                 mn+0.1*(mx-mn),
                 mn+0.1*(mx-mn) -2.*okjwst,
                 facecolor='limegreen', alpha=0.25, label=wl)

ax1.legend(loc=1, ncol=4, framealpha=0.2)


ax1.plot(ephall['datetime_jd'],ephall['magRaw'],label="Mag", c=mycol, linestyle=":", alpha=0.6)
ax1.plot(ephall['datetime_jd'],ephall['magPh'], c=mycol)


ax1.set_ylabel("Mag",color=mycol)
####ax1.text(limityears[1].jd,  mn, "M110={:.1f}  ".format(m110), ha='right')
ax1.set_ylim(mx+0.05*(mx-mn), mn-0.4*(mx-mn))

_ = ax1.set_yticks(range(int(mn-0.4*(mx-mn)-1), int(mx+0.05*(mx-mn)+1)),  minor=True)
_ = ax1.tick_params('y', colors=mycol)
_ = ax1.grid(axis='y',color='k',linewidth=0.1,alpha=0.1,linestyle='-', visible=True)

doXticks()




ithisplot += 1
#=PLOT 6===================== visibility
ax1 = axAll[ithisplot]
ax1.cla()

ax1.set_ylabel("Visibility [h]",color='k')
ax1.tick_params('y', colors='k')
_ = ax1.set_ylim(0.,15.)
_ = ax1.set_yticks(np.arange(0.,13.,2.))
_ = ax1.set_yticks(np.arange(0.,13.,1.),  minor=True)
_ = ax1.grid(axis='y',color='k',linewidth=0.1,alpha=0.1,linestyle='-')


doXticks()



##== PaO
LatOBSr = Angle(-24.67*u.deg)
obsPAO, TPAO, nightLenPAO = computeObs()
ax1.plot(ephall['datetime_jd'], 2.*TPAO.deg/15., c='g', linestyle=':')
ax1.plot(ephall['datetime_jd'],obsPAO.deg/15., c='g',label='PaO')

###== LSO
LatOBSr = Angle(-29.25*u.deg)
obsLSO, TLSO, nightLenLSO = computeObs()
ax1.plot(ephall['datetime_jd'], 2.*TLSO.deg/15., c='k', linestyle=':')
ax1.plot(ephall['datetime_jd'],obsLSO.deg/15., c='k', label="LSO")


###== CAN
LatOBSr = Angle(28.3*u.deg)
obsCAN, TCAN, nightLenCAN = computeObs()
ax1.plot(ephall['datetime_jd'], 2.*TCAN.deg/15., c='c', linestyle=':')
ax1.plot(ephall['datetime_jd'],obsCAN.deg/15., c='c', label="CAN")


##== MKO
LatOBSr = Angle(19.8*u.deg)
obsMKO, TMKO, nightLenMKO = computeObs()
ax1.plot(ephall['datetime_jd'], 2.*TMKO.deg/15., c='b', linestyle=':')
ax1.plot(ephall['datetime_jd'],obsMKO.deg/15., c='b', label="MKO")

#==
ax1.legend(loc=2, ncol=4, framealpha=0.)


ithisplot += 1
#=PLOT 7===================== Summary
ax1 = axAll[ithisplot]
ax1.cla()




y = 8
moond = ephall['lunar_elong'] > 90.
fli = ephall['lunar_illum'] < 25.
moondorfli = ((moond*1) + (fli*1)) > 0
moonok = moondorfli*1 + y
moonnotok = (1-moondorfli)*1 + y
ax1.fill_between(ephall['datetime_jd'],y, moonok,   linewidth=0, facecolor='g', alpha=0.25)
ax1.fill_between(ephall['datetime_jd'],y, moonnotok,linewidth=0, facecolor='r', alpha=0.25)

y -= 1.5
ephall['GlxDens'] = densGlx( ephall['GlxLon'],ephall['GlxLat'] )
ax1.fill_between(ephall['datetime_jd'], y, y+(ephall['GlxDens'] > 2), facecolor='r', alpha=0.5)
ax1.fill_between(ephall['datetime_jd'], y,
                 y+((ephall['GlxDens'] > 1)&(ephall['GlxDens'] < 2)),
                 facecolor='orange', alpha=0.5)
ax1.fill_between(ephall['datetime_jd'], y,
                 y+((ephall['GlxDens'] > 0.5)&(ephall['GlxDens'] < 1)),
                 facecolor='yellow', alpha=0.5)
ax1.fill_between(ephall['datetime_jd'], y,
                 y+(ephall['GlxDens'] < 0.5), facecolor='g', alpha=0.5)
ax1.text(limityears[1].jd,y,'Galaxy   ', ha='right', fontsize=7)


y -= 1.5
obsactive = active*1 + y
ax1.fill_between(ephall['datetime_jd'],y, obsactive, linewidth=0, facecolor='c', alpha=0.35)
ax1.text(limityears[1].jd,y,'Active   ', ha='right', fontsize=7)

y -= 1.5
obsjwst = (okjwst * elongJwstOk) +y
ax1.fill_between(ephall['datetime_jd'],y, obsjwst, linewidth=0, facecolor='limegreen', alpha=0.5)
ax1.text(limityears[1].jd,y,'JWST   ', ha='right', fontsize=7)


y -= 1.5
obsspectr = (okspectr*elongOk * (obsPAO.deg/15)>1.)*1 +y
obsphot   = (okphot  *elongOk * (obsPAO.deg/15)>1.)*1 +y
obsastrom = (okastrom*elongOk * (obsPAO.deg/15)>1.)*1 +y
ax1.fill_between(ephall['datetime_jd'],y, obsspectr, linewidth=0, facecolor='b', alpha=0.5)
ax1.fill_between(ephall['datetime_jd'],y, obsphot  , linewidth=0, facecolor='g', alpha=0.5)
ax1.fill_between(ephall['datetime_jd'],y, obsastrom, linewidth=0, facecolor='y', alpha=0.5)
ax1.text(limityears[1].jd,y,'PaO   ', ha='right', fontsize=7)

y -= 1.5
obsspectr = (okspectr*elongOk * (obsMKO.deg/15)>1.)*1 +y
obsphot   = (okphot  *elongOk * (obsMKO.deg/15)>1.)*1 +y
obsastrom = (okastrom*elongOk * (obsMKO.deg/15)>1.)*1 +y
ax1.fill_between(ephall['datetime_jd'],y, obsspectr, linewidth=0, facecolor='b', alpha=0.5)
ax1.fill_between(ephall['datetime_jd'],y, obsphot  , linewidth=0, facecolor='g', alpha=0.5)
ax1.fill_between(ephall['datetime_jd'],y, obsastrom, linewidth=0, facecolor='y', alpha=0.5)
ax1.text(limityears[1].jd,y,'MKO   ', ha='right', fontsize=7)




ax1.set_ylim(0.,10)
_ = ax1.set_yticks([])

# this is for year labels
_ = ax1.set_xticks(tickyears.jd)
_ = ax1.set_xticklabels(yearlabels)
_ = ax1.set_xticks(tickmonths.jd, minor=True)
_ = ax1.tick_params(axis='x', which='both', direction='in', bottom=True, top=True)

# and for month labels
if tickMonthFlag:
    _ = ax1.set_xticks(tickmonths.jd)
    _ = ax1.set_xticklabels(monthlabel)
_ = ax1.set_xlim(limityears.jd)
_ = ax1.grid(axis='x',color='k',linewidth=0.1,alpha=0.1,linestyle='-', visible=True)



#=CLOSE

fig.subplots_adjust(hspace=0, wspace=0)
plt.savefig(my_dir+"/"+my_file+".pdf")
print( "V: output plot in "+my_dir+"/"+my_file+".pdf")
