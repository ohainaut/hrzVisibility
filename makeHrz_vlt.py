#!/usr/bin/python3
# makeHrz_vlt.py

import numpy as np
import argparse
from datetime import datetime
from astropy.time import Time
from astropy.coordinates import Angle
import astropy.units as u
from astroquery.jplhorizons import Horizons


def reqSeeing(mag):
    if mag > 26.:
        return 0.6
    elif mag > 25.:
        return 0.8
    elif mag > 24.:
        return 1.0
    elif mag > 22:
        return 1.2
    else:
        return 1.4

expTmag25s10 = 240. #s for mag=25 snr=10
pointing = 375. #s
readOut = 23. #s

#------------------------------------------------------------------------------

parser = argparse.ArgumentParser(description='Generate a visibility plot for a solar system object')
parser.add_argument('-f','--outFile', default="HRZ",
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
parser.add_argument('-t','--step', default='10m', 
                        help='ephem step, in JPL format')
parser.add_argument('-l','--location', default=309,
                        help='JPL location code')


myargs = parser.parse_args()



# Epochs required:
Ts = Time(myargs.start)
Te = Time(myargs.end)
epochs = {'start': Ts.value,
    'stop' : Te.value,
    'step' : myargs.step }


# READ DATA
ephall = Horizons( id=myargs.object, location=myargs.location, epochs=epochs ).ephemerides()
print(f'Ephemerides in from Horizon; {len(ephall)} lines for {ephall["targetname"][0]}')

lEph = len(ephall[  ( ephall['EL'] > 27 ) & ( ephall["solar_presence"] != "*")      ])
if lEph == 0:
    raise ValueError(f'0 ephemeride line with el>27 during night')
print(f'                              {lEph} lines observable')

# output
f = open(myargs.outFile+'v.eph', 'w') 
f.write(
f'''PAF.HDR.START;                                             # Start of PAF Header
PAF.TYPE                  "Instrument Setup";              # Type of PAF
PAF.ID                    "";                              # ID for PAF
PAF.NAME                  "{myargs.outFile}";                   # Name of PAF
PAF.DESC                  "Target body name: {ephall["targetname"][0]}     "
PAF.DESC                  "Center body name: Earth (399)  "
PAF.DESC                  "Center-site name: {myargs.location} "
PAF.DESC                  "Start time      : {Ts.isot}            "
PAF.DESC                  "Stop  time      : {Te.isot}            "
PAF.DESC                  "Step-size       : {myargs.step}        "
PAF.DESC                  "Atmos refraction: NO (AIRLESS)                                             "
PAF.CRTE.NAME             "makeHrz_vlt";                # Name of creator
PAF.CRTE.DAYTIM           "{datetime.now().isoformat()}";           # Civil Time for creation
PAF.LCHG.NAME             "";                              # Name of person/appl. changing
PAF.LCHG.DAYTIM           "";                              # Timestamp of last change
PAF.CHCK.NAME             "";                              # Name of appl. checking
PAF.HDR.END;                                               # End of PAF Header

#------------------------------------------------------------------------------

TPL.FILE.DIRNAME          "$INS_ROOT/SYSTEM/DETDATA";      # Storage filename

#------------------------------------------------------------------------------
''')

ra0 = -99999.
de0 = -99999.
ddmax = -99
ddcount = 0
el0 = 9999.
jd0 =  int( ephall["datetime_jd"][0] ) 
fresh = True # this day has not meas yet
l0 = ephall[0]

print(f'DateTime          \tmag \t"/h \tSg "'+
                  f'\ttMx s \tDIT s \tNDIT \texpTs \t'+
                  f'Tel m \tsnr1 \tsnrT \tstep " '+
                  f'\tGlxLt \tSMAA"@Th \tRA+Dec')


for il in np.arange(len(ephall)):

    l = ephall[il]

    if l["EL"] > 27. and l["solar_presence"] != '*': #filter high airmasses and day

        # some conversions

        t = Time(l["datetime_jd"], format='jd').isot
        ra = Angle(l["RA"], 'degree').hms
        if  f'{ra[2]:09.6f}' == "60.000000": # catch rounding too close
            ra[2] = 59.9999

        de = Angle(l["DEC"], 'degree').signed_dms
        if  f'{de[3]:08.5f}' == "60.00000": # catch rounding too close
            de[3] = 59.9999

        des = "+" if de[0] >0 else "-"
        dra = l["RA_rate"]/3600.
        dde = l["DEC_rate"]/3600.


        # write VLT ephem

        f.write(f'INS.EPHEM.RECORD          "{t}, {l["datetime_jd"]:17.9f}, '+
                f'{int(ra[0]):02d} {int(ra[1]):02d} {ra[2]:09.6f}, '+
                f'{des}{int(de[1]):02d} {int(de[2]):02d} {de[3]:08.5f},'+
                f' {dra:+9.6f}, {dde:+9.6f}, , "\n')
        
        # check steps

        step = np.sqrt(  (l["RA"]-ra0)**2 + (l["DEC"] - de0)**2 )*3600.
        if step < 4e4 and step > 30.:
            ddcount += 1
            ddmax = max(ddmax, step)
        
        # info for OB

        if int( l["datetime_jd"]) != jd0:
            fresh = True # new night

        if fresh and  (  l['EL'] < el0     or   int( l["datetime_jd"]) != jd0)  : # close to transit, new day
            jd0 = int( l["datetime_jd"]) 
            fresh = False

            mag = l0["V"]
            speed = np.sqrt( l0["RA_rate"]**2  + l0["DEC_rate"]**2) # arcsec/h
            seeing = reqSeeing( mag )
            ditMax = 3600./speed * seeing
            dit = min( max( int( ditMax/5. )*5, 5) , 60 ) # max 60s


            expTs10 = expTmag25s10 *10.**( 0.8* (mag - 25.) ) # time for 10sigma
            snrDit = np.sqrt(dit / expTs10)*10. # snr for 1DIT

            nDit = max( 5, int( expTs10 / dit + 1)) 
            expTtot = dit * nDit
            telTtot = (pointing + nDit * ( readOut + dit )) / 60. # min
            snrTot = np.sqrt(expTtot / expTs10)*10.

            print(f'{l0["datetime_str"]} \t{mag:.1f} \t{speed:.1f} \t{seeing} '+
                  f'\t{ditMax:.1f} \t{dit} \t{nDit} \t{expTtot} '+
                  f'\t{telTtot:.2f} \t{snrDit:.1f} \t{snrTot:.1f} \t{step:.1f} '+
                  f'\t{l0["GlxLat"]:.1f} \t{l0["SMAA_3sigma"]:.1f}"@{l0["Theta_3sigma"]:.1f}'+
                  f'\t{int(ra[0]):02d}:{int(ra[1]):02d}{des}{int(de[1]):02d} ')

        l0 = l # preserve valid line for print if needed.

    ra0 = l["RA"] # current line;  valid or not
    de0 = l["DEC"]
    el0 = l["EL"]

if ddcount > 0:
    print(f'step too large {myargs.step} on {ddcount} epochs, max {ddmax:.2f}"')

f.write("\n")
f.close()