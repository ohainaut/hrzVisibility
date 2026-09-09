#!/usr/bin/python3
# makeHrz_vlt.py

# generate an ephemerides in VLT format,
# and a summary table

import numpy as np
import argparse
from datetime import datetime
from astropy.table import Table
from astropy.time import Time
from astropy.coordinates import Angle
import astropy.units as u
from astroquery.jplhorizons import Horizons

#----------------------------------------------------------------------
def reqSeeing(mag):
    '''return a seeing value [arcsec] appropriate to observe the input magnitude'''
    
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



#==============================================================================
#==============================================================================
#==============================================================================


expTmag25s10 = 240. #s for mag=25 snr=10
pointing = 375. #s, pointing overhead
readOut  = 23.  #s, readout overhead per image

elMin = 25. # deg, minimum elevation

print('VLT hrz')

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

print('get ephem...')

try:
    ephall = Table.read(  myargs.outFile+'v.ecsv' )
    print('read ephemerides from local file '+myargs.outFile+'v.ecsv')
except FileNotFoundError:
    ephall = Horizons( id=myargs.object, location=myargs.location, epochs=epochs ).ephemerides()

    # to select same night, work with pseudo MJD centred on Paranal Noon
    ephall['intMJD'] = (ephall['datetime_jd']-2400000.25).astype(int) # corresponding to noon PaO
    ephall.write(myargs.outFile+'v.ecsv', overwrite=True)
    print('wrote ephemerides to local file '+myargs.outFile+'v.ecsv')
print('        ...got ephem')





lEph = len(ephall[  ( ephall['EL'] > 27 ) & ( ephall["solar_presence"] != "*")      ])
if lEph == 0:
    raise ValueError(f'0 ephemeride line with elev>27 during night')
print(f'                              {lEph} lines observable')

# output
eph_file    = open(myargs.outFile+'v.eph', 'w') 
readme_file = open(myargs.outFile+'v.txt', 'w') 



readme_header = (f'Ephemerides in from Horizon; {len(ephall)} lines for {ephall["targetname"][0]}')
print(readme_header)
readme_file.write(readme_header+"\n\n")



eph_file.write(
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
fresh = True # this day has not meas yet
l0 = ephall[0]

          
jd0 =  ephall["intMJD"][0] 

# overall qualification of the moon, from 0 to 10k.
# myMoon<1000 is acceptable
ephall['myMoon'] = ((ephall["lunar_presence"] != "") * (180.-ephall["lunar_elong"])/1.8 * ephall["lunar_illum"]).astype(int) 



readme_header = (f'DateTime     \tmag \t"/h \tSg "'+
                  '\ttMx s \tDIT s \tNDIT \texpTs \t'+
                  'Tel m \tsnr1 \tsnrT \tstep " '+
                  '\tGlxLt \tSMAA"@Th \tRA+Dec   '+
                  '\tFLI@Elon \tObs[h] [FromTo]')
print(readme_header)
readme_file.write(readme_header+"\n")


for il in np.arange(len(ephall)):

    l = ephall[il]
    if l["EL"] > 23. and l["solar_presence"] != '*': #filter high airmasses and day

        # some conversions
        t = Time(l["datetime_jd"], format='jd').isot
        ra = Angle(l["RA"], 'degree').hms
        ras =  f'{ra[2]:09.6f}'
        if ras ==  "60.000000" : ras = "59.999999"

        de = Angle(l["DEC"], 'degree').signed_dms
        des =  f'{de[3]:08.5f}'
        if des == "60.00000": des = "59.999999"

        desig = "+" if de[0] >0 else "-"
        dra = l["RA_rate"]/3600.
        dde = l["DEC_rate"]/3600.


        # write VLT ephem
        eph_file.write(f'INS.EPHEM.RECORD          "{t}, {l["datetime_jd"]:17.9f}, '+
                f'{int(ra[0]):02d} {int(ra[1]):02d} {ras}, '+
                f'{desig}{int(de[1]):02d} {int(de[2]):02d} {des},'+
                f' {dra:+9.6f}, {dde:+9.6f}, , "\n')
        
        # check steps
        step = np.sqrt(  (l["RA"]-ra0)**2 + (l["DEC"] - de0)**2 )*3600.
        if step < 4e4 and step > 30.:
            ddcount += 1
            ddmax = max(ddmax, step)
        
        # info for OB
        if  l["intMJD"] != jd0  : #  new day
            jd0 = l["intMJD"]
            fresh = False

            # get all the lines for the same Night
            ephMyJD = ephall[(ephall["EL"] > elMin ) & 
                             (ephall["intMJD"] == jd0 ) &   
                             (ephall['solar_presence'] != "N" ) &  (ephall['solar_presence'] != "C" )  &  (ephall['solar_presence'] != "*" )  &
                             (ephall["myMoon"] < 1000 )]

            try:
                mag = l0["V"]
            except:
                mag = l0["RA_rate"]*0.+20.
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



            outLine =  (f'{l0["datetime_str"][:11]} \t{mag:.1f} \t{speed:.1f} \t{seeing} '
              + f'\t{ditMax:.1f} \t{dit} \t{nDit} \t{expTtot:<5.0f} '
              + f'\t{telTtot:<6.1f} \t{snrDit:.1f} \t{snrTot:.1f} \t{step:.1f} '
              + f'\t{l0["GlxLat"]:.1f} \t{l0["SMAA_3sigma"]:.1f}"@{l0["Theta_3sigma"]:.1f}'
              + f'\t{int(ra[0]):02d}:{int(ra[1]):02d}{desig}{int(de[1]):02d} ')
            if len(ephMyJD) == 0:
                outLine += '\t-NO-'
            else:
                outLine += f'\t{l0["lunar_illum"]/100:.2f}@{l0["lunar_elong"]:.0f}d' + \
                      f'\t{(ephMyJD["datetime_jd"][-1] - ephMyJD["datetime_jd"][0])*24.:.1f}h' + \
                      f' ({ephMyJD["datetime_str"][0][9:17]}-{ephMyJD["datetime_str"][-1][12:17]})'  

            print(outLine)
            readme_file.write(outLine+"\n") 


        l0 = l # preserve valid line for print if needed.

    ra0 = l["RA"] # current line;  valid or not
    de0 = l["DEC"]
    el0 = l["EL"]

if ddcount > 0:
    print(f'step too large {myargs.step} on {ddcount} epochs, max {ddmax:.2f}"')

eph_file.write("\n")
eph_file.close()
