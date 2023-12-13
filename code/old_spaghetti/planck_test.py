import numpy as np
import matplotlib.pyplot as plt
import pixell
import healpy as hp
import astropy
import pandas as pd

# Header from file
"""
SIMPLE  =                    T /Written by IDL:  Tue May 19 09:42:53 2015       BITPIX  =                    8 / 
NAXIS   =                    0 /                                                EXTEND  =                    T /File contains extensions         
DATE    = '2014-11-08'         /Date file created:                              INSTRUME= 'LFI     '           /Low Frequency Instrument                       
VERSION = 'PCCS2   '           /Version of PCCS (PCCS2/PCCS2_E)                 ORIGIN  = 'LFI-DPC '           /Name of organisation responsible for data       
TELESCOP= 'PLANCK  '           /Planck                                          CREATOR = 'OCT12-v1.1'         /Pipeline version                                
DATE-OBS= '2009-08-13'         /Start time of the survey                        DATE-END= '2013-08-03'         /End time of the survey                          
FWHM    = '32.293  '           /From an elliptical Gaussian fit to the beam     OMEGA_B = '1190.06 '           /Area of the main beam                           
FWHM_EFF= '32.408  '           /Computed from OMEGA_B assuming a Gaussian beam  OMEGA_B1= '1117.30 '           /Beam area within 1xFWHM_EFF radius              
OMEGA_B2= '1188.93 '           /Beam area within 2xFWHM_EFF radius              END           
XTENSION= 'BINTABLE'           /Written by IDL:  Tue May 19 09:42:53 2015       BITPIX  =                    8 /                                                
NAXIS   =                    2 /Binary table                                    NAXIS1  =                  208 /Number of bytes per row                         
NAXIS2  =                 1560 /Number of rows                                  PCOUNT  =                    0 /Random parameter count                          
GCOUNT  =                    1 /Group count                                     TFIELDS =                   35 /Number of columns                               
DATE    = '2015-05-19'         /Creation date                                   INSTRUME= 'LFI     '           /Low Frequency Instrument                        
VERSION = 'PCCS2   '           /Version of PCCS (PCCS2/PCCS2_E)                 ORIGIN  = 'LFI-DPC '           /Name of organisation responsible for data       
TELESCOP= 'PLANCK  '           /Planck                                          CREATOR = 'OCT12-v1.1'         /Pipeline version                                
DATE-OBS= '2009-08-13'         /Start time of the survey                        DATE-END= '2013-08-03'         /End time of the survey                          
FWHM    = '32.293  '           /From an elliptical Gaussian fit to the beam     OMEGA_B = '1190.06 '           /Area of the main beam                           
FWHM_EFF= '32.408  '           /Computed from OMEGA_B assuming a Gaussian beam  OMEGA_B1= '1117.30 '           /Beam area within 1xFWHM_EFF radius              
OMEGA_B2= '1188.93 '           /Beam area within 2xFWHM_EFF radius              EXTNAME = 'PCCS2_f030'         /Extension name                                  
TFORM1  = '23A     '           /Character string                                TTYPE1  = 'NAME    '           /Label for column 1                              
TUNIT1  = 'None    '           /Units of column 1                               TFORM2  = '1D      '           /Real*8 (double precision)                       
TTYPE2  = 'GLON    '           /Label for column 2                              TUNIT2  = 'degrees '           /Units of column 2                               
TFORM3  = '1D      '           /Real*8 (double precision)                       TTYPE3  = 'GLAT    '           /Label for column 3                              
TUNIT3  = 'degrees '           /Units of column 3                               TFORM4  = '1D      '           /Real*8 (double precision)                       
TTYPE4  = 'RA      '           /Label for column 4                              TUNIT4  = 'degrees '           /Units of column 4                               
TFORM5  = '1D      '           /Real*8 (double precision)                       TTYPE5  = 'DEC     '           /Label for column 5                              
TUNIT5  = 'degrees '           /Units of column 5                               TFORM6  = '1E      '           /Real*4 (floating point)                         
TTYPE6  = 'DETFLUX '           /Label for column 6                              TUNIT6  = 'mJy     '           /Units of column 6                               
TFORM7  = '1E      '           /Real*4 (floating point)                         TTYPE7  = 'DETFLUX_ERR'        /Label for column 7                              
TUNIT7  = 'mJy     '           /Units of column 7                               TFORM8  = '1E      '           /Real*4 (floating point)                         
TTYPE8  = 'APERFLUX'           /Label for column 8                              TUNIT8  = 'mJy     '           /Units of column 8                               
TFORM9  = '1E      '           /Real*4 (floating point)                         TTYPE9  = 'APERFLUX_ERR'       /Label for column 9                              
TUNIT9  = 'mJy     '           /Units of column 9                               TFORM10 = '1E      '           /Real*4 (floating point)                         
TTYPE10 = 'PSFFLUX '           /Label for column 10                             TUNIT10 = 'mJy     '           /Units of column 10                              
TFORM11 = '1E      '           /Real*4 (floating point)                         TTYPE11 = 'PSFFLUX_ERR'        /Label for column 11                             
TUNIT11 = 'mJy     '           /Units of column 11                              TFORM12 = '1E      '           /Real*4 (floating point)                         
TTYPE12 = 'GAUFLUX '           /Label for column 12                             TUNIT12 = 'mJy     '           /Units of column 12                              
TFORM13 = '1E      '           /Real*4 (floating point)                         TTYPE13 = 'GAUFLUX_ERR'        /Label for column 13                             
TUNIT13 = 'mJy     '           /Units of column 13                              TFORM14 = '1E      '           /Real*4 (floating point)                         
TTYPE14 = 'GAU_SEMI1'          /Label for column 14                             TUNIT14 = 'arcmin  '           /Units of column 14                              
TFORM15 = '1E      '           /Real*4 (floating point)                         TTYPE15 = 'GAU_SEMI1_ERR'      /Label for column 15                             
TUNIT15 = 'arcmin  '           /Units of column 15                              TFORM16 = '1E      '           /Real*4 (floating point)                         
TTYPE16 = 'GAU_SEMI2'          /Label for column 16                             TUNIT16 = 'arcmin  '           /Units of column 16                              
TFORM17 = '1E      '           /Real*4 (floating point)                         TTYPE17 = 'GAU_SEMI2_ERR'      /Label for column 17                             
TUNIT17 = 'arcmin  '           /Units of column 17                              TFORM18 = '1E      '           /Real*4 (floating point)                         
TTYPE18 = 'GAU_THETA'          /Label for column 18                             TUNIT18 = 'degrees '           /Units of column 18 
TFORM19 = '1E      '           /Real*4 (floating point)                         TTYPE19 = 'GAU_THETA_ERR'      /Label for column 19                             
TUNIT19 = 'degrees '           /Units of column 19                              TFORM20 = '1E      '           /Real*4 (floating point)                         
TTYPE20 = 'GAU_FWHM_EFF'       /Label for column 20                             TUNIT20 = 'arcmin  '           /Units of column 20                              
TFORM21 = '1E      '           /Real*4 (floating point)                         TTYPE21 = 'P       '           /Label for column 21                             
TUNIT21 = 'mJy     '           /Units of column 21                              TFORM22 = '1E      '           /Real*4 (floating point)                         
TTYPE22 = 'P_ERR   '           /Label for column 22                             TUNIT22 = 'mJy     '           /Units of column 22                              
TFORM23 = '1E      '           /Real*4 (floating point)                         TTYPE23 = 'ANGLE_P '           /Label for column 23                             
TUNIT23 = 'degrees '           /Units of column 23                              TFORM24 = '1E      '           /Real*4 (floating point)                         
TTYPE24 = 'ANGLE_P_ERR'        /Label for column 24                             TUNIT24 = 'degrees '           /Units of column 24                              
TFORM25 = '1E      '           /Real*4 (floating point)                         TTYPE25 = 'APER_P  '           /Label for column 25                             
TUNIT25 = 'mJy     '           /Units of column 25                              TFORM26 = '1E      '           /Real*4 (floating point)                         
TTYPE26 = 'APER_P_ERR'         /Label for column 26                             TUNIT26 = 'mJy     '           /Units of column 26                              
TFORM27 = '1E      '           /Real*4 (floating point)                         TTYPE27 = 'APER_ANGLE_P'       /Label for column 27                             
TUNIT27 = 'degrees '           /Units of column 27                              TFORM28 = '1E      '           /Real*4 (floating point)                         
TTYPE28 = 'APER_ANGLE_P_ERR'   /Label for column 28                             TUNIT28 = 'degrees '           /Units of column 28                              
TFORM29 = '1E      '           /Real*4 (floating point)                         TTYPE29 = 'P_UPPER_LIMIT'      /Label for column 29                             
TUNIT29 = 'mJy     '           /Units of column 29                              TFORM30 = '1E      '           /Real*4 (floating point)                         
TTYPE30 = 'APER_P_UPPER_LIMIT' /Label for column 30                             TUNIT30 = 'mJy     '           /Units of column 30                              
TFORM31 = '1I      '           /Integer*2 (short integer)                       TTYPE31 = 'EXTENDED'           /Label for column 31                             
TUNIT31 = '0/1     '           /Units of column 31                              TNULL31 =                   -1 /Null value for column 31                        
TFORM32 = '1I      '           /Integer*2 (short integer)                       TTYPE32 = 'EXT_VAL '           /Label for column 32                             
TUNIT32 = '0/1/2/3 '           /Units of column 32                              TNULL32 =                   -1 /Null value for column 32                        
TFORM33 = '24A     '           /Character string                                TTYPE33 = 'ERCSC   '           /Label for column 33                             
TUNIT33 = 'None    '           /Units of column 33                              TFORM34 = '23A     '           /Character string                               
TTYPE34 = 'PCCS    '           /Label for column 34                             TUNIT34 = 'None    '           /Units of column 34                              
TFORM35 = '1I      '           /Integer*2 (short integer)                       TTYPE35 = 'HIGHEST_RELIABILITY_CAT' /Label for column 35                        
TUNIT35 = '        '           /Units of column 35                              TNULL35 =                   -1 /Null value for column 35
"""

wmap_map_I = hp.read_cl("../../data/raw_planck_data/COM_PCCS_030_R2.04.fits")[:, 0]
print(wmap_map_I)

#hp.mollview(
#    wmap_map_I,
#    coord=["G", "E"],
#    title="Histogram equalized Ecliptic",
#    unit="mK",
#    norm="hist"
#)
#
plt.show()
