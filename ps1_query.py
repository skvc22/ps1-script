#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
script to query the Pan-STARRS API for a given RA, Dec, and radius, with options to filter the output 
for only stars or only galaxies. adapted from code by Drs. Ryan Ridden-Harper and Edward Schlafly.
"""

import pandas as pd
import argparse
from pdastro import pdastroclass, AnotB, AandB, AorB
import numpy as np
import os

class PS1_query_class:
    def __init__(self):
        self.ra = None
        self.dec = None
        self.radius = None

        self.filename = None
        
        self.verbose = 0
        self.debug = False

        self.filterlist = []
        self.outputtype = 'skinny'
        self.outputcols = None

        self.allcols = []

        self.maglimlist = [] 

        self.catalog = pdastroclass() # self.catalog.t to access the catalog dataframe

    def add_options(self, parser=None, usage=None):
        if parser == None:
            parser = argparse.ArgumentParser(usage=usage, conflict_handler="resolve")

        parser.add_argument('RA', help='right ascention coordinate')
        parser.add_argument('Dec', help='declination coordinate')
        parser.add_argument('radius', help='the radius in degrees around those coordinates to retireve data for')
        parser.add_argument('-os', '--onlystars', default=False, action="store_true", help='filter out all non-star objects')
        parser.add_argument('-og', '--onlygalaxies', default=False, action="store_true", help='filter out all non-galaxy objects')
        parser.add_argument('-v', '--verbose', action="count", dest="verbose",default=0)
        parser.add_argument('-d','--debug', default=False, action="store_true", help='Debugging output')
        parser.add_argument('-f', '--filename', type=str, help='file name of output file')
        parser.add_argument('--version', type=str, default='dr2', choices=['dr1','dr2'],help='the api version to use. default is dr2')
        parser.add_argument('--outputtype', type=str, default='skinny', choices=['skinny', 'medium','full'],
                            help='Skinny: return ra, dec, mag, magErr for each chosen filter. Medium: also include ra and dec error, more mag info, indicate if a row is a star. Full: return full catalog')
        parser.add_argument('--sexagesimal', default=False, action='store_true', help='return ra and dec in sexagesimal as well as decimal degrees')
        parser.add_argument('--maglim', nargs=3, default=None, action='append', help='upper and lower magnitude limits for each band. can pass multiple times to place limits on multiple bands')
        parser.add_argument('-g',  default=False, action="store_true", help='g band')
        parser.add_argument('-r',  default=False, action="store_true", help='r band')
        parser.add_argument('-i',  default=False, action="store_true", help='i band')
        parser.add_argument('-z',  default=False, action="store_true", help='z band')
        parser.add_argument('-y',  default=False, action="store_true", help='y band')
                
        return(parser)
    
    # read in command line arguments to determine desired grizy bands; if none are passed, return data for all of them
    def get_filter_list(self):
        if self.options.g: self.filterlist.append('g')
        if self.options.r: self.filterlist.append('r')
        if self.options.i: self.filterlist.append('i')
        if self.options.z: self.filterlist.append('z')
        if self.options.y: self.filterlist.append('y')

        if self.filterlist == []:        
            self.filterlist = ['g','r','i','z','y']

        if self.options.verbose >= 1:
            print("Filter list: ", self.filterlist)

    # add a column 'star' to the dataframe, which contains a 1 for every star and a 0 for every galaxy
    def mark_stars(self, Qf_lim=0.85,psfkron_diff=0.05):
        # first find all entries that have a QfPerfect for the given filterlist
        ixs_good_QfPerfect = self.catalog.getindices()
        for filt in self.filterlist:
            ixs_good_QfPerfect = self.catalog.ix_inrange(f'{filt}QfPerfect',Qf_lim,None,indices=ixs_good_QfPerfect)

        # cut on PSF-Kron mag, usiong only the entries that also passed the QfPerfect cut
        self.catalog.t['rPSF_rKron'] = self.catalog.t['rMeanPSFMag'] - self.catalog.t['rMeanKronMag']
        ixs_stars = self.catalog.ix_inrange('rPSF_rKron',None,psfkron_diff,indices=ixs_good_QfPerfect)

        # set 'star'=1 for stars, 0 for the rest
        self.catalog.t['star']=0
        self.catalog.t.loc[ixs_stars,'star']=1
        ixs_galaxies = AnotB(self.catalog.getindices(),ixs_stars)

        return(ixs_stars,ixs_galaxies)
    
    """
    establish 3 types of output, with varying amounts of information.
    skinny is the bare minimum (coordinates, magnitude, and magnitude error)
    medium includes object names and more information on both coordinates and magnitude 
    full includes every column the API outputs 
    """
    def get_outputcols(self):
        if self.options.outputtype == 'skinny':
            if self.options.sexagesimal:
                self.outputcols = ['raMean','decMean','raMeanSex','decMeanSex']
            else:
                self.outputcols = ['raMean', 'decMean']
            for f in self.filterlist:
                filtcols = [f'{f}MeanPSFMag',f'{f}MeanPSFMagErr']
                for i in filtcols:
                    self.outputcols.append(i)

        if self.options.outputtype == 'medium':
            if self.options.sexagesimal:
                self.outputcols = ['objName','raMean','decMean','raMeanErr','decMeanErr','raMeanSex','decMeanSex','epochMean',]
            else: 
                self.outputcols = ['objName','raMean','decMean','raMeanErr','decMeanErr','epochMean',]
            for f in self.filterlist:
                filtcols = [f'{f}MeanPSFMag',f'{f}MeanPSFMagErr',f'{f}MeanPSFMagStd',f'{f}MeanPSFMagNpt',f'{f}MeanPSFMagMin',f'{f}MeanPSFMagMax']
                for i in filtcols:
                    self.outputcols.append(i)

        if self.options.outputtype == 'full':  
            self.outputcols = self.allcols[:self.allcols.index(f'{self.filterlist[0]}MeanPSFMag')]

            if self.options.sexagesimal:
                firstcoord = self.allcols.index('raStack')
                lastcoord = self.allcols.index('epochMean')
                sexlist = ['raStack','decStack','raStackErr','decStackErr','raMean','decMean','raMeanErr','decMeanErr','raMeanSex','decMeanSex']
                self.outputcols[firstcoord:lastcoord] = sexlist

            for f in self.filterlist:
                filtcols = [f'{f}MeanPSFMag',f'{f}MeanPSFMagErr',f'{f}MeanPSFMagStd',f'{f}MeanPSFMagNpt',f'{f}MeanPSFMagMin',f'{f}MeanPSFMagMax']
                for i in filtcols:
                    self.outputcols.append(i)

            self.outputcols.append('distance')

        self.outputcols.append('star')

    """
    sort a passed list of lists (maglimlist) alphabetically according to grizy, using the first value in each sublist.
    require that there are only grizy bands in the final, sorted list
    """
    def custom_sort(self, key):
        grizy = ['g','r','i','z','y']

        if key[0] not in grizy:
            raise ValueError("Filter must be one of " + str(grizy) + ". You passed " + str(key[0]))
        else:
            return grizy.index(key[0])

    # return a list of all indices where the magnitude in a given band is outside of the provided bounds for that band 
    def magnitude_cut(self):
        cut_ixs = []
        if self.options.maglim is not None:
            sorted_maglims = sorted(self.options.maglim, key=self.custom_sort)
            for i in sorted_maglims:
                colname = f'{i[0]}MeanPSFMag'
                if sorted_maglims.index(i) == 0:
                    cut_ixs = self.catalog.ix_outrange(colname, lowlim=float(i[1]), uplim=float(i[2]))
                else:
                    current_band_cut = self.catalog.ix_outrange(colname, lowlim=float(i[1]), uplim=float(i[2]))
                    cut_ixs = AorB(cut_ixs, current_band_cut)
        return cut_ixs

    # Query the API for the provided coordinates and radius. Adapted from code by Dr. Ridden-Harper.
    def query_ps1(self, ra, dec, radius, 
                  version='dr2'):

        # Determine which version of the API to use. Default is dr2
        # if (version.lower() != 'dr2') & (version.lower() != 'dr1'):
        #     m = 'Version must be dr2, or dr1'
        #     raise ValueError(m)
        
        # query the API and read the response into a dataframe: self.catalog.t
        str = f'https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/{version.lower()}/mean?ra={ra}&dec={dec}&radius={radius}&nDetections.gte=1&pagesize=-1&format=csv'
        try:
            self.catalog.t = pd.read_csv(str)
        except pd.errors.EmptyDataError:
            print('No detections')
            self.catalog.t = []
            return(1)
        
        self.allcols = self.catalog.t.columns.tolist()

        return(0)

    # call query_ps1(). finalize the outputted file. perform any requested cuts based on magitude or type of object. convert coordinates to sexagesimal if desired.  
    def get_catalog(self, ra, dec, radius, 
                    version='dr2', 
                    onlystars=False, onlygalaxies=False, 
                    Qf_lim=0.85,psfkron_diff=0.05,
                    outputfile=None):
        
        self.query_ps1(ra, dec, radius, 
                       version=version)

        # determine which rows contain stars and which contain galaxies. if onlystars or onlygalaxies, discard all rows not containing the selected body, 
        ixs = self.catalog.getindices()
        (ixs_stars,ixs_galaxies) = self.mark_stars(Qf_lim=Qf_lim,psfkron_diff=psfkron_diff)
        if onlystars: 
            ixs = ixs_stars
        elif onlygalaxies: 
            ixs = ixs_galaxies

        # determine which type of output (skinny, medium, or full) to return
        self.get_outputcols()
        
        # discard all rows containing magnitudes outside of the passed magnitude ranges 
        cut_ixs = self.magnitude_cut()
        all_ixs = AnotB(ixs, cut_ixs)

        # add columns to display coordinates in sexagesimal format as well as decimal
        if self.options.sexagesimal:
            self.catalog.assert_radec_cols_sexagesimal('raMean', 'decMean','raMeanSex', 'decMeanSex')
        
        # print the output (or part of it) to the terminal, along with saving it to a file
        if self.options.debug:            
            self.catalog.write(columns=['raMean', 'decMean','rMeanPSFMag','rMeanPSFMagErr','star'],indices=all_ixs)
        elif self.options.verbose>2:
            print('The first 10 entries of the catalog')
            self.catalog.write(columns=self.outputcols,indices=all_ixs[:10])
            
        # save the catalog to the passed file, or indicate that no filename was passed 
        if outputfile is None:
            print('WARNING: no output filename specified, not saving the catalog!')
        else:
            self.catalog.write(outputfile,columns=self.outputcols,indices=all_ixs,verbose=self.verbose)

if __name__=='__main__':
    query=PS1_query_class()
    parser = query.add_options(usage='getPS1cat.py RA Dec radius -o [filename] -[desired grizy bands]')

    query.options = parser.parse_args()
    
    # determine which filters are wanted in the output file
    query.get_filter_list()

    query.get_catalog(query.options.RA,query.options.Dec,query.options.radius,
                      version=query.options.version,
                      onlystars=query.options.onlystars,
                      onlygalaxies=query.options.onlygalaxies,
                      outputfile=query.options.filename)
    
    print('Success! Output file: %s' % os.path.abspath(os.getcwd()) + "/" + query.options.filename)


    
              
