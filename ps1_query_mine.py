"""
script to query the Pan-STARRS API for a given RA, Dec, and radius, with options to filter the output 
for only stars or only galaxies. adapted from code by Drs. Ryan Ridden-Harper and Edward Schlafly.
"""

import pandas as pd
import sys, os,re,types,string,math,random
import argparse
from pdastro import pdastroclass, AorB, AnotB, AandB

class PS1_query_class:
    def __init__(self):
        self.ra = None
        self.dec = None
        self.radius = None

        self.filename = None

        self.filterlist=[]

        self.catalog = pdastroclass() # self.catalog.t to access the catalog dataframe

    def add_options(self, parser=None, usage=None):
        if parser == None:
            parser = argparse.ArgumentParser(usage=usage, conflict_handler="resolve")

        # fix defaults to remove elif stuff in main
        parser.add_argument('RA', help='right ascention coordinate')
        parser.add_argument('Dec', help='declination coordinate')
        parser.add_argument('radius', help='the radius in degrees around those coordinates to retireve data for')
        parser.add_argument('-os', '--onlystars', default=False, action="store_true", help='filter out all non-star objects')
        parser.add_argument('-og', '--onlygalaxies', default=False, action="store_true", help='filter out all non-galaxy objects')
        parser.add_argument('-v', '--verbose', action="count", dest="verbose",default=0)
        parser.add_argument('-d','--debug', default=False, action="store_true", help='Debugging output')
        parser.add_argument('-ps', '--pagesize', default=1, type=int, help='positve integer number of results to return. -1 to return all data')
        parser.add_argument('-f', '--filename', type=str, help='file name of output file')
        parser.add_argument('--version', type=str, default='dr2', help='the api version to use. default is dr2')
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

        if self.options.verbose:
            print("Filter list: ", self.filterlist)

    # Determine which indices contain stars and which contain galaxies, and filter for only 
    # one type of course if the boolean arguments `only_stars` or `only_galaxies` are set to 
    # true. Create column(s) called star or galaxy which contain a 1 if a row is a star or 
    # galaxy, respectively, and a 0 if not. Label which a row is even if both booleans are False.
    # Adapted from code by Dr. Ridden-Harper.
    def isolate_stars_galaxies(self, cat,Qf_lim=0.85,psfkron_diff=0.05,only_stars=False, only_galaxies=False):
        qf_ind = ((cat.t.gQfPerfect.values > Qf_lim) & (cat.t.rQfPerfect.values > Qf_lim) & (cat.t.iQfPerfect.values > Qf_lim) & (cat.t.zQfPerfect.values > Qf_lim))
        kron_ind = (cat.t.rMeanPSFMag.values - cat.t.rMeanKronMag.values) < psfkron_diff
        ind = qf_ind & kron_ind # indices which contain stars
        
        if only_stars:
            if self.options.verbose: 
                print("Discarding all non-star rows...")
            cat.t = cat.t.iloc[ind]
            cat.t.loc[:,'star'] = 1
        elif only_galaxies:
            if self.options.verbose: 
                print("Discarding all non-galaxy rows...")
            cat.t = cat.t.iloc[AnotB(cat.t.index.to_numpy(), ind)]
            cat.t.loc[:,'galaxy'] = 1
        else:
            cat.t.loc[:,'star'] = 0
            cat.t.loc[ind,'star'] = 1

            cat.t.loc[:,'galaxy'] = 1
            cat.t.loc[ind,'galaxy'] = 0
        return cat.t 
    
    # Query the API for the provided coordinates and radius. Call isolate_stars_galaxies(). Adapted from code by Dr. Ridden-Harper.
    def query_ps1(self, ra,dec,radius,only_stars=False, only_galaxies=False, version='dr2', page_size=1, outputfile=None):
        if outputfile is None:
            raise RuntimeError('No output file name provided')

        self.get_filter_list() # determine which columns are wanted in the output file

        # Determine which version of the API to use. Default is dr2
        if (version.lower() != 'dr2') & (version.lower() != 'dr1'):
            m = 'Version must be dr2, or dr1'
            raise ValueError(m)
        
        # query the API and read the response into a dataframe: self.catalog.t
        str = f'https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/{version.lower()}/mean?ra={ra}&dec={dec}&radius={radius}&nDetections.gte=1&pagesize={page_size}&format=csv'
        try:
            self.catalog.t = pd.read_csv(str)
        except pd.errors.EmptyDataError:
            print('No detections')
            self.catalog.t = []

        self.catalog.t = self.isolate_stars_galaxies(self.catalog,only_stars=only_stars, only_galaxies=only_galaxies)

        # Drop excess columns
        keepcols = ['objName', 'raMean', 'decMean', 'raMeanErr', 'decMeanErr', 'epochMean',]
        # Add to the list of kept columns the PSF magnitude values for any grizy bands passed through argparse
        for f in self.filterlist:
            filtcols = [f'{f}MeanPSFMag',f'{f}MeanPSFMagErr',f'{f}MeanPSFMagStd',f'{f}MeanPSFMagNpt',f'{f}MeanPSFMagMin',f'{f}MeanPSFMagMax']
            for i in filtcols:
                keepcols.append(i)
        # keep the star and/or galaxy columns as needed
        if(only_stars): 
            keepcols.append('star')
        elif(only_galaxies): 
            keepcols.append('galaxy')
        else: 
            keepcols.append('star')
            keepcols.append('galaxy')
        cutcat = self.catalog.t[keepcols] # copy the kept columns into a new dataframe 
        
        # save the dataframe to the passed output file; drop index numbers
        cutcat.to_csv(outputfile, index=False) 

if __name__=='__main__':
    query=PS1_query_class()
    parser = query.add_options(usage='getPS1cat.py RA Dec radius -o [filename] -[desired grizy bands]')

    query.options = parser.parse_args()

    query.query_ps1(query.options.RA, query.options.Dec, query.options.radius, only_stars=query.options.onlystars, only_galaxies=query.options.onlygalaxies, 
                    version=query.options.version, page_size=query.options.pagesize, outputfile=query.options.filename)