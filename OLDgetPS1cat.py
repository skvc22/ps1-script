#!/usr/bin/env python

from astropy.io import fits
import sys, os,re,types,string,math,random
import optparse
import numpy as np
import requests
# from PS1toINST import PS1toINSTclass
def makepath(path,raiseError=1):
    if path == '':
        return(0)
    if not os.path.isdir(path):
        os.makedirs(path)
        if not os.path.isdir(path):
            if raiseError == 1:
                raise RuntimeError('ERROR: Cannot create directory %s' % path)
            else:
                return(1)

    return(0)


def makepath4file(filename,raiseError=1):
    path = os.path.dirname(filename)
    if not os.path.isdir(path):
        return(makepath(path,raiseError=raiseError))
    else:
        return(0)

def deg2sex(degrees, ra=False, outputformatRA='%02d:%02d:%06.3f',outputformatDEC='%1s%02d:%02d:%05.2f'):
    if type(degrees) is  bytes:
        degrees=float(degrees)
    if ra:
        # a.k.a. indeg and outhours
        if degrees < 0.0:
            while degrees<0:degrees+=360.
        if degrees > 360.0:
            while degrees>360.0:degrees-=360.
        degrees /= 15.

    if degrees < 0:
        sign = '-'
    else:
        sign = '+'

    degrees = abs(degrees)

    d1  = (degrees - (degrees % 1))
    rem = (degrees % 1) * 60
    d2  = rem - (rem % 1)
    srem = (rem % 1) * 60
    d3 = srem

    if ra:
      return outputformatRA % (d1, d2, d3)
    else:
      return outputformatDEC % (sign, d1, d2, d3)

### Converts sexigesimal notation to decimal degrees or decimal hours (if option 'ra=True' invoked)
def sex2deg(sexigecimal, ra=False):
    ### 2005/12/02 - AR: make sure it is in sexagesimal format!
    # is it a string? if not check if it is None
    if not (type(sexigecimal) is bytes):
        if type(sexigecimal) == None:
            raise RuntimeError("ERROR: sex2deg cannot handle 'None'!")
        return sexigecimal
    # Does it have ':' or ' '? If not, it must be a float in string format, just convert it and return
    if re.search('\:|\s',sexigecimal) == None:
        return(string.atof(sexigecimal))

    s1, s2, s3 = list(map(string.atof, re.split('[DHMSdhms:\s]', sexigecimal.strip())[0:3]))

    # Get the sign
    if re.search('-', sexigecimal):
        sign = -1
    else:
        sign = 1

    deg = abs(s1) + s2 / 60. + s3 / 3600.

    deg *= sign

    if ra:
        deg *= 15.

    return deg

# Returns the passed in RA in decimal degrees
# input RA can be in 'HH:MM:SS.ss', 'HH MM SS.ss' or in decimal degrees 
def RaInDeg(Ra):
    import types
    if type(Ra)==bytes:
        if re.search('[DHMShms: ]',Ra.strip()):
            return(sex2deg(Ra,ra=True))
    return(float(Ra))

# Returns the passed in Dec in decimal degrees
# input Dec can be in 'DD:MM:SS.ss', 'DD MM SS.ss' or in decimal degrees 
def DecInDeg(Dec):
    import types
    if type(Dec)==bytes:
        if re.search('[DHMShms: ]',Dec.strip()):
            return(sex2deg(Dec,ra=False))
    return(float(Dec))

def rmfile(filename,raiseError=1,gzip=False):
    " if file exists, remove it "
    if os.path.lexists(filename):
        os.remove(filename)
        if os.path.isfile(filename):
            if raiseError == 1:
                raise RuntimeError('ERROR: Cannot remove %s' % filename)
            else:
                return(1)
    if gzip and os.path.lexists(filename+'.gz'):
        os.remove(filename+'.gz')
        if os.path.isfile(filename+'.gz'):
            if raiseError == 1:
                raise RuntimeError('ERROR: Cannot remove %s' % filename+'.gz')
            else:
                return(2)
    return(0)

def rmfiles(filenames,raiseError=1,gzip=False):
    if not (type(filenames) is list):
        raise RuntimeError("List type expected as input to rmfiles")
    errorflag = 0
    for filename in filenames:
        errorflag |= rmfile(filename,raiseError=raiseError,gzip=gzip)
    return(errorflag)

class getPS1catclass:
    def __init__(self):
        self.ra = None
        self.dec = None

        self.alldata={}

        self.filename = None

        self.filterlist=[]
        self.colorfilterlist=[]
        self.errorcolsFlag= True
        self.Nonestring = 'nan'

        # Eddies columns are named 'median(0)' and 'err(0)' for g, etc
        self.filt2index = {'g':0,'r':1,'i':2,'z':3,'y':4}
        self.supercal={'g':0.020,'r':0.033,'i':0.024,'z':0.028,'y':0.011}

        self.inst=''


    def autooutfilename(self):
        if self.ra==None or self.dec==None:
            raise RuntimeError('No RA/Dec defined to be used for auto filename!')

        filename = '%011.7f_%011.7f.PS1.' % (self.ra,self.dec)
        filename += ''.join(self.filterlist)
        filename += '.cat'
        
        return(filename)

    def add_options(self, parser=None, usage=None):
        if parser == None:
            parser = optparse.OptionParser(usage=usage, conflict_handler="resolve")

        parser.add_option('-v', '--verbose', action="count", dest="verbose",default=0)
        parser.add_option('-d','--debug', default=False, action="store_true",
                          help='Debugging output')
        parser.add_option('--cathtml'  , default='http://faun.rc.fas.harvard.edu/eschlafly/wise/fitscat_qy', type="string",
                          help='root html for Eddies cats (default=%default)')
        parser.add_option('--tmpdir'  , default='/tmp', type="string",
                          help='Directory for temporary files (default=%default)')
        parser.add_option('--cat_dir', default='/projects/caps/uiucsn/PS1SKYMAPPER_fitscat_qy/fitscat_qy',
                          help="Gautham's local mirror of Eddie's cats (default=%default)")
        parser.add_option('-B',  default=False, action="store_true",
                          help='B band')
        parser.add_option('-V',  default=False, action="store_true",
                          help='V band')
        parser.add_option('-R',  default=False, action="store_true",
                          help='R band')
        parser.add_option('-I',  default=False, action="store_true",
                          help='I band')
        parser.add_option('-u',  default=False, action="store_true",
                          help='u band')
        parser.add_option('-g',  default=False, action="store_true",
                          help='g band')
        parser.add_option('-r',  default=False, action="store_true",
                          help='r band')
        parser.add_option('-i',  default=False, action="store_true",
                          help='i band')
        parser.add_option('-z',  default=False, action="store_true",
                          help='z band')
        parser.add_option('-y',  default=False, action="store_true",
                          help='y band')
        parser.add_option('--sexagesimal',  default=False, action="store_true",
                          help='output RA/Dec are in sexagesimal format')
        parser.add_option('-s','--show',  default=False, action="store_true",
                          help='print catalog to screen')
        parser.add_option('--keepnans',  default=False, action="store_true",
                          help='keep entries with NaN values')
        parser.add_option('--requiregriz',  default=False, action="store_true",
                          help='only keep objects that have griz measurements')
        parser.add_option('-o','--outfile'  , default=None , type="string",
                          help='file name of output file. If not specified, then automatic filename depending on RA,Dec and filters (default=%default)')
        parser.add_option('-c','--clobber',  default=False, action="store_true",
                          help='overwrite file if it already exists')
        parser.add_option('-s','--size'  , default=None , type="string",
                          help='sidelength of box in degree. If format AxB, then rectangle with sidelength A,B along Ra, Dec, respectively(default=%default)')
        parser.add_option('--skipsave',  default=False, action="store_true",
                          help='Don\'t save the catalog')
        parser.add_option('--Mmax'  , default=None , type="float",
                          help='cut all objects with M<Mmin for specified filters (default=%default)')
        parser.add_option('--Mmin'  , default=None , type="float",
                          help='cut all objects with M>Mmax for specified filters (default=%default)')
        parser.add_option('--dMmax'  , default=None , type="float",
                          help='cut all objects with dM<dMmin for specified filters (default=%default)')
        parser.add_option('--dMmin'  , default=None , type="float",
                          help='cut all objects with dM>dMmax for specified filters (default=%default)')
        parser.add_option('--decam'  , default=False , action="store_true",
                          help='convert the PS1 mag to DECAM(DES)')
        parser.add_option('--megacam'  , default=False , action="store_true",
                          help='convert the PS1 mag to MEGACAM(SNLS)')
        parser.add_option('--swope'  , default=False , action="store_true",
                          help='convert the PS1 mag to SWOPE(CSP)')
        parser.add_option('--transfdir'  , default=False , type="string",
                          help='path to the PS1transformations dir (default=%default)')
        parser.add_option('--skip_sc'  , default=False ,action="store_true",
                          help='Skip the supercal correction to the PS1 mag (default=%default)')
        parser.add_option('--spline'  , default=False ,action="store_true",
                          help='Use spline fit (default=linear)')
        

        return(parser)

    def check4TRANSFdir(self):

        if self.options.transfdir: 
            print("check4 1st if true")
            self.path2name=self.options.transfdir
        else:
            if 'PIPE_PS1TRANSFDIR' in os.environ: 
                print("thing is in wnviron")
                if os.path.isdir(os.getenv('PIPE_PS1TRANSFDIR')): self.path2name=os.getenv('PIPE_PS1TRANSFDIR')+'/PS1trans_linear.txt'
                else:
                    print('WARNING: %s does not exist, looking into %s' %(os.getenv('PIPE_PS1TRANSFDIR'),os.getcwd()+'/PS1transformations'))
                    if os.path.isdir('./PS1transformations'): 
                        self.path2name='./PS1transformations/PS1trans_linear.txt' 
                    else: sys.exit('ERROR: Unable to find PS1transformations dir using the current path: %s. Try again with --transfdir option specifing the new path for the dir' % os.getcwd())
            else:
                print("thing not there")
                print('WARNING: setenv PIPE_PS1TRANSFDIR does not exist, looking into %s' %os.getcwd()+'/PS1transformations')
                if os.path.isdir('./PS1transformations'): 
                    print("its okay tho") # but it breaks when i go to argparse? i think? figure out tomorrow
                    self.path2name='./PS1transformations/PS1trans_linear.txt' 
                else: sys.exit('ERROR:  Unable to find PS1transformations dir using the current path: %s. Try again  with --transfdir option specifing the new path for the dir' % os.getcwd())


    def getfilterlist_from_options(self):

        if self.options.decam:  
            self.inst='DES0'
        if self.options.megacam:  
            self.inst='SNLS0'
            self.options.spline=True
            self.options.skip_sc=True
        if self.options.swope:
            self.inst='CSP1'

        if self.options.u and self.inst == None: 
            print('> PS1 has no u band. Use --decam or --decam option.')
            sys.exit(0)
       
        if not self.inst:
            if self.options.u: self.filterlist.append('u')
            if self.options.g: self.filterlist.append('g')
            if self.options.r: self.filterlist.append('r')
            if self.options.i: self.filterlist.append('i')
            if self.options.z: self.filterlist.append('z')
            if self.options.y: self.filterlist.append('y')
 

        if self.filterlist == []:
           
            self.filterlist = ['g','r','i','z','y']
            
       

    def add2list(self,alldata,col1,newdata,col2):
        Nrowsnew = len(newdata[col2])
        if col1 in alldata:
            Nrows1 = len(alldata[col1])
            if Nrowsnew>0:
                tmp=np.zeros((Nrows1+Nrowsnew))
                tmp[:Nrows1]=alldata[col1]
                tmp[Nrows1:]=newdata[col2]
                alldata[col1]=tmp
        else:
            if Nrowsnew>0:
                alldata[col1]=newdata[col2]
        return(0)




    def fluxcolname(self,filt):
        return('median(%d)' % self.filt2index[filt])

    def errcolname(self,filt):
        return('err(%d)' % self.filt2index[filt])    
 
    def getdataforRADEC(self,ra,dec,alldata=None,urldonelist=None,ramin=None,ramax=None,decmin=None,decmax=None,tmpfiledir='.'):
        tmpfiledir = re.sub('\/$','',os.path.abspath(os.path.expanduser(tmpfiledir)))

        if alldata==None:
            alldata=self.alldata

        if ra<0.0:ra+=360.0
        if ra>=360.0:ra-=360.0

        if ramin!=None:
            if (ra-ramin)<-180:
                ramin-=360.0
                ramax-=360.0
            elif (ra-ramin)>180:
                ramin+=360.0
                ramax+=360.0


        # root http
        http_addr=re.sub('\/$','',self.options.cathtml)
        local_dir = self.options.cat_dir
        # just in case that ra,dec are integers, add some wiggle room...
        epsil = 0.000001
        # name of the catalog
        web_name='cat.l='+repr(int(math.floor(ra)))+'_'+repr(int(math.ceil(ra+epsil)))+'.b='+repr(int(math.floor(dec)))+'_'+repr(int(math.ceil(dec+epsil)))+'.fits'

        if self.options.verbose>1:
            if ramin!=None:
                print('RA=%.6f Dec=%.6f (ramin=%.6f ramax=%.6f): %s' % (ra,dec,ramin,ramax,web_name))
            else:
                print('RA=%.6f Dec=%.6f: %s' % (ra,dec,web_name))
                
        # put it all together into a url 
        url = '%s/%s.gz' % (http_addr,web_name) # keep this in here for now just in case
        
        ### If we are on the Illinois Campus cluster, just use the local mirror
        ### Else, just use the remote one
        if os.path.exists(local_dir):
            template_file = '%s/%s.gz' % (local_dir,web_name)
            tmpfile = template_file
        else:
            # Did we already add this url to alldata?
            if (type(urldonelist) is list) and (url in urldonelist):
                if self.options.verbose:
                    print('%s already downloaded and added to cat, skipping...' % (web_name+'.gz'))
                return(0,url)

            
            tmpdir = re.sub('\/$','',self.options.tmpdir)
            if self.options.debug:
                tmpfile = '%s/%s.delme.gz' % (tmpfiledir,web_name) 
            else:
                tmpfile = '%s/%s.%d.delme.gz' % (tmpfiledir,web_name,random.randint(0,1000000)) 
            
            if self.options.verbose>1:
                print('gz file temporarely saved into %s' % tmpfile)

            makepath4file(tmpfile)

            if self.options.debug and os.path.isfile(tmpfile):
                print('DEBUG: skipping downloading %s into %s, it already exists!' % (url,tmpfile))
            else:
                rmfile(tmpfile)
                makepath4file(tmpfile)
                print('Downloading ',url)
                r = requests.get(url)
                print('VVVVVVVVV',r)
                if r.status_code!=200:
                    print('status message:',r.text)
                    raise RuntimeError('ERROR: could not get to url %s, status code %d' % (url,r.status_code))
                with open(tmpfile, 'wb') as f:
                    for chunk in r.iter_content(8096):
                        f.write(chunk)
                    f.close()
            

        hdulist = fits.open(tmpfile) #open file
        tbdata = hdulist[1].data

        # calculate mags
        if self.options.verbose>1:
            print('calculating mags...')

        if ramin!=None:
            mask = ((tbdata['ra']<ramax) & (tbdata['ra']>ramin) & (tbdata['dec']<decmax) & (tbdata['dec']>decmin))
            data2keep = tbdata[mask]
        else:
            data2keep = tbdata
        if self.options.requiregriz:
            mask =  ((data2keep[self.fluxcolname('g')] !=0) & (data2keep[self.fluxcolname('r')] !=0) & (data2keep[self.fluxcolname('i')] !=0) & (data2keep[self.fluxcolname('z')] !=0)) 
            data2keep = data2keep[mask]
            for elno in range(len(data2keep)): print((data2keep[elno][0],data2keep[elno][1],data2keep[elno][2],data2keep[elno][3],data2keep[elno][4],data2keep[elno][5]))
            
        # get the mags and dmags, save them temporarely in the same space than flux and fluxerror
        #print self.filterlist
       

        #color2keep=data2keep
        if not self.options.skip_sc: print('>>>>>> Adding supercal correction')  
        for filt in self.filterlist:
           
            data2keep[self.errcolname(filt)] =  1.06*data2keep[self.errcolname(filt)]/data2keep[self.fluxcolname(filt)] #convert errors from flux to info
            data2keep[self.fluxcolname(filt)] = -2.5*np.log10(data2keep[self.fluxcolname(filt)]) # convert flux to mag of observed values
            
            if not self.options.skip_sc:  data2keep[self.fluxcolname(filt)] = data2keep[self.fluxcolname(filt)] + self.supercal[filt] # add supercal correction
              
                
           
        # check if there are limits on the mags and uncertainties
        for filt in self.filterlist:
            if self.options.Mmin!=None:
                mask =  (data2keep[self.fluxcolname(filt)]>=self.options.Mmin)
                data2keep = data2keep[mask]
            if self.options.Mmax!=None:
                mask =  (data2keep[self.fluxcolname(filt)]<=self.options.Mmax)
                data2keep = data2keep[mask]
            if self.options.dMmin!=None:
                mask =  (data2keep[self.errcolname(filt)]>=self.options.dMmin)
                data2keep = data2keep[mask]
            if self.options.dMmax!=None:
                mask =  (data2keep[self.errcolname(filt)]<=self.options.dMmax)
                data2keep = data2keep[mask]
        
         # now add the data to self.alldata 
        self.add2list(alldata,'ra',data2keep,'ra')
        self.add2list(alldata,'dec',data2keep,'dec')

        for filt in self.filterlist:
            self.add2list(alldata,filt    ,data2keep,self.fluxcolname(filt))
            self.add2list(alldata,'d'+filt,data2keep,self.errcolname(filt))

 


        hdulist.close()

        # Just like above, if using the local mirror, don't delete the file
        # But do delete it if using the remote one

        #===========================================commented out bc it was deleting the files==============================================================

        # if os.path.exists('/projects/caps/uiucsn/PS1SKYMAPPER_fitscat_qy/fitscat_qy'):
        #     pass
        # else:
        #     if not self.options.debug:
        #         rmfile(tmpfile)

        del tbdata,data2keep#,color2keep
    
        return(0,url)

    def getcatalog(self, ra, dec,switch, tmpfiledir = '.'):
        # keep track of RA,Dec
        self.ra  = RaInDeg(ra)
        self.dec = DecInDeg(dec)
        if self.options.verbose>1:
            print('RA: %.7f, Dec:%.7f' % (self.ra,self.dec))
        RAboxsize = DECboxsize = None
        if self.options.size!=None:
            # get the boxsize, and ra/dec limits
            if re.search('\w+x\w+',self.options.size):
                m=re.search('(^\S+)x(\S+$)',self.options.size)
                if m!=None:
                    RAboxsize = float(m.groups()[0])
                    DECboxsize = float(m.groups()[1])
                    print('box sizes: ',RAboxsize,DECboxsize)
                else: 
                    raise RuntimeError('Could not parse size %s for' % self.options.size)
            else:
                RAboxsize = DECboxsize = float(self.options.size)

            # get the maximum 1.0/cos(DEC) term: used for RA cut
            minDec = self.dec-0.5*DECboxsize
            if minDec<=-90.0:minDec=-89.9
            maxDec = self.dec+0.5*DECboxsize
            if maxDec>=90.0:maxDec=89.9

            invcosdec = max(1.0/math.cos(self.dec*math.pi/180.0),
                            1.0/math.cos(minDec  *math.pi/180.0),
                            1.0/math.cos(maxDec  *math.pi/180.0))

            # get the (conservative) boxlimits
            ramin = self.ra-0.5*RAboxsize*invcosdec
            ramax = self.ra+0.5*RAboxsize*invcosdec
            decmin = self.dec-0.5*DECboxsize
            decmax = self.dec+0.5*DECboxsize
            # check for url for center and all 4 corners...
            radeclist =[(self.ra,self.dec),
                        (ramin,decmin), 
                        (ramin,decmax), 
                        (ramax,decmin), 
                        (ramax,decmax)] 
            
            
            if ramax-ramin > 1 and decmax-decmin > 1:
                radeclist=[]
                for n in range(int(ramax-ramin)+1):
                    app_ra=ramin+n
                    for l in range(int(decmax-decmin)+1):
                        app_dec=decmin+l
                        radeclist.extend([(app_ra,app_dec)])
                    
            elif ramax-ramin > 1 and decmax-decmin < 1:
                radeclist=[]
                for n in range(int(ramax-ramin)+1):
                    app_ra=ramin+n
                    for dec in [decmin,decmax]:
                        app_dec=dec
                        radeclist.extend([(app_ra,app_dec)])

            elif ramax-ramin < 1 and decmax-decmin > 1:
                radeclist=[]
                for ra in [ramin,ramax]:
                    app_ra=ra
                    for l in range(int(decmax-decmin)+1):
                        app_dec=decmin+l
                        radeclist.extend([(app_ra,app_dec)])


        else:
            ramin = ramax = decmin = decmax = None
            radeclist = [(self.ra,self.dec)]

        # Which filters?
        self.getfilterlist_from_options()
        if self.options.verbose>1:
            print('Filters:',self.filterlist)


        # list of urls for which the data is already parsed. so they don't have to be done again...
        urldonelist=[]
        print(radeclist)
        for (ra4url,dec4url) in radeclist:
            print('### checking %f %f'% (ra4url,dec4url))
            (errorflag,url) = self.getdataforRADEC(ra4url,dec4url,urldonelist=urldonelist,ramin=ramin,ramax=ramax,decmin=decmin,decmax=decmax,tmpfiledir=tmpfiledir)
            urldonelist.append(url)

        return(0)

    def get_output_format_info(self,list_of_filters):
        # get cols ...
        cols=['ra','dec']
        colsformat = ['%11.7f','%11.7f']
        header = '#%11s %11s' % ('ra','dec')
        for filt in list_of_filters:
            colsformat.append('%7.4f')
            cols.append(filt)
            header += ' %7s' % (filt)
            if self.errorcolsFlag:
                colsformat.append('%7.4f')
                cols.append('d'+filt)
                header += ' %7s' % ('d'+filt)

        return(cols,colsformat,header)

    def getformattedcatalog(self,cols,colsformat, alldata=None, indices = None, Nonestring=None, addreturnFlag=False):
        if alldata==None:
            alldata=self.alldata

        if indices == None:
            indices = range(len(alldata['ra']))

        if Nonestring==None:
            Nonestring = self.Nonestring

        cs = range(len(cols))
        
        lines=[]
        for i in indices:
            line = ''
            nanflag = False
            for c in cs:                    
                #print len(alldata[cols[c]]),i
                if math.isnan(alldata[cols[c]][i]) or math.isinf(alldata[cols[c]][i]):
                    nanflag=True
                    line += ' %7s' %  Nonestring
                else:
                    if cols[c] in ['ra','dec'] and self.options.sexagesimal:
                        line += ' '+deg2sex(alldata[cols[c]][i],ra=(cols[c]=='ra'))
                    else:
                        line += ' '+colsformat[c] % alldata[cols[c]][i] 
                    
            # skip bad lines if not --keepnans
            if nanflag and (not self.options.keepnans):
                continue

            if addreturnFlag:
                line += '\n'

            lines.append(line)

        return(lines)
            

    def savecatalog(self, outfilename, keys = None, cols = None, errorcolsFlag=True, Nonestring='NaN',clobber=False):
        makepath4file(outfilename)

        # file already exists?
        if os.path.isfile(outfilename):
            if clobber:
                if self.options.verbose>0:
                    print('file %s exists, clobbering it...' % outfilename)
                rmfile(outfilename)
            else:
                print('file %s already exists, stopping (if you want to overwrite file, use --clobber)' % outfilename)
                return(1)

        if self.options.verbose>1:
            print('Getting formatted catalog...')

        (cols,colsformat,header) = self.get_output_format_info(self.filterlist)

        (lines) = self.getformattedcatalog(cols,colsformat,Nonestring=Nonestring,addreturnFlag=True)

        if self.options.verbose:
            print('Saving %d entries into %s' % (len(lines),outfilename))
        f=open(outfilename,'w')
        f.writelines([header+'\n'])
        f.writelines(lines)
        f.close()

        # keep track of filename
        self.filename = outfilename

        if not os.path.isfile(outfilename):
            raise RuntimeError('Could not write to %s' % outfilename)

        return(0)

    def showcatalog(self, keys = None, cols = None, errorcolsFlag=True, Nonestring='NaN'):
        (cols,colsformat,header) = self.get_output_format_info()
        (lines) = self.getformattedcatalog(cols,colsformat,Nonestring=Nonestring)
        print(header)
        for line in lines: print(line)
        return(0)
        
    def converter(self,outfilename):
        #added by gio 01/04/2016 convertion from PS1 to DECAM
        #surv_name=self.inst
        selected_band=[]
        #print self.filterlist
        if self.options.B: 
            selected_band.append('B')
            #print selected_band
        elif self.options.V: 
            selected_band.append('V')
            #print selected_band
        if self.options.R:
            selected_band.append('R')
        if self.options.I:
            selected_band.append('I')
        elif self.options.u: 
            selected_band.append('u')
            #print selected_band
        elif self.options.g: 
            selected_band.append('g')
            selected_band.append('r')
            #print selected_band
        elif self.options.i: 
            selected_band.append('i')
            #print selected_band
        elif self.options.z: 
            selected_band.append('z')
            #print selected_band
        elif self.options.y: 
            selected_band.append('y')
            #print selected_band
        else:
            selected_band.append('u')
            selected_band.extend(self.filterlist)


        
        # PS1toINST=PS1toINSTclass()
        # PS1toINST.convert(self.path2name,self.inst,outfilename,selected_band,self.options.sexagesimal,self.options.spline)

if __name__=='__main__':

    getPS1cat=getPS1catclass()
    parser = getPS1cat.add_options(usage='getPS1cat.py RA Dec')

    
    if len(sys.argv)<3:
        options,  args = parser.parse_args(args=['-h'])        
        sys.exit(0)

    
    # this is a hack to allow negative numbers as arguments (The one bug in optparse). 
    # this means that all options must be after the arguments.
    #getPS1cat.options,  args = parser.parse_args()
    getPS1cat.options,  dummy = parser.parse_args(args=sys.argv[3:])
    args = sys.argv[1:3]
    print("options: ", getPS1cat.options)
    print("args: ", args)
    print("dummy: ", dummy)
    if len(dummy)!=0:
        print('ERROR: extra arguments!',dummy)
        sys.exit(0)

    (ra,dec)=args[:2]

    getPS1cat.check4TRANSFdir()

    # outfilename
    if getPS1cat.options.outfile!=None:
        outfilename = getPS1cat.options.outfile
    else:
        outfilename = getPS1cat.autooutfilename()

    # directory for temporary files
    tmpfiledir = os.path.dirname(os.path.abspath(os.path.expanduser(outfilename)))

    # get the catalog from the web
    getPS1cat.getcatalog(ra,dec,switch=0,tmpfiledir=tmpfiledir)    

    # show the catalog?
    if getPS1cat.options.show:
        getPS1cat.showcatalog()

    # save the catalog?
    if not getPS1cat.options.skipsave:
        getPS1cat.savecatalog(outfilename,clobber=getPS1cat.options.clobber)
        

    
        
    #move PS1cat to Surv?       
    if getPS1cat.options.decam or getPS1cat.options.megacam or getPS1cat.options.swope:
        getPS1cat.converter(outfilename)
        
        
                
    print('SUCCESS %s' % os.path.basename(sys.argv[0]))
