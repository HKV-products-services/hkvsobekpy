import os
import multiprocessing as mp
try:
    from pathlib import Path
except:    
    from pathlib2 import Path
import argparse
import sys
import itertools
import geopandas as gpd
from datetime import datetime, timedelta
import copy
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
import numpy as np
import pandas as pd
from scipy import optimize
import fire
import warnings
import configparser
from hkvsobekpy.core.utils import *
try:
    shell = get_ipython().__class__.__name__
    if shell == 'ZMQInteractiveShell':
        from tqdm import tqdm_notebook as tqdm # Jupyter notebook or qtconsole        
    elif shell == 'TerminalInteractiveShell':
        from tqdm import tqdm # Terminal running IPython
    else:
        from tqdm import tqdm # Other type (?)
except NameError:
    from tqdm import tqdm # Probably standard Python interpreter
    
#for debugging within notebook    

try:
    from IPython.core.debugger import set_trace    
except:
    pass

# for python 2/3 compatibility
try:
    xrange
except NameError:
    xrange = range    
try:
    get_ipython().magic('matplotlib inline')    
except:
    pass

class __his_class(object):
    class __errors__(object):
        """
        error class met verschillende foutmeldingen
        """
        @staticmethod
        def fileNotFound(path_file):
            raise IOError("File doesn't exist. Your path '{}' is OK?".format(path_file))
        @staticmethod
        def metadataNotSet():
            raise AttributeError('Metadata is not set. Use hkvsobekpy.ReadMetadata(his_file)')
        @staticmethod
        def metadataError():
            raise AttributeError('Metadata could not be set using function ReadMetadata()')
        @staticmethod
        def administratieError():            
            raise AttributeError('Administratieblok van his-file kan niet worden uitgelezen')
        @staticmethod
        def locationNotFound(loc):            
            raise AttributeError("Could not find location '{}'. Is it an existing location?".format(loc))
        @staticmethod
        def longLocationNotFound(error):
            raise AttributeError('There is an error reading the locations in the hia-file. Full error is:\n{}'.format(error))
        @staticmethod
        def longParameterNotFound(error):
            raise AttributeError('There is an error reading the parameters in the hia-file. Full error is:\n{}'.format(error))
        @staticmethod
        def unexpectedT0Error(timestampInfo):            
            raise AttributeError("T0 could not be read. Timestamp info is '{}'.".format(timestampInfo))            
        @staticmethod
        def parameterNotFound(param):            
            raise AttributeError("Could not find parameter '{}'. Is it an existing parameter?".format(param))
        @staticmethod
        def annualmaxError():
            raise ValueError('Voor de MultiWaardenArray functie kan alleen `year` of `none` gebruikt worden als voor de jaarmax_as parameter')
        @staticmethod
        def gewogenGemiddeldeError():
            raise ValueError('Kan gewogen gemiddelde niet bepalen. Vereiste is minimaal 2 punten aan de linker en rechterzijde van de gekozen T')
        def gebeurtenissenError():
            raise ValueError('Meer gebeurtenissen dan N.')
    
    def __init__(self):
        """
        SobekResultaten class. Class met methodes om een .his bestand 
        uit te lezen.
        WaterstandStatistiek class. Class met methodes om een GumbleFit te bepalen
        """
        # initatie van lege immutable objects
        # zie: https://stackoverflow.com/a/22526544/2459096
        # python2.7: obj = argparse.Namespace()
        # python3.x: obj = types.SimpleNamespace()
        self.hisFile = argparse.Namespace()
        self.hisFile.timestampInfo = argparse.Namespace()
        self.hisFile.parameterInfo = argparse.Namespace()
        self.hisFile.locationInfo = argparse.Namespace()       
        self.hiaFile = argparse.Namespace()
        self.hiaFile.parameterInfo = argparse.Namespace()
        self.hiaFile.locationInfo = argparse.Namespace()       
                
    def _leesHeader(self, f):
        """
        Lees de header van het bestand

        Parameters
        ----------
        f : BufferedReader    
            Leesbaar ruw IO object

        Returns
        -------
        header : array
            De offset van de location in bytes    
        """
        header = []
        for i in range(4):    
            header.append(f.read(40).decode('windows-1252'))  
        return header           
        
    def _infoFromHeader(self,timestampInfo):
        """
        Lees informatie uit de header

        Parameters
        ----------
        timestampInfo : str
            De location in string format

        Returns
        -------
        beginDate : datetime 
            De begindatum T0, waarop de andere datetime gebasseerd zijn
        timeStepInterval : str
            De eenheid waarin de stapgrote van de datetime objecten berekend zijn, 
            dit kan zijn s [seconds], m [minutes], h[hours], d[days]
        timeStepFactor : int
            Factor welke op de stapgrote toegepast moet worden
        """

        if len(timestampInfo) == 40:
            year = timestampInfo[4:8]
            month = timestampInfo[9:11]
            day = timestampInfo[12:14]
            hours = timestampInfo[15:17]
            minutes = timestampInfo[18:20]
            seconds = timestampInfo[21:23]

            beginDate = datetime(int(year), int(month), int(day), int(hours), int(minutes), int(seconds))

            timeStepInterval = timestampInfo[-2:-1]
            timeStepFactor = int(timestampInfo[30:-2])
        else:
            self.__errors__.unexpectedT0Error(timestampInfo)

        return beginDate, timeStepInterval, timeStepFactor    

    def _leesHia(self, f):
        """
        Lees hia bestand
        """
                
        try:
            p = f.resolve()
            self.hiaFile.hiaFileName = str(p)            

            config = configparser.ConfigParser()
            config.read(self.hiaFile.hiaFileName)
            try:
                # parse hia locations to list of dict
                locations = [{'index':int(val)-1, 'long name':config['Long Locations'][val]} for val in config['Long Locations']]  
                self.hiaFile.locationInfo.locations = locations
            except KeyError as e:
                self.__errors__.longLocationNotFound()
            
            try:
                # parse hia parameters to list of dict
                parameters = [{'index':int(val)-1, 'long name':config['Long Parameters'][val]} for val in config['Long Parameters']]
                self.hiaFile.parameterInfo.parameters = parameters
            except KeyError as e:
                self.__errors__.longParametersNotFound()

        except:
            # hia doesn't exist return empty objects
            self.hiaFile.locationInfo.locations = []
            self.hiaFile.parameterInfo.parameters = []
    
    def _leesAdmin(self, f):
        """
        Lees het administratie blok
        """
        self.hisFile.header = self._leesHeader(f)

        # lees info van header
        self.hisFile.timestampInfo.beginDate, self.hisFile.timestampInfo.timeStepInterval, self.hisFile.timestampInfo.timeStepFactor = self._infoFromHeader(self.hisFile.header[3])        

        # lees aantal parameters
        self.hisFile.parameterInfo.numPar = np.fromfile(f,np.int,1)[0]    

        # lees aantal locations
        self.hisFile.locationInfo.numLoc = np.fromfile(f,np.int,1)[0]    

        # lambda function to split string
        split_string = lambda x, n: [x[i:i+n] for i in range(0, len(x), n)]
        
        # lees parameters
        self.hisFile.parameterInfo.parameters = []
        param_string = f.read(self.hisFile.parameterInfo.numPar * 20).decode('windows-1252')
        self.hisFile.parameterInfo.parameters = split_string(param_string, 20)

#         for i in range(self.hisFile.parameterInfo.numPar):
#             self.hisFile.parameterInfo.parameters.append(f.read(20).decode('windows-1252'))    

        # Lees locations
        self.hisFile.locationInfo.id = []
        self.hisFile.locationInfo.locations = []
        
        #self.hisFile.locationInfo.locations = np.fromfile(f, np.dtype("i2, V20" ), 2)#, U20"), 1)
        for i in range(self.hisFile.locationInfo.numLoc):            
            self.hisFile.locationInfo.id.append(np.fromfile(f,np.int,1)[0])
            self.hisFile.locationInfo.locations.append(f.read(20).decode('windows-1252').rstrip())
            

        # lees timestamps
        self.hisFile.headerSize = f.tell()
        self.hisFile.locationInfo.numTime = int((self.hisFile.hisFileSize - self.hisFile.headerSize) /
                      (self.hisFile.locationInfo.numLoc * self.hisFile.parameterInfo.numPar + 1) / 4)

        aantalBytesPerStap = int((self.hisFile.hisFileSize - self.hisFile.headerSize) / self.hisFile.locationInfo.numTime)

        byteNr = self.hisFile.headerSize

        self.hisFile.timestampInfo.moments = []
        
        self.hisFile.timestampInfo.offset = []
        for i in range(self.hisFile.locationInfo.numTime):
            f.seek(byteNr)
            moment = np.fromfile(f,np.int,1)[0]

            if self.hisFile.timestampInfo.timeStepInterval == 's':
                self.hisFile.timestampInfo.moments.append(self.hisFile.timestampInfo.beginDate + timedelta(seconds = (float(moment) * self.hisFile.timestampInfo.timeStepFactor)))
            elif self.hisFile.timestampInfo.timeStepInterval == 'm':
                self.hisFile.timestampInfo.moments.append(self.hisFile.timestampInfo.beginDate + timedelta(minutes = (float(moment) * self.hisFile.timestampInfo.timeStepFactor)))
            elif self.hisFile.timestampInfo.timeStepInterval == 'h':
                self.hisFile.timestampInfo.moments.append(self.hisFile.timestampInfo.beginDate + timedelta(hours= (float(moment) * self.hisFile.timestampInfo.timeStepFactor)))
            elif self.hisFile.timestampInfo.timeStepInterval == 'd':
                self.hisFile.timestampInfo.moments.append(self.hisFile.timestampInfo.beginDate + timedelta(days= (float(moment) * self.hisFile.timestampInfo.timeStepFactor)))
            self.hisFile.timestampInfo.offset.append(f.tell())
            byteNr += aantalBytesPerStap
        self.hisFile.timestampInfo.offset = np.array(self.hisFile.timestampInfo.offset, dtype=np.int64)
        self.hisFile.timestampInfo.N = np.array(self.hisFile.timestampInfo.moments).max().year - np.array(self.hisFile.timestampInfo.moments).min().year + 1
        return True

    def _locOffset(self, loc, param):
        """
        Krijg de offset van de parameter

        Parameters
        ----------
        loc : string    
            De location in kwestie
        param : string
            De parameter in kwestie        

        Returns
        -------
        locOffset : int
            De offset van de location in bytes
        """
        #  Zoek de index van de parameter
        try:
            parFound = self.hisFile.parameterInfo.parameters.index(param)
        except :
            self.__errors__.parameterNotFound(param)            

        # Zoek de index/id van de location
        try:
            # removed dictionary
            # locFound = next((item for item in self.hisFile.locationInfo.locations if item['location'] == loc))
            locFound = self.hisFile.locationInfo.locations.index(loc)
        except :
            self.__errors__.locationNotFound(loc)

        locOffset = locFound * self.hisFile.parameterInfo.numPar * 4 + parFound * 4

        return locOffset   
    
    def _mergeHiaHisLocations(self):
        """
        Merge locations his en hia
        """
        df_his_locs = pd.DataFrame(self.hisFile.locationInfo.locations, columns=['his name'])
        
        df_hia_locs = pd.DataFrame(self.hiaFile.locationInfo.locations)
        # only merge if df_hia_locs contains items
        if not df_hia_locs.empty:
            df_hia_locs.set_index('index', inplace=True)
        
            df_locs = pd.merge(df_his_locs, df_hia_locs, left_index=True, right_index=True, how='outer')
            df_locs.loc[df_locs.isnull().any(axis=1), 'long name'] = df_locs[df_locs.isnull().any(axis=1)]['his name']
            
            self.hisFile.locationInfo.locations = df_locs['long name'].tolist()
        
    def _mergeHiaHisParameters(self):
        """
        Merge parameters his en hia
        """
        df_his_pars = pd.DataFrame(self.hisFile.parameterInfo.parameters, columns=['his name'])
        
        df_hia_pars = pd.DataFrame(self.hiaFile.parameterInfo.parameters)
        # only merge if df_hia_pars contains items
        if not df_hia_pars.empty:        
            df_hia_pars.set_index('index', inplace=True)
                    
            df_pars = pd.merge(df_his_pars, df_hia_pars, left_index=True, right_index=True, how='outer')
            df_pars.loc[df_pars.isnull().any(axis=1), 'long name'] = df_pars[df_pars.isnull().any(axis=1)]['his name']
            
            self.hisFile.parameterInfo.parameters = df_pars['long name'].tolist()
            
    def _process_wrapper(self, chunk, values):
        with open(self.hisFile.hisFileName, "rb") as f:
            for offset in chunk:
                seek = f.seek(offset)
                values.append({'offset_val': offset, 'value':np.fromfile(f,np.float32,1)[0]}) 
                
    def _chunkify(self, array, no_chunks):
        chunks = np.array_split(array, no_chunks)
        chunks = [x for x in chunks if x.size > 0]
        return chunks                
        
    def GetLocations(self):
        """
        Get the available locations from the his-file. 
            
        Returns
        -------
        locations : list
            The locations available within the his-file
        """        
        locations = self.hisFile.locationInfo.locations
        return locations          
        
    def KrijgLokaties(self):
        """
        Krijg de locations beschikbaar in de his-file. 
        
        Parameters
        ----------
        his_file : str 
            pad naar his_file
        
        Returns
        -------
        locations : list
            De locations welke bekend zijn binnen het his-bestand
        """
        warnings.warn("this function will deprecate in the future, use function hkvsobekpy.GetLocations")
        locations = self.hisFile.locationInfo.locations
        return locations
    
    def GetParameters(self):
        """
        Get the available parameters from the his-file. 
            
        Returns
        -------
        parameters : list
            The parameters available within the his-file
        """        
        parameters = self.hisFile.parameterInfo.parameters
        return parameters    
    
    def KrijgParameters(self):#, his_file):
        """
        Krijg de parameters beschikbaar in de his-file. 
        
        Parameters
        ----------
        his_file : str 
            pad naar his_file        
        
        Returns
        -------
        parameters : list
            De parameters welke bekend zijn binnen het his-bestand
        """
        warnings.warn("this function will deprecate in the future, use function hkvsobekpy.GetParameters")
        parameters = self.hisFile.parameterInfo.parameters
        return parameters

    def GetTimestamps(self):
        """
        Get the available timestamps from the his-file. 
            
        Returns
        -------
        timestamps : list of datetime objects
            The timestamps available within the his-file
        """
        timestamps = self.hisFile.timestampInfo.moments
        return timestamps

    def KrijgTijdstappen(self):#, his_file):
        """
        Krijg de timestamps beschikbaar in de his-file. 
        
        Parameters
        ----------
        his_file : str 
            pad naar his_file
            
        Returns
        -------
        timestamps : list
            De timestamps welke bekend zijn binnen het his-bestand
        """
        warnings.warn("this function will deprecate in the future, use function hkvsobekpy.GetTimestamps")
        #self.LeesMetadata(his_file)
        timestamps = self.hisFile.timestampInfo.moments
        return timestamps    
    
    def ReadMetadata(self, his_file, hia_file='auto'):
        """
        Initiate file by reading headers. This function can be followed by the function DataFrame()

        Parameters
        ----------
        his_file : str
            Full path to the file
        hia_file : anyOf(path_to_file, 'auto', 'none')
            Choose from 'auto', 'none' or full path of file.
            Default is 'auto', it will try to use his_file.with_suffix('.hia')
            When 'none' is set, it won't use any .hia-file
            Other values will be treated as full path to the .hia-file
        """
        self.hisFile.metaDataIngelezen = False
        myHisFile = Path(his_file)
            
        try:
            p = myHisFile.resolve()
        except:
            # doesn't exist
            self.__errors__.fileNotFound(myHisFile)

        # exists    
        self.hisFile.hisFileName = str(p)            

        # set length filesize (bytes)
        statInfo = os.stat(his_file)
        self.hisFile.hisFileSize = statInfo.st_size

        # read his file meta data
        with open(self.hisFile.hisFileName, "rb") as f:
            try:
                self._leesAdmin(f)

            except:
                self.__errors__.administratieError()

        # continue with hia file, if not 'none'
        if not hia_file=='none':
            # try to parse a hia file
            # on auto-path
            if hia_file == 'auto':
                myHiaFile = Path(his_file).with_suffix('.hia')
            # on full-path
            else:
                myHiaFile = Path(hia_file)
                try:
                    myHiaFile.resolve()
                except:
                    # doesn't exist
                    self.__errors__.fileNotFound(myHiaFile)
            # file exist, read files 
            # merge with locations and parameters where possible
            self._leesHia(myHiaFile)
            self._mergeHiaHisLocations()
            self._mergeHiaHisParameters()                
            
        self.hisFile.metaDataIngelezen = True
        return self
        

    def LeesMetadata(self,his_file, hia_file='auto'):
        """
        Open het bestand

        Parameters
        ----------
        his_file : str
            Volledige pad van het bestand
        hia_file : anyOf(path_to_file, 'auto', 'none')
            Kies uit 'auto', 'none' of het volledig pad naar het bestand.
            Default is 'auto', het zal dan proberen om his_file.with_suffix('.hia') te gebruiken
            Wanneer 'none' gezet wordt zal er geen link met een .hia bestand gelegd worden.
            Anders is het het volledige pad van het bestand

        Returns
        -------
        Str : 'Metadata Ingelezen'
        
        """
        warnings.warn("this function will deprecate in the future, use function hkvsobekpy.ReadMetadata instead")
        self.hisFile.metaDataIngelezen = False
        myHisFile = Path(his_file)


            
        try:
            p = myHisFile.resolve()
        except:
            # doesn't exist
            self.__errors__.fileNotFound(myHisFile)

        # exists    
        self.hisFile.hisFileName = str(p)            

        # set length filesize (bytes)
        statInfo = os.stat(his_file)
        self.hisFile.hisFileSize = statInfo.st_size

        # read his file meta data
        with open(self.hisFile.hisFileName, "rb") as f:
            try:
                self._leesAdmin(f)

            except:
                self.__errors__.administratieError()

        # continue with hia file, if not 'none'
        if not hia_file=='none':
            # try to parse a hia file
            # on auto-path
            if hia_file == 'auto':
                myHiaFile = Path(his_file).with_suffix('.hia')
            # on full-path
            else:
                myHiaFile = Path(hia_file)
                try:
                    myHiaFile.resolve()
                except:
                    # doesn't exist
                    self.__errors__.fileNotFound(myHiaFile)
            # file exist, read files 
            # merge with locations and parameters where possible
            self._leesHia(myHiaFile)
            self._mergeHiaHisLocations()
            self._mergeHiaHisParameters()                
            
        self.hisFile.metaDataIngelezen = True
        return self
    
    def SelectPeriodeWaardenArray(self, df, startMMdd=(1,1), endMMdd=(12,31), jaarmax_as='none'):
        """
        Selecteer op basis van een DataFrame de gebeurtenissen binnen een bepaalde periode (bijv. een groeiseizoen).
        Waarbij er de mogelijkheid is om de gebeurtenissen binnen deze periode te groeperen.
        Pas op: in de `EnkeleWaardenArray` en `MultiWaardenArray` wordt deze functie intern ook aangeroepen.
        
        Parameters
        ----------
        df : pandas.DataFrame
            DataFrame met datetime als index
        startMMdd : tuple
            Tuple in het formaat (M,d). Bepaling van de start datum van de periode. Genoemde datum is inclusief
        endMMdd : tuple
            Tuple in het formaat (M,d). Bepaling van de eind datum van de periode. Genoemde datum is inclusief
        jaarmax_as : str
            Mogelijkheden om de DataFrame te groeperen op basis van jaar om het jaarmaxima te bepalen. 
            Het jaarmaxima wordt bepaald nadat de slice van de jaarlijkse periode is toegepast. 
            Keuze bestaat uit:            
            'date' - bepaalt de jaarlijkse maxima en geeft de maximale waarde terug met de exacte datum van deze gebeurtenis
            'year' - bepaalt de jaarlijkse maxima en geeft de maximale waarde terug met het jaar van de gebeurtenis
            'none' - retourneert alle gebeurtenissen in elk jaar
            
        Returns
        -------
        df : pandas.DataFrame
            DataFrame met datetime als index
            
        Examples
        --------
        Voor bepaling van het groeiseizoen:
        Bijvoorbeeld het jaarlijks groeseizoen loopt van 15 april tot en met 11 oktober
        startMMdd = (4,15)
        endMMdd = (10,11)
        
        betekent, slice elk jaar in april na de 15e inclusief tot aan oktober voor de 11e inclusief
        
        Voor bepaling van het jaarmaxima:
        groupby = 'date'
                    A
        2012-10-06  1501
        2013-04-22  1534
        2014-04-18  1591

        groupby = 'year'
                    A
        2012        1501
        2013        1534
        2014        1591
        """
        self.startMM, self.startdd = startMMdd
        self.endMM, self.enddd = endMMdd    
        
        # https://stackoverflow.com/a/45996897/2459096
        # maak een month_day dataframe van de MultiYear DataFrame
        month_day = pd.concat([
                        df.index.to_series().dt.month, 
                        df.index.to_series().dt.day
                    ], axis=1).apply(tuple, axis=1)        
        
        # selecteer alleen de gebeurtenissen binnen een jaarlijks terugkerende periode 
        df = df[(month_day >= (self.startMM, self.startdd)) & (month_day <= (self.endMM, self.enddd))]

        if jaarmax_as=='year':
            # groupby year en selecteer de max. 
            # hetvolgende retourneert alleen jaren + max, willen van elk jaar de datum + max
            df = df.groupby(df.index.year).max()

        elif jaarmax_as=='date':
            # groupy year en selecteer de max. Returned volle datums waar de max van dat jaar was 
            key = df.columns.levels[1][0]
            level = df.columns.names[1]
            slice_col = df.columns.levels[0][0]

            # next get the year maxima of the remaining gebeurtenissen
            idx = df.xs(key=key, level=level, axis=1).groupby([df.index.year])[slice_col].transform(max) == df.xs(key=key, level=level, axis=1)[slice_col]
            df = df[idx]            
            
            # if there are multiple events in a single year with same value take the first
            df_unique = pd.DataFrame(df.index, columns=['date'])
            df_unique['index'] = df_unique['date'].apply(lambda x:x.year)
            df_unique = df_unique.groupby('index').first()            
            #slice_unique_values = df_unique['date'].values
            #df = df.loc[slice_unique_values]
            ix2 = pd.np.isin(df.index.values, df_unique.date.values)
            df = df[ix2]

        
        elif jaarmax_as=='none':
            #doe niets
            pass
    
        return df

    def EnkeleWaardenArray(self, location, parameter, startMMdd=(1,1), endMMdd=(12,31), jaarmax_as='none'):
        """
        Lees de waarden van een enkele parameter op een enkele location

        Parameters
        ----------
        location : string    
            De location in kwestie (alleen de eerste 21 karakters, anders cutoff)
        parameter : string
            De parameter in kwestie (alleen de eerste 21 karakters, anders cutoff)
        startMMdd : tuple
            Tuple in het formaat (M,d). Bepaling van de start datum van de periode. Genoemde datum is inclusief
        endMMdd : tuple
            Tuple in het formaat (M,d). Bepaling van de eind datum van de periode. Genoemde datum is inclusief
        jaarmax_as : str
            Mogelijkheden om de DataFrame te groeperen op basis van jaar om het jaarmaxima te bepalen. 
            Het jaarmaxima wordt bepaald nadat de slice van de jaarlijkse periode is toegepast. 
            Keuze bestaat uit:            
            'date' - bepaalt de jaarlijkse maxima en geeft de maximale waarde terug met de exacte datum van deze gebeurtenis
            'year' - bepaalt de jaarlijkse maxima en geeft de maximale waarde terug met het jaar van de gebeurtenis
            'none' - retourneert alle gebeurtenissen in elk jaar            

        Returns
        -------
        df : DataFrame
            DataFrame met een multicolumn van parameter en locations met datetime index en 
            de bijbehorende waarden
            
        Examples
        -----------
        location = '1'
        parameter = 'Waterlevel maximum (mNAP)'
        
        Intern worden alle inputs afgekort op [0:20] en ziet de tool het als:
        location = '1'
        parameter = 'Waterlevel maximum ('

        Voor bepaling van het groeiseizoen:
        Bijvoorbeeld het jaarlijks groeseizoen loopt van 15 april tot en met 11 oktober
        startMMdd = (4,15)
        endMMdd = (10,11)
        
        betekent, slice elk jaar in april na de 15e inclusief tot aan oktober voor de 11e inclusief
        
        Voor bepaling van het jaarmaxima:
        groupby = 'date'
                    A
        2012-10-06  1501
        2013-04-22  1534
        2014-04-18  1591

        groupby = 'year'
                    A
        2012        1501
        2013        1534
        2014        1591
        """

        if not hasattr(self.hisFile, 'metaDataIngelezen'):
            self.__errors__.metadataNotSet()
        
        if self.hisFile.metaDataIngelezen == False:
            self.__errors__.metadataError()
        
        loc = location
        param = parameter
        
        paramLocOffset = self._locOffset(loc, param)
        offsets = self.hisFile.timestampInfo.offset + paramLocOffset
        values= []
        
        with open(self.hisFile.hisFileName, "rb") as f: 

            for offset in offsets:               
                seek = f.seek(offset)
                values.append(np.fromfile(f,np.float32,1)[0])

        # maak dataframe
        df = pd.DataFrame(data=values, index=self.hisFile.timestampInfo.moments, columns=[(loc,param)])
        df.columns = pd.MultiIndex.from_tuples(df.columns, names=['location','parameter'])
        
        df = self.SelectPeriodeWaardenArray(df, startMMdd=startMMdd, endMMdd=endMMdd, jaarmax_as=jaarmax_as)

        return df#.T.squeeze()
    
    def MultiWaardenArray(self, locations, parameters, startMMdd=(1,1), endMMdd=(12,31), 
                          jaarmax_as='none', drop_lege_jaren=True):
        """
        Lees de waarden van meerdere parameters op meerdere locations

        Parameters
        ----------
        locations : list    
            Een list van de locations (type: str) in kwestie (alleen de eerste 21 karakters, anders cutoff) 
        parameters : list
            Een list van de parameters in kwestie (alleen de eerste 21 karakters, anders cutoff) 
        startMMdd : tuple
            Tuple in het formaat (M,d). Bepaling van de start datum van de periode. Genoemde datum is inclusief
        endMMdd : tuple
            Tuple in het formaat (M,d). Bepaling van de eind datum van de periode. Genoemde datum is inclusief
        jaarmax_as : str
            Mogelijkheden om de DataFrame te groeperen op basis van jaar om het jaarmaxima te bepalen. 
            Het jaarmaxima wordt bepaald nadat de slice van de jaarlijkse periode is toegepast. 
            Keuze bestaat uit:            
            'year' - bepaalt de jaarlijkse maxima en geeft de maximale waarde terug met het jaar van de gebeurtenis
            'none' - retourneert alle gebeurtenissen in elk jaar  
        drop_lege_jaren : boolean
            Mogelijkheid om de jaren te verwijderen welke geen gebeurtenissen bevatten voor de gegeven selectie. 
            Wanneer dit als `True` is geselecteerd zal het jaar alleen gedropt worden als deze niet bestaat voor 
            alle locations/parameters in de selectie

        Returns
        -------
        df : DataFrame
            DataFrame met een multicolumn van parameter en locations met datetime index en 
            de bijbehorende waarden
            
        Examples
        -----------
        locations = ['1','6']
        parameters = ['Waterlevel maximum (mNAP)']
        
        Intern worden alle inputs afgekort op [0:20] en ziet de tool het als:
        locations = ['1','6']
        parameters = ['Waterlevel maximum (']
        """
        # error checking
        if jaarmax_as not in ('year', 'none'):
            self.__errors__.annualmaxError()
            
        # define empty dataframe
        clmns = pd.MultiIndex.from_product([parameters,locations], names=['parameters','locations'])
        
        # voor een MultiWaardenArray wordt alleen een jaar meegenomen als index
        # het hangt van een seizoen filter af of er gebeurtenissen voor dat jaar zijn. Bepaal dit eerst
        loc0 = self.hisFile.locationInfo.locations[0]
        par0 = self.hisFile.parameterInfo.parameters[0]
        idx = self.EnkeleWaardenArray(loc0, par0, jaarmax_as=jaarmax_as).index
        
        self.hisFile.timestampInfo.N
        
        df_full = pd.DataFrame(index=idx, columns=clmns)
        df_full.index.name='date'
        for ix in clmns:
            df = self.EnkeleWaardenArray(location=ix[1], parameter=ix[0], 
                                         startMMdd=startMMdd, endMMdd=endMMdd, jaarmax_as=jaarmax_as)
            df_full.loc[:,(ix[0],ix[1])] = df.T.squeeze()
        
        # bepaling of lege jaren gedropt moeten worden
        if drop_lege_jaren:
            df_full.dropna(axis=0, how='all', inplace=True)
        
        return df_full
        
    def read_series(self, his_file, location, parameter, his_folder, normalize_by_unicode=True, 
                include_simularity=True, sequence_simularity=0.82, return_matching_parameter=False):
        """
        Extract timeseries from his-file
        
        Parameters
        ----------
        his_file : str
            name of his-file to query (e.g: 'reachseg', 'calcpnt', 'struc')
        location : str
            location to query in his-file
        parameter : str
            parameter to query in his-file
        his_folder : str
            path to folder containing the his_file
        normalize_by_unicode : boolean (default True)
            inlcude this option to include NFKD unicode compatibility decomposition.
            see: http://unicode.org/reports/tr15/
        include_simluratiy : boolean (default True)
            include this option to include Ratcliff/Obershelp pattern recognition
        sequence_simularity : float (default 0.82)
            number between 0.0 and 1.0, function as threshold, where only a simularity above 
            this value is mapped            
        
        Returns
        -------
        df : pandas.DataFrame
            timeseries containing all available timesteps in his-file
        """
        # get metadata from appropriate his-file
        hisfile = self.LeesMetadata(os.path.join(his_folder, '{0}.his'.format(his_file).upper()))
        his_parameters = hisfile.KrijgParameters()
        his_locations = hisfile.KrijgLokaties()
        
        # normalize unicode and/or check simularity of parameters 
        parameter = compare_df_parameter_his_parameter(parameter.ljust(20), his_parameters, normalize_by_unicode, 
                                               include_simularity, sequence_simularity)
                                               
        # get dataframe for single location/parameter combination
        df = self.EnkeleWaardenArray(
            location=location, 
            parameter=parameter, 
            jaarmax_as='none')
        if return_matching_parameter == True:
            return df, parameter
        else:
            return df
        
    def DataFrame(self):
        """
        Extract all timeseries from his-file as pandas.DataFrame. This is quick.
        
        Returns
        -------
        df : pandas.DataFrame
            containing all available timesteps for all locations and parameters in his-file
        """
        # get metadata from appropriate his-file
        no_moments = len(self.hisFile.timestampInfo.moments)
        no_params = self.hisFile.parameterInfo.numPar
        no_locs = self.hisFile.locationInfo.numLoc

        # the byte range goes like:
        # for t in times: for l in locs: for p in params: print (t, l, p)
        with open(self.hisFile.hisFileName, "rb") as f: 
            start = self.hisFile.timestampInfo.offset[0]
            no_values = (no_moments * no_params * no_locs) + no_moments
            seek = f.seek(start-4)    
            values = np.fromfile(f, np.float32, no_values) # np.float32 == np.dtype('f4')

        # reshape on number of timestamps
        values = values.reshape(no_moments, -1)
        # filter by slice as first column contains 4 non-assigned bytes
        values = np.ravel(values[:,1::])
        # reshape to prepare for dataframe
        values = values.reshape((no_params, no_locs, no_moments), order='F')
        values = values.T.reshape(no_moments, -1)              
        
        # create multi-index for the columns
        lable_locs = self.hisFile.locationInfo.locations
        lable_params = self.hisFile.parameterInfo.parameters
        cols = pd.MultiIndex.from_product([lable_locs, lable_params])

        # parse into dataframe
        df = pd.DataFrame(values, index=self.hisFile.timestampInfo.moments, columns=cols)    
        df = df.swaplevel(i=0, j=1, axis=1)
        df.sort_index(axis=1,inplace=True)
        return df