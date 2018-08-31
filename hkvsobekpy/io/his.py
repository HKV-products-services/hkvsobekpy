import os
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
        def fileNotFound():
            raise IOError('HIS-bestand bestaat niet. Is je pad goed?')        
        @staticmethod
        def metadataNotSet():
            raise AttributeError('Metadata is niet bekend. Zet met sobek.LeesMetadata(his_file)')
        @staticmethod
        def metadataError():
            raise AttributeError('Metadata kon tijdens de functie LeesMetadata() niet bepaald worden')
        @staticmethod
        def administratieError():            
            raise AttributeError('Administratieblok van his-file kan niet worden uitgelezen')
        @staticmethod
        def locationNotFound():            
            raise AttributeError('location niet gevonden. Is het een bestaande location?')
        @staticmethod
        def longLocationNotFound():
            raise AttributeError('Long location niet gevonden in hia.')
        @staticmethod
        def longParameterNotFound():
            raise AttributeError('Long Parameter niet gevonden in hia.')
        @staticmethod
        def unexpectedT0Error():            
            raise AttributeError('T0 kon niet uitgelezen worden. Contact HKV')            
        @staticmethod
        def variabeleNotFound():            
            raise AttributeError('Parameter niet gevonden. Is het een bestaande parameter?')
        @staticmethod
        def jaarmaxError():
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
        self.hisFile.tijdstapInfo = argparse.Namespace()
        self.hisFile.variabeleInfo = argparse.Namespace()
        self.hisFile.locationInfo = argparse.Namespace()       
        self.hiaFile = argparse.Namespace()
        self.hiaFile.variabeleInfo = argparse.Namespace()
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
        
    def _infoFromHeader(self,tijdstapInfo):
        """
        Lees informatie uit de header

        Parameters
        ----------
        tijdstapInfo : str
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

        if len(tijdstapInfo) == 40:
            year = tijdstapInfo[4:8]
            month = tijdstapInfo[9:11]
            day = tijdstapInfo[12:14]
            hours = tijdstapInfo[15:17]
            minutes = tijdstapInfo[18:20]
            seconds = tijdstapInfo[21:23]

            beginDate = datetime(int(year), int(month), int(day), int(hours), int(minutes), int(seconds))

            timeStepInterval = tijdstapInfo[-2:-1]
            timeStepFactor = int(tijdstapInfo[30:-2])
        else:
            self.__errors__.unexpectedT0Error()

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
                self.hiaFile.variabeleInfo.variabelen = parameters
            except KeyError as e:
                self.__errors__.longParametersNotFound()

        except:
            # hia doesn't exist return empty objects
            self.hiaFile.locationInfo.locations = []
            self.hiaFile.variabeleInfo.variabelen = []
    
    def _leesAdmin(self, f):
        """
        Lees het administratie blok
        """
        self.hisFile.header = self._leesHeader(f)

        # lees info van header
        self.hisFile.tijdstapInfo.beginDate, self.hisFile.tijdstapInfo.timeStepInterval, self.hisFile.tijdstapInfo.timeStepFactor = self._infoFromHeader(self.hisFile.header[3])        

        # lees aantal variablelen
        self.hisFile.variabeleInfo.numVar = np.fromfile(f,np.int,1)[0]    

        # lees aantal locations
        self.hisFile.locationInfo.numLoc = np.fromfile(f,np.int,1)[0]    

        # lees variabelen
        self.hisFile.variabeleInfo.variabelen = []
        for i in range(self.hisFile.variabeleInfo.numVar):
            self.hisFile.variabeleInfo.variabelen.append(f.read(20).decode('windows-1252'))    

        # Lees locations
        self.hisFile.locationInfo.id = []
        self.hisFile.locationInfo.locations = []
        for i in range(self.hisFile.locationInfo.numLoc):            
            self.hisFile.locationInfo.id.append(np.fromfile(f,np.int,1)[0])
            self.hisFile.locationInfo.locations.append(f.read(20).decode('windows-1252').rstrip())
            #self.hisFile.locationInfo.locations.append(locationInfo)

        # lees tijdstappen
        self.hisFile.headerSize = f.tell()
        self.hisFile.locationInfo.numTime = int((self.hisFile.hisFileSize - self.hisFile.headerSize) /
                      (self.hisFile.locationInfo.numLoc * self.hisFile.variabeleInfo.numVar + 1) / 4)

        aantalBytesPerStap = int((self.hisFile.hisFileSize - self.hisFile.headerSize) / self.hisFile.locationInfo.numTime)

        byteNr = self.hisFile.headerSize

        self.hisFile.tijdstapInfo.moments = []
        
        self.hisFile.tijdstapInfo.offset = []
        for i in range(self.hisFile.locationInfo.numTime):
            f.seek(byteNr)
            moment = np.fromfile(f,np.int,1)[0]

            if self.hisFile.tijdstapInfo.timeStepInterval == 's':
                self.hisFile.tijdstapInfo.moments.append(self.hisFile.tijdstapInfo.beginDate + timedelta(seconds = (float(moment) * self.hisFile.tijdstapInfo.timeStepFactor)))
            elif self.hisFile.tijdstapInfo.timeStepInterval == 'm':
                self.hisFile.tijdstapInfo.moments.append(self.hisFile.tijdstapInfo.beginDate + timedelta(minutes = (float(moment) * self.hisFile.tijdstapInfo.timeStepFactor)))
            elif self.hisFile.tijdstapInfo.timeStepInterval == 'h':
                self.hisFile.tijdstapInfo.moments.append(self.hisFile.tijdstapInfo.beginDate + timedelta(hours= (float(moment) * self.hisFile.tijdstapInfo.timeStepFactor)))
            elif self.hisFile.tijdstapInfo.timeStepInterval == 'd':
                self.hisFile.tijdstapInfo.moments.append(self.hisFile.tijdstapInfo.beginDate + timedelta(days= (float(moment) * self.hisFile.tijdstapInfo.timeStepFactor)))
            self.hisFile.tijdstapInfo.offset.append(f.tell())
            byteNr += aantalBytesPerStap
        self.hisFile.tijdstapInfo.N = np.array(self.hisFile.tijdstapInfo.moments).max().year - np.array(self.hisFile.tijdstapInfo.moments).min().year + 1
        return True

    def _locOffset(self, loc, var):
        """
        Krijg de offset van de variabele

        Parameters
        ----------
        loc : string    
            De location in kwestie
        var : string
            De variabele in kwestie        

        Returns
        -------
        locOffset : int
            De offset van de location in bytes
        """
        #  Zoek de index van de variabele
        try:
            varFound = self.hisFile.variabeleInfo.variabelen.index(var)
        except :
            self.__errors__.variabeleNotFound()            

        # Zoek de index/id van de location
        try:
            # removed dictionary
            # locFound = next((item for item in self.hisFile.locationInfo.locations if item['location'] == loc))
            locFound = self.hisFile.locationInfo.locations.index(loc)
        except :
            self.__errors__.locationNotFound()

        locOffset = locFound * self.hisFile.variabeleInfo.numVar * 4 + varFound * 4

        return locOffset   
    
    def _mergeHiaHisLocations(self):
        """
        Merge locations his en hia
        """
        df_his_locs = pd.DataFrame(self.hisFile.locationInfo.locations, columns=['his name'])
        
        df_hia_locs = pd.DataFrame(self.hiaFile.locationInfo.locations)
        df_hia_locs.set_index('index', inplace=True)
        
        df_locs = pd.merge(df_his_locs, df_hia_locs, left_index=True, right_index=True, how='outer')
        df_locs.loc[df_locs.isnull().any(axis=1), 'long name'] = df_locs[df_locs.isnull().any(axis=1)]['his name']
        
        self.hisFile.locationInfo.locations = df_locs['long name'].tolist()
        
    def _mergeHiaHisParameters(self):
        """
        Merge parameters his en hia
        """
        df_his_pars = pd.DataFrame(self.hisFile.variabeleInfo.variabelen, columns=['his name'])
        
        df_hia_pars = pd.DataFrame(self.hiaFile.variabeleInfo.variabelen)
        df_hia_pars.set_index('index', inplace=True)
                
        df_pars = pd.merge(df_his_pars, df_hia_pars, left_index=True, right_index=True, how='outer')
        df_pars.loc[df_pars.isnull().any(axis=1), 'long name'] = df_pars[df_pars.isnull().any(axis=1)]['his name']
        
        self.hisFile.variabeleInfo.variabelen = df_pars['long name'].tolist()
        
    def KrijgLokaties(self):#, his_file):
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
        #self.LeesMetadata(his_file)
        locations = self.hisFile.locationInfo.locations
        return locations
    
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
        #self.LeesMetadata(his_file)
        parameters = self.hisFile.variabeleInfo.variabelen
        return parameters

    def KrijgTijdstappen(self):#, his_file):
        """
        Krijg de tijdstappen beschikbaar in de his-file. 
        
        Parameters
        ----------
        his_file : str 
            pad naar his_file
            
        Returns
        -------
        tijdstappen : list
            De tijdstappen welke bekend zijn binnen het his-bestand
        """
        #self.LeesMetadata(his_file)
        tijdstappen = self.hisFile.tijdstapInfo.moments
        return tijdstappen    

    def LeesMetadata(self,his_file, hia_file='auto'):
        """
        Open het bestand

        Parameters
        ----------
        his_file : str
            Volledige pad van het bestand

        Returns
        -------
        Str : 'Metadata Ingelezen'
        
        """
        self.hisFile.metaDataIngelezen = False
        myHisFile = Path(his_file)
        
        # Try to parse a hia file
        if hia_file == 'auto':
            myHiaFile = Path(his_file).with_suffix('.hia')
        else:
            myHiaFile = Path(hia_file)
            try:
                myHiaFile.resolve()
            except:
                # doesn't exist
                self.__errors__.fileNotFound()
        self._leesHia(myHiaFile)
            
        try:
            p = myHisFile.resolve()
        except:
            # doesn't exist
            self.__errors__.fileNotFound()
        else:
            # exists    
            self.hisFile.hisFileName = str(p)            
            # set length filesize (bytes)
            statInfo = os.stat(his_file)
            self.hisFile.hisFileSize = statInfo.st_size
            #f = open(str(p), "rb")
            with open(self.hisFile.hisFileName, "rb") as f:
                try:
                    self._leesAdmin(f)
                    
                    self._mergeHiaHisLocations()
                    self._mergeHiaHisParameters()
                except:
                    self.__errors__.administratieError()
            
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
            slice_unique_values = df_unique['date'].values
            df = df.loc[slice_unique_values]

        
        elif jaarmax_as=='none':
            #doe niets
            pass
    
        return df

    def EnkeleWaardenArray(self, location, parameter, startMMdd=(1,1), endMMdd=(12,31), jaarmax_as='none'):
        """
        Lees de waarden van een enkele variabele op een enkele location

        Parameters
        ----------
        location : string    
            De location in kwestie (alleen de eerste 21 karakters, anders cutoff)
        parameter : string
            De variabele in kwestie (alleen de eerste 21 karakters, anders cutoff)
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
            DataFrame met een multicolumn van variabele en locations met datetime index en 
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
        var = parameter
        
        with open(self.hisFile.hisFileName, "rb") as f: 
            varLocOffset = self._locOffset(loc, var)
            values= []

            for i in range(self.hisFile.locationInfo.numTime):
                offset = self.hisFile.tijdstapInfo.offset[i] + varLocOffset                
                seek = f.seek(offset)
                values.append(np.fromfile(f,np.float32,1)[0])

        # maak dataframe
        df = pd.DataFrame(data=values, index=self.hisFile.tijdstapInfo.moments, columns=[(loc,var)])
        df.columns = pd.MultiIndex.from_tuples(df.columns, names=['location','parameter'])
        
        df = self.SelectPeriodeWaardenArray(df, startMMdd=startMMdd, endMMdd=endMMdd, jaarmax_as=jaarmax_as)

        return df#.T.squeeze()
    
    def MultiWaardenArray(self, locations, parameters, startMMdd=(1,1), endMMdd=(12,31), 
                          jaarmax_as='none', drop_lege_jaren=True):
        """
        Lees de waarden van meerdere variabelen op meerdere locations

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
            DataFrame met een multicolumn van variabele en locations met datetime index en 
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
            self.__errors__.jaarmaxError()
            
        # define empty dataframe
        clmns = pd.MultiIndex.from_product([parameters,locations], names=['parameters','locations'])
        
        # voor een MultiWaardenArray wordt alleen een jaar meegenomen als index
        # het hangt van een seizoen filter af of er gebeurtenissen voor dat jaar zijn. Bepaal dit eerst
        loc0 = self.hisFile.locationInfo.locations[0]
        par0 = self.hisFile.variabeleInfo.variabelen[0]
        idx = self.EnkeleWaardenArray(loc0, par0, jaarmax_as=jaarmax_as).index
        
        self.hisFile.tijdstapInfo.N
        
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