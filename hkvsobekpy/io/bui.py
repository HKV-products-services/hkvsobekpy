import argparse
import pandas as pd
import datetime
import numpy as np
# import logging
# from importlib import reload
# # initiatie logging
# reload(logging)
# logging.basicConfig(
    # #filename="{0}/{1}.log".format(logPath, fileName),
    # format='%(asctime)s [%(threadName)-12.12s] [%(levelname)-5.5s]  %(message)s', 
    # level=logging.INFO)
# # to close file:
# #logging.shutdown()


class __bui_class(object):
    def __init__(self):
        return

    def read_bui(self,filename):
        """
        parse een bui-file naar een dataframe
        Parameters
        ----------
        filename : str
            pad naar het .bui bestand
        
        Returns
        -------
        df : pandas.DataFrame
            dataframe met in de kolommen de neerslag in mm voor de stations
            en in de rijen de waarnemingstijdstappen
        """
        # initatie van lege immutable objects
        # zie: https://stackoverflow.com/a/22526544/2459096
        # python2.7: obj = argparse.Namespace()
        # python3.x: obj = types.SimpleNamespace()        
        self.buiFile = argparse.Namespace()
        with open(filename, 'r') as infile:
            f = infile.readlines()
            for i in range(len(f)):
                #print ('regel {0}'.format(str(i)))
                line = f[i]
                if line[0] == '*':
                    # logging.info('regel {0} is comment'.format(str(i)))
                    
                    if 'GEBRUIK DE DEFAULT DATASET' in line.upper():
                        self.buiFile.default_dataset = int(f[i+1])
                    if 'AANTAL STATION' in line.upper():
                        self.buiFile.aantal_stations = int(f[i+1])
                    if 'NAMEN VAN STATION' in line.upper():
                        stations = f[i+1:i+1+self.buiFile.aantal_stations]
                        self.buiFile.stations = [station.rstrip().replace("'","") for station in stations]
                    if 'AANTAL SECONDEN PER WAARNEMINGSSTAP' in line.upper():
                        geb_sec = f[i+1].split()
                        self.buiFile.aantal_gebeurtenissen = int(geb_sec[0])
                        self.buiFile.aantal_seconden = int(geb_sec[1])
                    if 'HET FORMAT IS: YYYY' in line.upper():
                        T0_raw = f[i+2]#.split()                
                        idx_block = i+1
            # get data block
            # voor elk station de neerslag in mm per tijdstap
            s = pd.Series(f[idx_block+2:len(f)-1])

        # parse datablock to DataFrame
        df = s.str.rstrip().str.lstrip().str.split(' +', expand=True)
        df = df.apply(pd.to_numeric)            

        # startdatum en -tijd
        T0_str = T0_raw.rstrip().lstrip().split()
        # Het format is: yyyymmdd:hhmmss
        self.buiFile.start_datum        = datetime.datetime(*list(map(int,T0_str[0:5]))) #seconden is eruit gehaald anders str[0:6]
        # lengte van de gebeurtenis in dd hh mm ss
        lengte_gebeurtenis = list(map(int,T0_str[-4::]))
        self.buiFile.delta_dag          = lengte_gebeurtenis[0]
        self.buiFile.delta_uur          = lengte_gebeurtenis[1]
        self.buiFile.delta_minuut       = lengte_gebeurtenis[2]
        self.buiFile.delta_seconde      = lengte_gebeurtenis[3]

        # date_range gebeurtenis
        gebeurtenis_date_range = pd.date_range(
            start=self.buiFile.start_datum, 
            end=None, 
            periods=df.shape[0], 
            freq='{0}S'.format(self.buiFile.aantal_seconden)
        )

        df.set_index(gebeurtenis_date_range, inplace=True)
        df.columns = self.buiFile.stations    

        return df
    
    @staticmethod
    def seconds_toperiods(total_seconds):
        """
        Ontleedt een getal in seconden naar lijst met vaste periodes: days, hours, minutes, seconds
        waarbij day    =          24 * 60 * 60 =     86,400 seconden 
                hour   =               60 * 60 =      3,600 seconden 
                minute =                    60 =         60 seconden
        
        Parameters
        --------------
        total_seconds : float or int 
            bijvoorbeeld Timedelta.total_seconds() van Pandas dataframe bv index[1]-index[0]
        
        Returns
        -------------
        list : list of integer values [days, hours, minutes, seconds]
        """
        
        #years, remainder = divmod(total_seconds,365.25*3600*24) #total seconds per year
        days, remainder = divmod(total_seconds, 24*3600)
        hours, remainder = divmod (remainder, 3600)
        minutes, seconds = divmod(remainder, 60)
        
        return [int(days), int(hours), int(minutes), int(seconds)]
        
    def write_bui(self,df,filename,dataset="bui",comment=None):
        """ 
        Maak van een pandas dataframe een sobek buifile
        
        Parameters
        ----------
        df : pandas dataframe
            index df is datumtijdas. aanname is uniforme tijdstap delta time
            in kolommen neerslag in mm/tijdstap
            kolomnaam = stationnaam
        
        filename: str
            path + naam buifile "{str:8}.bui" lengte buinaam is 8 tekens
        
        dataset: str [default = "bui"] options ["bui", "reeks"]
            type dataset dit is voorgedefinieerd door sobek.
            
        comment: str
            extra opmerkingenregels die in de header van de buifile worden geplaatst
        
        
        RETURNS
        -----------
        sobekbuifile : gebruiksklare buifile die door sobek2 kan worden ingelezen.
        
        Standaard format sobekbuifile: 
        * user defined comment
        * Gebruik de default dataset voor overige invoer (altijd 1 bij bui, 0 bij reeks)
        1
        * Aantal stations
        2
        * Namen van Stations
        'Station1'
        'Station2'
        * Aantal gebeurtenissen en het aantal seconden per waarnemingsstap
        1 3600
        * Elke commentaarregel wordt begonnen met een * (asteriks).
        * Meteo data: neerslag stations; voor elk station: neerslag intensiteit in mm
        * Eerste record bevat start datum & tijd en lengte van de gebeurtenis in dd hh mm ss
        * Het format is: yyyy mm dd hh mm ss
        * Daarna voor elk station de neerslag in mm per tijdstap.
        2010 11 07 00 00 00 14 00 00 00
        ... space seperated values in [mm/timestep] ... 
        """
        with open(filename, 'w') as bui:
            if not comment==None: 
                bui.write("* "+comment+"\n")
            
            bui.write("* Deze buifile is automatisch aangemaakt door de buifile generator van HKV Lijn in Water"+"\n")
            bui.write("* Gebruik de default dataset voor overige invoer (altijd 1 bij bui, 0 bij reeks)"+"\n")
            if dataset.upper() == "BUI":
                bui.write("1"+"\n")
            if dataset.upper() == "REEKS":
                #throw error
                print("not supported")
            
            bui.write("* Aantal stations"+"\n")
            bui.write(str(df.shape[1])+"\n")
            
            
            bui.write("* Namen van stations"+"\n")
            for station in df.columns:
                bui.write("'{}'".format(str(station))+"\n")
            
            bui.write("* Aantal gebeurtenissen en het aantal seconden per waarnemingsstap"+"\n")
            timedelta = int(datetime.datetime.timestamp(df.index[1]) - datetime.datetime.timestamp(df.index[0]))
            if dataset.upper() == "BUI": 
                bui.write("1 "+str(timedelta)+"\n")
        
            bui.write("* Elke commentaarregel wordt begonnen met een * (asteriks)."+"\n")
            bui.write("* Meteo data: neerslag stations; voor elk station: neerslag intensiteit in mm."+"\n")
            bui.write("* Eerste record bevat start datum & tijd en lengte van de gebeurtenis in dd hh mm ss"+"\n")
            bui.write("* Het format is: yyyy mm dd hh mm ss"+"\n")
            bui.write("* Daarna voor elk station de neerslag in mm per tijdstap."+"\n")
            
            lengte_bui = self.seconds_toperiods((df.index[-1] - df.index[0]).total_seconds()) 
            bui.write(df.index[0].strftime('%Y %m %d %H %M %S ')+"{:0^2} {:0^2} {:0^2} {:0^2}".format(lengte_bui[0],lengte_bui[1],lengte_bui[2],lengte_bui[3])+"\n")
            
            df.to_csv(bui, sep = " ", index = False, header=False)
        
            bui.close()
        
        
        
        
        
        