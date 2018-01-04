import os
try:
    from pathlib import Path
except:    
    from pathlib2 import Path
import types
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
from hkvsobekpy.io.his import __his_class
his_reader = __his_class()

class __waterlevelstat_class(object):

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
    def _initNamespaces(self):
        try:
            self.stats = types.SimpleNamespace()
            self.stats.plotposities = types.SimpleNamespace()
            self.stats.stap1 = types.SimpleNamespace()
            self.stats.stap2 = types.SimpleNamespace()
            self.stats.stap3 = types.SimpleNamespace()
            self.stats.stap4 = types.SimpleNamespace()     
        except Exception as e:
            print(e)
            self.stats = argparse.Namespace()
            self.stats.plotposities = argparse.Namespace()
            self.stats.stap1 = argparse.Namespace()
            self.stats.stap2 = argparse.Namespace()
            self.stats.stap3 = argparse.Namespace()
            self.stats.stap4 = argparse.Namespace()
        return self

    class __errors__(object):
        """
        error class met verschillende foutmeldingen
        """
        @staticmethod
        def gewogenGemiddeldeError():
            raise ValueError('Kan gewogen gemiddelde niet bepalen. Vereiste is minimaal 2 punten aan de linker en rechterzijde van de gekozen T')
        def gebeurtenissenError():
            raise ValueError('Meer gebeurtenissen dan N. Zorg ervoor dat je tijdreeks slechts 1 gebeurtenis per jaar heeft')
        def schrijfTabelError():
            raise ValueError('Kan deze tabel alleen wegschrijven met T10,25,50,100 voor Gumbel en T10 voor gewogen gemiddelde')

    def _best_fit_slope_and_intercept(self,xs,ys):
        try:
            warnings.filterwarnings("error")
            m = (((xs.mean()*ys.mean()) - (xs*ys).mean()) / ((xs.mean()*xs.mean()) - (xs*xs).mean()))
        except RuntimeWarning:
            #set_trace()
            m = 0#np.nan
        b= ys.mean() - m*xs.mean()
        return m, b

    def _squared_error(self,ys_orig,ys_line):
        return sum((ys_line - ys_orig) * (ys_line - ys_orig))

    def _coefficient_of_determination(self,ys_orig,ys_line):
        y_mean_line = [ys_orig.mean() for y in ys_orig]
        squared_error_regr = self._squared_error(ys_orig, ys_line)
        squared_error_y_mean = self._squared_error(ys_orig, y_mean_line)
        try:
            warnings.filterwarnings("error")
            _r_squared = 1 - (squared_error_regr/squared_error_y_mean)
        except RuntimeWarning:
            _r_squared = 0
        return _r_squared

    def _r_squared(self,xs,ys):
        """
        Bepaalt de r**2 waarde op basis van twee input arrays

        Parameters
        ----------
        xs : numpy.array
            array met orginele waarden
        ys : numpy.array
            array met afgeleide waarden

        Returns
        -------
        r_squared : float
            de afgeleide r**2 waarde
        """
        m, b = self._best_fit_slope_and_intercept(xs,ys)
        regression_line = [(m*x)+b for x in xs]

        r_squared = self._coefficient_of_determination(ys,regression_line)
        return r_squared        
        
    def _igumbelFunc(self, arr, Td, t):
        """
        Gumbel Functie
        
        Parameters
        ----------
        arr : numpy.array
            array met gesorteerde waterstanden
        Td : 
        
        Returns
        -------
        gumbelFuncFunc : numpy.array
            Gumbel Functie
        """
        a,d = arr
        self.stats.plotposities.gumbelFunc = (float(-1.0)* float(a)*np.log(1.0/t * float(Td))+float(d))
        return self.stats.plotposities.gumbelFunc
    
    def _igumbelFuncFit(self, arr, h, Td, t):
        """
        Gumbel Functie
        
        Parameters
        ----------
        arr : numpy.array
            array met gesorteerde waterstanden
        h : int/float
        
        Td : int/float
        
        T : int/float
        
        Returns
        -------
        gumbelFuncFit : numpy.array
            Fitted Gumbel Functie
        """        
        a,d = arr
        self.stats.plotposities.gumbelFuncFit = h-(float(-1.0)* float(a)*np.log(1.0/t * float(Td))+float(d))        
        return self.stats.plotposities.gumbelFuncFit
    
    def _gewogenGemiddelde4TOI(self, Tarray, WSarray, TOI):
        """
        Functie om het gewogen gemiddelde te krijgen op basis van de vier dichtsbijzijnde waarden rondom de T of interest

        Parameters
        ----------
        Tarray : numpy.array
            Array met de terugkeertijden, 
        WSarray : numpy.array
            Array met de parameterwaarden, bijvoorbeeld waterstanden
        TOI : int
            Waarde met de terugkeertijd of interest, bijvoorbeeld 10 voor T=10

        Returns
        -------

        """
        # get the index of the closest value to T10 on its left side, 
        # T is sorted in descending order so reverse array first to get ascending order
        cl = self._find_closest(Tarray[::-1], TOI)    

        # get the indices of the four closest points
        # [..]  *     *      *      *   |     *     *  [..]
        #           cl-2   cl-1    cl  T10  cl+1  cl+2
        #                  ----    --       ----  ----
        try:
            indices = np.array([cl-1,cl,cl+1,cl+2])

            WSarray_sel  = np.take(WSarray[::-1], indices)
            Tarray_sel   = np.take(Tarray[::-1], indices) 
        except:
            self.__errors__.gewogenGemiddeldeError()

        # inverse distance weighting
        weights = 1 / abs(Tarray_sel - TOI)
        weights /= weights.sum()    

        # weight_cl-1 * ws_cl-1  +  weight_cl * ws_cl  +  weight_cl+1 * ws_cl+1  +  weight_cl+2 * ws_cl+2
        WSarray_TOI = np.dot(weights, WSarray_sel)
        return WSarray_TOI, WSarray_sel, Tarray_sel

    def _find_closest(self, Tarray, TOI):
        # T must be sorted (in ascending order)
        idx = Tarray.searchsorted(TOI)
        idx = np.clip(idx, 1, len(Tarray)-1)
        left = Tarray[idx-1]
        right = Tarray[idx]
        idx -= TOI - left < right - TOI
        return idx    
        
    def _calcFi_T(self, len_arr,N, Ggi=0.44, GgN=0.12):
        """
        Bepaal plotpositie jaarmaximum als kans        

        Parameters
        ----------
        len_arr : int
            aantal gebeurtenissen
        N : int
            aantal jaren waarover de Gumble fit bepaald moet worden
        Ggi : float
            Gringgorten plotposities i (standaard: 0.44)
        GgN : float
            Gringgorten plotposities N (standaard: 0.12)

        Returns
        -------
        T : np.array
            array van terugkeertijden
        Fijaar : np.array
            array van plotposities jaarmaxima
            
        Literature
        ----------
        Voor bepaling Gringorten coefficienten zie bijv: http://glossary.ametsoc.org/wiki/Gringorten_plotting_position
        """
        Fi = []        
        T = []
        for i in range(len_arr):
            i += 1
            f = float(1 - ((i - Ggi) / (N+GgN)))
            Fi.append(f)
        try:
            warnings.filterwarnings("error")
            Fijaar = (np.log(Fi) )* -1
        except RuntimeWarning:
            #set_trace()
            pass
        T = 1 / (Fijaar)
        # set set self
        self.stats.plotposities.Fijaar = Fijaar
        self.stats.plotposities.T = T
        return self.stats.plotposities.T, self.stats.plotposities.Fijaar
    
    def PlotFiguur(self,stat_object,out_folder='none'):
        """
        plot figuur
        
        Parameters
        ----------
        ID : str
            locationnaam
        
        ## STEP1 Gumbel Fit voor T25,T50,T100 voor alle buien
        S1_vAmin : int
            venster array minimum
        S1_vAmax : int
            venster array maximum
        S1_T : ndarray
            Terugkeertijden van jaarmaxima gebeurtenissen genomen over gehele jaar
        S1_ws_srt : ndarray
            Waterstanden welke horen bij de terugkeertijden van jaarmaxima gebeurtenissen genomen over gehele jaar
        S1_GumbelT : ndarray
            Terugkeertijden welke horen bij de Gumbel fit (eg, 25,50,100)
        S1_GumbelWS : ndarray
            Waterstanden welke horen bij de terugkeertijden welke horen bij de Gumbel fit (eg, 25,50,100)
        S1_GumbelWS_line : ndarray
            Waterstanden welke horen bij de terugkeertijden welke horen bij de Gumbel fit vallend binnend het venster
        S1_r_squared : float
            r**2 tussen de Gumbel fit en gebeurtenenissed binnen je venster
        
        ## STEP2 Gumbel Fit voor T10 voor alleen zomerbuien
        S2_T : ndarray
            Terugkeertijden van jaarmaxima gebeurtenissen genomen over het groeiseizoen
        S2_ws_srt : ndarray
            Waterstanden welke horen bij de terugkeertijden van jaarmaxima gebeurtenissen genomen over het groeiseizoen
        
        ## STEP4 Gewogen gemiddelde voor T10 voor alleen de zomerbuien
        S4_Tarray_sel_jaar : ndarray
            Terugkeertijden van jaarmaxima gebeurtenissen genomen over het groeiseizoen welke meegenomen zijn voor bepaling gewogen gemiddeld
        S4_WSarray_sel_jaar : ndarray
            Waterstanden welke horen bij terugkeertijden van jaarmaxima gebeurtenissen welke meegenomen zijn voor bepaling gewogen gemiddeld
        TOI : int
            Terugkeertijd of interest voor bepaling van het gewogen gemiddelde (normaal gesproken 10)
        S4_WSarray_TOI_jaar : float
            Waterstand welke hoort bij de TOI
        
        Returns
        -------
        Maakt figuur en slaat deze in de folder waar het script gerund wordt.
        """

        ID = stat_object.stats.stap1.ID
        S1_vAmin = stat_object.stats.stap1.vAmin
        S1_vAmax = stat_object.stats.stap1.vAmax
        S1_T = stat_object.stats.stap1.T
        S1_ws_srt = stat_object.stats.stap1.ws_srt
        S1_GumbelT = stat_object.stats.stap1.GumbelT
        S1_GumbelWS = stat_object.stats.stap1.GumbelWS
        S1_GumbelWS_line = stat_object.stats.stap1.GumbelWS_line
        S1_r_squared = stat_object.stats.stap1.r_squared
        
        S2_T = stat_object.stats.stap2.T
        S2_ws_srt = stat_object.stats.stap2.ws_srt
        
        S3_Tarray_sel_jaar = stat_object.stats.stap3.Tarray_sel_jaar
        S3_WSarray_sel_jaar = stat_object.stats.stap3.WSarray_sel_jaar
        TOI3 = stat_object.stats.stap3.TOI
        S3_WSarray_TOI_jaar = stat_object.stats.stap3.WSarray_TOI_jaar          
        
        S4_Tarray_sel_jaar = stat_object.stats.stap4.Tarray_sel_jaar
        S4_WSarray_sel_jaar = stat_object.stats.stap4.WSarray_sel_jaar
        TOI = stat_object.stats.stap4.TOI
        S4_WSarray_TOI_jaar = stat_object.stats.stap4.WSarray_TOI_jaar  
        
        
        # Voor het figuur mag de TOI niet geplot worden in de gumbel fit. 
        # Kijk of TOI binnen de GumbelT bestaan en verwijder dit punt als zo.
        if TOI in S1_GumbelT:            
            idx = S1_GumbelT.index(TOI)
            S1_GumbelT = np.delete(S1_GumbelT, idx)
            S1_GumbelWS = np.delete(S1_GumbelWS, idx)       
        
        fig=plt.figure(figsize=(8,6))
        ax = fig.add_subplot(111)
        ax.set_xscale('log')

        # plot jaarmaxima gebeurtenissen over gehele jaar ## OPEN CIRKELS
        ax.plot(S1_T, S1_ws_srt, fillstyle='none', label='waterstanden (jaarmaxima)',
                color='#2E589F', linestyle='', marker='o', markeredgewidth=2, markersize=10, markerfacecoloralt='gray')

        # plot maxima gebeurtenissen over de zomerbuien ## DICHTE CIRKELS
        ax.plot(S2_T, S2_ws_srt, fillstyle='full', label='waterstanden (maxima groeiseizoen)',
               color='cornflowerblue', linestyle='', marker='o', markersize=8, markerfacecoloralt='gray', alpha=0.75)

        # plot gumble fit voor T25, T50, T100 in lijn en punt
        ax.plot(S1_T[S1_vAmin:S1_vAmax],S1_GumbelWS_line, '-', color='#3BBB75', label='Gumbel fit (venster)')
        ax.scatter(S1_GumbelT, S1_GumbelWS,s=100, marker='s', facecolors='#3BBB75', zorder=10,label='waterstand op T25, T50 en T100')

        # plot gewogen gemiddelde voor T10 op basis van geheel jaar in punt
        ax.scatter(TOI3, S3_WSarray_TOI_jaar,s=100, marker='s',facecolors='#FF5AC3', zorder=10,label='waterstand op T10 (stedelijk)') # geheel jaar

        # plot gewogen gemiddelde voor T10 op basis van geheel jaar in punt ## HALF-GEVULDE CIRKELS        
        ax.plot(S3_Tarray_sel_jaar,S3_WSarray_sel_jaar, fillstyle='right', label='gebeurtenissen voor bepaling T'+str(TOI3)+' (stedelijk)',
               color='#865FC5', linestyle='none', marker='o', markersize=8, markerfacecoloralt='white')        
        
        # plot gewogen gemiddelde voor T10 op basis van alleen zomerbuien in punt
        ax.scatter(TOI, S4_WSarray_TOI_jaar,s=100, marker='s',facecolors='#CE0002', zorder=10,label='waterstand op T10 (landelijk)') # zomerbuien               
        
        # plot gewogen gemiddelde voor T10 op basis van zomerbuien in punt ## HALF-GEVULDE CIRKELS        
        ax.plot(S4_Tarray_sel_jaar,S4_WSarray_sel_jaar, fillstyle='right', label='gebeurtenissen voor bepaling T'+str(TOI)+ ' (landelijk)',
               color='#F5AC1B', linestyle='none', marker='o', markersize=8, markerfacecoloralt='white')
                 
        # plot r**2 rechtsbovenin
        ax.text(0.975, 0.975, '$r^2$: '+str(np.round(S1_r_squared,2)),
                horizontalalignment='right',
                verticalalignment='top',
                transform=ax.transAxes)
        
        # axes settings
        ax.set_xlim(1,300)
        Ymin = min(S1_ws_srt) - 0.05
        Ymax = max(S1_ws_srt) + 0.25
        ax.set_ylim(Ymin,Ymax)
        ax.set_ylabel('Waterstand (m+NAP)')
        ax.set_xlabel('Terugkeertijd (jaren)')
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
        ax.grid(True,which='both', axis='both', linestyle='-', color= '0.75', zorder=0)

        # legend settings
        legend = plt.legend(loc='upper left',prop={'size':10})
        legend.get_frame().set_facecolor('white')
        legend.get_frame().set_alpha(1)
        legend.get_frame().set_linewidth(1)
        legend.get_frame().set_edgecolor('#8B8B8B')

        # table settings
        col_labels=['Terugkeertijd','Waterstand (m+NAP)']        
        table_vals = pd.DataFrame([['T'+str(i) for i in S1_GumbelT],S1_GumbelWS.round(decimals=2).tolist()]).T.values.tolist()        
        table_vals = [['T'+str(TOI)+' (landelijk)', str(round(S4_WSarray_TOI_jaar,2))]] + table_vals
        table_vals = [['T'+str(TOI3)+' (stedelijk)', str(round(S3_WSarray_TOI_jaar,2))]] + table_vals

        # the rectangle is where I want to place the table
        the_table = ax.table(cellText=table_vals,
                          colWidths = [0.4]*2,
                          colLabels=col_labels,
                          loc='center right',bbox=[1.01, 0.0, 0.5, 0.5])                

        plt.title('location: '+ID)
        plt.tight_layout()
        if out_folder=='none':
            plt.show()
        else:
            out_folder_png = os.path.join(out_folder,'png')
            out_folder_png_file = os.path.join(out_folder_png, 'waterstandstatistiek_loc_'+ID+'.png')
            plt.savefig(out_folder_png_file, bbox_inches='tight', dpi=90, pad_inches=0.2) 
            plt.clf()
            plt.close(fig)
    
    def _gewogenGemiddelde(self, df_enkel, N=109,TOI=10):
        """
        """
        # lengte van array van waterstanden
        Nws = df_enkel.size        
        T, Fijaar = self._calcFi_T(Nws,N)
        
        # kijk eerst of de input dataframe een multi-column heeft en reduceerd indien ja
        try:
            df_enkel = df_enkel.iloc[:,0]
        except:
            pass

        if len(df_enkel.name) == 2:
            enkel_location = df_enkel.name[0]
            enkel_parameter = df_enkel.name[1]
        else: 
            enkel_location = enkel_parameter = 'niet gedefinieerd'          
        
        df_enkel_sorted = df_enkel.sort_values(ascending=False)        
        
        # bepaal gewogen gemiddelde voor parameter waarde voor terugkeertijd of interest op basis van 4 dichtsbijzijnde waarden
        WSarray_TOI, WSarray_sel, Tarray_sel = self._gewogenGemiddelde4TOI(T, df_enkel_sorted.values, TOI)
        return WSarray_TOI, WSarray_sel, Tarray_sel
    
    def _enkeleGumbelFit(self, df_enkel, N, vensterArray=[0,9], GumbelT=[25,50,100.0], Ggi=0.44, GgN=0.12):
        """
        TOI : int
            Terugkeertijd of interest
            
        N : int
            Pas op: de bepaling voor over hoeveel jaren het aantal plotposities berekend moet worden, 
            moet gedaan worden op basis van het aantal mogelijke jaren in de gebeurtenissen reeks. 
            Bepaal dit op een gebeurtenissenreeks alvorens een periode filter is toegepast. 
        """        

        # lengte van array van waterstanden
        Nws = df_enkel.size  
        if Nws > N:
            #set_trace()
            self.__errors__.gebeurtenissenError()        
        T, Fijaar = self._calcFi_T(Nws, N, Ggi, GgN)        
        

        # kijk eerst of de input dataframe een multi-column heeft en reduceerd indien ja
        try:
            df_enkel = df_enkel.iloc[:,0]
        except:
            pass

        if len(df_enkel.name) == 2:
            enkel_location = df_enkel.name[0]
            enkel_parameter = df_enkel.name[1]
        else: 
            enkel_location = enkel_parameter = 'niet gedefinieerd'   
            
        # get min/max van venster array
        vAmin = vensterArray[0]
        vAmax = vensterArray[1]
        
        df_enkel_sorted = df_enkel.sort_values(ascending=False)                
        
        # bepaal de terugkeertijden waarden en afgeleiden op basis van venster array voor Gumbel fit
        dT = np.log(T[vAmin]) - np.log(T[vAmax])
        dh = df_enkel_sorted.iloc[vAmin] - df_enkel_sorted.iloc[vAmax]
        dTpunt = -np.log(T[vAmax])
        dhpunt = (dh/dT)*dTpunt
        a0 = [dh/dT, df_enkel_sorted.iloc[vAmax]+dhpunt]
        
        # apply de gumbel fit functie
        a_out, cov_x, infodict, msg, flag = optimize.leastsq(self._igumbelFuncFit, a0, args = (df_enkel_sorted.iloc[vAmin:vAmax], 1, T[vAmin:vAmax]), full_output=True)
        # a_out holds the paramters of the fitted line [a and b]
        a = a_out[0]
        b = a_out[1]
        # get the locations on the fit corresponding to the required return periods
        GumbelWS = self._igumbelFunc(a_out, 1, np.array(GumbelT))
        # get the locations on the fit corresponding to the slice on all Ts using the vensterArray
        GumbelWS_lineVensterArray = self._igumbelFunc(a_out, 1, T[vAmin:vAmax])
        #set_trace()
        
        xs = df_enkel_sorted.iloc[vAmin:vAmax]
        ys = GumbelWS_lineVensterArray
        r_squared = self._r_squared(xs,ys)
        
        # schrijf weg naar gumbel dataframe
        T_columns = ['T'+str(int(item)) for item in GumbelT]
        df_gumbel = pd.DataFrame(index=[enkel_location], columns=T_columns+['solution'])
        df_gumbel.index.name='location'

        df_gumbel.loc[enkel_location][T_columns] = GumbelWS
        if flag > 4:
            df_gumbel.loc[enkel_location]['solution'] = 'no solution'
        else:
            df_gumbel.loc[enkel_location]['solution'] = 'optimization succesful'  
            
        return (vAmin,vAmax, df_enkel_sorted.values, T, GumbelT, GumbelWS, enkel_location, a,b, GumbelWS_lineVensterArray,r_squared), df_gumbel

    def AfleidingParameters(self, df_enkel, N, vensterArray, GumbelT, TOI, startMMdd=(1,1), endMMdd=(12,31), 
                            jaarmax_as='date', Ggi=0.44, GgN=0.12):
        """
        
        """
        self = self._initNamespaces()
        # STEP 1 Gumbel Fit voor geselecteerde T's [eg. T10,T25,T50,T100] voor alle buien
        S1_param4fig, S1_df_gumbel = self._enkeleGumbelFit(df_enkel, N=N, vensterArray=vensterArray, GumbelT=GumbelT,Ggi=Ggi,GgN=GgN)
        self.stats.stap1.vAmin = S1_param4fig[0]
        self.stats.stap1.vAmax = S1_param4fig[1]
        self.stats.stap1.ws_srt = S1_param4fig[2]
        self.stats.stap1.T = S1_param4fig[3]
        self.stats.stap1.GumbelT = S1_param4fig[4]
        self.stats.stap1.GumbelWS = S1_param4fig[5]
        self.stats.stap1.ID = S1_param4fig[6]
        self.stats.stap1.a = S1_param4fig[7]
        self.stats.stap1.b = S1_param4fig[8]
        self.stats.stap1.GumbelWS_line = S1_param4fig[9]
        self.stats.stap1.r_squared = S1_param4fig[10]
        
        # STEP 2 Gumbel Fit voor T10 voor alleen zomerbuien
        df_enkel_groei = his_reader.SelectPeriodeWaardenArray(df_enkel, startMMdd=startMMdd, endMMdd=endMMdd, jaarmax_as=jaarmax_as)
        S2_param4fig, S2_df_gumbel = self._enkeleGumbelFit(df_enkel_groei, N=N, GumbelT=[TOI], vensterArray = vensterArray,Ggi=Ggi,GgN=GgN)
        self.stats.stap2.ws_srt = S2_param4fig[2]
        self.stats.stap2.T = S2_param4fig[3]
        self.stats.stap2.GumbelT = S2_param4fig[4]
        self.stats.stap2.GumbelWS = S2_param4fig[5]
        self.stats.stap2.startMMdd = startMMdd
        self.stats.stap2.endMMdd = endMMdd
        
        # STEP 3 Gewogen gemiddelde voor T10 voor alle buien
        self.stats.stap3.WSarray_TOI_jaar, self.stats.stap3.WSarray_sel_jaar, self.stats.stap3.Tarray_sel_jaar = self._gewogenGemiddelde(df_enkel, N=N, TOI=TOI)
        self.stats.stap3.TOI = TOI
        self.stats.stap3.startMMdd = startMMdd
        self.stats.stap3.endMMdd = endMMdd        

        # STEP 4 Gewogen gemiddelde voor T10 voor alleen de zomerbuien
        self.stats.stap4.WSarray_TOI_jaar, self.stats.stap4.WSarray_sel_jaar, self.stats.stap4.Tarray_sel_jaar = self._gewogenGemiddelde(df_enkel_groei, N=N, TOI=TOI)
        self.stats.stap4.TOI = TOI
        self.stats.stap4.startMMdd = startMMdd
        self.stats.stap4.endMMdd = endMMdd
        
        return self
        
    def EnsembleRunner(self, out_folder, his_file, shp_file, shp_key='nodeID', parameter='auto',
                       startMMdd=(5,15), endMMdd=(10,15), vensterArray=[0,10], GumbelT=[10,25,50,100], 
                       TOI=10, draw_plot=True, write_table=True):
        """
        Lees SOBEK weggeschreven variablene uit en schrijf de afgeleide waterstanden 
        toebehorend aan gebruikers gedefineerde terugkeertijden naar zowel csv/dbf 
        en png. Voor het afleiden van de waterstanden welke horen bij de verschillende
        terugkeertijden kan gebruik gemaakt worden van plotposities/Gumbel fit en een 
        gewogen gemiddelde.

        Input
        -----
        out_folder : path
            een folder waarin de csv en shp/dbf bestand worden weggeschreven.
            voor de figuren moet er ook een folder binnen de out_folder bestaan 
            getiteld 'png'
        his_file : path
            absoluut pad naar een his-file welke ingelezen moet gaan worden
        shp_file : path
            absoluut pad naar een shp-file van waaruit de locations gelinkt moeten 
            worden
        shp_key : str
            string met de naam van de kolom in de shp-file waarin de locations staan
        parameter : str
            naam van de parameter in het his-file welke gebruikt moet gaan worden. 
            Default is 'auto'. In dit geval zal de eerste parameter gebruikt worden 
            vanuit het his-file
        startMMdd : tuple
            Tuple in het formaat (maand,dag). Bepaling van de start datum van de 
            periode. Genoemde datum is inclusief. Voorbeeld (5,15) staat voor 
            maand 5, dag 15.
        endMMdd : tuple
            Tuple in het formaat (maand,dag). Bepaling van de eind datum van de 
            periode. Genoemde datum is inclusief. Voorbeeld (10,15) staat voor 
            maand 10, dag 15.
        vensterArray : array
            het venster welke gebruikt wordt als filter om de gebeurtenissen mee te 
            nemen voor de bepaling van de Gumbel fit. Het venster is een array van 
            twee waarden, vaak wordt [0,10] gekozen, waar 0 overeenkomt met de meest 
            extreme waarde en 10 de op 10 na meest extreme waarde. In ander woorden, 
            in dit geval is het venster de 10 meest extreme waarden.
        GumbelT : array
            een array met terugkeertijden welke meegenomen moet worden als locations 
            waarover een de waterstand moet bepaald worden aan de hand van de 
            Gumbel fit. Een vaak gebruikte array is [10,25,50,100] wat overeenkomt met 
            de T10, T25, T50 en T100
        TOI : int
            Waarde met de terugkeertijd of interest voor het bepalen van het gewogen 
            gemiddelde. Dit wordt gebruikt voor een terugkeertijd welke tenminste 2 
            gebeurtenissen en 2 gebeurtenissen na zicht heeft. Gewoonlijk kan dit 
            gebruikt worden voor het bepalen van de T10, welke soms nog in de 
            'knik' ligt.
        draw_plot : boolean
            Aangeven of je de figuren ook wilt wegschrijven naar png's. Indien zo, 
            moet er een 'png' folder bestaan binnen de out_folder. Voor het figuur 
            zal er gekeken worden of de TOI in de GumbelT bevindt. Indien dit zo is, 
            zal alleen de terugkeertijd van de TOI ingetekend worden.
        write_table : boolean
            Aangeven of je de waarden ook wilt wegschrijven naar tabellen. De tabellen 
            welke worden weggeschreven zijn 1 csv bestand en 1 shp/dbf bestand.

        Return
        ------
        geeft geen return binnen python
        """
        print ('read his-file')
        # create frame table
        his_object = his_reader.LeesMetadata(his_file)
        stat_init = self._initNamespaces()
        stat_table1 = pd.DataFrame()#
        stat_table2 = pd.DataFrame()#self._initTabellen(stat_init,shp_key)
        print ('read shp-file')
        locations_df = gpd.read_file(shp_file)
        
        if parameter=='auto':
            parameter = his_object.hisFile.variabeleInfo.variabelen[0]
        print ('create png-file')
        # iterate locationlist and create item in table
        pbar = tqdm(locations_df[shp_key].tolist())
        for location in pbar:
            pbar.set_description("location %s" % location.ljust(21))
            df_enkel = his_object.EnkeleWaardenArray(location, parameter, startMMdd=(1, 1), endMMdd=(12, 31), jaarmax_as='date')
            if df_enkel.size > his_object.hisFile.tijdstapInfo.N:
                #set_trace()
                self.__errors__.gebeurtenissenError()
            stat_object = self.AfleidingParameters(df_enkel, his_object.hisFile.tijdstapInfo.N, vensterArray, GumbelT, TOI, startMMdd, endMMdd, jaarmax_as='date')
            if draw_plot == True:
                self.PlotFiguur(stat_object, out_folder)
            if write_table ==True:
                stat_tables = self.SchrijfTabellen(stat_object, shp_key=shp_key, init_tabel=False)
                #set_trace()
                stat_table1 = stat_table1.append(stat_tables.stats.df_table1)
                stat_table2 = stat_table2.append(stat_tables.stats.df_table2)
        print ('save csv-file')
        # join on input shapefile and write tables
        # save table 1
        stat_table1.to_csv(os.path.join(out_folder,'waterstandstatistiek_output_1.csv'),index=False)
        print ('save shp-file')
        # prepare table 2 and save
        # set_trace()
        locations_df = locations_df.set_index(shp_key).join(stat_table2.set_index(shp_key)).reset_index()
        locations_df.to_file(os.path.join(out_folder,'waterstandstatistiek_output_2.shp'), driver='ESRI Shapefile')
        stat_table2 = locations_df
        print ('done')

    def _initTabellen(self,stat_object,shp_key):
        """
        Initieer tabelstructuur
        
        Parameters
        ----------
        shp_key : str
            key kolom vanuit de shp/dbf bestand
        """
        stat_object.stats.df_table1 = pd.DataFrame(columns=[shp_key,'H','T','seizoen','methode'])
        stat_object.stats.df_table2 = pd.DataFrame(columns=[shp_key,'T10_LANDELIJK','T10_STEDELIJK',
                                                     'T25','T50','T100','r^2','a','b'])
        return stat_object
        
      
    def SchrijfTabellen(self,stat_object,shp_key,decimals=3, init_tabel=True):
        """
        Schrijf afgeleide parameters naar tabel
        
        Parameters
        ----------
        shp_key : str
            key kolom vanuit de shp/dbf bestand
        decimals : int
            aantal decimalen om af te ronden
        init_tabel : boolean
            parameter om aan te geven of de _initTabellen aangeroepen moet worden. Moet op True staan wanneer functie 
            buiten de SingleRunner functie aangeroepen wordt.
        """
        stat_object._initTabellen(stat_object,shp_key)
        # if init_tabel==True:
            # stat_object._initTabellen(stat_object,shp_key)
        # else:
            # stat_empty_init = self._initNamespaces()
            # stat_empty_table = self._initTabellen(stat_empty_init, shp_key)
            # stat_object.stats.df_table1 = types.SimpleNamespace()
            # stat_object.stats.df_table2 = types.SimpleNamespace()
            # try:
                # stat_object.stats.df_table1 = init_tabel.df_table1
                # stat_object.stats.df_table2 = init_tabel.df_table2
            # except Exception as e:
                # print(e)
                # set_trace()            
        
        table1_list = []
        table2_list = []
        # gewogen gemiddelde voor jaarmaxima alle buien (# STEP 3)
        table1_list.append({shp_key:stat_object.stats.stap1.ID,
                   'H':round(stat_object.stats.stap3.WSarray_TOI_jaar,decimals),
                   'T':stat_object.stats.stap3.TOI,
                   'seizoen':'geheel jaar',
                   'methode':'gewogen gemiddelde'})

        # gewogen gemiddelde voor maxima zomerbuien (# STEP 4)
        table1_list.append({shp_key:stat_object.stats.stap1.ID,
                   'H':round(stat_object.stats.stap4.WSarray_TOI_jaar,decimals),
                   'T':stat_object.stats.stap4.TOI,
                   'seizoen':'alleen zomer',
                   'methode':'gewogen gemiddelde'})

        # Gumbel fit voor enkel maxima zomerbuien (# STEP 2)
        table1_list.append({shp_key:stat_object.stats.stap1.ID,
                   'H':round(stat_object.stats.stap2.GumbelWS[0],decimals),
                   'T':stat_object.stats.stap2.GumbelT[0],
                   'seizoen':'alleen zomer',
                   'methode':'Gumbel fit'})

        # Gumbel fit voor jaarmaxima alle buien (# STEP 1)
        for idx, T in enumerate(stat_object.stats.stap1.GumbelT):
            item = {shp_key:stat_object.stats.stap1.ID,
                    'H':round(stat_object.stats.stap1.GumbelWS[idx],decimals),
                    'T':T,
                    'seizoen':'geheel jaar',
                    'methode':'Gumbel fit'}    
            table1_list.append(item)

        
        # append to dataframe 1
        stat_object.stats.df_table1 = stat_object.stats.df_table1.append(table1_list, ignore_index=True)
        
        try:
            table2_list.append({shp_key:stat_object.stats.stap1.ID,
                                'T10_LANDELIJK':stat_object.stats.df_table1[(stat_object.stats.df_table1['T'] == 10) & 
                                                (stat_object.stats.df_table1['methode'] == 'gewogen gemiddelde') & 
                                                (stat_object.stats.df_table1['seizoen'] == 'alleen zomer') &
                                                (stat_object.stats.df_table1[shp_key] == stat_object.stats.stap1.ID)]['H'].values[0],
                                'T10_STEDELIJK':stat_object.stats.df_table1[(stat_object.stats.df_table1['T'] == 10) & 
                                                (stat_object.stats.df_table1['methode'] == 'gewogen gemiddelde') & 
                                                (stat_object.stats.df_table1['seizoen'] == 'geheel jaar') &
                                                (stat_object.stats.df_table1[shp_key] == stat_object.stats.stap1.ID)]['H'].values[0],                                            
                                'T25':stat_object.stats.df_table1[(stat_object.stats.df_table1['T'] == 25) & 
                                                (stat_object.stats.df_table1[shp_key] == stat_object.stats.stap1.ID)]['H'].values[0],
                                'T50':stat_object.stats.df_table1[(stat_object.stats.df_table1['T'] == 50) & 
                                                (stat_object.stats.df_table1[shp_key] == stat_object.stats.stap1.ID)]['H'].values[0],
                                'T100':stat_object.stats.df_table1[(stat_object.stats.df_table1['T'] == 100) & 
                                                (stat_object.stats.df_table1[shp_key] == stat_object.stats.stap1.ID)]['H'].values[0],
                                'r^2':round(stat_object.stats.stap1.r_squared,decimals),
                                'a':round(stat_object.stats.stap1.a,decimals),
                                'b':round(stat_object.stats.stap1.b,decimals)})    
        except (ValueError,IndexError):
            self.__errors__.schrijfTabelError() 
        # append to dataframe 2
        stat_object.stats.df_table2 = stat_object.stats.df_table2.append(table2_list, ignore_index=True)
        if init_tabel==True:
            return stat_object.stats.df_table1, stat_object.stats.df_table2
        else:        
            return stat_object