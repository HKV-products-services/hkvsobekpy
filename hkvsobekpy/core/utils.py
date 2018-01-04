from difflib import SequenceMatcher
import unicodedata

def similar(a, b):
    return SequenceMatcher(None, a, b).ratio()

def normalize_caseless(text):
    return unicodedata.normalize("NFKD", text.casefold())

def caseless_equal(left, right):
    return normalize_caseless(left) == normalize_caseless(right)
    
def compare_df_column_his_list(df, df_column_key, his_parameters, normalize_by_unicode=True, 
                                include_simularity=False, sequence_simularity=0.82):
    """
    function to apply unicode normalization and similarity checking for two list of columns
    
    Parameters
    ----------
    df : geopandas.GeoDataFrame
        dataframe containing the input columns, normally is the source a shp or dbf file
    df_column_key : str
        name of column for usage as key
    his_parameters : list
        list of strings of locations/parameters of his-file to compare against the column 
        in the geodataframe
    normalize_by_unicode : boolean
        inlcude this option to include NFKD unicode compatibility decomposition.
        see: http://unicode.org/reports/tr15/
    include_simluratiy : boolean
        include this option to include Ratcliff/Obershelp pattern recognition
    sequence_simularity : float
        number between 0.0 and 1.0, function as threshold, where only a simularity above 
        this value is mapped
    
    Returns
    df : geopandas.GeoDataFrame
        dataframe where column matching the column key is updated with matching values 
    
    """

    for df_idx, df_parameter in enumerate(df[df_column_key]):

        for his_parameter in his_parameters:
            # option to have only a unicode normalization
            if normalize_by_unicode == True and include_simularity == False:
                if caseless_equal(his_parameter, df_parameter) == True:
                    df.loc[df_idx, df_column_key] = his_parameter
                    print('{0} changed into {1}'.format(
                        df_parameter, his_parameter))

            # option to apply only parameter similarity functionality
            if include_simularity == True and normalize_by_unicode == False:
                smty = similar(his_parameter, df_parameter)
                # print(smty)
                if sequence_simularity <= smty < 1.0:
                    df.loc[df_idx, df_column_key] = his_parameter
                    print('{0} changed into {1}'.format(
                        df_parameter, his_parameter))

            # option to apply both unicode normalization and similarity functionality
            if include_simularity == True and normalize_by_unicode == True:
                # first check unicode normalization
                if caseless_equal(his_parameter, df_parameter) == True:
                    df.loc[df_idx, df_column_key] = his_parameter
                    print('{0} changed into {1}'.format(
                        df_parameter, his_parameter))
                # if not succeed try similarity
                else:
                    
                    smty = similar(normalize_caseless(his_parameter), normalize_caseless(df_parameter))
                    #print(smty)
                    if sequence_simularity <= smty < 1.0:
                        df.loc[df_idx, df_column_key] = his_parameter
                        print('{0} changed into {1}'.format(
                            df_parameter, his_parameter))
    return df
    
def compare_df_parameter_his_parameter(df_parameter, his_parameters, normalize_by_unicode=True, 
                                       include_simularity=False, sequence_simularity=0.82):
    """
    function to apply unicode normalization and similarity checking for two list of columns
    
    Parameters
    ----------
    df_parameter : str
        string of locations/parameter from a dbf of shp file
    his_parameters : list
        list of strings of locations/parameters of his-file to compare against the df value
    normalize_by_unicode : boolean
        inlcude this option to include NFKD unicode compatibility decomposition.
        see: http://unicode.org/reports/tr15/
    include_simluratiy : boolean
        include this option to include Ratcliff/Obershelp pattern recognition
    sequence_simularity : float
        number between 0.0 and 1.0, function as threshold, where only a simularity above 
        this value is mapped
    
    Returns
    df_parameter : str
        df parameter updated with matching values of his parmater/location if any 
    
    """
    for his_parameter in his_parameters:
        if normalize_by_unicode == True and include_simularity == False:
            if caseless_equal(his_parameter, df_parameter) == True:
                print('Normalized unicode {0} matches {1}'.format(
                    df_parameter, his_parameter))            
                df_parameter = his_parameter


        # option to apply only parameter similarity functionality
        if include_simularity == True and normalize_by_unicode == False:
            smty = similar(his_parameter, df_parameter)
            # print(smty)
            if sequence_simularity <= smty < 1.0:
                print('{0} changed into {1}, since similarity is {2}'.format(
                    df_parameter, his_parameter, round(smty,2)))            
                df_parameter = his_parameter

        # option to apply both unicode normalization and similarity functionality
        if include_simularity == True and normalize_by_unicode == True:
            # first check unicode normalization
            if caseless_equal(his_parameter, df_parameter) == True:
                print('Normalized unicode {0} matches {1}'.format(
                    df_parameter, his_parameter))            
                df_parameter = his_parameter

            # if not succeed try similarity
            else:
                smty = similar(normalize_caseless(his_parameter), normalize_caseless(df_parameter))
                #print(smty)
                if sequence_simularity <= smty < 1.0:
                    print('{0} changed into {1}, since similarity is {2}'.format(
                        df_parameter, his_parameter, round(smty,2)))                
                    df_parameter = his_parameter

    return df_parameter    