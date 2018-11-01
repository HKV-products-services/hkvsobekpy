# hkvsobekpy
package for water-statistics and plausibility checker using his- and bui-files. 
A his-file is a binary file object used for read/write in- and output in SOBEK. The file contains parameters, locations and timesteps. Since it is a binary format it always has been difficult to use this files directly in explorative analysis or automated workflows. This package provide options to work with these his-files within Python.

A bui-file is a text file object used for storing precipitation events for multiple locations. The file is also used for input in SOBEK models. This function provide options to parse this bui-file into a Pandas DataFrame.

Next to reading the his and bui-files, this package also contains some modules for applying statistics to obtain the return period for T10, T25, T50 and T100. This information is useful to compare against actual waterlevels.

# installation
install using pypip:

`pip install hkvsobekpy`

# dependencies
hkvsobekpy will install `tqdm` and `fire` as dependencies. The other required packgages are not included since Windows cannot compile these packages from source. These are:
- numpy
- scipy
- pandas
- GDAL
- shapely
- pyproj
- Fiona
- geopandas
- fire
- tqdm

Use `pip install <package>` to install or go to https://www.lfd.uci.edu/~gohlke/pythonlibs to download these packages on Windows (and use `pip install path/to/package.whl` to install the package).

# usage package
Import the package and define path to his-file

    import hkvsobekpy as hkv
    his_file = r'input_data\waterstand-statistiek\CALCPNT.HIS'
    
Metadata of the his-file is read first and using this metadata block subequent functions can be applied.

    calcpnt = hkv.read_his.ReadMetadata(his_file)
    
Such as the functions to get the locations, timesteps and parameters:
    
    locaties = calcpnt.GetLocations()
    tijdstappen = calcpnt.GetTimestamps()
    parameters = calcpnt.GetParameters()

    print("""\
    first 5 locations:     {0}
    last 2 timesteps:      {1}
    all parameters:        {2}""".format(locaties[0:5],
                                         tijdstappen[-2:],
                                         parameters))

    first 5 locations:     ['1', '126', '11', '8', '14']
    last 2 timesteps:      [datetime.datetime(2014, 10, 18, 15, 20), datetime.datetime(2014, 12, 12, 10, 20)]
    all parameters:        ['Waterlevel max. (m A', 'Waterdepth max. (m) ']
    
To read the whole his-file into a pandas Dataframe use 

    df = calcpnt.DataFrame()

To read a single timeseries use:

    df = calcpnt.EnkeleWaardenArray(locaties[0],
        parameters[0],
        startMMdd=(1, 1),
        endMMdd=(12, 31),
        jaarmax_as='date'))
    df.plot(legend=True)

![alt text](https://github.com/HKV-products-services/hkvsobekpy/blob/master/img/waterlevel.png "single timeseries using location and parameter")

For a notebook with the example as kentioned above, click here http://nbviewer.jupyter.org/github/HKV-products-services/hkvsobekpy/blob/master/notebook/revisit%20lezen%20HIS%20file.ipynb

See the jupyter notebook 'waterstand statistiek.ipynb' in the notebook folder for more usage examples. For example how to get the return periods T25, T50 and T100 using a Gumbel function fit,  while the T10 is computed using a weigted average of the four nearest events

![alt text](https://github.com/HKV-products-services/hkvsobekpy/blob/master/img/stats.png "waterlevel statistics")

# contact
We at HKV provide expert advice and research in the field of water and safety. Using `hkvsobekpy` we optimize our research and provide options to derive threshold levels to be used in custom-build operational apps and dashboards for river, coasts and deltas providing early-warnings and forecasts for risk and disaster management.

Interested? Head to https://www.hkv.nl/en/ or drop me a line at m.vanhoek [at] hkv.nl
