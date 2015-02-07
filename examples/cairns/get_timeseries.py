"""
Generate time series of nominated "gauges" read from project.gauge_filename. This 
is done by first running sww2csv_gauges on two different directories to make 
'csv' files. Then running csv2timeseries_graphs detailing the two directories 
containing the csv file and produces one set of graphs in the 'output_dir' containing
the details at the gauges for both these sww files.

Note, this script will only work if pylab is installed on the platform
"""

from os import sep
import project
import anuga

try:
    anuga.sww2csv_gauges('cairns_slide.sww',
                project.gauge_filename,
                quantities=['stage','speed','depth','elevation'],
                verbose=True)
except:
    print 'Failed to process cairns_slide'

try:                
    anuga.sww2csv_gauges('cairns_fixed_wave.sww',
               project.gauge_filename,
               quantities=['stage', 'speed','depth','elevation'],
               verbose=True)
except:
    print 'Failed to process cairns_fixed_wave'

try: 
    import pylab
    anuga.csv2timeseries_graphs(directories_dic={'fixed_wave'+sep: ['Fixed Wave',0,0]},
                          output_dir='fixed_wave'+sep,
                          base_name='gauge_',
                          plot_numbers='',
                          quantities=['stage','speed','depth'],
                          extra_plot_name='',
                          assess_all_csv_files=True,                            
                          create_latex=False,
                          verbose=True)
except ImportError:
    #ANUGA does not rely on pylab to work 
    print 'must have pylab installed to generate plots'


try: 
    import pylab
    anuga.csv2timeseries_graphs(directories_dic={'slide'+sep: ['Slide',0, 0]},
                          output_dir='slide'+sep,
                          base_name='gauge_',
                          plot_numbers='',
                          quantities=['stage','speed','depth'],
                          extra_plot_name='',
                          assess_all_csv_files=True,                            
                          create_latex=False,
                          verbose=True)
except ImportError:
    #ANUGA does not rely on pylab to work 
    print 'must have pylab installed to generate plots'
