
from builtins import zip
import sys
import os
import re
import anuga.utilities.log as log
import csv

defaultOutputFile ='timing.csv'
timingDelimiter = log.TimingDelimiter

def analyse_log(path, output_file, log_file='anuga.log'):
    """
    path - the directory to look for log files in.
    output_file - the file produced, with meta log info
    log_file - the name of the log files
     
    Read log files and create a meta log file
    
    WARNING: if a 
    """
    log_pairs = build_log_info(path, log_file)
    if log_pairs is not None:
        write_meta_log(log_pairs, output_file)

def build_log_info(path, log_file):
    log_pairs = []
    for (path, dirs, files) in os.walk(path):    
        for file in files:
            if log_file in file: 
                dictResults = {}
                file_handle = open(os.path.join(path,file))
                for line in file_handle:
                    if line.find(timingDelimiter)>-1:
                        key_value = line.split(timingDelimiter)[1]
                        # FIXME remove the magic comma
                        key_value_list = key_value.split(',')
                        key = key_value_list[0].strip()
                        value =  key_value_list[1].strip()
                        if key in dictResults:
                            print("WARNING. Multiple values ignored.")
                        dictResults[key] = value
                        #print "key", key 
                        #print "value", value
                log_pairs.append(dictResults)
                file_handle.close()
    return log_pairs
  
def write_meta_log(log_pairs, output_file):
    """Write the info from the log files to a file"""
    
    all_keys = {} # values aren't needed, but are there
    for log_p in log_pairs:
        all_keys.update(log_p)
    sorted_all_keys = sorted(all_keys.keys()) # sort the keys alphabetacally
    
    han = open(output_file, 'w')
    writer = csv.DictWriter(han, delimiter=',', fieldnames=sorted_all_keys,
                        extrasaction='ignore')
    # Title 
    writer.writerow(dict(list(zip(sorted_all_keys, sorted_all_keys))))
    
    for pair in log_pairs: # Write the main body
        writer.writerow(pair)
        
    han.close()

####################################################
if __name__ == '__main__':
    
    path = sys.argv[1]

    if len(sys.argv) < 3:
        outputFile = open(defaultOutputFile, "a")
    else:
        outputFile = open(sys.argv[2], "a")

    analyse_log(path, outputFile)

    
 
