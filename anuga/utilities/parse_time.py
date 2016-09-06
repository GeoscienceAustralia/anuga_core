
def parse_time(time = None, verbose=False, debug=False):
    """
    Time: seconds since epoch  or 
    string of form '20120229'  '20120229_1210' '20120229 1210' '201202291210'
    """
    
    if time is None: return None
    
    if not isinstance(time, basestring):
        
        try:
            time = float(time)
            return time
        except ValueError:
            pass
    
    year, month, day, hour, minute  = 1970, 1, 1, 0, 0
    
    try:
        year = int(time[0:4])
    except:
        year, month, day = 1970,1,1
        
        
    #month = int(time[4:6])
    
    try:
        month = int(time[4:6])
    except:
        month = 1
        
    try:
        day = int(time[6:8])
    except:
        day = 1
            
    try:
        dash = time[8:9]
        assert dash == '_' or dash == ':' or dash == '/' or dash == ' '
    except:
        dash = None

            
    try:
        if dash is None:
            hour = int(time[8:10])
        else:
            hour = int(time[9:11])
    except:
        hour = 0
            
    try:
        if dash is None:
            minute = int(time[10:12])
        else:
            minute = int(time[11:13])
    except:
        minute = 0
            
    if debug:
        print year, month, day, hour, minute
        print 'Convert to epoch'

            
                
    import datetime
    time = int((datetime.datetime(year,month,day,hour,minute) - datetime.datetime(1970,1,1)).total_seconds())

    if debug: print time
    
    return float(time)