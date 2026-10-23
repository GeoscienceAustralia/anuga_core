def parse_time(time = None, verbose=False, debug=False):
    """
    Time: seconds since epoch  or
    string of form '20120229'  '20120229_1210' '20120229 1210' '201202291210'
    or '20120229_121000' '20120229 121000' '20120229121000'
    """

    if time is None: return None

    if not isinstance(time, str):

        try:
            time = float(time)
            return time
        except ValueError:
            pass

    year, month, day, hour, minute, second  = 1970, 1, 1, 0, 0, 0

    try:
        year = int(time[0:4])
    except (ValueError, TypeError):
        year, month, day = 1970,1,1


    #month = int(time[4:6])

    try:
        month = int(time[4:6])
    except (ValueError, TypeError):
        month = 1

    try:
        day = int(time[6:8])
    except (ValueError, TypeError):
        day = 1

    try:
        dash = time[8:9]
        assert dash == '_' or dash == ':' or dash == '/' or dash == ' '
    except (AssertionError, IndexError):
        dash = None


    try:
        if dash is None:
            hour = int(time[8:10])
        else:
            hour = int(time[9:11])
    except (ValueError, TypeError):
        hour = 0

    try:
        if dash is None:
            minute = int(time[10:12])
        else:
            minute = int(time[11:13])
    except (ValueError, TypeError):
        minute = 0

    try:
        if dash is None:
            second = int(time[12:14])
        else:
            second = int(time[13:15])
    except (ValueError, TypeError):
        second = 0

    if debug:
        print(year, month, day, hour, minute, second)
        print('Convert to epoch')



    import datetime
    time = int((datetime.datetime(year,month,day,hour,minute,second) - datetime.datetime(1970,1,1)).total_seconds())

    if debug: print(time)

    return float(time)

def seconds_to_hhmmss(seconds):
    hours = seconds // 3600
    minutes = (seconds % 3600) // 60
    secs = seconds % 60

    parts = []
    if hours:
        parts.append(f"{hours}h")
    if minutes:
        parts.append(f"{minutes}m")
    if secs or not parts:
        parts.append(f"{secs}s")

    return ':'.join(parts)
