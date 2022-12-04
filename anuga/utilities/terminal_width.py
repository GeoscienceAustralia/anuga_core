#!/usr/bin/env python
'''Function to get terminal width, Windows or Linux.'''


######
# This code was found at http://code.activestate.com/recipes/440694/
# 
# Any errors here are mine, as the code below changed the code found above.
######

def terminal_width():
    """Get the current terminal width.

    Returns the terminal width in characters.

    If the width cannot be found, return 80 as a default.
    """

    # First, try Windows.
    try:
        # fails if not Windows
        from ctypes import windll, create_string_buffer

        # stdin handle is -10
        # stdout handle is -11
        # stderr handle is -12
        h = windll.kernel32.GetStdHandle(-12)
        csbi = create_string_buffer(22)
        res = windll.kernel32.GetConsoleScreenBufferInfo(h, csbi)

        if res:
            import struct
            (bufx, bufy, curx, cury, wattr, left,
             top, right, bottom, maxx, maxy) = \
                    struct.unpack("hhhhHhhhhhh", csbi.raw)
            width = right - left + 1
        else:
            width = 80      # can't determine size - return default values
    # No, try Linux
    except:
        width = 0
        try:
            import struct, fcntl, termios

            s = struct.pack('HHHH', 0, 0, 0, 0)
            x = fcntl.ioctl(1, termios.TIOCGWINSZ, s)
            width = struct.unpack('HHHH', x)[1]
        except (IOError, ImportError):
            pass
        if width <= 0:
            try:
                width = int(os.environ['COLUMNS'])
            except:
                pass
        if width <= 0:
            width = 80      # can't determine size - return default values

    return width

if __name__ == '__main__':
    print('terminal width=%d' % terminal_width())
