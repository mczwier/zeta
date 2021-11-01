'''
Created on Oct 30, 2021

@author: mzwier
'''

from collections import deque

def ffloat(text):
    '''Convert a Fortran-formatted float string to a float. Also trims
    whitespace. Not exactly fast when called a few million times.'''
    return float(text.strip().replace('D', 'E').replace('d', 'e'))
    
class ParseError(RuntimeError):
    pass

class FailedMatchError(ParseError):
    pass

class TextFileParser:
    '''A line-buffering text file parser. Stream-compatible. Lookaheads cache
    lines. Last line read and last regexp matched are stored for convenience.
    '''
    def __init__(self, textfile):
        self.textfile = textfile
        
        self.last_match = None # last match object
        self.last_line = None  # last line read 
        self._buffer = deque()
        
    def _match(self, regexp, string):
        m = regexp.search(string)
        self.last_match = m
        return m            
        
    def readline(self, strip=False):
        '''Read a line from the file, optionally stripping whitespace'''
        if self._buffer:
            line = self._buffer.popleft()
        else:
            line = self.textfile.readline()
            
        if strip is True:
            line = line.strip()
        elif strip is False:
            pass
        else:
            line = strip(line)
            
        self.last_line = line
        return line
        
    def peek(self):
        line = self.readline()
        self._buffer.append(line)
        return line
    
    def discard_buffer(self):
        self._buffer = deque()
        
    def matches(self, regexp):
        line = self.peek()
        return self._match(regexp, line)

    def read_and_match(self, regexp):
        line = self.readline()
        m = self._match(regexp, line)
        if not m:
            raise FailedMatchError('expected match to {!r}'.format(regexp))
        return m

    def discard_to_match(self, regexp):
        line = self.readline()
        while line:
            m = self._match(regexp, line)
            if m:
                return m
            line = self.readline()
                        
    def test_lookahead(self, regexp, nlines):
        '''Returns True if regexp matches within the next nline lines'''
        n = 1
        while n <= nlines:
            line = self.readline()
            if self._match(regexp, line):
                return True
            n += 1
        return False
            
    def skip_blanks(self):
        line = self.readline()
        if not line: return #eof
        while line.strip() == '':
            line = self.readline()
        self._buffer.append(line)
                
    def skip_to_blank(self):
        line = self.readline()
        if line.strip() != '':
            line = self.readline()
        self._buffer.append(line)
                
    def scan_and_dispatch(self, dispatch_table):
        '''Given a list of (regexp, fn) pairs, scan until at least one regexp
        matches, then dispatch to function as fn(match object, line)'''
        line = self.readline()
        while line:
            for regexp, fn in dispatch_table:
                if self._match(regexp, line):
                    return fn(self.last_match, line)
            line = self.readline()
            