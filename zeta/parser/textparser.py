'''
Created on Oct 30, 2021

@author: mzwier
'''

from collections import deque
import re

from . import ParseError

def ffloat(text):
    '''Convert a Fortran-formatted float string to a float. Also trims
    whitespace. Not exactly fast when called a few million times.'''
    return float(text.strip().replace('D', 'E').replace('d', 'e'))

class FailedMatchError(ParseError):
    def __init__(self, text, predicate, result):        
        self.text = text
        self.predicate = predicate
        self.result = result
    

class TextPredicateResult:
    '''A record of a predicate test. This must at least evaluate to True/False
    in a boolean context, and may carry additional information about the match
    (such as a regexp match result, line and/or character position of the
    match, etc.)'''
    def __init__(self, predicate, text, matched, result, context):
        
        # The predicate used to test
        self.predicate = predicate

        # The text scanned
        self.text = text

        # Whether the predicate evaluated to true
        self.matched = matched

        # Additional result information
        self.result = result

        # Additional context        
        self.context = context
                
    def __bool__(self):
        return self.matched

class RegexpPredicateResult(TextPredicateResult):
    def __init__(self, predicate, text, result, context):
        super().__init__(predicate, text, bool(result), result, context)
        
    def __getitem__(self, key):
        return self.result[key]
    
    def __getattr__(self, name):
        '''Forward everything not otherwise covered to the regexp match object'''
        return getattr(self.result, name)

class TextPredicate:   
    def __call__(self, text, context = None):
        '''Test this predicate on `text`. Returns an instance of 
        `record_type`, a subclass of `TextPredicateResult`, 
        tracking the success/failure of the match and any associated 
        context (like filename or line number) or match data (like 
        position of match or regexp match object).'''
        
        raise NotImplementedError         

        
class ContainsText(TextPredicate):
    def __init__(self, match_text, case_sensitive=True):
        super().__init__()
        self.match_text = match_text
        self.case_sensitive = case_sensitive
        if case_sensitive:
            self._norm_match_text = match_text 
        else:
            # Case insensitive match
            self._norm_match_text = match_text.casefold()
        
    def __call__(self, text, context=None):
        if self.case_sensitive:
            norm_text = text
        else:
            norm_text = text.casefold()
        loc = norm_text.find(self._norm_match_text)
        return TextPredicateResult(self, text, (loc > -1), loc, context)
        
class RegexpMatch(TextPredicate):

    def __init__(self, pattern_or_text):
        super().__init__()
        if isinstance(pattern_or_text, str):
            self.regexp = re.compile(pattern_or_text)
        else:
            self.regexp = pattern_or_text
        
    def __call__(self, text, context=None):
        m = self.regexp.search(text)
        return RegexpPredicateResult(self, text, m, context)

class InvertPredicate(TextPredicate):
    def __init__(self, predicate):
        super().__init__()
        self.predicate = predicate
    
    def __call__(self, text, context=None):
        return not self.predicate(text, context)

# A few common predicates    
whitespace_only = RegexpMatch('^\s*$')
not_whitespace_only = InvertPredicate(whitespace_only)
empty_line = RegexpMatch('^$')
nonempty_line = InvertPredicate(empty_line)

class TextFileParser:
    '''A line-buffering text file parser. Stream-compatible. Lookaheads cache
    lines. The last line read and the last predicate result (if applicable) are
    available.
    
    Principles: read-and-then-act, not act-then-read, except where obvious
    '''
    def __init__(self, textfile):
        self.textfile = textfile
        
        # Buffer for lookaheads
        # Deque of (linenum, line) pairs
        self._buffer = deque()
        
        # The last line read from the stream/buffer
        self.line = None 
        
        # The current (1-based) line number
        self.linenum = 0
        
        # The predicate result for current_line, if any
        self.presult = None 
        
        
        # Implementation details
        self._stream_last_read = None # line most recently read from the stream
        self._stream_linenum = 0 # actual lines read from the stream
         
        # Vomit lots of information
        self._debug_firehose = False
        
    
    def _stream_readline(self):
        line = self.textfile.readline()
        self._stream_last_read = line
        self._stream_linenum += 1 # one past the end, at the end
        
        if self._debug_firehose:
            print('<<< line {:d} {!r}'.format(self._stream_linenum, self._stream_last_read))
    
        
    def _apply_predicate(self, predicate):
        result = predicate(self.line, context={'stream': self.textfile,
                                               'linenum': self.linenum})
        self.presult = result
        
    def buffer_is_empty(self):
        return not bool(self._buffer)
    
    def fast_forward_linenum(self, nlines):
        '''Manually advance the line number by `nlines`, for cases where
        another routine reads directly from `textfile`.'''
        assert self.buffer_is_empty()
        self._stream_last_read = None
        self._stream_linenum += nlines
        self.line = None
        self.linenum = self._stream_linenum

    def nextline(self):
        '''Read the next line, which comes from the buffer if any lines are
        buffered, otherwise from the underlying stream. Returns True if a 
        a line is read successfully, False otherwise (e.g. end-of-file).'''
        
        if self._buffer:
            linenum, line = self._buffer.popleft()
            if self._debug_firehose:
                print('<< unbuffered line {:d} {!r}'.format(linenum, line))
                print('>> (nextline) current buffer {}'.format(self._buffer))
        else:
            self._stream_readline()
            line = self._stream_last_read
            linenum = self._stream_linenum
        
        self.line = line
        self.linenum = linenum
                
        return bool(line)
    
    def readline(self):
        '''Compatibility function for line-oriented text streams. Reads and
        returns one line, or the empty string at end-of-file. Not sure this
        should even exist.'''
        self.nextline()
        return self.line
        
    def testp(self, predicate):
        self._apply_predicate(predicate)
        return bool(self.presult)
            
    def assertp(self, predicate):
        '''Assert that the current line matches `predicate`, raising 
        `FailedMatchError` if not. Does not read an additional line.'''
        
        if self.linenum == 0:
            raise ParseError('no lines have been read')
        self._apply_predicate(predicate)
        if not self.presult:
            raise FailedMatchError(self.line, predicate, self.presult)

    def testp_within_next(self, predicate, nlines):
        '''Test whether the given predicate matches within the next `nlines`
        lines'''
        buffer = []
        try:
            for _ in range(nlines):
                self.nextline()
                buffer.append((self.linenum, self.line))
                if self.testp(predicate):
                    return True
            return False
        finally:
            self._buffer.extendleft(reversed(buffer))
            if self._debug_firehose:
                print('>> (testp_within_next) current buffer: {}'.format(self._buffer))                
            
    def assertp_within_next(self, predicate, nlines):
        if not self.testp_within_next(predicate, nlines):
            raise ParseError('could not match within {:d} lines'.format(nlines))

    def read_and_assertp(self, predicate):
        '''Read one line and assert that it matches `predicate`, raising
        `FailedMatchError` if it does not.'''
        self.nextline()
        if not self.testp(predicate):
            raise FailedMatchError(self.line, predicate, self.presult)

    
    def read_until_match(self, predicate):
        '''A loop predicate that reads one line, tests against `predicate`, and
        returns False when that predicate matches. This supports loops of the form
            while read_until_match(predicate):
                do_something_with_the_nonmatching_line(s)
        '''
        self.nextline()
        if not self.line:
            raise ParseError('end of file before match')
        return (not self.testp(predicate))

    def read_while_match(self, predicate):
        return self.read_until_match(InvertPredicate(predicate))
            
        
    def discard_until_match(self, predicate, maxlines=None):
        '''Read and discard lines until `predicate` is matched, after which the
        current line is the line that matched the predicate. If `maxlines` is
        given, the match must occur within `maxlines` lines or else an error
        is raised.'''
        
        n = 0
        while self.read_until_match(predicate):
            n += 1
            if maxlines and n > maxlines:
                raise ParseError('did not match within {:d} lines'.format(maxlines))
            elif self.presult:
                return
        
    def discard_while_match(self, predicate, maxlines=None):
        '''Read and discard lines while `predicate` matches.'''
        self.discard_until_match(InvertPredicate(predicate), maxlines=maxlines)
        
    def discardn(self, nlines):
        '''Read and discard `nlines` lines'''
        for _ in range(nlines):
            self.nextline()
                    
    def scan_and_dispatch(self, dispatch_table):
        '''Given a list of (regexp, fn) pairs, scan until at least one regexp
        matches, then dispatch to function as fn(self)'''
        self.nextline()
        while self.line:
            for predicate, fn in dispatch_table:
                if self.testp(predicate):
                    return fn(self)
            self.nextline()
            