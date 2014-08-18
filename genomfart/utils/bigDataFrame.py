import gzip
import re
from cachetools import LRUCache
from bisect import bisect_left

## Compile regular expressions
int_re = re.compile(r'^-*[\d]+$')
float_re = re.compile(r'(^-*[\d]+\.[\d]+$|^-*[\d]{1}([\.]?[\d])*(E|e)-[\d]+$)')

class BigDataFrame(object):
    """ Class acting as a DataFrame for large datasets

    Not all of the dataset is loaded into memory. Rather, a cache of rows is
    maintained along with a handle and a position in the file
    """
    def __init__(self, filename, header = True, sep='\t',
                 detect_type = True, assume_uniform_types = True,
                 ignore_chars = None, maxsize=10000,
                 byte_record_interval=10000):
        """ Instantiates the big data frame

        Parameters
        ----------
        filename : str
            The file containing delimited data
        header : boolean, optional
            Whether the file has a header
        rowLabs : int, optional
            0-based column containing headers for the rows
        sep : str, optional
            The separateor of columns in the file (Can be a compiled regular expression)
        detect_type : boolean, optional
            Whether type should be guessed on loading. If false, everythin returned
            as str
        assume_uniform_types : boolean, optional
            Whether it should be assumed that the types in the first row apply to all rows.
            Only matters if detect_type is True
        ignore_chars : Set
            Characters that signal a line in the file should be skipped if it
            starts with one of those characters
        maxsize : int, optional
            The maximum number of rows to cache
        byte_record_interval : int, optional
            How often to record the starting byte of a row
        """
        self.byte_record_interval = byte_record_interval
        self.ignore_chars = ignore_chars if ignore_chars else set()
        ## Dictionary with headers for each of the columns, where values are integer
        # values
        self.colLabels = {}
        ## Dictionary with headers for each of the rows, where values are integer
        # values
        self.rowLabels = {}
        ## LRU cache holding rows of data, with row indices as keys
        self.data = LRUCache(maxsize = maxsize)
        ## Function used to convert to correct type
        self.type_func = lambda x: x[1]
        self.type_dict = {}
        if detect_type:
            if assume_uniform_types:
                self.type_func = lambda x: self.type_dict[x[0]](x[1])
            else:
                def type_func(x):
                    if float_re.match(x[1]): return float(x[1])
                    elif int_re.match(x[1]): return int(x[1])
                    else:
                        return str(x[1])
                self.type_func = type_func
        ## Function used to split lines
        if type(sep) == type(int_re):
            self.split_func = lambda x: re.split(sep, x)
        elif type(sep) == str:
            self.split_func = lambda x: x.split(sep)
        else:
            raise TypeError("Sep must be compiled regular expression or string")
        ## Get the file handle
        if filename.endswith('.gz'):
            self.handle = gzip.open(filename)
        else:
            self.handle = open(filename)
        ## Skip anything that needs to be skipped
        current_byte = self.handle.tell()
        while 1:
            current_byte = self.handle.tell()
            line = self.handle.readline().strip()
            if line[0] in self.ignore_chars: continue
            else:
                break
        self.handle.seek(current_byte)
        ## Incorporate header if necessary
        if header:
            line = self.handle.readline().strip()
            if ignore_chars:
                while 1:
                    if line[0] in self.ignore_chars:
                        line = self.handle.readline().strip()
                    else:
                        break
            line = self.split_func(line)
            for i,colname in enumerate(line):
                self.colLabels[colname] = i
        ## Create dictionary of line index -> byte
        self.row_byte_dict = {0: self.handle.tell()}                
        ## Get the types
        line = self.handle.readline().strip()
        line = self.split_func(line)
        for i, val in enumerate(line):
            if float_re.match(val):
                self.type_dict[i] = float
            elif int_re.match(val):
                self.type_dict[i] = int
            else:
                self.type_dict[i] = str
        self.handle.seek(self.row_byte_dict[0])
        ## Ordered list of the line indices that have been recorded
        self.ordered_row_inds = [0]
        ## The current line index
        self.current_line_ind = 0
    def __getitem__(self, inds):
        """ Gets an entry or a row from the data frame

        Parameters
        ----------
        inds : int or tuple
            Either the zero-based index of a row to get or a tuple of
            two numbers, where the first is the row and the second is the column,
            or a tuple of a number and a string, where the string is the name
            of the column to get

        Returns
        -------
        Either a single entry or a dictionary for the row
        """
        if type(inds) == tuple:
            inds = list(inds)
            if len(inds) > 2:
                raise TypeError("Two indices maximum")
            elif type(inds[1]) == str:
                # Get the index associated with the column
                try:
                    inds[1] = self.colLabels[inds[1]]
                except KeyError:
                    raise KeyError("%s is not a column label" % inds[1])
            if type(inds[0]) != int:
                raise TypeError("Row index must be integer")
            # Get from the cache
            row = self._fetch_from_cache(inds[0])
            return row[inds[1]]
        elif type(inds) == int:
            row = self._fetch_from_cache(inds)
            # Return the row as a dictionary
            if len(self.colLabels) > 0:
                return dict((k,row[v]) for k,v in self.colLabels.items())
            else:
                return dict((i,x) for i,x in enumerate(row))
        else:
            raise TypeError("Row index must be integer")
    def __iter__(self):
        """ Iterates through the rows of the data frame

        Returns
        -------
        Generator of dictionaries for each row
        """
        self.current_line_ind = self._goto_closest_line(0)
        header = True if len(self.colLabels) > 0 else False
        for line in self.handle:
            row = tuple(map(self.type_func, enumerate(self.split_func(line.strip()))))
            self.data[self.current_line_ind] = row
            if header:
                yield dict((k,row[v]) for k,v in self.colLabels.items())
            else:
                yield dict((i,x) for i,x in enumerate(row))            
            self.current_line_ind += 1
    def _fetch_from_cache(self, row_ind):
        """ Gets a row from the cache if present. Otherwise, puts it in the cache
        and fetches it

        Parameters
        ----------
        row_ind : int
            The index of the row

        Returns
        -------
        The row as a tuple
        """
        if row_ind in self.data:
            return self.data[row_ind]
        elif self.current_line_ind == row_ind:
            line = self.handle.readline()
            if (self.current_line_ind + 1) % self.byte_record_interval == 0:
                # Record the current byte if necessary
                self.row_byte_dict[self.current_line_ind+1] = self.handle.tell()
            if self.current_line_ind == row_ind:
                ## If we're on the line we want
                # Put the line in the cache after processing
                self.data[self.current_line_ind] = tuple(map(self.type_func,
                                                enumerate(self.split_func(line.strip()))))
                return self.data[self.current_line_ind]           
        else:
            # Send the handle to the closest line
            self.current_line_ind = self._goto_closest_line(row_ind)
            for line in self.handle:
                if (self.current_line_ind + 1) % self.byte_record_interval == 0:
                    # Record the current byte if necessary
                    self.row_byte_dict[self.current_line_ind+1] = self.handle.tell()
                if self.current_line_ind == row_ind:
                    ## If we're on the line we want
                    # Put the line in the cache after processing
                    self.data[self.current_line_ind] = tuple(map(self.type_func,
                                                    enumerate(self.split_func(line.strip()))))
                    return self.data[self.current_line_ind]
                self.current_line_ind += 1
        # If nothing else, there's an index error
        raise IndexError("Line %d is not in the file" % row_ind)
            
    def _goto_closest_line(self, line_ind):
        """ Sends the handle to the line closest to (but before) the one you
        want that is already recorded in the row_byte_dict

        Parameters
        ----------
        line_ind : int
            The line you want to go to

        Returns
        -------
        The 0-based index of the row where the handle currently sits
        """
        closest_ind_ind = max(bisect_left(self.ordered_row_inds, line_ind)-1, 0)
        # Send the handle to the location
        self.handle.seek(self.row_byte_dict[closest_ind_ind])
        self.current_line_ind = self.ordered_row_inds[closest_ind_ind]
        return self.ordered_row_inds[closest_ind_ind]
    def get_current_row_ind(self):
        """ Gets the index of the current row

        Returns
        -------
        The zero-based index of the current row
        """
        return self.current_line_ind
        
