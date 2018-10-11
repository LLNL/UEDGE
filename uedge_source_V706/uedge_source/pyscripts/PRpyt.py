"""
# Copyright (c) 1996, 1997, The Regents of the University of California.
# All rights reserved.  See Legal.htm for full text and disclaimer.

HDF basic writer class PR by David Grote, LLNL
Heavily modified from PR.py originally written by Paul Dubois, LLNL, to use
HDF files.
$Id: PRpyt.py,v 7.0 2018/02/28 18:43:48 meyer8 Exp $

This provides a nice wrapper over the pytables for reading general
self-describing data sets written out with PWpyt.PW.

"""

import pprint
import tables
import sys
import re
import cPickle
import string

_version = '0.4'

class PRError(Exception):
    pass

class PR:
    """HDF file read-access class.
Access to data in the file is available via attributes. For example, to get
the value of x from a file, do the following...

ff = PR(filename)
x = ff.x

There is also the read method
x = ff.read('x')
which is useful for names that are not usable as python attributes.
    """
    file_type = "HDF"

    """
Extra documentation:

The _cache is a dictionary for convenience. When a name is first looked up,
it is added to the cache so the next time, that value in the cache is used
instead of having to read in the data from the file again. This is meant to
save time. When the file is opened, all of the variables in the _ints and
_floats are added to the cache, as well as everything that was in the pickle
dict.  

The _names is a list of all variables in the top level of the file, i.e. 
group='/'. This includes the names in the _ints and _floats tables, the names
in the pickle dict, and all of the other names that were directly written out.

When the file is opened, a scan is done looking for hashed names, names that
were too long to include directly in the _ints and _floats tables. The
translation to the original name is stored in _hashednames. It is only used
when the tables are read into the cache.

The delimiter can be used to write complicated names and is used heavily by
the Forthon package. pytables issues complaints if '@' is used, so an option
is given in PW to specify it. PW writes the delimiter is uses to the file so
that PR can read it in.

    """

    def __del__(self):
        self.close()

    def __getattr__(self, name):
        # --- Read in the variable from the file (or get it from the cache).
        return self.read(name)

    def __init__(self, filename, group='/', verbose = 1):
        """PR(filename='', group='/',verbose=1) opens filename"""

        self._file = None
        self._cache = {}
        self._names = []
        self._hashednames = {}
        self.set_verbosity(verbose)
        self.set_group(group)
        self._fixdelimiter = re.compile('@')
        if filename:
            self.open(filename, group,mode='r')

    def __repr__(self):
        if self.is_open():
            current_mode = 'opened for reading'
        else:
            current_mode = '(PR object not open on any file)'
        return 'HDF file %s %s.' % (self.inquire_filename(), current_mode)

    __str__ = __repr__

    def check_open(self):
        "check_open(): raise exception if not open for read."
        if not self.is_open():
            raise PRError, 'PR object not open for read.'

    def close(self):
        "close(): close the file."
        if self.is_open():
            if self.inquire_verbosity():
                print "Closing HDF file ",self.inquire_filename()
            self._file.close()

            self._file = None
            self._cache = {}
            self._names = []
            self._hashednames = {}

    def inquire_filename(self):
        "inquire_filename() = name of this file."
        if self.is_open ():
            return self._file.filename
        else:
            return ''

    def inquire_ls(self,group=None):
        """inquire_ls(group=None) = list of those names 
        in the file which represent objects in the given group.
        The default groups gets all."""
        self.check_open()
        if group is None or group == '/':
          # --- Return the list of all names, including names from the ints
          # --- and floats tables, and from the pickledict.
          return self._names
        else:
          # --- Get only names in the specified group
          ll = self._file.listNodes(group,classname='Leaf')
          ll = map(lambda l:l.name,ll)
          # --- Fix the delimiter in all of the names if needed.
          if self._delimiter != '@':
              ll = map(lambda l:self._fixdelimiter.sub(self._delimiter,l),ll)
          return ll

    def inquire_mode(self):
        "inquire_mode() = mode ('r', 'w', or 'a') of this file."
        self.check_open()
        return self._file.mode

    def inquire_names(self):
        "inquire_names() = sorted list of registered names"
        return self.inquire_ls()

    def inquire_group(self):
        "inquire_group() = present HDF group"
        return self._group

    def inquire_verbosity(self):
        "inquire_verbosity() = current value of verbose flag."
        return self._verbose

    def is_open(self):
        "is_open() = true if file is open for read"
        return not (self._file == None)

    def print_names(self, file=sys.stdout):
        """print_list(file=sys.stdout): pretty print the list of
        registered names to file."""
        pprint.pprint(self.inquire_names(), file)

    def open(self, filename, group=None,mode='r'):
        "open(name, group=None,mode='r'): open file, optionally starting at group"
        self.close()
        if filename:
            if group is None: group = self.inquire_group()
            self._file = tables.openFile(filename,mode=mode,rootUEP=group)
            self.set_group(group)
            self._read_delimiter()
            self._scan_Nodes()
            self._load_cache_with_pickledict()
            self._ints = self._file.getNode(group,'_ints')
            self._readtable(self._ints)
            self._floats = self._file.getNode(group,'_floats')
            self._readtable(self._floats)

    def read(self,name):
        "read(name) = the value of name as a Python object."
        self.check_open()
        try:
            return self._cache[name]
        except KeyError:
            # --- Replace the delimiter in the input name with the one that
            # --- was used when writing the file.
            if self._delimiter != '@':
                name = self._fixdelimiter.sub(self._delimiter,name)

            # --- Read in the data.
            node = self._file.getNode(self.inquire_group(),name)
            data = node.read()

            # --- Handle special cases
            if node.title == "Pickled":
              data = cPickle.loads(data)
            if node.title == 'ZeroLength':
              data = numpy.zeros(data[:-1],dtype=string.ascii_letters[data[-1]])

            # --- Add the data to the cache
            self._cache[name] = data

            # --- And finally return it.
            return data

    def _read_delimiter(self):
        "read in the delimiter"
        node = self._file.getNode(self.inquire_group(),'_delimiter_')
        if node.title != 'Delimiter':
            raise PRError,"Error: the delimiter was overwritten"
        self._delimiter = node.read()

    def _scan_Nodes(self):
        self.check_open()
        group = self.inquire_group()
        for node in self._file.iterNodes(group,classname='Leaf'):
            
            if node.title == 'HashedName':
                # --- Find hashed names in the tables
                hashedname = node.name
                origname = node.read()
                self._hashednames[hashedname] = origname
            elif node.name in ['_ints','_floats','_pickledict']:
                # --- These names are special
                pass
            else:
                # --- Save the rest of the names in the names list
                # --- Fix the delimiter in all of the names if needed.
                if self._delimiter != '@':
                    name = self._fixdelimiter.sub(self._delimiter,name)
                self._names.append(node.name)

    def _load_cache_with_pickledict(self):
        # --- Read in the pickled dictionary
        node = self._file.getNode(self.inquire_group()+'_pickledict')
        if node.title != 'PickleDict':
            print "Warning: the pickle dictionary was written over, the data may not be reliable"
        globalpickledict = cPickle.loads(node.read())

        # --- Add the unpickled data to the cache
        self._cache.update(globalpickledict)

        # --- Add all of the names to the list
        self._names += globalpickledict.keys()

    def _readtable(self,table):
        for i in table.iterrows():
            name = i['name']
            try: name = self._hashednames[name]
            except KeyError: pass
            self._cache[name] = i['value']
            self._names.append(name)

    def set_group(self, name):
        """set_group(name)
        -- change HDF group to name, return status"""
        if name[0] == '/':
            group = name
        elif len(self.inquire_group()) == 1:
            group = '/' + name
        else:
            group = self.inquire_group() + '/' + name
        if group[-1] != '/': group = group + '/'
        self._group = group

    def set_verbosity(self, flag):
        """verbose(flag) sets verbosity level to flag.
        0 for quiet operation, 
        1 to report closing only, 
        2 to report access to data."""
        if 0 <= flag <= 2:
            self._verbose = flag
        else:
            raise PRError, 'Illegal value for verbosity: '+`flag`

