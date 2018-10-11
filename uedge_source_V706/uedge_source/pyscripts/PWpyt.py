"""
# Copyright (c) 1996, 1997, The Regents of the University of California.
# All rights reserved.  See LEGAL.LLNL for full text and disclaimer.

HDF basic writer class PW by David Grote, LLNL
Heavily modified from PW.py originally written by Paul Dubois, LLNL, to use
HDF files.
$Id: PWpyt.py,v 7.0 2018/02/28 18:43:48 meyer8 Exp $

This writes all scalars into two tables, one for ints and one for floats.
Everything else is written out using createArray.
Note that writing out scalars to tables is faster than writing then out
individually using createArray.
For scalars, the setup requires a maximum name length. Any names longer are
hashed.
Arrays are written out directly. Everything else is put into a dictionary
which is written out on close.

"""
import tables
from tables.exceptions import NaturalNameWarning
import cPickle
import re
import string
import numpy
import warnings
import types

_version = '0.4'

class PWError(Exception):
    pass

_max_name_len = 200
_transtable = (10*string.ascii_lowercase)[:256]

class IntScalar(tables.IsDescription):
    name = tables.StringCol(_max_name_len)
    value = tables.IntCol()
class FloatScalar(tables.IsDescription):
    name = tables.StringCol(_max_name_len)
    value = tables.FloatCol()

class PW:
    "HDF file writer class."
    file_type = "HDF"

    def __del__(self):
        "Close any file open when this object disappears."
        self.close()

    def __init__(self, filename=None, mode="w", group='/',
                 verbose = 1, delimiter='@'):
        "PW(filename='', verbose=1) creates filename if given" 
        self.__dict__['_delimiter'] = delimiter
        self.__dict__['_file'] = None
        self.set_verbosity(verbose)
        self.set_group(group)
        if filename is not None:
            self.open(filename,group,mode)

        self.__dict__['_globalpickledict'] = {}

        # --- This allows replacement of '@' with a user specified
        # --- delimiter.
        self.__dict__['_fixdelimiter'] = re.compile('@')

        # --- Turn off the NaturalNameWarning in this case since the
        # --- delimiter cannot be used in a python name.
        if self._delimiter == '@':
            warnings.simplefilter('ignore',NaturalNameWarning)

    def __setattr__(self, name, value):
        self.write(name, value)

    def __repr__(self):
        if self.is_open():
            current_mode = 'opened for writing'
            return 'HDF file %s %s.' % \
               (self.inquire_filename(), current_mode)
        else:
            return PWError,'(PW object not open on any file)'

    __str__ = __repr__

    def check_open(self):
        "check_open(): raise exception if not open for write."
        if not self.is_open():
            raise PWError, 'PW object not open for write.'

    def close(self):
        "close(): close the file."
        h = self.inquire_file()
        if h is not None:
            if self.inquire_verbosity():
                print "Closing HDF file being written:",self.inquire_filename()

            # --- Write out the list of things that have to be pickled,
            # --- pickling them all at once.
            q = cPickle.dumps(self._globalpickledict,-1)
            h.createArray(self.inquire_group(),'_pickledict',q,
                          title='PickleDict')

            # --- Write out the delimiter used.
            h.createArray(self.inquire_group(),
                          '_delimiter_',self._delimiter,
                          title='Delimiter')

            h.flush()
            h.close()

        # --- Reset various quantities
        self.__dict__['_file'] = None
        self.__dict__['_ints'] = None
        self.__dict__['_floats'] = None
        self.__dict__['_globalpickledict'] = {}

        if self._delimiter == '@':
            warnings.simplefilter('default',NaturalNameWarning)

    def inquire_filename(self):
        "inquire_filename() = name of this file."
        if self.is_open():
            return self.inquire_file().filename
        else:
            return ''

    def inquire_file(self):
        "inquire_file() = object open on this file."
        return self._file

    def inquire_mode(self):
        "inquire_mode() = mode('w', or 'a') of this file."
        self.check_open()
        return self._mode

    def inquire_group(self):
        "inquire_group() = present HDF group"
        return self.__dict__['_group']

    def inquire_verbosity(self):
        "inquire_verbosity() = current value of verbose flag."
        return self._verbose_flag

    def is_open(self):
        "is_open() = true if file is open"
        if self.inquire_file() is None: return 0
        return self.inquire_file().isopen

    def open(self, filename,group=None,mode = "w"):
        "open(filename, 'w')"
        self.close()
        assert mode in ['w','a'],Exception("Improper mode: " + mode)
        if group is None: group = self.inquire_group()
        self.set_group(group)
        self.__dict__['_mode'] = mode
        self.__dict__['_file'] = tables.openFile(filename,mode=mode,
                                                 rootUEP=self.inquire_group())
        if mode == 'a':
            # --- If mode is append, check if the ints and floats tables
            # --- have been created. If not, catch the error.
            try:
                self.__dict__['_ints'] = self.inquire_file().root._ints
                self.__dict__['_floats'] = self.inquire_file().root._floats
            except LookupError:
                pass

        # --- If mode is 'w' or if the mode is 'a' but the tables have
        # --- not yet been written, then create them.
        if self.__dict__['_ints'] is None:
            self.__dict__['_ints'] = self.inquire_file().createTable(
                                          self.inquire_group(),
                                          '_ints',IntScalar,"Scalar Ints")
        if self.__dict__['_floats'] is None:
            self.__dict__['_floats'] = self.inquire_file().createTable(
                                          self.inquire_group(),
                                          '_floats',FloatScalar,"Scalar Floats")

    def make_group(self, name, title=''):
        """make_group(name) 
        -- create a new HDF group, return status"""
        self.check_open()
        self.inquire_file().createGroup(self.inquire_group(),name,title=title)
        return 1

    def set_group(self, name):
        """set_group(name) 
        -- change HDF group to name, return status"""
        if name[0] == '/':
            group = name
        elif len(self.inquire_group()) == 1:
            group = '/' + name
        else:
            group = self.inquire_group() + '/' + name
        self.__dict__['_group'] = group

    def set_verbosity(self, flag):
        """set_verbosity(flag) sets verbosity level to flag.
        0 for quiet operation, 
        1 to report closing only, 
        2 to report access to data."""
        if 0 <= flag <= 2:
            self.__dict__['_verbose_flag'] = flag
        else:
            self.__dict__['_verbose_flag'] = 2

    def fixlongnames(self,name):
      """
Handle names longer than the maximum length allowed by IntScalar and
FloatScalar.
      """
      origname = name
      if len(name) > _max_name_len:
          hash = string.translate(md5.new(origname).digest(),_transtable)
          hashedname = hash[:_max_name_len]
          h.createArray(self.inquire_group(),
                        hashedname,
                        origname,
                        title='HashedName')
      else:
          hashedname = name
      return origname,hashedname

    def write(self, name, quantity, title=''):
        """Write quantity to file as 'name'"""
        self.check_open()
        if self.inquire_verbosity() > 1: print "PW::write writing", name
        h = self.inquire_file()

        if self._delimiter != '@':
            name = self._fixdelimiter.sub(self._delimiter,name)

        if isinstance(quantity,types.IntType):
            # --- integers are put into the ints table
            origname,name = self.fixlongnames(name)
            self._ints.row['name'] = name
            self._ints.row['value'] = quantity
            self._ints.row.append()

        elif isinstance(quantity,types.FloatType):
            # --- floats are put into the floats table
            origname,name = self.fixlongnames(name)
            self._floats.row['name'] = name
            self._floats.row['value'] = quantity
            self._floats.row.append()

        elif isinstance(quantity,numpy.ndarray):
            # --- Arrays are written out directly
            if min(quantity.shape) == 0:
                # --- Arrays with a dimension of length zero cannot be
                # --- direcrly written out. Instead, write out its shape
                # --- and type so it can be later reconstructed.
                arrayinfo = (list(quantity.shape) +
                             [string.ascii_letters.find(quantity.dtype.char)])
                h.createArray(self.inquire_group(),name,arrayinfo,"ZeroLength")
            else:
                h.createArray(self.inquire_group(),name,quantity,title=title)

        else:
            # --- Everything else is either skipped or put into the global
            # --- pickle dict

            try:
                try:
                    if quantity.__class__.__module__ == '__main__':
                        if self.inquire_verbosity():
                            print name,' is being skipped since it is an instance of a class defined in main and therefore cannot be unpickled'
                        return
                except:
                    pass
                try:
                    if quantity.__module__ == '__main__':
                        if self.inquire_verbosity():
                            print name,' is being skipped since it is a class defined in main and therefore cannot be unpickled'
                        return
                except:
                    pass
                # --- Check if an object can be pickled. This will be inefficient
                # --- for large objects since they will be pickled twice, but this
                # --- is the only safe way.
                try:
                    q = cPickle.dumps(quantity,-1)
                    del q
                except:
                    if self.inquire_verbosity():
                      print name,' is being skipped since it could not be written or pickled'
                    return
                # --- Things that need to be pickled will all be written out at
                # --- once in the same pickle when the file is being closed so that
                # --- multiple references to the same object are handled correctly.
                self._globalpickledict[name] = quantity
                return
            except:
                pass

            raise PWError,"Could not write the variable %s"%name

if __name__ == "__main__":
    f=PW("foo.hdf")
    a = 1
    b = 2.0
    c = "Hello world"
    from multiarray import *
    d = array([1.,2., 3.])
    e = array([1,2,3])
    g = array(["hello", "world", "array"])
    h = array([[1.,2.,3.], [4.,5.,6]])
    k = 3
    f.a = a
    f.b = b
    f.c = c
    f.d = d
    f.e = e
    f.g = g
    f.h = h
    f.close()
    f.open("foo.hdf", "a")
    f.k = k
    f.close()
# read-back test
    from PRpyt import PR
    f = PR('foo.hdf')
    for x in f.inquire_names():
        print x, "is", eval(x), ", in file it is", eval('f.'+x)
    f.close()
# record-writing
    g = PW('goo.hdf')
    g.set_verbosity(1)
    xh = array([0.]*4)
    for i in range(len(xh)):
        x = i / 10.
        g.write('xh', x, i + 1)
        xh [i] = x
    g.close()
    g = PR('goo.hdf')
    print "xh is", xh, ", file it is ", g.xh
    g.close()






