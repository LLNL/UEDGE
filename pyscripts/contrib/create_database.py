# Holm10 Nov 5 2019, created from scratch
'''
Module defining functions for reading a UEDGE hdf5 dump
'''
from h5py import File
import uedge
from uedge.hdf5 import hdf5_dump
from uedge.contrib.utils import readcase


def natsort(l): 
    from re import split
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)


class SETUP():
    ''' Class containing a list of CASE objects '''
    def __init__(self,caselist):
        ''' Stores the list sorted by midplane separatrix electron density '''
        self.cases=caselist
        self.sort_mp('bbb.ne')
        self.ev=1.602e-19


    '''==========================================
    Check species and index
    =========================================='''
    def check_species(self,var,s):
        ''' Boolean that checks whether the parameter requested has a species index '''
        # See if the variable has several species
        if len(self.cases[0].get(var).shape)>2:
            # See if a specific species is requested
            if s is None:   # If not, warn user that standard index is used
                print('WARNING! Multi-species array requested. Using species index 0/{}.'.format(self.cases[0].get(var).shape[-1]-1))
                return 0
            else:           # If requested, don't issue warning
                return s
        else:   # Do nothing if single-species array
            return False

    def itind(self,var):
        ''' Function that returns the index of the outer target paramet of var - different for flux parameters and volumetric quantities '''
        
        # Lookup table listing quatities given at, or over, cell faces. To be amended in the future as more variables are identified
        fluxvar=[   'fqp','fq2','fqx','fqxb','fdiaxlb','fdiaxrb','floxebgt','floxibgt','fqy','fqyb',
                    'fqyn','fqym','fqya','fqydt','fqyao','fqyd','fqygp','fq2d','fqypneo','fq2pneo',
                    'fqyqneo','fq2qneo','fnix','fnixcb','fniy','fniy4ord','fniycb','flnix','flniy',
                    'fmix','fmiy','fmity','fmixy','feex','feey','feexy','feey4ord','feix','feiy',
                    'qipar','fniycbo','feiycbo','feeycbo','feixy','feiy4ord','fngx','flngx','fngxs',
                    'fngy','flngy','fngxy','flngxy','fngyx','uup','up','upi']
        if var in fluxvar:
            return 0
        else:
            return 1


    '''==========================================
    Handle the case list
    =========================================='''

    def sort_core(self,var,s=None):
        ''' Sorts list by value of var at core midplane
        sort_core(var,**keys)

        Variables:
        var:        String of variable to be used in sorting, e.g. 'bbb.ne'

        Keyword parameters:
        s[=0]:      Species index to be used, defaults to 0

        Returns:
        Void
        '''
        # Check whether the requested parameter has a species index
        suse=self.check_species(var,s) 
        if suse is False:
            self.cases.sort(key=lambda case: case.get(var)[case.ixmp(),0])
        else: 
            self.cases.sort(key=lambda case: case.get(var)[case.ixmp(),0,suse])


    def sort_mp(self,var,s=None):
        ''' Sorts list by value of var at  midplane immediately outside the separatrix
        sort_mp(var,**keys)

        Variables:
        var:        String of variable to be used in sorting, e.g. 'bbb.ne'

        Keyword parameters:
        s[=0]:      Species index to be used, defaults to 0

        Returns:
        Void
        '''
        # Check whether the requested parameter has a species index
        suse=self.check_species(var,s) 
        if suse is False:
            self.cases.sort(key=lambda case: case.get(var)[case.ixmp(),case.iysptrx()+1])
        else: 
            self.cases.sort(key=lambda case: case.get(var)[case.ixmp(),case.iysptrx()+1,suse])



    def sort_ind(self,var,ix,iy,s=None):
        ''' Sorts list by value of var at index (ix,iy) 
        sort_mp(var,ix,iy,**keys)

        Variables:
        var:        String of variable to be used in sorting, e.g. 'bbb.ne'
        ix:         Polidal index for sorting
        iy:         Radial index for sorting

        Keyword parameters:
        s[=0]:      Species index to be used, defaults to 0

        Returns:
        Void
        '''
        # Check whether the requested parameter has a species index
        suse=self.check_species(var,s) 
        if suse is False:
            self.cases.sort(key=lambda case: case.get(var)[ix,iy])
        else: 
            self.cases.sort(key=lambda case: case.get(var)[ix,iy,suse])


    def sort_rowmax(self,var,row,s=None):
        ''' Sorts list by value of var at index (ix,iy) 
        sort_mp(var,row,**keys)

        Variables:
        var:        String of variable to be used in sorting, e.g. 'bbb.ne'
        Row:        Poloidal row along which the max value is found, and used to sort the cases

        Keyword parameters:
        s[=0]:      Species index to be used, defaults to 0

        Returns:
        Void
        '''
        # Check and warn if row is in guard cells
        if row in [0,-1,self.cases[0].nx()+1]:
            print('WARNING! Requested row is guard cell row')

        # Check whether the requested parameter has a species index
        suse=self.check_species(var,s) 
        if suse is False:
            self.cases.sort(key=lambda case: max(case.get(var)[row,:]))
        else: 
            self.cases.sort(key=lambda case: max(case.get(var)[row,:,suse]))


    def get_closest_mp(self,var,val,s=None):
        ''' Get the index of the case with var closest to val at MP-sep'''
        return abs(self.mp(var,s=s)-val).argmin()

    def get_closest_core(self,var,val,s=None):
        ''' Get the index of the case with var closest to val at core MP'''
        return abs(self.core(var,s=s)-val).argmin()

    def get_closest_index(self,var,val,ix,iy,s=None):
        ''' Get the index of the case with var closest to val at location specified by index'''
        return abs(self.index(var,ix,iy,s=s)-val).argmin()


    '''==========================================
    Get parameter
    =========================================='''
    def get_name(self):
        ''' Returns a list of case names '''
        return [case.get('bbb.label')[0].decode('utf-8').strip() for case in self.cases]


    def get(self,var,s=None):
        ''' Returns a list of arrays containing var
        get(var,**keys)

        Variables:
        var:        String of variable to be returned, e.g. 'bbb.ne'

        Keyword parameters:
        s[=0]:      Species index to be used, defaults to 0
        '''
        from numpy import asarray
        suse=self.check_species(var,s) 
        if suse is False:
           return asarray([case.get(var) for case in self.cases])
        else: 
           return asarray([case.get(var)[:,:,suse] for case in self.cases])


    '''==========================================
    Get locations
    =========================================='''

    def row(self,var,row,s=None):
        ''' Returns a list of 1D arrays containing var for row
        row(var,row,**keys)

        Variables:
        var:        String of variable to be returned, e.g. 'bbb.ne'
        row:        Row along which to return variable

        Keyword parameters:
        s[=0]:      Species index to be used, defaults to 0
        '''
        from numpy import asarray

        # Check and warn if row is in guard cells
        if row in [0,-1,self.cases[0].nx()+1]:
            print('WARNING! Requested row is guard cell row')

        suse=self.check_species(var,s) 
        if suse is False:
           return asarray([case.get_row(var,row) for case in self.cases])
        else: 
           return asarray([case.get_row(var,row,s=suse) for case in self.cases])


    def row_mp(self,var,s=None):
        ''' Returns a list of 1D arrays of var at the midplane
        row_mp(var,**keys)

        Variables:
        var:        String of variable to be returned, e.g. 'bbb.ne'

        Keyword parameters:
        s[=0]:      Species index to be used, defaults to 0
        '''
        return self.row(var,self.cases[0].ixmp(),s)

    def row_it(self,var,s=None):
        ''' Returns a list of 1D arrays of var at the inner target
        row_it(var,row,**keys)

        Variables:
        var:        String of variable to be returned, e.g. 'bbb.ne'

        Keyword parameters:
        s[=0]:      Species index to be used, defaults to 0
        '''
        return self.row(var,self.itind(var),s)


    def row_ot(self,var,s=None):
        ''' Returns a list of 1D arrays of var at the outer target
        row_ot(var,row,**keys)

        Variables:
        var:        String of variable to be returned, e.g. 'bbb.ne'

        Keyword parameters:
        s[=0]:      Species index to be used, defaults to 0
        '''
        return self.row(var,-2,s)


    def ft(self,var,ft,s=None):
        ''' Returns a list of 1D arrays containing var along flux-tube ft
        ft(var,ft,**keys)

        Variables:
        var:        String of variable to be returned, e.g. 'bbb.ne'
        ft:         Flux tube along which to return variable

        Keyword parameters:
        s[=0]:      Species index to be used, defaults to 0
        '''
        from numpy import asarray

        # Check and warn if row is in guard cells
        if ft in [0,-1,self.cases[0].ny()+1]:
            print('WARNING! Requested flux tube is guard cell flux tube')

        suse=self.check_species(var,s) 
        if suse is False:
           return asarray([case.get_ft(var,ft) for case in self.cases])
        else: 
           return asarray([case.get_ft(var,ft,s=suse) for case in self.cases])

    def sep(self,var,s=None):
        ''' Returns a list of 1D arrays containing var along (first radial index outside) the separatrix
        sep(var,**keys)

        Variables:
        var:        String of variable to be returned, e.g. 'bbb.ne'

        Keyword parameters:
        s[=0]:      Species index to be used, defaults to 0
        '''
        return self.ft(var,self.cases[0].iysptrx()+1)

    def index(self,var,ix,iy,s=None):
        ''' Returns a list values of var at location (ix,iy)
        index(var,ix,iy,**keys)

        Variables:
        var:        String of variable to be returned, e.g. 'bbb.ne'
        ix:         Polidal index to use
        iy:         Radial index to use

        Keyword parameters:
        s[=0]:      Species index to be used, defaults to 0
        '''
        if ix in [0,-1,self.cases[0].nx()+1]:
            print('WARNING! Requested poloidal index is guard cell index')
        if iy in [0,-1,self.cases[0].ny()+1]:
            print('WARNING! Requested radial index is guard cell index')
        from numpy import asarray

        suse=self.check_species(var,s) 
        if suse is False:
           return asarray([case.get_ind(var,ix,iy) for case in self.cases])
        else: 
           return asarray([case.get_ind(var,ix,iy,s=suse) for case in self.cases])


    def mp(self,var,s=None):
        ''' Returns a list values of var at the separatrix (immediately outside) midplane
        mp(var,**keys)

        Variables:
        var:        String of variable to be returned, e.g. 'bbb.ne'

        Keyword parameters:
        s[=0]:      Species index to be used, defaults to 0
        '''
        return self.index(var,self.cases[0].ixmp(),self.cases[0].iysptrx()+1,s)

    def core(self,var,s=None):
        ''' Returns a list values of var at the separatrix (immediately outside) midplane
        mp(var,**keys)

        Variables:
        var:        String of variable to be returned, e.g. 'bbb.ne'

        Keyword parameters:
        s[=0]:      Species index to be used, defaults to 0
        '''
        from numpy import asarray
        return asarray([case[self.cases[0].ixmp(),0] for case in self.get(var,s)])

    '''==========================================
    Get min/max
    =========================================='''

    def row_min(self,var,row,s=None):
        ''' Returns a list values of the minimum value of var the specified row
        row_min(var,row,**keys)

        Variables:
        var:        String of variable to be returned, e.g. 'bbb.ne'
        row:        Row along which to find the minimum

        Keyword parameters:
        s[=0]:      Species index to be used, defaults to 0
        '''
        from numpy import asarray

        # Check and warn if row is in guard cells
        if row in [0,-1,self.cases[0].nx()+1]:
            print('WARNING! Requested row is guard cell row')

        suse=self.check_species(var,s) 
        if suse is False:
           return asarray([min(case.get_row(var,row)) for case in self.cases])
        else: 
           return asarray([min(case.get_row(var,row,s=suse)) for case in self.cases])

    def ot_min(self,var,s=None):
        ''' Returns a list values of the minimum value of var the outer target
        ot_min(var,**keys)

        Variables:
        var:        String of variable to be returned, e.g. 'bbb.ne'

        Keyword parameters:
        s[=0]:      Species index to be used, defaults to 0
        '''
        return self.row_min(var,-2,s)

    def it_min(self,var,s=None):
        ''' Returns a list values of the minimum value of var the inner target
        it_min(var,**keys)

        Variables:
        var:        String of variable to be returned, e.g. 'bbb.ne'

        Keyword parameters:
        s[=0]:      Species index to be used, defaults to 0
        '''
        return self.row_min(var,self.itind(var),s)


    def mp_min(self,var,s=None):
        ''' Returns a list values of the minimum value of var the midplane
        mp_min(var,**keys)

        Variables:
        var:        String of variable to be returned, e.g. 'bbb.ne'

        Keyword parameters:
        s[=0]:      Species index to be used, defaults to 0
        '''
        return self.row_min(var,self.cases[0].ixmp(),s)


    def row_max(self,var,row,s=None):
        ''' Returns a list values of the maximum value of var the specified row
        row_max(var,row,**keys)

        Variables:
        var:        String of variable to be returned, e.g. 'bbb.ne'
        row:        Row along which to find the minimum

        Keyword parameters:
        s[=0]:      Species index to be used, defaults to 0
        '''
        from numpy import asarray

        # Check and warn if row is in guard cells
        if row in [0,-1,self.cases[0].nx()+1]:
            print('WARNING! Requested row is guard cell row')

        suse=self.check_species(var,s) 
        if suse is False:
           return asarray([max(case.get_row(var,row)) for case in self.cases])
        else: 
           return asarray([max(case.get_row(var,row,s=suse)) for case in self.cases])


    def ot_max(self,var,s=None):
        ''' Returns a list values of the maximum value of var the outer target
        ot_max(var,**keys)

        Variables:
        var:        String of variable to be returned, e.g. 'bbb.ne'

        Keyword parameters:
        s[=0]:      Species index to be used, defaults to 0
        '''
        return self.row_max(var,-2,s)

    def it_max(self,var,s=None):
        ''' Returns a list values of the maximum value of var the inner target
        it_max(var,**keys)

        Variables:
        var:        String of variable to be returned, e.g. 'bbb.ne'

        Keyword parameters:
        s[=0]:      Species index to be used, defaults to 0
        '''
        return self.row_max(var,self.itind(var),s)


    def mp_max(self,var,s=None):
        ''' Returns a list values of the maximum value of var the midplane
        mp_max(var,**keys)

        Variables:
        var:        String of variable to be returned, e.g. 'bbb.ne'

        Keyword parameters:
        s[=0]:      Species index to be used, defaults to 0
        '''
        return self.row_max(var,self.cases[0].ixmp(),s)





        
    def ft_min(self,var,ft,s=None):
        ''' Returns a list values of the minimum value of var the specified flux tube
        ft_min(var,ft,**keys)

        Variables:
        var:        String of variable to be returned, e.g. 'bbb.ne'
        ft:         Flux tube along which to find the minimum

        Keyword parameters:
        s[=0]:      Species index to be used, defaults to 0
        '''
        from numpy import asarray

        # Check and warn if row is in guard cells
        if ft in [0,-1,self.cases[0].nx()+1]:
            print('WARNING! Requested flux tube is guard cell flux tube')

        suse=self.check_species(var,s) 
        if suse is False:
           return asarray([min(case.get_ft(var,ft)) for case in self.cases])
        else: 
           return asarray([min(case.get_ft(var,ft,s=suse)) for case in self.cases])


        
    def sep_min(self,var,s=None):
        ''' Returns a list values of the minimum value of var along (immediately radially outside) the separatrix
        sep_min(var,**keys)

        Variables:
        var:        String of variable to be returned, e.g. 'bbb.ne'

        Keyword parameters:
        s[=0]:      Species index to be used, defaults to 0
        '''
        return self.ft_min(var,self.cases[0].iysptrx()+1,s)



    def ft_max(self,var,ft,s=None):
        ''' Returns a list values of the maximum value of var the specified flux tube
        ft_max(var,ft,**keys)

        Variables:
        var:        String of variable to be returned, e.g. 'bbb.ne'
        ft:         Flux tube along which to find the minimum

        Keyword parameters:
        s[=0]:      Species index to be used, defaults to 0
        '''
        from numpy import asarray

        # Check and warn if row is in guard cells
        if ft in [0,-1,self.cases[0].nx()+1]:
            print('WARNING! Requested flux tube is guard cell flux tube')

        suse=self.check_species(var,s) 
        if suse is False:
           return asarray([max(case.get_ft(var,ft)) for case in self.cases])
        else: 
           return asarray([max(case.get_ft(var,ft,s=suse)) for case in self.cases])

 
    def sep_max(self,var,s=None):
        ''' Returns a list values of the maximum value of var along (immediately radially outside) the separatrix
        sep_max(var,**keys)

        Variables:
        var:        String of variable to be returned, e.g. 'bbb.ne'

        Keyword parameters:
        s[=0]:      Species index to be used, defaults to 0
        '''
        return self.ft_max(var,self.cases[0].iysptrx()+1,s)





    '''==========================================
    Get locations
    =========================================='''

    def get_maxlocation(self, var,s=None,xind=(1,-1),yind=(1,-1)):
        ''' Returns a list of (R,Z) coordinates of max of var in xind/yind index interval 
        get_maxlocation(var,**keys)

        Variables:
        var:            String of variable to be returned, e.g. 'bbb.ne'

        Keyword parameters:
        s[=0]:          Species index to be used, defaults to 0
        xind[=(1,-1)]   Poloidal index space to search
        yind[=(1,-1)]   Radial index space to search
        ''' 
        from numpy import asarray
        return asarray([case.get_maxlocation(var,s,xind,yind) for case in self.cases])

    def get_minlocation(self, var,s=None,xind=(1,-1),yind=(1,-1)):
        ''' Returns a list of (R,Z) coordinates of min of var in xind/yind index interval 
        get_maxlocation(var,**keys)

        Variables:
        var:            String of variable to be returned, e.g. 'bbb.ne'

        Keyword parameters:
        s[=0]:          Species index to be used, defaults to 0
        xind[=(1,-1)]   Poloidal index space to search
        yind[=(1,-1)]   Radial index space to search
        ''' 
        from numpy import asarray
        return asarray([case.get_minlocation(var,s,xind,yind) for case in self.cases])
    
    '''==========================================
    Get special
    =========================================='''
        
    def get_nusole(self):
        from numpy import asarray
        ''' Returns the electron SOL collisionality defined at the Sep-MP'''
        # TODO: implement this one: requires calculation of the connection length
        return asarray([(1e-16*case.get_mp('bbb.ne')*case.get_rightL()[case.ixmp(),case.iysptrx()+1])/
                ((case.get_mp('bbb.te')/self.ev)**2) for case in self.cases])

    def get_nusoli(self):
        from numpy import asarray
        ''' Returns the electron SOL collisionality defined at the Sep-MP'''
        # TODO: implement this one: requires calculation of the connection length
        return asarray([(1e-16*case.get_mp('bbb.ni',s=0)*case.get_rightL()[case.ixmp(),case.iysptrx()+1])/
                ((case.get_mp('bbb.ti')/self.ev)**2) for case in self.cases])

    def get_rowsum(self):
        ''' To be implemented '''
        # TODO: implement
        return None

    def get_ftsum(self):
        ''' To be implemented '''
        # TODO implement
        return None

    def get_allsum(self,var):
        ''' To be implemented '''
        # TODO implement sum in bot spatial dimensions, not in index
        return None

    def get_rangesum(self):
        ''' To be implemented '''
        # TODO implement sum in bot spatial dimensions within range, not in index
        return None

    '''==========================================
    Get case indices
    =========================================='''
    def get_closest(self,val,var,ix,iy,s=None):
        ''' Get index of  case with closest value to var at (ix,iy) 
        get_closest(vel,var,ix,iy,**keys)

        Variables:
        val:            Value to match
        var:            String of variable to be returned, e.g. 'bbb.ne'
        ix:             Poloidal location index
        iy:             Radial location index

        Keyword parameters:
        s[=0]:          Species index to be used, defaults to 0
        '''
        ind=abs(self.index(var,ix,iy,s)-val).argmin()
        print('Closest value is: '+str(self.index(var,ix,iy,s)[ind]))
        return ind
        
    
    def get_closest_mp(self,val,var,s=None):
        ''' Get index of  case with closest value to var at midplane separatrix (immediately radially outside of) 
        get_closest_mp(vel,var,**keys)

        Variables:
        val:            Value to match
        var:            String of variable to be returned, e.g. 'bbb.ne'

        Keyword parameters:
        s[=0]:          Species index to be used, defaults to 0
        '''
        ''' Get index of  case with closest value to var at OMP '''
        return self.get_closest(val,var,self.cases[0].ixmp(),self.cases[0].iysptrx()+1,s)

    
    def get_closest_core(self,val,var,s=None):
        ''' Get index of  case with closest value to var at core midplane 
        get_closest_core(vel,var,**keys)

        Variables:
        val:            Value to match
        var:            String of variable to be returned, e.g. 'bbb.ne'

        Keyword parameters:
        s[=0]:          Species index to be used, defaults to 0
        '''
        from numpy import asarray
        return self.get_closest(val,var,self.cases[0].ixmp(),1,s)




class CASE():
    ''' Class containing UEDGE run data '''
    def __init__(self,variables=None):
        from uedge import com,aph,api,bbb,flx,grd,svr,wdf
        from numpy import copy
        ''' Reads 'variables' from 'packages' into dictionary '''
        self.data=dict()
        # Read variables
        for var in variables:
            try:
                exec('self.data["'+var+'"]=copy('+var+')') in globals(),locals()
            except:
                print('Warning! '+var+' not found')


    '''==========================================
    Helper functions
    =========================================='''

    def assert_shape(self,arr,s):
        ''' Returns a 2D array of species s, or the array itself if already a 2D array '''
        if len(arr.shape)>2:
            return arr[:,:,s]
        else:
            return arr

    '''==========================================
    Get functions
    =========================================='''

    def get(self,var):
        ''' Returns the requested parameter '''
        try:
            return self.data[var]
        except:
            print('Error! '+var+' not found in data!')
            return

    def get_ft(self,var,ft,s=None):
        ''' Get the requested FT '''
        return self.assert_shape(self.get(var),s)[1:-1,ft]
        
    def get_row(self,var,row,s=None):
        ''' Get the requested row '''
        return self.assert_shape(self.get(var),s)[row,1:-1]

    def get_ind(self,var,ix,iy,s=None):
        ''' Get the requested index '''
        return self.assert_shape(self.get(var),s)[ix,iy]

    def get_core(self,var,s=None):
        ''' Gets the parameter at the core mp '''
        return self.assert_shape(self.get(var),s)[self.ixmp(),0]

    def get_mp(self,var,s=None):
        ''' Gets the parameter immediately outside the sep at the mp '''
        return self.assert_shape(self.get(var),s)[self.ixmp(),self.iysptrx()+1]

    '''==========================================
    Get min/max
    =========================================='''

    def get_ftmax(self,var,ft,s=None):
        ''' Get the max value in the poloidal direction for flux-tube 'ft' '''
        return max(self.assert_shape(self.get(var),s)[1:-1,ft])

    def get_ftmin(self,var,ft,s=None):
        ''' Get the max value in the poloidal direction for flux-tube 'ft' '''
        return min(self.assert_shape(self.get(var),s)[1:-1,ft])
            
    def get_rowmax(self,var,row,s=None):
        ''' Get the max value in the radial direction for row 'row' '''
        return max(self.assert_shape(self.get(var),s)[row,1:-1])
    
    def get_rowmin(self,var,row,s=None):
        ''' Get the min value in the radial direction for row 'row' '''
        return min(self.assert_shape(self.get(var),s)[row,1:-1])
    
    '''==========================================
    Get indices
    =========================================='''

    def ixpt1(self):
        ''' Returns the left X-point poloidal index '''
        return self.get('com.ixpt1')[0]    

    def ixpt2(self):
        ''' Returns the right X-point poloidal index '''
        return self.get('com.ixpt2')[0]    

    def iysptrx(self):
        ''' Returns the radial index of the last core FT '''
        return self.get('com.iysptrx')

    def ixmp(self):
        ''' Returns the poloidal midplane index '''
        return self.get('bbb.ixmp')

    def nx(self):
        ''' Returns the poloidal midplane index '''
        return self.get('com.nx')

    def ny(self):
        ''' Returns the poloidal midplane index '''
        return self.get('com.ny')

    def get_coordinates(self,pos,node=0):
        ''' Returns the (R,Z) coordinates of index pos=(x,y) at node(=0 default)'''
        (x,y)=pos
        return (self.get('com.zm')[x,y,node], self.get('com.rm')[x,y,node])

    
    def get_maxindlocation(self,var,s=None,xind=(1,-1),yind=(1,-1)):
        ''' Returns the (ix,iy) location of the maximum value in range '''
        from numpy import unravel_index
        (xmin,xmax),(ymin,ymax)=xind,yind
        ret=unravel_index(  self.assert_shape(self.get(var),s)[xmin:xmax,ymin:ymax].argmax(axis=None),
                            self.assert_shape(self.get(var),s)[xmin:xmax,ymin:ymax].shape)
        return (ret[0]+xmin, ret[1]+ymin)   # Increase the indices to get the actual locations in the full array

    
    def get_minindlocation(self,var,s=None,xind=(1,-1),yind=(1,-1)):
        ''' Returns the (ix,iy) location of the maximum value in range '''
        from numpy import unravel_index
        (xmin,xmax),(ymin,ymax)=xind,yind
        ret=unravel_index(  self.assert_shape(self.get(var),s)[xmin:xmax,ymin:ymax].argmin(axis=None),
                            self.assert_shape(self.get(var),s)[xmin:xmax,ymin:ymax].shape)
        return (ret[0]+xmin, ret[1]+ymin)   # Increase the indices to get the actual locations in the full array


    '''==========================================
    Get locations
    =========================================='''

    def get_maxlocation(self,var,s=None,xind=(1,-1),yind=(1,-1)):
        ''' Returns the (R,Z) location of the maximum value in range '''
        return self.get_coordinates(self.get_maxindlocation(var,s=s,xind=xind,yind=yind))
         
    
    def get_minlocation(self,var,s=None,xind=(1,-1),yind=(1,-1)):
        ''' Returns the (R,Z) location of the maximum value in range '''
        return self.get_coordinates(self.get_minindlocation(var,s=s,xind=xind,yind=yind))


    '''==========================================
    Get special
    =========================================='''

    def get_rightL(self):
        ''' Returns an array containing the connection lengths to the right target.
            Returns zeros for closed flux surfaces.
        '''
        from numpy import zeros_like
        ret=zeros_like(self.get('bbb.ne'))
        for ix in range(self.nx()+1,0,-1):
            ret[ix,self.iysptrx()+1:]=((self.get('com.bphi')**2+self.get('com.bpol')**2)**0.5)[ix,self.iysptrx()+1:,0]*self.get('com.dx')[ix,self.iysptrx()+1:]
        for ix in range(self.nx()+1,self.ixpt2(),-1):
            ret[ix,:self.iysptrx():]=((self.get('com.bphi')**2+self.get('com.bpol')**2)**0.5)[ix,:self.iysptrx():,0]*self.get('com.dx')[ix,:self.iysptrx():]
        return ret



def createdump(name,packages=['bbb','com']):
    ''' Creates a dump 
        Paremeters:
            name        Name of the dump (str)
            packages    List of strings of packages to include
    '''
    from uedge.uedge_lists import list_packages
    # Store the requested UEDGE objects to list
    packages=[list_packages(objects=1)[list_packages().index(x)] for x in packages]
    hdf5_dump(name,packages=packages)
    
def create_dict(dump,variables=None):
    ''' Creates a dictionary of all variables in 'dump'
        Parameters:
            dump        Path to dump to be dictionarized
            variables   List of variables to be saved (default is all) - no packages
    '''
    from numpy import array
    ret=dict()          # Dictionary to be returned
    f=File(dump,'r')    # Open readable dump as f
    for pack in f.keys():   # Loop through the list of packages
        p=f.get(pack)       # Get the package object p
        for var in p.keys():    # Loop through all variables in the package
            if variables is None:   # By default, save all variables to dict
                ret[pack+'.'+var]=array(p.get(var))    # Store the variable to the dictionary as array
            elif var in variables:  # If variables listed, save only listed variables
                ret[pack+'.'+var]=array(p.get(var))    # Store the variable to the dictionary as array
    return ret

def create_database(savename=None,sortlocation='mp',outpath='.',path='.',subpath='data',commands=[],ret=True,variables=None):
    ''' Creates a database
        Parameters:
            savename        If set, saves dict as pickle named 'savename'
            sortlocation    Location for sorting by te: 'core' or 'mp'
            outpath         Path to location where pickle is saved
            path            Path to parent directory
            subpath         Path to input.py within child directories of path: 
                            path/*/supath/input.py
            commands        List of commands to be executed before restoring solution
            variables       List of all variable names, including package, to be stored
            ret             Boolean whether to return dict or not
    '''
    from os import getcwd,chdir,remove,walk
    from os.path import abspath  
    from uedge.uexec import uexec 
    from pickle import dump
    from uedge import bbb
    from importlib import reload

    outpath=abspath(outpath)    # Get absolute path of out directory
    chdir(path)                 # Go to the parent directory
    parent=getcwd()               # Get path to parent
    # Get list of subdirectories in path
    dirs=natsort(next(walk(path))[1])
    # Omit supporting file directories
    try:
        dirs.remove('grid')
    except:    
        pass
    try:
        dirs.remove('rates')
    except:
        pass
    try:
        dirs.remove('ignore')
    except:
        pass
    
    if len(dirs)==0:
        return 'No directories found! Aborting...'


    # Empty list to return
    retl=[]
    
    if variables is None:
        variables=default_variables()


    for child in dirs:      # Loop through all directories
        print('******************************')
        print('*** Directory: '+child+' ***')
        print('******************************')
        readcase(child,subpath) # Restore the case


        # Execute any commands before executing
        for cmd in commands:
            exec(cmd) in globals(),locals()
        # Read and repopulate all arrays
        bbb.issfon=0;bbb.ftol=1e20;bbb.exmain()


        retl.append(CASE(variables))

        chdir(parent)
    
    
    lst=SETUP(retl) 
    # Get the sep and xpt locations
    if sortlocation=='core':
        lst.sort_core('bbb.ne')
    elif sortlocation=='mp':
        lst.sort_mp('bbb.ne')
    else:
        print('Error! Unknown sortin location "'+sortlocation+'". Terminating...')
    chdir(outpath)
    # Check if save requested
    if savename is not None:
        with open(savename,'wb') as f:
            dump(lst,f)
    if ret:
        return lst


def restore_database(name):    
    ''' Restores pickled case '''
    from pickle import load
        
    with open(name,'rb') as f:
        ret=load(f)

    try:
        with open(name,'rb') as f:
            ret=load(f)
    except:
        print('ERROR! Could not find file "'+name+'"! Aborting...')
        ret=[]
    return ret


def default_variables():
    ''' Returns a list of commonly requested parameters '''
    return [    'bbb.feex', 'bbb.feey',
                'bbb.feix', 'bbb.feiy',
                'bbb.fegx', 'bbb.fegy',
                'bbb.fnix', 'bbb.fniy',
                'bbb.fngx', 'bbb.fngy',
                'bbb.te', 'bbb.ti', 'bbb.tg',
                'bbb.ne', 'bbb.ni', 'bbb.ng',
                'bbb.up', 'bbb.uup', 'bbb.vy',
                'bbb.vex', 'bbb.upe', 'bbb.vey',
                'com.xcs', 'com.yyc',
                'com.rm', 'com.zm',
                'com.gxf', 'com.gx',
                'com.gyf', 'com.gy',
                'com.dx', 'com.dy',
                'bbb.eqpg', 'bbb.eqp',
                'com.vol','com.nx','com.ny',
                'bbb.ixmp', 'com.ixpt1', 'com.ixpt2', 'com.iysptrx',
                'com.sx', 'com.sy', 
                'com.bpol','com.bphi',
                'bbb.hcxg', 'bbb.hcyg','bbb.floxge','bbb.floyge','bbb.conxge','bbb.conyge',
                'bbb.kxg_use','bbb.kyg_use','bbb.psor','bbb.psordis','bbb.psorrgc',
                'bbb.pradiz','bbb.pradrc','bbb.pbinde','bbb.pbindrc','bbb.prdiss','bbb.pibirth',
                'bbb.pradc','bbb.pradz','bbb.pradzc','bbb.prad','bbb.pradht',
                'bbb.erliz','bbb.erlrc',
                'bbb.psorbgg','bbb.psorbgz',
                'bbb.ziin', 'bbb.minu', 'bbb.mi', 'bbb.mg', 'bbb.ziin',
                'bbb.label' ]




