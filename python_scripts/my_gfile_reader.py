#gfile file reader
#mfitz 2014
#requires https://pypi.python.org/pypi/fortranformat
# use python/3.5.1 on FREIA
# modified by James O to vectorise.

#Will write and read files according to http://w3.pppl.gov/ntcc/TORAY/G_EQDSK.pdf
#Will also attempt to read FLUSH incorrect eqdsk output (provided that blank lines are first removed by hand)

from scipy import interpolate
import matplotlib.pyplot as plt
import numpy as np

class Gfile:
    def __init__(self):
        self.A52=""
        self.nr=0
        self.nz=0
        self.xdim=0
        self.zdim=0
        self.rcentr=0
        self.rleft=0
        self.zmid=0
        self.rmaxis=0
        self.zmaxis=0
        self.psimag=0
        self.psibnd=0
        self.bcentr=0
        self.current=0
        self.ssimag=0
        self.ssibry=0
        self.f_dia=[]
        self.pressure=[]
        self.ffprime=[]
        self.pprime=[]
        self.psi_grid=[[]]
        self.q=[]
        self.nbnd=0
        self.limitr=0
        self.R_bnd=[]
        self.Z_bnd=[]
        self.rlim=[]
        self.zlim=[]
        self.psin=[]
        self.xdum=0      # just a dummy to mimic eqdsktohelena behavior
        self.PsiSpline=0 # Bivariate Spline of Psi Map
        self.FSpline=0   # Spline of F = RBphi
        self.QSpline=0   # Spline of Q
        self.R=[]        # R coordinates
        self.Z=[]        # Z coordinates
    def set(self,obj,value):  self.__dict__[obj]=value
    def get(self,obj): return self.__dict__[obj]

def getGfile(filename):
    import numpy as np
    from io import StringIO
    gfile=Gfile()
    import fortranformat as ff
    f2020=ff.FortranRecordReader('5e16.9')
    f2022=ff.FortranRecordReader('2i5')

    def read(handle,gfile,varList,File):
        line=handle.read(File.readline())
        for pel,el in enumerate(varList):
            gfile.set(el,line[pel])
        return gfile

    def get1DArray(File,npsi):
        res=[]
        while(len(res)<npsi):
            line=f2020.read(File.readline())
            for l in line:
                if(l is not None): #don't append Nones
                    res.append(l)
        return res

    File=open(filename,'r')
    try: #try strict FORTRAN formatting
        gfile.A52=File.read(52)
        readStatement=lambda handle,varList:read(handle,gfile,varList,File)
        gfile=readStatement(ff.FortranRecordReader('2i4'),['nr','nz'])

        gfile=readStatement(f2020,['xdim','zdim','rcentr','rleft','zmid'])
        gfile=readStatement(f2020,['rmaxis','zmaxis','psimag','psibnd','bcentr'])
        gfile=readStatement(f2020,['current','ssimag','xdum','rmaxis','xdum'])
        gfile=readStatement(f2020,['zmaxis','xdum','ssibry','xdum','xdum'])

        npsi=gfile.nr
        gfile.f_dia=get1DArray(File,npsi)
        gfile.pressure=get1DArray(File,npsi)
        gfile.ffprime=get1DArray(File,npsi)
        gfile.pprime=get1DArray(File,npsi)

        gfile.psi_grid=np.zeros((gfile.nr,gfile.nz))
        stack=[]
        for j in range(gfile.nz):
            for i in range(gfile.nr):
                if(len(stack)==0):
                    stack=f2020.read(File.readline())
                gfile.psi_grid[i][j]=stack.pop(0)

        gfile.q=get1DArray(File,npsi)
        gfile=readStatement(f2022,['n_bnd','limitr'])
        stack=[]
        while(len(gfile.R_bnd)<gfile.n_bnd):
            if(len(stack)==0):
                stack=f2020.read(File.readline())
            gfile.R_bnd.append(stack.pop(0))
            if(len(stack)==0):
                stack=f2020.read(File.readline())
            gfile.Z_bnd.append(stack.pop(0))

        stack=[]
        while(len(gfile.rlim)<gfile.limitr):
            if(len(stack)==0):
                stack=f2020.read(File.readline())
            gfile.rlim.append(stack.pop(0))
            if(len(stack)==0):
                stack=f2020.read(File.readline())
            gfile.zlim.append(stack.pop(0))

        for i in range(npsi):
            gfile.psin.append(float(i)/float(npsi-1))
    except:#failed strict formatting, trying FLUSH formatting but please remove
        #blank lines in file before trying.
        File.seek(0)
        def getInt(Buff):
            return int(Buff.pop(0))
        def getFloat(Buff):
            return float(Buff.pop(0))
        def getFloats(gfile,varList,Buff):
            for el in varList:
                gfile.set(el,getFloat(Buff))
            return gfile
        def get1DArray(Buff,n):
            return [getFloat(Buff) for i in range(n)]
        print('Formatting is dodgy')
        Buff=File.readline()
        gfile.A52=Buff[0:52]
        Buff=Buff[52:].split()
        gfile.nr,gfile.nz=getInt(Buff),getInt(Buff)
        Buff=File.readline().split()
        gfile=getFloats(gfile,['xdim','zdim','rcentr','rleft','zmid'],Buff)
        Buff=File.readline().split()
        gfile=getFloats(gfile,['rmaxis','zmaxis','psimag','psibnd','bcentr'],Buff)
        Buff=File.readline().split()
        gfile=getFloats(gfile,['current','ssimag','xdum','rmaxis','xdum'],Buff)
        Buff=File.readline().split()
        gfile=getFloats(gfile,['zmaxis','xdum','ssibry','xdum','xdum'],Buff)

        #now alas, have to manipulate list instead
        Buff=File.read().split()
        npsi=gfile.nr
        gfile.f_dia=get1DArray(Buff,npsi)
        gfile.pressure=get1DArray(Buff,npsi)
        gfile.ffprime=get1DArray(Buff,npsi)
        gfile.pprime=get1DArray(Buff,npsi)

        gfile.psi_grid=[[0 for i in range(gfile.nr)] for j in range(gfile.nz)]
        for j in range(gfile.nz-1,-1,-1):
            for i in range(gfile.nr):
                gfile.psi_grid[i][j]=getFloat(Buff)

        gfile.q=get1DArray(Buff,npsi)

        gfile.n_bnd=getInt(Buff)
        gfile.limitr=getInt(Buff)

        for i in range(gfile.n_bnd):
            gfile.R_bnd.append(getFloat(Buff))
            gfile.Z_bnd.append(getFloat(Buff))

        for i in range(gfile.limitr):
            gfile.rlim.append(getFloat(Buff))
            gfile.zlim.append(getFloat(Buff))

        for i in range(npsi):
            gfile.psin.append(float(i)/float(npsi-1))
        File.close()

    # Construct axes
    gfile.R = [ gfile.rleft + gfile.xdim * float(i) / float( gfile.nr - 1 ) for i in range( gfile.nr ) ]
    gfile.Z = [ gfile.zmid + gfile.zdim * ( float(i) / float( gfile.nz - 1 ) - 0.5 ) for i in range( gfile.nz ) ]

    # Perform Bivariate Spline of Psi Map to get Derivatives etc.
    gfile.PsiSpline = interpolate.RectBivariateSpline(gfile.R, gfile.Z, gfile.psi_grid, kx=3, ky=3, s=0)
    gfile.FSpline   = interpolate.interp1d(gfile.psin, gfile.f_dia, kind='cubic')
    gfile.QSpline   = interpolate.interp1d(gfile.psin, gfile.q, kind='cubic')

    return gfile

def getGfile_FIESTA(filename):

    import numpy as np
    from io import StringIO

    gfile=Gfile()
    import fortranformat as ff
    f2020=ff.FortranRecordReader('5e16.9')
    f2022=ff.FortranRecordReader('2i5')

    def read(handle,gfile,varList,File):
        line=handle.read(File.readline())
        for pel,el in enumerate(varList):
            gfile.set(el,line[pel])
        return gfile

    def get1DArray(File,npsi):
        res=[]
        while(len(res)<npsi):
            line=f2020.read(File.readline())
            for l in line:
                if(l is not None): #don't append Nones
                    res.append(l)
        return res

    File=open(filename,'r')

    try: #try strict FORTRAN formatting
        gfile.A52=File.read(52)
        readStatement=lambda handle,varList:read(handle,gfile,varList,File)
        gfile=readStatement(ff.FortranRecordReader('2i4'),['nr','nz'])

        gfile=readStatement(f2020,['xdim','zdim','rcentr','rleft','zmid'])
        gfile=readStatement(f2020,['rmaxis','zmaxis','psimag','psibnd','bcentr'])
        gfile=readStatement(f2020,['current','ssimag','xdum','rmaxis','xdum'])
        gfile=readStatement(f2020,['zmaxis','xdum','ssibry','xdum','xdum'])

        npsi=gfile.nr
        gfile.f_dia=get1DArray(File,npsi)
        gfile.pressure=get1DArray(File,npsi)
        gfile.ffprime=get1DArray(File,npsi)
        gfile.pprime=get1DArray(File,npsi)

        gfile.psi_grid=np.zeros((gfile.nr,gfile.nz))
        stack=[]

        def chunk(string,length):
            return(string[0+i:length+i] for i in range(0,len(string),length))

        for j in range(gfile.nz):
            for i in range(gfile.nr):

                if(len(stack)==0):
                    stack=File.readline().strip()
                    stack=list(chunk(stack,16))
                    stack=[ float(x) for x in stack ]

                gfile.psi_grid[i][j]=stack.pop(0)

        gfile.q=get1DArray(File,npsi)

        gfile=readStatement(f2022,['n_bnd','limitr'])

        f2024=ff.FortranRecordReader('2e16.9')

        stack=[]
        while(len(gfile.R_bnd)<gfile.n_bnd):
            if(len(stack)==0):
                stack=f2024.read(File.readline())
            gfile.R_bnd.append(stack.pop(0))
            if(len(stack)==0):
                stack=f2024.read(File.readline())
            gfile.Z_bnd.append(stack.pop(0))

        stack=[]
        while(len(gfile.rlim)<gfile.limitr):
            if(len(stack)==0):
                stack=f2024.read(File.readline())
            gfile.rlim.append(stack.pop(0))
            if(len(stack)==0):
                stack=f2024.read(File.readline())
            gfile.zlim.append(stack.pop(0))

        for i in range(npsi):
            gfile.psin.append(float(i)/float(npsi-1))

    except:

        exit('Error reading GFILE')

    # Construct axes
    gfile.R = [ gfile.rleft + gfile.xdim * float(i) / float( gfile.nr - 1 ) for i in range( gfile.nr ) ]
    gfile.Z = [ gfile.zmid + gfile.zdim * ( float(i) / float( gfile.nz - 1 ) - 0.5 ) for i in range( gfile.nz ) ]

    # Perform Bivariate Spline of Psi Map to get Derivatives etc.
    gfile.PsiSpline = interpolate.RectBivariateSpline(gfile.R, gfile.Z, gfile.psi_grid, kx=3, ky=3, s=0)
    gfile.FSpline   = interpolate.interp1d(gfile.psin, gfile.f_dia, kind='cubic')
    gfile.QSpline   = interpolate.interp1d(gfile.psin, gfile.q, kind='cubic')

    return gfile

def writeGfile(gfile,filename):
        File=open(filename,'w')
        import numpy as np
        from io import StringIO
        import fortranformat as ff
        f2020=ff.FortranRecordWriter('5e16.9')
        f2022=ff.FortranRecordWriter('2i5')

        def writeStatement(handle,varList):
            lst=[]
            for pel,el in enumerate(varList):
                lst.append(gfile.get(el))
            File.write(handle.write(lst))
            File.write("\n")

        def writeArray(handle,variable):
            File.write(handle.write(gfile.get(variable)))
            File.write("\n")

        def writeOrderedPairs(handle,var1,var2):
            longArrayOfPairs=[]
            v1=gfile.get(var1)
            v2=gfile.get(var2)
            for pv,_ in enumerate(v1):
                longArrayOfPairs.append(v1[pv])
                longArrayOfPairs.append(v2[pv])
            #and pretend it's an array
            File.write(handle.write(longArrayOfPairs))
            File.write("\n")

        File.write(gfile.A52[0:52])
        writeStatement(ff.FortranRecordWriter('2i4'),['nr','nz'])
        writeStatement(f2020,['xdim','zdim','rcentr','rleft','zmid'])
        writeStatement(f2020,['rmaxis','zmaxis','psimag','psibnd','bcentr'])
        writeStatement(f2020,['current','ssimag','xdum','rmaxis','xdum'])
        writeStatement(f2020,['zmaxis','xdum','ssibry','xdum','xdum'])

        writeArray(f2020,'f_dia')
        writeArray(f2020,'pressure')
        writeArray(f2020,'ffprime')
        writeArray(f2020,'pprime')

        ###2-D psi profile
        File.write(f2020.write(np.array(gfile.psi_grid).flatten(order='F')))
        File.write("\n")
        #####

        writeArray(f2020,'q')
        writeStatement(f2022,['n_bnd','limitr'])
        writeOrderedPairs(f2020,'R_bnd','Z_bnd')
        writeOrderedPairs(f2020,'rlim','zlim')

        File.close()

def GetPsiAndB(gfile, R, Z):

    Psi    = gfile.PsiSpline.__call__(R, Z, dx=0, dy=0, grid=False)
    dPsidR = gfile.PsiSpline.__call__(R, Z, dx=1, dy=0, grid=False)
    dPsidZ = gfile.PsiSpline.__call__(R, Z, dx=0, dy=1, grid=False)

    psin = ( Psi - gfile.ssimag ) / ( gfile.ssibry - gfile.ssimag )

    BR = -1.0 * dPsidZ / R
    BZ = dPsidR / R
    BT = gfile.bcentr*gfile.rcentr / R

    ### Need to check signs because people are stupid and don't conform to a single standard :/
    # If psi is decreasing then Wesson equations imply Bz < 0 on outer midplane => Bp clockwise and Ip > 0
    # If psi is increasing then Wesson equations imply Bz > 0 on outer midplane => Bp anti-clockwise and Ip < 0
    dpsi = ( gfile.psibnd - gfile.psimag )
    Ip   = gfile.current

    if ( dpsi < 0.0 and Ip < 0.0 ) or ( dpsi > 0.0 and Ip > 0.0 ):
        BR = -1.0 * BR
        BZ = -1.0 * BZ

    if isinstance(psin, float):
        if psin < 0.0: psin = 0.0
        if psin <= 1.0:
            BT = gfile.FSpline(psin) / R
    else:
        psin = np.array(psin)
        psin[psin<0.0] = 0.0
        PlasmaBools = (psin <= np.amax(gfile.psin)) & (psin >= np.amin(gfile.psin))
        BT[PlasmaBools] = gfile.FSpline(psin[PlasmaBools]) / R[PlasmaBools]

    return Psi, psin, BR, BT, BZ


def GetNPsiOnly(gfile, R, Z):

    Psi    = gfile.PsiSpline.__call__(R, Z, dx=0, dy=0, grid=False)
    psin = ( Psi - gfile.ssimag ) / ( gfile.ssibry - gfile.ssimag )

    if isinstance(psin, float):
        if psin < 0.0: psin = 0.0
    else:
        psin = np.array(psin)
        psin[psin<0.0] = 0.0

    return psin


def Make2DPlot( x, y, z, xlab, ylab, title, equal=True ):
    cMap = plt.get_cmap('plasma')
    fig = plt.figure()
    ax = plt.contourf( x, y, z, 100, cmap=cMap, vmin = np.floor(np.amin(z)), vmax = np.ceil(np.amax(z)) )
    plt.xlabel(xlab, fontsize=14)
    plt.ylabel(ylab, fontsize=14)
    cbar = plt.colorbar( cmap=cMap, orientation='vertical')
    cbar.ax.set_ylabel(title,size=14, labelpad=12, rotation=0)
    plt.tick_params(labelsize=14)
    cbar.ax.tick_params(labelsize=14)
    plt.axes().set_aspect('equal')
    plt.show()


def PlotFields(gfile):

    BR_grid=np.zeros((gfile.nr,gfile.nz))
    BT_grid=np.zeros((gfile.nr,gfile.nz))
    BZ_grid=np.zeros((gfile.nr,gfile.nz))

    for j in range(gfile.nz):
        for i in range(gfile.nr):
            Psi, Psin, BR_grid[i][j], BT_grid[i][j], BZ_grid[i][j] = GetPsiAndB(gfile,gfile.R[i],gfile.Z[j])

    Make2DPlot( gfile.R, gfile.Z, gfile.psi_grid.T, "R [m]", "Z [m]", "$\\Psi_p$" )
    Make2DPlot( gfile.R, gfile.Z, BR_grid.T, "R [m]", "Z [m]", "$B_R$ [T]" )
    Make2DPlot( gfile.R, gfile.Z, BT_grid.T, "R [m]", "Z [m]", "$B_\\phi$ [T]" )
    Make2DPlot( gfile.R, gfile.Z, BZ_grid.T, "R [m]", "Z [m]", "$B_Z$ [T]" )


def FindCharacteristicB(gfile):
    #Magnetic field strength on the magnetic axis (Normalised psi = 0.0)
    junk, junk, BR, junk, BZ = GetPsiAndB(gfile, gfile.rmaxis, gfile.zmaxis)
    B0 = np.sqrt((gfile.FSpline(0.0) / gfile.rmaxis)**2.0 + BR**2.0 + BZ**2.0)
    print('B0 = %.3f T'%B0)

    R, Z = np.array(gfile.R), np.array(gfile.Z)
    Z = Z[(R >= min(gfile.R_bnd))*(R <= max(gfile.R_bnd))]
    R = R[(R >= min(gfile.R_bnd))*(R <= max(gfile.R_bnd))]

    Btot = np.zeros([len(Z), len(R)])
    for j, z in enumerate(Z):
        Psi, Psin, BR, BT, BZ = GetPsiAndB(gfile, R, z)
        Btot[j,:] = np.sqrt(BR**2 + BT**2 + BZ**2)

    indMin = np.unravel_index(np.argmin(Btot, axis=None), Btot.shape)
    indMax = np.unravel_index(np.argmax(Btot, axis=None), Btot.shape)
    Bmin, Rmin, Zmin = Btot[indMin], gfile.R[indMin[1]], gfile.Z[indMin[0]]
    Bmax, Rmax, Zmax = Btot[indMax], gfile.R[indMax[1]], gfile.Z[indMax[0]]
    print('Bmin = %.3f T at (R, Z) = (%.3f, %.3f) m.'%(Bmin, Rmin, Zmin) )
    print('Bmax = %.3f T at (R, Z) = (%.3f, %.3f) m.'%(Bmax, Rmax, Zmax) )

    # Need psi(a) and the magnetic field strength on the LFS (B_LFS != Bmin)
    junk, junk, BR, BT, BZ = GetPsiAndB(gfile, np.amax(gfile.R_bnd), 0.0)
    B_LFS = np.sqrt(BR**2.0 + BT**2.0 + BZ**2.0)
    print('B(a_LFS) = %.3f T'%B_LFS)

    return B0, Bmin, Rmin, Zmin, Bmax, Rmax, Zmax, gfile.psibnd, B_LFS
#gfile file reader
#mfitz 2014 and James O from 2020

#Will write and read files according to http://w3.pppl.gov/ntcc/TORAY/G_EQDSK.pdf
#Will also attempt to read FLUSH incorrect eqdsk output (provided that blank lines are first removed by hand)

import numpy as np
import fortranformat as ff
from scipy.interpolate import RectBivariateSpline, interp1d, UnivariateSpline

class Gfile:
    def __init__(self):
        self.A52 = " "*48+"   0"   # = ""
        self.nr = 0
        self.nz = 0
        self.rdim = 0
        self.zdim = 0
        self.rcentr = 0
        self.rleft = 0
        self.zmid = 0
        self.rmaxis = 0
        self.zmaxis = 0
        self.psimag = 0
        self.psibnd = 0
        self.bcentr = 0
        self.current = 0
        self.simag = 0
        self.sibry = 0
        self.f_dia = []
        self.pressure = []
        self.ffprime = []
        self.pprime = []
        self.psi_grid = [[]]
        self.q = []
        self.n_bnd = 0
        self.limitr = 0
        self.R_bnd = []
        self.Z_bnd = []
        self.rlim = []
        self.zlim = []
        self.psin = []
        self.xdum = 0      # just a dummy to mimic eqdsktohelena behavior
        self.PsiSpline = 0 # Bivariate Spline of Psi Map
        self.FSpline = 0   # Spline of F  =  RBphi
        self.QSpline = 0   # Spline of Q 
        self.R = []        # R coordinates
        self.Z = []        # Z coordinates
        self.cocos = None
    def set(self, obj, value):  self.__dict__[obj] = value
    def get(self, obj): return self.__dict__[obj]


def getGfile(filename):
    from io import StringIO
    gfile = Gfile() 
    f2020 = ff.FortranRecordReader('5e16.9')
    f2022 = ff.FortranRecordReader('2i5')

    def read(handle, gfile, varList, File):
        line = handle.read(File.readline())
        for pel, el in enumerate(varList):
            gfile.set(el, line[pel])  
        return gfile

    def get1DArray(File, npsi):
        res = []
        while(len(res)<npsi):
            line = f2020.read(File.readline())
            for l in line:
                if(l is not None): #don't append Nones
                    res.append(l)
        return res
        
    File = open(filename, 'r')
    try: #try strict FORTRAN formatting
        gfile.A52 = File.read(52)
        readStatement = lambda handle, varList:read(handle, gfile, varList, File)
        gfile = readStatement(ff.FortranRecordReader('2i4'), ['nr', 'nz'])
        gfile = readStatement(f2020,['rdim','zdim','rcentr','rleft','zmid'])
        gfile = readStatement(f2020,['rmaxis','zmaxis','psimag','psibnd','bcentr'])
        gfile = readStatement(f2020,['current','simag','xdum','rmaxis','xdum']) 
        gfile = readStatement(f2020,['zmaxis','xdum','sibry','xdum','xdum'])  
    
        npsi = gfile.nr
        gfile.f_dia = get1DArray(File, npsi)
        gfile.pressure = get1DArray(File, npsi)
        gfile.ffprime = get1DArray(File, npsi)
        gfile.pprime = get1DArray(File, npsi)

        gfile.psi_grid = np.zeros((gfile.nr, gfile.nz))
        stack = []
        for j in range(gfile.nz):
            for i in range(gfile.nr):
                if(len(stack) == 0):
                    stack = f2020.read(File.readline())
                gfile.psi_grid[i][j] = stack.pop(0)

        gfile.q = get1DArray(File, npsi)
        gfile = readStatement(f2022,['n_bnd','limitr']) 
        stack = []
        while(len(gfile.R_bnd)<gfile.n_bnd):
            if(len(stack)==0):
                stack = f2020.read(File.readline())
            gfile.R_bnd.append(stack.pop(0))
            if(len(stack)==0):
                stack = f2020.read(File.readline())
            gfile.Z_bnd.append(stack.pop(0))
 
        stack = []
        while(len(gfile.rlim)<gfile.limitr):
            if(len(stack)==0):
                stack = f2020.read(File.readline())
            gfile.rlim.append(stack.pop(0))
            if(len(stack)==0):
                stack = f2020.read(File.readline())
            gfile.zlim.append(stack.pop(0))  
        
        for i in range(npsi):
            gfile.psin.append(float(i)/float(npsi-1))
    except:#failed strict formatting, trying FLUSH formatting but please remove
        #blank lines in file before trying.
        File.seek(0)
        def getInt(Buff):
            return int(Buff.pop(0))
        def getFloat(Buff):
            return float(Buff.pop(0))
        def getFloats(gfile, varList, Buff):
            for el in varList:
                gfile.set(el, getFloat(Buff))
            return gfile
        def get1DArray(Buff, n):
            return [getFloat(Buff) for i in range(n)]
        print('Formatting is dodgy')
        Buff = File.readline()
        gfile.A52 = Buff[0:52]
        Buff = Buff[52:].split()
        gfile.nr, gfile.nz = getInt(Buff), getInt(Buff)
        Buff = File.readline().split()
        gfile = getFloats(gfile,['rdim','zdim','rcentr','rleft','zmid'],Buff)
        Buff = File.readline().split()
        gfile = getFloats(gfile,['rmaxis','zmaxis','psimag','psibnd','bcentr'],Buff)
        Buff = File.readline().split()
        gfile = getFloats(gfile,['current','simag','xdum','rmaxis','xdum'],Buff) 
        Buff = File.readline().split()
        gfile = getFloats(gfile,['zmaxis','xdum','sibry','xdum','xdum'],Buff)  

        #now alas, have to manipulate list instead
        Buff = File.read().split()
        npsi = gfile.nr
        gfile.f_dia = get1DArray(Buff, npsi)
        gfile.pressure = get1DArray(Buff, npsi)
        gfile.ffprime = get1DArray(Buff, npsi)
        gfile.pprime = get1DArray(Buff, npsi)
        
        gfile.psi_grid = [[0 for i in range(gfile.nr)] for j in range(gfile.nz)]
        for j in range(gfile.nz-1, -1, -1):
            for i in range(gfile.nr):
                gfile.psi_grid[i][j] = getFloat(Buff)

        gfile.q = get1DArray(Buff, npsi)

        gfile.n_bnd = getInt(Buff)
        gfile.limitr = getInt(Buff)

        for i in range(gfile.n_bnd):
            gfile.R_bnd.append(getFloat(Buff))
            gfile.Z_bnd.append(getFloat(Buff))
        
        for i in range(gfile.limitr):
            gfile.rlim.append(getFloat(Buff))
            gfile.zlim.append(getFloat(Buff))
         
        for i in range(npsi):
            gfile.psin.append(float(i)/float(npsi-1))
    File.close()

    # Construct axes
    gfile.R = gfile.rleft+gfile.rdim*np.arange(gfile.nr)/float(gfile.nr-1)
    gfile.Z = gfile.zmid+gfile.zdim*(np.arange(gfile.nz)/float(gfile.nz-1)-0.5)

    # Perform Bivariate Spline of Psi Map to get Derivatives etc. 
    gfile.PsiSpline = RectBivariateSpline(gfile.R, gfile.Z, gfile.psi_grid, kx=3, ky=3, s=0)
    gfile.FSpline  = interp1d(gfile.psin, gfile.f_dia, kind='cubic')
    gfile.QSpline  = interp1d(gfile.psin, gfile.q, kind='cubic')
    
    return gfile



def getGfile_FIESTA(filename):

    from io import StringIO

    gfile = Gfile() 
    f2020 = ff.FortranRecordReader('5e16.9')
    f2022 = ff.FortranRecordReader('2i5')

    def read(handle, gfile, varList, File):
        line = handle.read(File.readline())
        for pel, el in enumerate(varList):
            gfile.set(el, line[pel])  
        return gfile

    def get1DArray(File, npsi):
        res = []
        while(len(res)<npsi):
            line = f2020.read(File.readline())
            for l in line:
                if(l is not None): #don't append Nones
                    res.append(l)
        return res
        
    File = open(filename, 'r')

    try: #try strict FORTRAN formatting
        gfile.A52 = File.read(52)
        readStatement = lambda handle, varList:read(handle, gfile, varList, File)
        gfile = readStatement(ff.FortranRecordReader('2i4'),['nr','nz'])

        gfile = readStatement(f2020,['rdim','zdim','rcentr','rleft','zmid'])
        gfile = readStatement(f2020,['rmaxis','zmaxis','psimag','psibnd','bcentr'])
        gfile = readStatement(f2020,['current','simag','xdum','rmaxis','xdum']) 
        gfile = readStatement(f2020,['zmaxis','xdum','sibry','xdum','xdum'])  
    
        npsi = gfile.nr
        gfile.f_dia = get1DArray(File, npsi)
        gfile.pressure = get1DArray(File, npsi)
        gfile.ffprime = get1DArray(File, npsi)
        gfile.pprime = get1DArray(File, npsi)

        gfile.psi_grid = np.zeros((gfile.nr, gfile.nz))
        stack = []

        def chunk(string, length):
            return(string[0+i:length+i] for i in range(0, len(string), length))

        for j in range(gfile.nz):
            for i in range(gfile.nr):

                if(len(stack)==0):
                    stack = File.readline().strip()
                    stack = list(chunk(stack, 16))
                    stack = [ float(x) for x in stack ]

                gfile.psi_grid[i][j] = stack.pop(0)

        gfile.q = get1DArray(File, npsi)

        gfile = readStatement(f2022,['n_bnd','limitr']) 

        f2024 = ff.FortranRecordReader('2e16.9')

        stack = []
        while(len(gfile.R_bnd)<gfile.n_bnd):
            if(len(stack)==0):
                stack = f2024.read(File.readline())
            gfile.R_bnd.append(stack.pop(0))
            if(len(stack)==0):
                stack = f2024.read(File.readline())
            gfile.Z_bnd.append(stack.pop(0))
 
        stack = []
        while(len(gfile.rlim)<gfile.limitr):
            if(len(stack)==0):
                stack = f2024.read(File.readline())
            gfile.rlim.append(stack.pop(0))
            if(len(stack)==0):
                stack = f2024.read(File.readline())
            gfile.zlim.append(stack.pop(0))  

        for i in range(npsi):
            gfile.psin.append(float(i)/float(npsi-1))
    
    except:

        exit('Error reading GFILE')

    # Construct axes
    gfile.R = gfile.rleft + gfile.rdim * np.arange(gfile.nr)/float(gfile.nr - 1)
    gfile.Z = gfile.zmid  + gfile.zdim * (np.arange(gfile.nz)/float(gfile.nz - 1) - 0.5)

    # Perform Bivariate Spline of Psi Map to get Derivatives etc. 
    gfile.PsiSpline = RectBivariateSpline(gfile.R, gfile.Z, gfile.psi_grid, kx=3, ky=3, s=0)
    gfile.FSpline   = interp1d(gfile.psin, gfile.f_dia, kind='cubic')
    gfile.QSpline   = interp1d(gfile.psin, gfile.q, kind='cubic')

    return gfile



def writeGfile(gfile, filename):
    File = open(filename, 'w')
    from io import StringIO
    f2020 = ff.FortranRecordWriter('5e16.9')
    f2022 = ff.FortranRecordWriter('2i5')

    def writeStatement(handle, varList):
        lst = []
        for pel, el in enumerate(varList):
            lst.append(gfile.get(el))
        File.write(handle.write(lst))
        File.write("\n")

    def writeArray(handle, variable):
        File.write(handle.write(gfile.get(variable)))
        File.write("\n")

    def writeOrderedPairs(handle, var1, var2):
        longArrayOfPairs = []
        v1 = gfile.get(var1)
        v2 = gfile.get(var2)
        for pv, _ in enumerate(v1):
            longArrayOfPairs.append(v1[pv])
            longArrayOfPairs.append(v2[pv])
        #and pretend it's an array
        File.write(handle.write(longArrayOfPairs))
        File.write("\n")

    File.write(gfile.A52[0:52])
    writeStatement(ff.FortranRecordWriter('2i4'),['nr','nz'])
    writeStatement(f2020,['rdim','zdim','rcentr','rleft','zmid'])
    writeStatement(f2020,['rmaxis','zmaxis','simag','sibry','bcentr'])
    writeStatement(f2020,['current','simag','xdum','rmaxis','xdum']) 
    writeStatement(f2020,['zmaxis','xdum','sibry','xdum','xdum'])  

    writeArray(f2020,'f_dia')
    writeArray(f2020,'pressure')
    writeArray(f2020,'ffprime')
    writeArray(f2020,'pprime')

    ###2-D psi profile
    File.write(f2020.write(np.array(gfile.psi_grid).flatten(order = 'F')))
    File.write("\n")
    #####

    writeArray(f2020, 'q')
    writeStatement(f2022,['n_bnd','limitr'])
    writeOrderedPairs(f2020,'R_bnd','Z_bnd')
    writeOrderedPairs(f2020,'rlim','zlim')

    File.close()



def GetPsiAndB(gfile, R, Z):

    Psi    = gfile.PsiSpline.__call__(R, Z, dx=0, dy=0, grid=False)
    dPsidR = gfile.PsiSpline.__call__(R, Z, dx=1, dy=0, grid=False)
    dPsidZ = gfile.PsiSpline.__call__(R, Z, dx=0, dy=1, grid=False)

    psin = ( Psi - gfile.simag ) / ( gfile.sibry - gfile.simag )

    BR = -1.0 * dPsidZ / R
    BZ = dPsidR / R
    BT = gfile.bcentr*gfile.rcentr / R

    ### Need to check signs because people are stupid and don't conform to a single standard :/
    # If psi is decreasing then Wesson equations imply Bz < 0 on outer midplane => Bp clockwise and Ip > 0
    # If psi is increasing then Wesson equations imply Bz > 0 on outer midplane => Bp anti-clockwise and Ip < 0
    dpsi = ( gfile.psibnd - gfile.psimag ) 
    Ip   = gfile.current 

    if ( dpsi < 0.0 and Ip < 0.0 ) or ( dpsi > 0.0 and Ip > 0.0 ):
        BR = -1.0 * BR
        BZ = -1.0 * BZ

    if isinstance(psin, float):
        if psin < 0.0: psin = 0.0
        if psin <= 1.0:
            BT = gfile.FSpline(psin) / R
    else:
        psin = np.array(psin)
        psin[psin<0.0] = 0.0
        PlasmaBools = (psin <= np.amax(gfile.psin)) & (psin >= np.amin(gfile.psin))
        BT[PlasmaBools] = gfile.FSpline(psin[PlasmaBools]) / R[PlasmaBools]

    return Psi, psin, BR, BT, BZ


def GetNPsiOnly(gfile, R, Z):

    Psi  = gfile.PsiSpline.__call__(R, Z, dx=0, dy=0, grid=False)
    psin = ( Psi - gfile.simag ) / ( gfile.sibry - gfile.simag )

    if isinstance(psin, float):
        if psin < 0.0: psin = 0.0
    else:
        psin = np.array(psin)
        psin[psin<0.0] = 0.0

    return psin


def FindCharacteristicB(gfile, noPos = True):
    #Magnetic field strength on the magnetic axis (Normalised psi = 0.0)
    Psi0, junk, BR, BT, BZ = GetPsiAndB(gfile, gfile.rmaxis, gfile.zmaxis)
    B0 = np.sqrt(BT**2.0 + BR**2.0 + BZ**2.0)

    Z = gfile.Z[(gfile.R >= min(gfile.R_bnd))*(gfile.R <= max(gfile.R_bnd))]
    R = gfile.R[(gfile.R >= min(gfile.R_bnd))*(gfile.R <= max(gfile.R_bnd))]

    Btot = np.zeros([len(Z), len(R)])
    for j, z in enumerate(Z):
        Psi, Psin, BR, BT, BZ = GetPsiAndB(gfile, R, z)
        Btot[j,:] = np.sqrt(BR**2 + BT**2 + BZ**2)

    indMin = np.unravel_index(np.argmin(Btot, axis=None), Btot.shape)
    indMax = np.unravel_index(np.argmax(Btot, axis=None), Btot.shape)
    Bmin, Rmin, Zmin = Btot[indMin], gfile.R[indMin[1]], gfile.Z[indMin[0]]
    Bmax, Rmax, Zmax = Btot[indMax], gfile.R[indMax[1]], gfile.Z[indMax[0]]

    # Need psi(a) and the magnetic field strength on the LFS (B_LFS != Bmin)
    junk, junk, BR, BT, BZ = GetPsiAndB(gfile, max(gfile.R_bnd), gfile.zmaxis)
    B_LFS = np.sqrt(BR**2.0 + BT**2.0 + BZ**2.0)
    
    if noPos:
        return B0, Bmin, Rmin, Zmin, Bmax, Rmax, Zmax, gfile.psibnd, B_LFS
    else:
        return B0, Bmin, Bmax, gfile.psibnd, B_LFS


##############################################################################
############################# from gfileReader ###############################
##############################################################################

def twogrid(x, y):
    nx, ny = len(x), len(y)
    px, py = np.mgrid[0:nx, 0:ny]
    X, Y = np.zeros(px.shape), np.zeros(py.shape)
    X.flat, Y.flat = x[px.flatten()], y[py.flatten()]
    return X, Y


def threegrid(x, y, z):
    nx, ny, nz = len(x), len(y), len(z)
    px, py, pz = np.mgrid[0:nx, 0:ny, 0:nz]
    X, Y, Z = np.zeros(px.shape), np.zeros(py.shape), np.zeros(pz.shape)
    X.flat, Y.flat, Z.flat = x[px.flatten()], y[py.flatten()], z[pz.flatten()]
    return X, Y, Z


def COCOS(sign_Bp,sign_RφZ, sign_ρθφ):
#Sauter O and Medvedev S Y 2013 Tokamak coordinate conventions: COCOS Comput. Phys. Commun. 184 293302
#Table 1
# sign_RφZ   {1:1,-1:2},{1:3,-1:4},{1:5,-1:6},{1:7,-1:8}
# sign_Bp    1:{1:5,-1:6},{1:1,-1:2}  -1:{1:3,-1:4},{1:7,-1:8}
# sign_Bp,sign_ρθφ,sign_RφZ   {1:{-1:{1:5,-1:6},1:{1:1,-1:2}} , -1:{-1:{1:3,-1:4},1:{1:7,-1:8}}}
    c={1:{-1:{1:5,-1:6},1:{1:1,-1:2}} , -1:{-1:{1:3,-1:4},1:{1:7,-1:8}}}
    try:
        return  c[sign_Bp][sign_ρθφ][sign_RφZ]
    except:
        return 0


def determineCOCOS(gfile, sign_RφZ = 1):
    #cocos is defined by 3 signs - flux direction, toroidal direction, poloidal direction

    #toroidal direction is defined by eqdsk standard as (R φ Z) allowing only cocos 1 3 5 7
    # φ is counter-clockwise
    
    #sign of Bz is physically determined by direction of Ip
    #clockwise (negative) current gives positive outboard Bz
    sign_Bz = -int(sign_RφZ*np.sign(gfile.current))

    # But Bz also depends on flux
    # Bz = -sign_Bp*sign_RφZ*(1/R(dψ/dR))
    sign_Bp = sign_Bz*sign_RφZ*(-np.sign(gfile.sibry-gfile.simag))  

    #now some information on the direction of poloidal angle is required.
    #this is only contained in q(edge) and its relationship to B0 and Ip
    #see discussion in COCOS paper below fig 1. Sign of q, for Ip and B0 going
    #same direction, along with toroidal direction choice, determines q
    sign_ρθφ = sign_RφZ*np.sign(gfile.current)*np.sign(gfile.bcentr)*np.sign(gfile.q[-1])

    return COCOS(sign_Bp, sign_RφZ, sign_ρθφ)


def what_sign_Bp(cocos):
    for sign_Bp in [-1, 1]:
        for sign_RφZ in [-1, 1]:
            for sign_ρθφ in [-1, 1]:
                if(COCOS(sign_Bp, sign_RφZ, sign_ρθφ)==cocos):
                    return sign_Bp
    print("invalid COCOS")
    raise ValueError


def what_sign_RPhiZ(cocos):
    for sign_Bp in [-1, 1]:
        for sign_RφZ in [-1, 1]:
            for sign_ρθφ in [-1, 1]:
                if(COCOS(sign_Bp, sign_RφZ, sign_ρθφ)==cocos):
                    return sign_RφZ
    print("invalid COCOS")
    raise ValueError


def what_sign_rhoThetaPhi(cocos):
    for sign_Bp in [-1, 1]:
        for sign_RφZ in [-1, 1]:
            for sign_ρθφ in [-1, 1]:
                if(COCOS(sign_Bp, sign_RφZ, sign_ρθφ)==cocos):
                    return sign_ρθφ 
    print("invalid COCOS")
    raise ValueError


def flip_sign_Bp(gfile):
    gfile.psibnd   = -np.array(gfile.psibnd)
    gfile.psimag   = -np.array(gfile.psimag)
    gfile.sibry   = -np.array(gfile.sibry)
    gfile.simag   = -np.array(gfile.simag)
    gfile.psi_grid = -np.array(gfile.psi_grid)
    gfile.pprime   = -np.array(gfile.pprime)
    gfile.ffprime  = -np.array(gfile.ffprime)
    return gfile
    

def flip_sign_rhoThetaPhi(gfile):
    gfile.q = -np.array(gfile.q)
    return gfile


def flip_sign_RPhiZ(gfile):
    gfile.psibnd   = -np.array(gfile.psibnd)
    gfile.psimag   = -np.array(gfile.psimag)
    gfile.sibry    = -np.array(gfile.sibry)
    gfile.simag    = -np.array(gfile.simag)
    gfile.psi_grid = -np.array(gfile.psi_grid)
    gfile.pprime   = -np.array(gfile.pprime)
    gfile.ffprime  = -np.array(gfile.ffprime)
    gfile.q        = -np.array(gfile.q)
    gfile.current  = -np.array(gfile.current)
    gfile.bcentr   = -np.array(gfile.bcentr)
    gfile.f_dia    = -np.array(gfile.f_dia)
    return gfile


def convertCOCOS(gfile, target, origin_RφZ = 1):
    origin = determineCOCOS(gfile, sign_RφZ = origin_RφZ)

    origin_Bp  = what_sign_Bp(origin)
    origin_ρθφ = what_sign_rhoThetaPhi(origin)
    target_Bp  = what_sign_Bp(target)
    target_ρθφ = what_sign_rhoThetaPhi(target)
    target_RφZ = what_sign_RPhiZ(target)
    
    if(origin_Bp != target_Bp):
        gfile = flip_sign_Bp(gfile)
    if(origin_ρθφ != target_ρθφ):
        gfile = flip_sign_rhoThetaPhi(gfile)
    if(origin_RφZ != target_RφZ):
        gfile = flip_sign_RPhiZ(gfile)
    return gfile


def printCOCOStable():
    cs = [1,2,3,4,5,6,7,8]
    print("COCOS \t sign_Bp \t sign_RφZ \t sign_ρθφ")
    for cocos in cs:
        sign_Bp = what_sign_Bp(cocos)
        sign_RφZ = what_sign_RPhiZ(cocos)
        sign_ρθφ = what_sign_rhoThetaPhi(cocos)
        print("{} \t\t {} \t\t {} \t\t {}".format(cocos, sign_Bp, sign_RφZ, sign_ρθφ))

def makeSplinesFromGfile(gfile):
    global __psi_spline
    global __F_spline
    global __psi_edge
    global __psi_axis
    global __q_spline
    global __p_spline
    global __Zmag
    global __path
    global __B_R_spline
    global __B_Z_spline
    global __B_Phi_spline
    global __b_R_spline
    global __b_Z_spline
    global __b_Phi_spline
    global __RB_Phi_spline
    global __mu0
    global __sign_Bp
    global __sign_RφZ
    global __sign_ρθφ
    global __modB_spline

    __mu0 = 4*np.pi*1e-7

    if(not gfile.cocos):
        gfile.cocos = determineCOCOS(gfile)

    if(not gfile.cocos):
        print("WARNING - cannot determine COCOS, signs will be wrong")
        gfile.cocos = 3

    __sign_Bp  = what_sign_Bp(gfile.cocos)
    __sign_RφZ = what_sign_RPhiZ(gfile.cocos)
    __sign_ρθφ = what_sign_rhoThetaPhi(gfile.cocos)

    __psi_spline = RectBivariateSpline(gfile.R, gfile.Z, gfile.psi_grid, s = 0)
    __F_spline = UnivariateSpline(gfile.psin, gfile.f_dia, s = 0, ext = 3) 
    __q_spline = UnivariateSpline(gfile.psin, gfile.q, s = 0, ext = 3) 
    __p_spline = UnivariateSpline(gfile.psin, gfile.pressure, s = 0, ext = 1) 

    __psi_axis = gfile.simag
    __psi_edge = gfile.sibry
    __Zmag = gfile.zmaxis
    
    #inside outside detection
    import matplotlib.path as mpath
    bound = np.column_stack((gfile.R_bnd, gfile.Z_bnd))
    __path = mpath.Path(bound)
    
    RR, ZZ = twogrid(gfile.R, gfile.Z)
    B_R, B_Z, B_Phi = [ np.zeros(RR.shape) for i in range(3)]

    B_R.flat = B_R_RZ(RR.flatten(), ZZ.flatten())
    B_Z.flat = B_Z_RZ(RR.flatten(), ZZ.flatten())
    B_Phi.flat = B_Phi_RZ(RR.flatten(), ZZ.flatten())

    #B component splines for derivatives
    __B_R_spline = RectBivariateSpline(gfile.R, gfile.Z, B_R, s = 0)
    __B_Z_spline = RectBivariateSpline(gfile.R, gfile.Z, B_Z, s = 0)
    __B_Phi_spline = RectBivariateSpline(gfile.R, gfile.Z, B_Phi, s = 0)
    __RB_Phi_spline = RectBivariateSpline(gfile.R, gfile.Z, RR*B_Phi, s = 0)

    #drift relevant splines
    modB = np.zeros(RR.shape)
    modB.flat = modB_RZ(RR.flatten(), ZZ.flatten())
    __modB_spline = RectBivariateSpline(gfile.R, gfile.Z, modB, s = 0)
    
    b_R, b_Z, b_Phi = [np.zeros(RR.shape) for i in range(3)]
    b_R.flat = b_R_RZ(RR.flatten(), ZZ.flatten())
    b_Z.flat = b_Z_RZ(RR.flatten(), ZZ.flatten())
    b_Phi.flat = b_Phi_RZ(RR.flatten(), ZZ.flatten())

    __b_R_spline = RectBivariateSpline(gfile.R, gfile.Z, b_R, s = 0)
    __b_Z_spline = RectBivariateSpline(gfile.R, gfile.Z, b_Z, s = 0)
    __b_Phi_spline = RectBivariateSpline(gfile.R, gfile.Z, b_Phi, s = 0)
    

def F_psin(psin):
    return __F_spline(psin)


def P_psin(psin):
    return __p_spline(psin)


def P_RZ(R, Z):
    return P_psin(psin_RZ(R, Z))


def psi_RZ(R, Z):
    return __psi_spline(R, Z, grid = False)


def psin_RZ(R, Z):
    return (__psi_spline(R, Z, grid = False)-__psi_axis)/(__psi_edge-__psi_axis)


def gradPsi_R_RZ(R, Z):
    return __psi_spline(R, Z, dx = 1, dy = 0, grid = False)


def gradPsi_Z_RZ(R, Z):
    return __psi_spline(R, Z, dx = 0, dy = 1, grid = False)


def gradB_R_RZ(R, Z):
    return __modB_spline(R, Z, dx = 1, dy = 0, grid = False)


def gradB_Z_RZ(R, Z):
    return __modB_spline(R, Z, dx = 0, dy = 1, grid = False)


def gradPsi_Z_RZ(R, Z):
    return __psi_spline(R, Z, dx = 0, dy = 1, grid = False)


def B_R_RZ(R, Z):
    return __sign_Bp*__sign_RφZ*(1.0/R)*gradPsi_Z_RZ(R, Z)


def B_Z_RZ(R, Z):
    return -__sign_Bp*__sign_RφZ*(1.0/R)*gradPsi_R_RZ(R, Z)


def B_Phi_RZ(R, Z):
    return F_psin(psin_RZ(R, Z))/R


def modB_RZ(R, Z):
    return (B_R_RZ(R, Z)**2.0+B_Z_RZ(R, Z)**2.0+B_Phi_RZ(R, Z)**2.0)**0.5


def q_R(R):
    return __q_spline(psin_RZ(R, __Zmag))


def isInside(R, Z):
    points = np.column_stack((R, Z))
    return __path.contains_points(points)


def j_R_RZ(R, Z):
    #(curl B)_R
    return -__B_Phi_spline(R, Z, dx = 0, dy = 1, grid = False)/__mu0


def j_Phi_RZ(R, Z):
    #(curl B)_Phi
    return (__B_R_spline(R, Z, dx = 0, dy = 1, grid = False)-__B_Z_spline(R, Z, dx = 1, dy = 0, grid = False))/__mu0


def j_Z_RZ(R, Z):
    #(curl B)_Z
    return (1/R)*__RB_Phi_spline(R, Z, dx = 1, dy = 0, grid = False)/__mu0


def vgradB_R_RZ(R, Z):
    return 0.5*B_Phi_RZ(R, Z)*gradB_Z_RZ(R, Z)/modB_RZ(R, Z)**2


def vgradB_Z_RZ(R, Z):
    return -0.5*B_Phi_RZ(R, Z)*gradB_R_RZ(R, Z)/modB_RZ(R, Z)**2


def vgradB_Phi_RZ(R, Z):
    return 0.5*(B_Z_RZ(R, Z)*gradB_R_RZ(R, Z)-B_R_RZ(R, Z)*gradB_Z_RZ(R, Z))/modB_RZ(R, Z)**2


def jxB_R_RZ(R, Z):
    return j_Phi_RZ(R, Z)*B_Z_RZ(R, Z)-j_Z_RZ(R, Z)*B_Phi_RZ(R, Z)


def jxB_Z_RZ(R, Z):
    return j_R_RZ(R, Z)*B_Phi_RZ(R, Z)-j_Phi_RZ(R, Z)*B_R_RZ(R, Z)


def jxB_Phi_RZ(R, Z):
    return j_Z_RZ(R, Z)*B_R_RZ(R, Z)-j_R_RZ(R, Z)*B_Z_RZ(R, Z)


def b_R_RZ(R, Z):
    return B_R_RZ(R, Z)/modB_RZ(R, Z)


def b_Z_RZ(R, Z):
    return B_Z_RZ(R, Z)/modB_RZ(R, Z)


def b_Phi_RZ(R, Z):
    return B_Phi_RZ(R, Z)/modB_RZ(R, Z)


def kappa_R_RZ(R, Z):
    #last terms are Christoffels!
    return b_R_RZ(R, Z)*__b_R_spline(R, Z, dx = 1, dy = 0, grid = False)+b_Z_RZ(R, Z)*__b_R_spline(R, Z, dx = 0, dy = 1, grid = False)-b_Phi_RZ(R, Z)**2/R


def kappa_Phi_RZ(R, Z):
    #last terms are Christoffels!
    return b_R_RZ(R, Z)*__b_Phi_spline(R, Z, dx = 1, dy = 0, grid = False)+b_Z_RZ(R, Z)*__b_Phi_spline(R, Z, dx = 0, dy = 1, grid = False)+b_Phi_RZ(R, Z)*b_R_RZ(R, Z)/R


def kappa_Z_RZ(R, Z):
    return b_R_RZ(R, Z)*__b_Z_spline(R, Z, dx = 1, dy = 0, grid = False)+b_Z_RZ(R, Z)*__b_Z_spline(R, Z, dx = 0, dy = 1, grid = False)


def vkappa_R_RZ(R, Z):
    return (B_Phi_RZ(R, Z)*kappa_Z_RZ(R, Z)-B_Z_RZ(R, Z)*kappa_Phi_RZ(R, Z))/modB_RZ(R, Z)


def vkappa_Phi_RZ(R, Z):
    return (B_Z_RZ(R, Z)*kappa_R_RZ(R, Z)-B_R_RZ(R, Z)*kappa_Z_RZ(R, Z))/modB_RZ(R, Z)


def vkappa_Z_RZ(R, Z):
    return (B_R_RZ(R, Z)*kappa_Phi_RZ(R, Z)-B_Phi_RZ(R, Z)*kappa_R_RZ(R, Z))/modB_RZ(R, Z)


def writeVtk(gfile, filename, nPhi = 6, gridType = "unstructured", nx = None, ny = None, nz = None):

    RR, ZZ, PP = threegrid(gfile.R, gfile.Z, np.linspace(0, np.pi*2, nPhi))
    
    if(gridType == "unstructured"):
        RR, ZZ, PP = threegrid(gfile.R, gfile.Z, np.linspace(0, np.pi*2, nPhi))
        nPoints = len(RR.flatten())
        XX = RR*np.cos(PP)
        YY = RR*np.sin(PP)
        pointHeader = "UNSTRUCTURED_GRID\nPOINTS {} double\n".format(nPoints)

    elif(gridType == "structured"):
        if(not nx):
            nx = 2*len(gfile.R)
        if(not ny):
            ny = 2*len(gfile.R)
        if(not nz):
            nz = len(gfile.Z) 

        X = np.linspace(-gfile.R[-1], gfile.R[-1], nx)
        Y = np.linspace(-gfile.R[-1], gfile.R[-1], ny)
        Z = np.linspace( gfile.Z[0],  gfile.Z[-1], nz)

        XX, YY, ZZ = threegrid(X, Y, Z)
        RR = np.sqrt(XX**2+YY**2)
        PP = np.arctan2(YY, XX)
        branchCut = np.where(PP.flatten()<0)
        PP.flat[branchCut] = PP.flat[branchCut]+2*np.pi
        nPoints = len(RR.flatten())
        pointHeader = "STRUCTURED_GRID\nDIMENSIONS {} {} {}\nPOINTS {} double\n".format(nx, ny, nz, nPoints)

    else:
        raise ValueError("unsupported grid type: {}".format(gridType))

    scalars = ["Pressure", "PsiPoloidal"]
    vectors = ["B","B_T","B_P","j","j_T","j_P","vgradB","jxB","kappa","vkappa"]
    items = {}
    

    def componentX(comp_Phi, comp_R):
        return -comp_Phi*np.sin(PP.flatten())+comp_R*np.cos(PP.flatten())
    
    def componentY(comp_Phi, comp_R):
        return comp_Phi*np.cos(PP.flatten()) + comp_R*np.sin(PP.flatten())
    

    items['PsiPoloidal'] = psi_RZ(RR.flatten(), ZZ.flatten())
    items['Pressure'] = P_RZ(RR.flatten(), ZZ.flatten())

    items['B_R'] = B_R_RZ(RR.flatten(), ZZ.flatten())
    items['B_Z'] = B_Z_RZ(RR.flatten(), ZZ.flatten())
    items['B_Phi'] = B_Phi_RZ(RR.flatten(), ZZ.flatten())

    items['j_R'] = j_R_RZ(RR.flatten(), ZZ.flatten())
    items['j_Z'] = j_Z_RZ(RR.flatten(), ZZ.flatten())
    items['j_Phi'] = j_Phi_RZ(RR.flatten(), ZZ.flatten())
    
    items['vgradB_R'] = vgradB_R_RZ(RR.flatten(), ZZ.flatten())
    items['vgradB_Z'] = vgradB_Z_RZ(RR.flatten(), ZZ.flatten())
    items['vgradB_Phi'] = vgradB_Phi_RZ(RR.flatten(), ZZ.flatten())

    items['jxB_R'] = jxB_R_RZ(RR.flatten(), ZZ.flatten())
    items['jxB_Z'] = jxB_Z_RZ(RR.flatten(), ZZ.flatten())
    items['jxB_Phi'] = jxB_Phi_RZ(RR.flatten(), ZZ.flatten())

    items['kappa_R'] = kappa_R_RZ(RR.flatten(), ZZ.flatten())
    items['kappa_Z'] = kappa_Z_RZ(RR.flatten(), ZZ.flatten())
    items['kappa_Phi'] = kappa_Phi_RZ(RR.flatten(), ZZ.flatten())

    items['vkappa_R'] = vkappa_R_RZ(RR.flatten(), ZZ.flatten())
    items['vkappa_Z'] = vkappa_Z_RZ(RR.flatten(), ZZ.flatten())
    items['vkappa_Phi'] = vkappa_Phi_RZ(RR.flatten(), ZZ.flatten())

    for vector in ["B", "j", "vgradB", "jxB", "kappa", "vkappa"]:
        for version in ["", "_T", "_P"]:
            for component in ["_X", "_Y", "_Z"]:
                name = vector+version+component
                poloidal = 1.0
                toroidal = 1.0
                if "_T"==version:
                    poloidal = 0.0
                if "_P"==version:
                    toroidal = 0.0
                if "_X"==component:
                    items[name] = componentX(toroidal*items[vector+"_Phi"], poloidal*items[vector+"_R"])
                if "_Y"==component:
                    items[name] = componentY(toroidal*items[vector+"_Phi"], poloidal*items[vector+"_R"])
                if "_Z"==component:
                    items[name] = poloidal*items[vector+"_Z"]
                    
                    
    def sanitize(data):
        data.flat[RR.flatten()<gfile.R[0]] = 0
        data.flat[RR.flatten()>gfile.R[-1]] = 0

    for scalar in scalars:
        item = items[scalar]
        sanitize(item)

    for vector in vectors:
        for component in ["X", "Y", "Z"]:
            name = vector+"_"+component
            item = items[name]
            sanitize(item)

    #include pressure only within LCFS
    items['Pressure'][np.logical_not(isInside(RR.flatten(), ZZ.flatten()))] = 0
    
    with open(filename, "w") as File:
        File.write("# vtk DataFile Version 2.0\n")
        File.write(gfile.A52+"\n")
        File.write("ASCII\n")
        File.write("DATASET "+pointHeader)
        x = XX.flatten()
        y = YY.flatten()
        z = ZZ.flatten()
        for p, _ in enumerate(x):            
            File.write("{} {} {} \n".format(x[p], y[p], z[p]))
        File.write("POINT_DATA {} \n".format(nPoints))
        
        for scalar in scalars:
            item = items[scalar]
            File.write("SCALARS {} double 1\n".format(scalar))
            File.write("LOOKUP_TABLE default\n")
            for p, _ in enumerate(x): 
                File.write("{} \n".format(item[p]))

        for vector in vectors:
            File.write("VECTORS {} double\n".format(vector))
            for p, _ in enumerate(x): 
                for component in ["X", "Y", "Z"]:
                    name = vector+"_"+component
                    item = items[name]
                    File.write("{} ".format(item[p]))
                File.write("\n")


        def replaceLastComma():
            File.seek(File.tell()-1)
            File.write("\n")
            
        scalars = ["PointID", "x", "y", "z", "Pressure", "PsiPoloidal"]
        items["PointID"] = np.arange(len(x))
        items["x"] = x
        items["y"] = y
        items["z"] = z

        with open(filename+".csv", "w") as File:
              
            #titles
            for scalar in scalars:
                File.write(scalar+", ")
            for vector in vectors:
                for component in ["X", "Y", "Z"]:
                    name = vector+"_"+component
                    File.write(name+", ")
   
            replaceLastComma()

            for p, _ in enumerate(x): #for each point
                #go through each scalar name in the local namespace
                for scalar in scalars:
                    item = items[scalar]
                    File.write(str(item[p]))
                    File.write(", ")
                #repeat for each components of the vectors
                for vector in vectors:
                    for component in ["X", "Y", "Z"]:
                        name = vector+"_"+component
                        item = items[name]
                        File.write(str(item[p]))
                        File.write(",")
                #end of line
                replaceLastComma()
