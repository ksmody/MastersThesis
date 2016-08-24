#! /usr/bin/env python
"""
This module provides some functions to automate at your disposal!
"""
import math
import os
from scipy import stats

#---Create a folder and change name------------------------------
def MKchDIR(FolderName):
    """Creates a new directory and changes to this
    MKchDIR("Folder Name")"""
    command = "mkdir %s" % (FolderName)
    os.system(command)
    pfad = os.getcwd() + "/" + FolderName
    os.chdir(pfad)
 
#---cd xy-----------------------------------------------------
def chDir(ABC):
    """Change Directory ABC! ABC is a string and must be present"""
    pfad = os.getcwd() + "/" + ABC
    os.chdir(pfad)

#---cd ..-------------------------------------------------------
def back():
	"""Goes up one directory"""
	back = os.getcwd() +"/.."
	os.chdir(back)

#---Degree of Turbulence-------------------------------------
def TuFun(Re):
	"""Calculates the degree of turbulence of the Reynolds number - Approximation
	Tu=TuFun(Re)"""
	exponent = -1.0/8
	Tun = 0.16*pow(Re, exponent)
	return Tun

#---Viscosity---------------------------------------
def nuFun(Re, r, U):
	"""Calculates the Viscosity for the given Reynolds Number, 
	characteristic length and velocity
	nu=nuFun(Re,r,U)"""
	nun = U*r/Re
	return nun

#---Reynolds Number--------------------------------------
def ReFun(nu, r, U):
	"""Calculates the Reynolds Number for the given Viscosity, Velocity
	and characteristic length
	Re=ReFun(nu,r,U)"""
	Ren = U*r/nu
	return Ren

#---k-Turbulent Kinetic Energy-------------------------------
def kFun(Tu, U):
	"""Calculate k from Tu and Velocity
	k=kFun(Tu,U)"""
	Udash = U*Tu
	kn = 1.5*pow(Udash,2)
	return kn

#---Omega---------------------------------------------
def omegaFun(k, LT):
	"""Calculate omega from k and LT
	omega=omegaFun(k,LT)"""
	rootK = math.sqrt(k)
	omegan = 1.8*rootK/LT
	return omegan

#---Epsilon----------------------------------------------------
def epsilonFun(k, LT):
	"""Calculate epsilon from k and LT
	epsilon=epsilonFun(k,LT)"""
	epsHoch = pow(k, 1.5)
	epsilonn = 0.16*(epsHoch/LT)
	return epsilonn

#---File Customization----------------------------------------------
def randBed(File, Patch, Typ, Val):
	"""randBed(File, Patch, Typ, Val)\n with 
	File eg. U, k, p, is the file to be customized\n
	Patch eg. Inlet, Wall, etc is a word in the File\n
	Typ eg. value or type\n
	Val eg. uniform (50 0 0) or zeroGradient\n 
	the string to be inserted is..."""
	string = "touch %snew" % File
	newFile = "%snew" % File
	os.system(string) #create new File
	file = open(File, "rw+")
	new = open(newFile, "rw+")
	zlr = 0
	pruef = 5
	for line in file:
		if zlr >= 1:
			zlr+=1
		if Patch in line:
			zlr = 1
		if ("}" in line) and (zlr >= 1):
			pruef = 0
			zlr = 0
		if (zlr >= 1) and (pruef != 0):
			if Typ in line:
				if Patch == Typ:
					line = "%s            %s;\n" % (Typ, Val)
					pruef = 0
				else:
					line = "        %s           %s;\n" % (Typ, Val)
		new.write(line)
	file.close()
	new.close()
	string = "rm -rf %s" % File
	os.system(string)
	string2 = "mv %s %s" % (newFile, File)
	os.system(string2)    

#---Pressure Loss Coefficient---------------------------------------
def zetaFun(Re, c_in, c_out, meu, p_in):
    zetan = ((1.0/c_out**2))*((2*p_in)+(c_in**2-c_out**2)-(2*meu))
    return zetan

#---ReynoldsNumber Function-------------------------------------
def Reynolds(ges, Re, c, d, nu):
    """Requires 5 transfer values: x = Reynolds (ges, Re, c, d, nu) \n
     with ges(string) is the sought value from ("Re", "c", "d" or "nu") \n
     and Re as Reynolds number, c as velocity, d as characteristic
     size and nu as kinematic viscosity. \ N
     Handed over this value as 0"""
    if ((ges=="Re") or (ges=="re") or (ges=="Reynolds")):
        return ((c*d)/nu)
    elif ((ges=="c") or (ges=="U") or (ges=="u")):
        return ((Re*nu)/d)
    elif ((ges=="d") or (ges=="D") or (ges=="r") or (ges=="R")):
        return ((Re*nu)/c)
    elif (ges=="nu"):
        return ((c*d)/Re)
    else:
        print "Error, %s not defined" % ges
        print "Possible entries for ges: Re, c, d, nu"

#---Evaluation Picks the Last Element of the File------------------------------------------------        
def evaluationLL(Pfad, File):
    """Evaluation of the entry in the last line of the last column
     file in path"""
    alterPfad = os.getcwd()
    os.chdir(alterPfad + "/" + Pfad)
    file = open(File, "r")
    for l in file:
        line = l
    line = line.strip()
    line = line.split("\t")
    os.chdir(alterPfad)
    return line[len(line)-1]
    
    
#---MESH zur Querschnittsaenderung------------------------------
def netzQuerschnitt(R,r,L,l,w,d1,d2,Xteiler,Yteiler, Pfad):
    """Entsprechend vorbereitetes blockMeshDict muss vorhanden sein\n
    R grosser Rad, r kleiner Rad, L grosse Laenge, l kleine L, w Winkel
    d1 u d2 Wandnahe Unterteilung, Xteiler u Yteiler Zahl der Unterteilungen
    Pfad Pfad zu blockMeshDict \n
    Alle Laengenangaben in mm!!!"""
    
    aktuellerPfad = os.getcwd()
    os.chdir(Pfad)
    
    wBogenHalb = (0.5*w)*(math.pi/180)
    sn = math.sin(wBogenHalb)
    cn = math.cos(wBogenHalb)
    
    x1 = cn*(r-d1-d2)
    y1 = sn*(r-d1-d2)
    x2 = cn*(r-d1)
    y2 = sn*(r-d1)
    x3 = cn*r
    y3 = sn*r
    x4 = cn*(R-d1-d2)
    y4 = sn*(R-d1-d2)
    x5 = cn*(R-d1)
    y5 = sn*(R-d1)
    x6 = cn*R
    y6 = sn*R
    
    z1 = l
    z2 = l+d1
    z3 = l+d1+d2
    z4 = l+L
    
    uxl = 2*L
    uxk = 2*l
    grb = Yteiler
    mtl = int(11*(R-r))
    fn1 = int(grb/3)
    fn2 = 20
    
    grdX1 = 0.2
    grdX2 = 1/grdX1
    grdY1 = 0.5
    grdY2 = 1/grdY1
    grdY3 = 0.5
    grdY4 = 1/grdY3
    
    os.system("touch temp")
    neu = open("temp", "rw+")
    datei = open("blockMeshDict", "rw+")
    for line in datei:
        line = line.replace("x1", str(x1))
        line = line.replace("x2", str(x2))
        line = line.replace("x3", str(x3))
        line = line.replace("x4", str(x4))
        line = line.replace("x5", str(x5))
        line = line.replace("x6", str(x6))
        line = line.replace("y1", str(y1))
        line = line.replace("y2", str(y2))
        line = line.replace("y3", str(y3))
        line = line.replace("y4", str(y4))
        line = line.replace("y5", str(y5))
        line = line.replace("y6", str(y6))
        line = line.replace("z1", str(z1))
        line = line.replace("z2", str(z2))
        line = line.replace("z3", str(z3))
        line = line.replace("z4", str(z4))
        line = line.replace("uxk", str(uxk))
        line = line.replace("uxl", str(uxl))
        line = line.replace("grb", str(grb))
        line = line.replace("mtl", str(mtl))
        line = line.replace("fn1", str(fn1))
        line = line.replace("fn2", str(fn2))
        line = line.replace("grdX1", str(grdX1))
        line = line.replace("grdY1", str(grdY1))
        line = line.replace("grdX2", str(grdX2))
        line = line.replace("grdY2", str(grdY2))
        line = line.replace("grdY3", str(grdY3))
        line = line.replace("grdY4", str(grdY4))
        neu.write(line)
    datei.close()
    neu.close()
    os.system("rm -rf blockMeshDict")
    os.system("mv temp blockMeshDict")
    
    os.chdir(aktuellerPfad)
    
#---Nullmatrix---------------------------------------
def NuMa3(a,b,c):
    return [[[0 for x in range(a)] for y in range(b)] for z in range(c)]  
    
def NuMa2(a,b):
    return [[0 for x in range(a)] for y in range(b)]
    
def NuMa4(a,b,c,d):
    return [NuMa3(a,b,c) for x in range(d)]
    
#---Linear Regression from Measuring Points on Inlet and Outlet---
def linReg(Coordinates, Val, ML, MR, minZ):
    """sL, sR, Rq = linReg(k,w,ml,mr,minz)\n
    linear regression for two areas that def in a matrix
     must be. [Left, Left, Left, ..., right, right, right] \n
     Is the slope on the left and right back \n
     k: Array \ Coordinates with the tag \n
     w: array with the values at the measuring points \n
     ML: Number of measuring points 1 range \n
     MR: Number of measuring points 2. range \n
     minZ: Should the minimum number of points, which will included .. \n
     SL res hill 1. range \n
     Sr: res slope 2. Area \n
     Rq: min Rquadrat Assertiveness"""
    Rql = 0
    Rqr = 0
    for MZ in range(minZ, ML, 1):
        s,a,r,p,e = stats.linregress(Coordinates[:MZ],Val[:MZ])
        if ((r**2) < Rql) and (MZ > minZ):
            break
        Rql = r**2
        sL = s
    for MZ in range(minZ, MR, 1):
        mzh = len(Coordinates) - MZ
        s,a,r,p,e = stats.linregress(Coordinates[mzh:], Val[mzh:])
        if ((r**2)< Rqr) and (MZ > minZ):
            break
        Rqr = r**2
        sR = s
    return sL, sR, min(Rqr, Rql)
    
#---Variable Measuring Point Coordinates------------------
def replace(name, val, pfad, file):
    """Replace in file in pathname name with val\n
    eg. replace("x1", "10", "constant/polymesh","blockMeshDict")"""
    alterPfad = os.getcwd()
    os.chdir(pfad)
    os.system("touch temp")
    neu = open("temp", "rw+")
    file = open(file, "rw+")
    for line in datei:
        line = line.replace(name, val)
        neu.write(line)
    file.close()
    neu.close()
    os.system("rm -rf %s") % file
    os.system("mv temp %s") % file
    os.chdir(alterPfad)

#---m4 File Customization Length----------------------------------------------
def m4SetL(File, Patch, Typ, Val):
	"""m4Set(File, Patch, Typ, Val)\n with 
	File eg. U, k, p, is the file to be customized\n
	Patch eg. Inlet, Wall, etc is a word in the File\n
	Typ eg. value or type\n
	Val eg. uniform (50 0 0) or zeroGradient\n 
	the string to be inserted is..."""
	string = "touch %snew" % File
	newFile = "%snew" % File
	os.system(string) #create new File
	file = open(File, "rw+")
	new = open(newFile, "rw+")
	zlr = 0
	pruef = 5
	for line in file:
		if zlr >= 1:
			zlr+=1
		if Patch in line:
			zlr = 1
		if ("}" in line) and (zlr >= 1):
			pruef = 0
			zlr = 0
		if (zlr >= 1) and (pruef != 0):
			if Typ in line:
				if Patch == Typ:
					line = "m4_define(%s, %s)\n" % (Typ, Val)
					pruef = 0
				else:
					line = "m4_define(%s, %s);\n" % (Typ, Val)
		new.write(line)
	file.close()
	new.close()
	string = "rm -rf %s" % File
	os.system(string)
	string2 = "mv %s %s" % (newFile, File)
	os.system(string2)

#---m4 File Customization Angle----------------------------------------------
def m4SetA(File, Patch, Typ, Val):
	"""m4Set(File, Patch, Typ, Val)\n with 
	File eg. U, k, p, is the file to be customized\n
	Patch eg. Inlet, Wall, etc is a word in the File\n
	Typ eg. value or type\n
	Val eg. uniform (50 0 0) or zeroGradient\n 
	the string to be inserted is..."""
	string = "touch %snew" % File
	newFile = "%snew" % File
	os.system(string) #create new File
	file = open(File, "rw+")
	new = open(newFile, "rw+")
	zlr = 0
	pruef = 5
	for line in file:
		if zlr >= 1:
			zlr+=1
		if Patch in line:
			zlr = 1
		if ("}" in line) and (zlr >= 1):
			pruef = 0
			zlr = 0
		if (zlr >= 1) and (pruef != 0):
			if Typ in line:
				if Patch == Typ:
					line = "m4_define(%s, rad(%s))\n" % (Typ, Val)
					pruef = 0
				else:
					line = "m4_define(%s, rad(%s));\n" % (Typ, Val)
		new.write(line)
	file.close()
	new.close()
	string = "rm -rf %s" % File
	os.system(string)
	string2 = "mv %s %s" % (newFile, File)
	os.system(string2)
