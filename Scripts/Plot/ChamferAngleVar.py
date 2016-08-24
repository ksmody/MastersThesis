# -*- coding: utf-8 -*-
"""
Created on Tue 28.07.2015 11:04:55

@author: modyk

Pressure Drop

"""

import utility as ut
import os
import matplotlib.pyplot as mpl
import math as mt
import numpy as np

def S(x):
    return np.sqrt(x)

Reynolds_range = [500,1000,1600,1900,2500,3500,5000,6000,7000] 
#Reynolds_Fine = [500,1000,1600,1900,2200,2500,3500,4200,5000,6000,7000] 
#Reynolds_VeryFine = [500,1000,1600,1900,2200,2500,3500,4200,5000,6000,7000]
PlotProfile = ['g-s','k-o','c-d','y-^','m-v','r-*']
#PlotProfileA = ['m','g','k','y','r','b']
#PlotProfileB = ['?','s','o','^','*','>']
#Angles = ['8','15','22','30','37','45','53','60','68','75','82','0']
Angles = ['0','8','15','22','30', '37']
#Files = ['0.5','0.75','1','1.25','1.5','2','3']

mpl.figure(1)
#mpl.set_cmap('spring')

for j in range(len(Angles)):
	ut.chDir(Angles[j])
	ut.chDir('1')

### General ###
	p1 = [0 for i in range(len(Reynolds_range))]
	dp = [0 for i in range(len(Reynolds_range))]
	alpha = [0 for i in range(len(Reynolds_range))]

	for i in range(len(Reynolds_range)):
		Re = Reynolds_range[i]
		ut.chDir("Re%s" % Re)
        	nu = ut.Reynolds("nu", Re, 95.5, 2e-3, 0)
    		
        	p1[i] = ut.evaluationLL("postProcessing/p_1/0", "faceSource.dat")
        	dp[i] = (float(p1[i]))
		alpha[i] = 95.5/mt.sqrt(2*dp[i])
		print "%s" %alpha[i]
		ut.back()

	mpl.plot(S(Reynolds_range), alpha, PlotProfile[j])
	ut.back()
	ut.back()
    
#plot 
#mpl.figure(1)
#mpl.plot(S(Reynolds_range), alphac, 'b->')
#mpl.plot(S(Reynolds_Fine), alphaf, 'g-<')
#mpl.plot(S(Reynolds_VeryFine), alphavf, 'r-^')
mpl.grid()
mpl.ylim([0.70,1.00])
mpl.title(r'Chamfer Angle Variation, $0^{\circ}$< $\Theta$ < $37^{\circ}$')
mpl.xlabel(r'Root of Reynolds Number, $\sqrt{Re_d}$ [-]')
mpl.ylabel(r'Discharge Coefficient, $C_D$ [-]')
mpl.legend(("$\Theta = 0^{\circ}$", "$\Theta = 8^{\circ}$", "$\Theta = 15^{\circ}$", "$\Theta = 22^{\circ}$", "$\Theta = 30^{\circ}$", "$\Theta = 37^{\circ}$"), loc='best', prop={'size':9})
mpl.savefig("ChamferAngle2", figsize=(16, 9), dpi=300)
# -*- coding: utf-8 -*-
"""
Created on Tue 28.07.2015 11:04:55

@author: modyk

Pressure Drop

"""

import utility as ut
import os
import matplotlib.pyplot as mpl
import math as mt
import numpy as np

def S(x):
    return np.sqrt(x)

Reynolds_range = [500,1000,1600,1900,2500,3500,5000,6000,7000]
PlotProfileA = ['r--v','g--s','b--o']
PlotProfileB = ['r-v','g-s','b-o']

#PlotProfileA = ['m','g','k','y','r','b']
#PlotProfileB = ['?','s','o','^','*','>']
#Angles = ['8','15','22','30','37','45','53','60','68','75','82','0']

#Angles = ['30','37','45','53','60']
Angles = ['45','53']
Files = ['0.5','1','2']
Filesaaa = ['0.25','0.5','1.0']
Chamfer = ['ChamferX2', 'ChamferX10']

for j in range(len(Angles)):
	ut.chDir(Chamfer[0])
	ut.chDir(Angles[j])

	mpl.figure(j)

	for l in range(len(Files)):
		ut.chDir(Files[l])


		p1 = [0 for i in range(len(Reynolds_range))]
		dp = [0 for i in range(len(Reynolds_range))]
		alpha = [0 for i in range(len(Reynolds_range))]

		for i in range(len(Reynolds_range)):
			Re = Reynolds_range[i]
			ut.chDir("Re%s" % Re)
		   	nu = ut.Reynolds("nu", Re, 95.5, 2e-3, 0)
  		
		  	p1[i] = ut.evaluationLL("postProcessing/p_1/0", "faceSource.dat")
		  	dp[i] = (float(p1[i]))
			alpha[i] = 95.5/mt.sqrt(2*dp[i])
			print "%s" %alpha[i]
			ut.back()

		label = "f/d = 0.1, l/d = %s" % Filesaaa[l]
		mpl.plot(S(Reynolds_range), alpha, PlotProfileA[l], label = label)
		ut.back()

	ut.back()

	ut.back()
	ut.chDir(Chamfer[1])

#for j in range(len(Angles)):
	ut.chDir(Angles[j])
	
	for l in range(len(Files)):
		ut.chDir(Files[l])

		p1x = [0 for i in range(len(Reynolds_range))]
		dpx = [0 for i in range(len(Reynolds_range))]
		alphax = [0 for i in range(len(Reynolds_range))]

		for i in range(len(Reynolds_range)):
			Re = Reynolds_range[i]
			ut.chDir("Re%s" % Re)
		   	nu = ut.Reynolds("nu", Re, 95.5, 2e-3, 0)
 		
		  	p1x[i] = ut.evaluationLL("postProcessing/p_1/0", "faceSource.dat")
		  	dpx[i] = (float(p1x[i]))
			alphax[i] = 95.5/mt.sqrt(2*dpx[i])
			print "%s" %alphax[i]
			ut.back()
		
		label = "f/d = 0.5, l/d = %s" % Filesaaa[l]
		mpl.plot(S(Reynolds_range), alphax, PlotProfileB[l], label=label)
		alphax = 0
		ut.back()

	ut.back()
	ut.back()

	mpl.grid()
	title = "Chamfer Length, $f = 0.2 mm$ vs $f = 1.0 mm$ at $\Theta = %s^{\circ}$" %Angles[j]
	mpl.title(title)
	mpl.ylim([0.80,1.0])
	mpl.legend(loc='best', prop={'size':9})
	mpl.xlabel(r'Root of Reynolds Number, $\sqrt{Re_d}$ [-]')
	mpl.ylabel(r'Discharge Coefficient, $C_D$ [-]')
#	mpl.legend(("l/d = 0.25", "l/d = 0.50", "l/d = 1.0"), loc='best', prop={'size':9})
	mpl.savefig(Angles[j], figsize=(16, 9), dpi=300)
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
ó
½	`Vc           @   sï   d  Z  d d l Z d d l Z d d l m Z d „  Z d „  Z d „  Z d „  Z d „  Z	 d	 „  Z
 d
 „  Z d „  Z d „  Z d „  Z d „  Z d „  Z d „  Z d „  Z d „  Z d „  Z d „  Z d „  Z d „  Z d „  Z d „  Z d S(   sC   
This module provides some functions to automate at your disposal!
iÿÿÿÿN(   t   statsc         C   s<   d |  } t  j | ƒ t  j ƒ  d |  } t  j | ƒ d S(   sF   Creates a new directory and changes to this
    MKchDIR("Folder Name")s   mkdir %st   /N(   t   ost   systemt   getcwdt   chdir(   t
   FolderNamet   commandt   pfad(    (    s:   /home/modyk/Desktop/Simulations/13Variations/newUtility.pyt   MKchDIR
   s    
c         C   s%   t  j ƒ  d |  } t  j | ƒ d S(   s9   Change Directory ABC! ABC is a string and must be presentR   N(   R   R   R   (   t   ABCR   (    (    s:   /home/modyk/Desktop/Simulations/13Variations/newUtility.pyt   chDir   s    c          C   s!   t  j ƒ  d }  t  j |  ƒ d S(   s   Goes up one directorys   /..N(   R   R   R   (   t   back(    (    s:   /home/modyk/Desktop/Simulations/13Variations/newUtility.pyR      s    c         C   s!   d d } d t  |  | ƒ } | S(   sX   Calculates the degree of turbulence of the Reynolds number - Approximation
	Tu=TuFun(Re)g      ð¿i   g{®GázÄ?(   t   pow(   t   Ret   exponentt   Tun(    (    s:   /home/modyk/Desktop/Simulations/13Variations/newUtility.pyt   TuFun   s    
c         C   s   | | |  } | S(   sn   Calculates the Viscosity for the given Reynolds Number, 
	characteristic length and velocity
	nu=nuFun(Re,r,U)(    (   R   t   rt   Ut   nun(    (    s:   /home/modyk/Desktop/Simulations/13Variations/newUtility.pyt   nuFun'   s    c         C   s   | | |  } | S(   sm   Calculates the Reynolds Number for the given Viscosity, Velocity
	and characteristic length
	Re=ReFun(nu,r,U)(    (   t   nuR   R   t   Ren(    (    s:   /home/modyk/Desktop/Simulations/13Variations/newUtility.pyt   ReFun/   s    c         C   s!   | |  } d t  | d ƒ } | S(   s.   Calculate k from Tu and Velocity
	k=kFun(Tu,U)g      ø?i   (   R   (   t   TuR   t   Udasht   kn(    (    s:   /home/modyk/Desktop/Simulations/13Variations/newUtility.pyt   kFun7   s    
c         C   s!   t  j |  ƒ } d | | } | S(   s3   Calculate omega from k and LT
	omega=omegaFun(k,LT)gÍÌÌÌÌÌü?(   t   matht   sqrt(   t   kt   LTt   rootKt   omegan(    (    s:   /home/modyk/Desktop/Simulations/13Variations/newUtility.pyt   omegaFun?   s    c         C   s!   t  |  d ƒ } d | | } | S(   s9   Calculate epsilon from k and LT
	epsilon=epsilonFun(k,LT)g      ø?g{®GázÄ?(   R   (   R   R    t   epsHocht   epsilonn(    (    s:   /home/modyk/Desktop/Simulations/13Variations/newUtility.pyt
   epsilonFunG   s    c         C   si  d |  } d |  } t  j | ƒ t |  d ƒ } t | d ƒ } d } d }	 xÏ | D]Ç }
 | d k rq | d 7} n  | |
 k r† d } n  d |
 k r­ | d k r­ d }	 d } n  | d k r|	 d k r| |
 k r| | k rö d | | f }
 d }	 q	d	 | | f }
 qn  | j |
 ƒ qR W| j ƒ  | j ƒ  d
 |  } t  j | ƒ d | |  f } t  j | ƒ d S(   sî   randBed(File, Patch, Typ, Val)
 with 
	File eg. U, k, p, is the file to be customized

	Patch eg. Inlet, Wall, etc is a word in the File

	Typ eg. value or type

	Val eg. uniform (50 0 0) or zeroGradient
 
	the string to be inserted is...s   touch %snews   %snews   rw+i    i   i   t   }s   %s            %s;
s           %s           %s;
s	   rm -rf %ss   mv %s %sN(   R   R   t   opent   writet   close(   t   Filet   Patcht   Typt   Valt   stringt   newFilet   filet   newt   zlrt   prueft   linet   string2(    (    s:   /home/modyk/Desktop/Simulations/13Variations/newUtility.pyt   randBedO   s8    

			


c         C   s2   d | d d | | d | d d | } | S(   Ng      ð?i   (    (   R   t   c_int   c_outt   meut   p_int   zetan(    (    s:   /home/modyk/Desktop/Simulations/13Variations/newUtility.pyt   zetaFunu   s    .c         C   sÆ   |  d k s$ |  d k s$ |  d k r0 | | | S|  d k sT |  d k sT |  d k r` | | | S|  d k s |  d k s |  d	 k s |  d
 k rœ | | | S|  d k r´ | | | Sd |  GHd GHd S(   s  Requires 5 transfer values: x = Reynolds (ges, Re, c, d, nu) 

     with ges(string) is the sought value from ("Re", "c", "d" or "nu") 

     and Re as Reynolds number, c as velocity, d as characteristic
     size and nu as kinematic viscosity. \ N
     Handed over this value as 0R   t   ret   Reynoldst   cR   t   ut   dt   DR   t   RR   s   Error, %s not defineds&   Possible entries for ges: Re, c, d, nuN(    (   t   gesR   R@   RB   R   (    (    s:   /home/modyk/Desktop/Simulations/13Variations/newUtility.pyR?   z   s    $$0	c         C   s   t  j ƒ  } t  j | d |  ƒ t | d ƒ } x | D] } | } q7 W| j ƒ  } | j d ƒ } t  j | ƒ | t | ƒ d S(   sM   Evaluation of the entry in the last line of the last column
     file in pathR   R   s   	i   (   R   R   R   R(   t   stript   splitt   len(   t   PfadR+   t	   alterPfadR1   t   lR5   (    (    s:   /home/modyk/Desktop/Simulations/13Variations/newUtility.pyt   evaluationLL   s    
c
   -      C   s¸  t  j ƒ  }
 t  j |	 ƒ d | t j d } t j | ƒ } t j | ƒ } | | | | } | | | | } | | | } | | | } | | } | | } | |  | | } | |  | | } | |  | } | |  | } | |  } | |  } | } | | } | | | } | | } d | } d | } | }  t d |  | ƒ }! t |  d ƒ }" d }# d }$ d |$ }% d }& d |& }' d }( d |( }) t  j d	 ƒ t	 d
 d ƒ }* t	 d d ƒ }+ x»|+ D]³}, |, j
 d t | ƒ ƒ }, |, j
 d t | ƒ ƒ }, |, j
 d t | ƒ ƒ }, |, j
 d t | ƒ ƒ }, |, j
 d t | ƒ ƒ }, |, j
 d t | ƒ ƒ }, |, j
 d t | ƒ ƒ }, |, j
 d t | ƒ ƒ }, |, j
 d t | ƒ ƒ }, |, j
 d t | ƒ ƒ }, |, j
 d t | ƒ ƒ }, |, j
 d t | ƒ ƒ }, |, j
 d t | ƒ ƒ }, |, j
 d t | ƒ ƒ }, |, j
 d t | ƒ ƒ }, |, j
 d t | ƒ ƒ }, |, j
 d t | ƒ ƒ }, |, j
 d t | ƒ ƒ }, |, j
 d t |  ƒ ƒ }, |, j
 d  t |! ƒ ƒ }, |, j
 d! t |" ƒ ƒ }, |, j
 d" t |# ƒ ƒ }, |, j
 d# t |$ ƒ ƒ }, |, j
 d$ t |& ƒ ƒ }, |, j
 d% t |% ƒ ƒ }, |, j
 d& t |' ƒ ƒ }, |, j
 d' t |( ƒ ƒ }, |, j
 d( t |) ƒ ƒ }, |* j |, ƒ qÂW|+ j ƒ  |* j ƒ  t  j d) ƒ t  j d* ƒ t  j |
 ƒ d+ S(,   s  Entsprechend vorbereitetes blockMeshDict muss vorhanden sein

    R grosser Rad, r kleiner Rad, L grosse Laenge, l kleine L, w Winkel
    d1 u d2 Wandnahe Unterteilung, Xteiler u Yteiler Zahl der Unterteilungen
    Pfad Pfad zu blockMeshDict 

    Alle Laengenangaben in mm!!!g      à?i´   i   i   i   i   gš™™™™™É?i   s
   touch tempt   temps   rw+t   blockMeshDictt   x1t   x2t   x3t   x4t   x5t   x6t   y1t   y2t   y3t   y4t   y5t   y6t   z1t   z2t   z3t   z4t   uxkt   uxlt   grbt   mtlt   fn1t   fn2t   grdX1t   grdY1t   grdX2t   grdY2t   grdY3t   grdY4s   rm -rf blockMeshDicts   mv temp blockMeshDictN(   R   R   R   R   t   pit   sint   cost   intR   R(   t   replacet   strR)   R*   (-   RD   R   t   LRK   t   wt   d1t   d2t   Xteilert   YteilerRI   t   aktuellerPfadt
   wBogenHalbt   snt   cnRO   RU   RP   RV   RQ   RW   RR   RX   RS   RY   RT   RZ   R[   R\   R]   R^   R`   R_   Ra   Rb   Rc   Rd   Re   Rg   Rf   Rh   Ri   Rj   t   neut   dateiR5   (    (    s:   /home/modyk/Desktop/Simulations/13Variations/newUtility.pyt   netzQuerschnittœ   sŽ    












c         C   sO   g  t  | ƒ D]> } g  t  | ƒ D]% } g  t  |  ƒ D] } d ^ q3 ^ q  ^ q S(   Ni    (   t   range(   t   at   bR@   t   zt   yt   x(    (    s:   /home/modyk/Desktop/Simulations/13Variations/newUtility.pyt   NuMa3ó   s    c         C   s6   g  t  | ƒ D]% } g  t  |  ƒ D] } d ^ q  ^ q S(   Ni    (   R~   (   R   R€   R‚   Rƒ   (    (    s:   /home/modyk/Desktop/Simulations/13Variations/newUtility.pyt   NuMa2ö   s    c         C   s)   g  t  | ƒ D] } t |  | | ƒ ^ q S(   N(   R~   R„   (   R   R€   R@   RB   Rƒ   (    (    s:   /home/modyk/Desktop/Simulations/13Variations/newUtility.pyt   NuMa4ù   s    c         C   s  d } d } xs t  | | d ƒ D]_ } t j |  |  | |  ƒ \ } }	 }
 } } |
 d | k  rn | | k rn Pn  |
 d } | } q Wxƒ t  | | d ƒ D]o } t |  ƒ | } t j |  | | | ƒ \ } }	 }
 } } |
 d | k  rô | | k rô Pn  |
 d } | } q• W| | t | | ƒ f S(   s/  sL, sR, Rq = linReg(k,w,ml,mr,minz)

    linear regression for two areas that def in a matrix
     must be. [Left, Left, Left, ..., right, right, right] 

     Is the slope on the left and right back 

     k: Array \ Coordinates with the tag 

     w: array with the values at the measuring points 

     ML: Number of measuring points 1 range 

     MR: Number of measuring points 2. range 

     minZ: Should the minimum number of points, which will included .. 

     SL res hill 1. range 

     Sr: res slope 2. Area 

     Rq: min Rquadrat Assertivenessi    i   i   (   R~   R    t
   linregressRH   t   min(   t   CoordinatesR.   t   MLt   MRt   minZt   Rqlt   Rqrt   MZt   sR   R   t   pt   et   sLt   mzht   sR(    (    s:   /home/modyk/Desktop/Simulations/13Variations/newUtility.pyt   linRegý   s     )

)

c         C   s»   t  j ƒ  } t  j | ƒ t  j d ƒ t d d ƒ } t | d ƒ } x- t D]% } | j |  | ƒ } | j | ƒ qK W| j ƒ  | j ƒ  t  j d ƒ | t  j d ƒ | t  j | ƒ d S(   sk   Replace in file in pathname name with val

    eg. replace("x1", "10", "constant/polymesh","blockMeshDict")s
   touch tempRM   s   rw+s	   rm -rf %ss
   mv temp %sN(	   R   R   R   R   R(   R|   Ro   R)   R*   (   t   namet   valR   R1   RJ   R{   R5   (    (    s:   /home/modyk/Desktop/Simulations/13Variations/newUtility.pyRo     s    

c         C   si  d |  } d |  } t  j | ƒ t |  d ƒ } t | d ƒ } d } d }	 xÏ | D]Ç }
 | d k rq | d 7} n  | |
 k r† d } n  d |
 k r­ | d k r­ d }	 d } n  | d k r|	 d k r| |
 k r| | k rö d | | f }
 d }	 q	d	 | | f }
 qn  | j |
 ƒ qR W| j ƒ  | j ƒ  d
 |  } t  j | ƒ d | |  f } t  j | ƒ d S(   sì   m4Set(File, Patch, Typ, Val)
 with 
	File eg. U, k, p, is the file to be customized

	Patch eg. Inlet, Wall, etc is a word in the File

	Typ eg. value or type

	Val eg. uniform (50 0 0) or zeroGradient
 
	the string to be inserted is...s   touch %snews   %snews   rw+i    i   i   R'   s   m4_define(%s, %s)
s   m4_define(%s, %s);
s	   rm -rf %ss   mv %s %sN(   R   R   R(   R)   R*   (   R+   R,   R-   R.   R/   R0   R1   R2   R3   R4   R5   R6   (    (    s:   /home/modyk/Desktop/Simulations/13Variations/newUtility.pyt   m4SetL.  s8    

			


c         C   si  d |  } d |  } t  j | ƒ t |  d ƒ } t | d ƒ } d } d }	 xÏ | D]Ç }
 | d k rq | d 7} n  | |
 k r† d } n  d |
 k r­ | d k r­ d }	 d } n  | d k r|	 d k r| |
 k r| | k rö d | | f }
 d }	 q	d	 | | f }
 qn  | j |
 ƒ qR W| j ƒ  | j ƒ  d
 |  } t  j | ƒ d | |  f } t  j | ƒ d S(   sì   m4Set(File, Patch, Typ, Val)
 with 
	File eg. U, k, p, is the file to be customized

	Patch eg. Inlet, Wall, etc is a word in the File

	Typ eg. value or type

	Val eg. uniform (50 0 0) or zeroGradient
 
	the string to be inserted is...s   touch %snews   %snews   rw+i    i   i   R'   s   m4_define(%s, rad(%s))
s   m4_define(%s, rad(%s));
s	   rm -rf %ss   mv %s %sN(   R   R   R(   R)   R*   (   R+   R,   R-   R.   R/   R0   R1   R2   R3   R4   R5   R6   (    (    s:   /home/modyk/Desktop/Simulations/13Variations/newUtility.pyt   m4SetAT  s8    

			


(   t   __doc__R   R   t   scipyR    R	   R   R   R   R   R   R   R#   R&   R7   R=   R?   RL   R}   R„   R…   R†   R–   Ro   R™   Rš   (    (    (    s:   /home/modyk/Desktop/Simulations/13Variations/newUtility.pyt   <module>   s0   											&				W						&# -*- coding: utf-8 -*-
"""
Created on Tue Aug 18 15:48:22 2015

@author:modyk

Run Simulation

"""

import newUtility as ut
import os
import matplotlib.pyplot as mpl
import math as mt
import subprocess

Chamfer_Angle = ['0','8','15','22','30']

Throat_Length = ['0.5','1','2','1.5','2.35','0.75']

Reynolds_range = [500,1000,1600,1900,2500,3500,5000,6000,7000]

for i in range(len(Chamfer_Angle)):

	Ang = Chamfer_Angle[i]
	bashCP = "mkdir %s" % Ang
	os.system(bashCP)
	ut.chDir("%s" % Ang)

	for j in range(len(Throat_Length)):

		Len = Throat_Length[j]
		bashCP = "mkdir %s" % Len
		os.system(bashCP)
		ut.chDir("%s" % Len)
		bashCP = "cp -r ../../TestCase%s TestCase" % Ang
		os.system(bashCP)
		bashCP = "cp ../../plotAlpha.py plotAlpha.py"
		os.system(bashCP)
		bashCP = "cp ../../newUtility.py utility.py"
		os.system(bashCP)
		ut.chDir("TestCase")
		ut.chDir("constant")
		ut.chDir("polyMesh")
		ut.m4SetL("ArcWedge_CX2_A45_L50.m4", "ThroatLength", "ThroatLength", Len)
		ut.m4SetA("ArcWedge_CX2_A45_L50.m4", "ChamferAngle", "ChamferAngle", Ang)
		bashCP = "m4 -P ArcWedge_CX2_A45_L50.m4 > blockMeshDict"
		os.system(bashCP)
		ut.back()
		ut.back()
#//Not Required		os.system(". $WM_PROJECT_DIR/etc/config/unset.sh")
#//Not Required		os.system(". /opt/openfoam240/etc/bashrc")
#//Not Required		os.system("blockMesh")
		sp = subprocess.Popen(["/bin/bash", "-i", "-c", "wmUNSET; OF240; blockMesh"])
		sp.communicate()
		ut.back()

		for k in range(len(Reynolds_range)):

			Re = Reynolds_range[k]
			bashCP = "cp -r TestCase Re%s" % Re
			os.system(bashCP)
			ut.chDir("Re%s" % Re)
			nu = ut.Reynolds("nu", Re, 95.5, 2e-3, 0)
			ut.chDir("constant")
			nuVal = "nu [ 0 2 -1 0 0 0 0 ] %s" % nu
			ut.randBed("transportProperties", "nu", "nu", nuVal)
			ut.back()

			ut.chDir("0")
			Tu = ut.TuFun(Re/3)
			LT = 0.1*6e-3
			k = ut.kFun(Tu, 10.61)
			omega = ut.omegaFun(k, LT)
			ut.randBed("omega", "internalField", "internalField", "uniform %s" % omega)
			ut.randBed("k", "internalField", "internalField", "uniform %s" % k)
			ut.back()

#//Not Required			os.system(". $WM_PROJECT_DIR/etc/config/unset.sh")
#//Not Required			os.system(". $HOME/OpenFOAM/OpenFOAM-2.3.0/etc/bashrc")
			os.system("decomposePar")
			os.system("mpirun -np 4 simpleFoam -parallel | tee log")
			os.system("reconstructPar")
   
			os.system("rm -r processor*")

			os.system("gnuplot plotP")
			os.system("gnuplot resPlotPng")

			ut.back()

		os.system("python plotAlpha.py")
		ut.back()
	ut.back()
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
ó
äœ¾Uc           @   sÝ   d  Z  d d l Z d d l Z d d l m Z d „  Z d „  Z d „  Z d „  Z d „  Z	 d	 „  Z
 d
 „  Z d „  Z d „  Z d „  Z d „  Z d „  Z d „  Z d „  Z d „  Z d „  Z d „  Z d „  Z d „  Z d S(   sC   
This module provides some functions to automate at your disposal!
iÿÿÿÿN(   t   statsc         C   s<   d |  } t  j | ƒ t  j ƒ  d |  } t  j | ƒ d S(   sF   Creates a new directory and changes to this
    MKchDIR("Folder Name")s   mkdir %st   /N(   t   ost   systemt   getcwdt   chdir(   t
   FolderNamet   commandt   pfad(    (    sA   /home/modyk/Desktop/Simulations/13Variations/ChamferX2/utility.pyt   MKchDIR
   s    
c         C   s%   t  j ƒ  d |  } t  j | ƒ d S(   s9   Change Directory ABC! ABC is a string and must be presentR   N(   R   R   R   (   t   ABCR   (    (    sA   /home/modyk/Desktop/Simulations/13Variations/ChamferX2/utility.pyt   chDir   s    c          C   s!   t  j ƒ  d }  t  j |  ƒ d S(   s   Goes up one directorys   /..N(   R   R   R   (   t   back(    (    sA   /home/modyk/Desktop/Simulations/13Variations/ChamferX2/utility.pyR      s    c         C   s!   d d } d t  |  | ƒ } | S(   sX   Calculates the degree of turbulence of the Reynolds number - Approximation
	Tu=TuFun(Re)g      ð¿i   g{®GázÄ?(   t   pow(   t   Ret   exponentt   Tun(    (    sA   /home/modyk/Desktop/Simulations/13Variations/ChamferX2/utility.pyt   TuFun   s    
c         C   s   | | |  } | S(   sn   Calculates the Viscosity for the given Reynolds Number, 
	characteristic length and velocity
	nu=nuFun(Re,r,U)(    (   R   t   rt   Ut   nun(    (    sA   /home/modyk/Desktop/Simulations/13Variations/ChamferX2/utility.pyt   nuFun'   s    c         C   s   | | |  } | S(   sm   Calculates the Reynolds Number for the given Viscosity, Velocity
	and characteristic length
	Re=ReFun(nu,r,U)(    (   t   nuR   R   t   Ren(    (    sA   /home/modyk/Desktop/Simulations/13Variations/ChamferX2/utility.pyt   ReFun/   s    c         C   s!   | |  } d t  | d ƒ } | S(   s.   Calculate k from Tu and Velocity
	k=kFun(Tu,U)g      ø?i   (   R   (   t   TuR   t   Udasht   kn(    (    sA   /home/modyk/Desktop/Simulations/13Variations/ChamferX2/utility.pyt   kFun7   s    
c         C   s!   t  j |  ƒ } d | | } | S(   s3   Calculate omega from k and LT
	omega=omegaFun(k,LT)gÍÌÌÌÌÌü?(   t   matht   sqrt(   t   kt   LTt   rootKt   omegan(    (    sA   /home/modyk/Desktop/Simulations/13Variations/ChamferX2/utility.pyt   omegaFun?   s    c         C   s!   t  |  d ƒ } d | | } | S(   s9   Calculate epsilon from k and LT
	epsilon=epsilonFun(k,LT)g      ø?g{®GázÄ?(   R   (   R   R    t   epsHocht   epsilonn(    (    sA   /home/modyk/Desktop/Simulations/13Variations/ChamferX2/utility.pyt
   epsilonFunG   s    c         C   si  d |  } d |  } t  j | ƒ t |  d ƒ } t | d ƒ } d } d }	 xÏ | D]Ç }
 | d k rq | d 7} n  | |
 k r† d } n  d |
 k r­ | d k r­ d }	 d } n  | d k r|	 d k r| |
 k r| | k rö d | | f }
 d }	 q	d	 | | f }
 qn  | j |
 ƒ qR W| j ƒ  | j ƒ  d
 |  } t  j | ƒ d | |  f } t  j | ƒ d S(   sî   randBed(File, Patch, Typ, Val)
 with 
	File eg. U, k, p, is the file to be customized

	Patch eg. Inlet, Wall, etc is a word in the File

	Typ eg. value or type

	Val eg. uniform (50 0 0) or zeroGradient
 
	the string to be inserted is...s   touch %snews   %snews   rw+i    i   i   t   }s   %s            %s;
s           %s           %s;
s	   rm -rf %ss   mv %s %sN(   R   R   t   opent   writet   close(   t   Filet   Patcht   Typt   Valt   stringt   newFilet   filet   newt   zlrt   prueft   linet   string2(    (    sA   /home/modyk/Desktop/Simulations/13Variations/ChamferX2/utility.pyt   randBedO   s8    

			


c         C   s2   d | d d | | d | d d | } | S(   Ng      ð?i   (    (   R   t   c_int   c_outt   meut   p_int   zetan(    (    sA   /home/modyk/Desktop/Simulations/13Variations/ChamferX2/utility.pyt   zetaFunu   s    .c         C   sÆ   |  d k s$ |  d k s$ |  d k r0 | | | S|  d k sT |  d k sT |  d k r` | | | S|  d k s |  d k s |  d	 k s |  d
 k rœ | | | S|  d k r´ | | | Sd |  GHd GHd S(   s  Requires 5 transfer values: x = Reynolds (ges, Re, c, d, nu) 

     with ges(string) is the sought value from ("Re", "c", "d" or "nu") 

     and Re as Reynolds number, c as velocity, d as characteristic
     size and nu as kinematic viscosity. \ N
     Handed over this value as 0R   t   ret   Reynoldst   cR   t   ut   dt   DR   t   RR   s   Error, %s not defineds&   Possible entries for ges: Re, c, d, nuN(    (   t   gesR   R@   RB   R   (    (    sA   /home/modyk/Desktop/Simulations/13Variations/ChamferX2/utility.pyR?   z   s    $$0	c         C   s   t  j ƒ  } t  j | d |  ƒ t | d ƒ } x | D] } | } q7 W| j ƒ  } | j d ƒ } t  j | ƒ | t | ƒ d S(   sM   Evaluation of the entry in the last line of the last column
     file in pathR   R   s   	i   (   R   R   R   R(   t   stript   splitt   len(   t   PfadR+   t	   alterPfadR1   t   lR5   (    (    sA   /home/modyk/Desktop/Simulations/13Variations/ChamferX2/utility.pyt   evaluationLL   s    
c
   -      C   s¸  t  j ƒ  }
 t  j |	 ƒ d | t j d } t j | ƒ } t j | ƒ } | | | | } | | | | } | | | } | | | } | | } | | } | |  | | } | |  | | } | |  | } | |  | } | |  } | |  } | } | | } | | | } | | } d | } d | } | }  t d |  | ƒ }! t |  d ƒ }" d }# d }$ d |$ }% d }& d |& }' d }( d |( }) t  j d	 ƒ t	 d
 d ƒ }* t	 d d ƒ }+ x»|+ D]³}, |, j
 d t | ƒ ƒ }, |, j
 d t | ƒ ƒ }, |, j
 d t | ƒ ƒ }, |, j
 d t | ƒ ƒ }, |, j
 d t | ƒ ƒ }, |, j
 d t | ƒ ƒ }, |, j
 d t | ƒ ƒ }, |, j
 d t | ƒ ƒ }, |, j
 d t | ƒ ƒ }, |, j
 d t | ƒ ƒ }, |, j
 d t | ƒ ƒ }, |, j
 d t | ƒ ƒ }, |, j
 d t | ƒ ƒ }, |, j
 d t | ƒ ƒ }, |, j
 d t | ƒ ƒ }, |, j
 d t | ƒ ƒ }, |, j
 d t | ƒ ƒ }, |, j
 d t | ƒ ƒ }, |, j
 d t |  ƒ ƒ }, |, j
 d  t |! ƒ ƒ }, |, j
 d! t |" ƒ ƒ }, |, j
 d" t |# ƒ ƒ }, |, j
 d# t |$ ƒ ƒ }, |, j
 d$ t |& ƒ ƒ }, |, j
 d% t |% ƒ ƒ }, |, j
 d& t |' ƒ ƒ }, |, j
 d' t |( ƒ ƒ }, |, j
 d( t |) ƒ ƒ }, |* j |, ƒ qÂW|+ j ƒ  |* j ƒ  t  j d) ƒ t  j d* ƒ t  j |
 ƒ d+ S(,   s  Entsprechend vorbereitetes blockMeshDict muss vorhanden sein

    R grosser Rad, r kleiner Rad, L grosse Laenge, l kleine L, w Winkel
    d1 u d2 Wandnahe Unterteilung, Xteiler u Yteiler Zahl der Unterteilungen
    Pfad Pfad zu blockMeshDict 

    Alle Laengenangaben in mm!!!g      à?i´   i   i   i   i   gš™™™™™É?i   s
   touch tempt   temps   rw+t   blockMeshDictt   x1t   x2t   x3t   x4t   x5t   x6t   y1t   y2t   y3t   y4t   y5t   y6t   z1t   z2t   z3t   z4t   uxkt   uxlt   grbt   mtlt   fn1t   fn2t   grdX1t   grdY1t   grdX2t   grdY2t   grdY3t   grdY4s   rm -rf blockMeshDicts   mv temp blockMeshDictN(   R   R   R   R   t   pit   sint   cost   intR   R(   t   replacet   strR)   R*   (-   RD   R   t   LRK   t   wt   d1t   d2t   Xteilert   YteilerRI   t   aktuellerPfadt
   wBogenHalbt   snt   cnRO   RU   RP   RV   RQ   RW   RR   RX   RS   RY   RT   RZ   R[   R\   R]   R^   R`   R_   Ra   Rb   Rc   Rd   Re   Rg   Rf   Rh   Ri   Rj   t   neut   dateiR5   (    (    sA   /home/modyk/Desktop/Simulations/13Variations/ChamferX2/utility.pyt   netzQuerschnittœ   sŽ    












c         C   sO   g  t  | ƒ D]> } g  t  | ƒ D]% } g  t  |  ƒ D] } d ^ q3 ^ q  ^ q S(   Ni    (   t   range(   t   at   bR@   t   zt   yt   x(    (    sA   /home/modyk/Desktop/Simulations/13Variations/ChamferX2/utility.pyt   NuMa3ó   s    c         C   s6   g  t  | ƒ D]% } g  t  |  ƒ D] } d ^ q  ^ q S(   Ni    (   R~   (   R   R€   R‚   Rƒ   (    (    sA   /home/modyk/Desktop/Simulations/13Variations/ChamferX2/utility.pyt   NuMa2ö   s    c         C   s)   g  t  | ƒ D] } t |  | | ƒ ^ q S(   N(   R~   R„   (   R   R€   R@   RB   Rƒ   (    (    sA   /home/modyk/Desktop/Simulations/13Variations/ChamferX2/utility.pyt   NuMa4ù   s    c         C   s  d } d } xs t  | | d ƒ D]_ } t j |  |  | |  ƒ \ } }	 }
 } } |
 d | k  rn | | k rn Pn  |
 d } | } q Wxƒ t  | | d ƒ D]o } t |  ƒ | } t j |  | | | ƒ \ } }	 }
 } } |
 d | k  rô | | k rô Pn  |
 d } | } q• W| | t | | ƒ f S(   s/  sL, sR, Rq = linReg(k,w,ml,mr,minz)

    linear regression for two areas that def in a matrix
     must be. [Left, Left, Left, ..., right, right, right] 

     Is the slope on the left and right back 

     k: Array \ Coordinates with the tag 

     w: array with the values at the measuring points 

     ML: Number of measuring points 1 range 

     MR: Number of measuring points 2. range 

     minZ: Should the minimum number of points, which will included .. 

     SL res hill 1. range 

     Sr: res slope 2. Area 

     Rq: min Rquadrat Assertivenessi    i   i   (   R~   R    t
   linregressRH   t   min(   t   CoordinatesR.   t   MLt   MRt   minZt   Rqlt   Rqrt   MZt   sR   R   t   pt   et   sLt   mzht   sR(    (    sA   /home/modyk/Desktop/Simulations/13Variations/ChamferX2/utility.pyt   linRegý   s     )

)

c         C   s»   t  j ƒ  } t  j | ƒ t  j d ƒ t d d ƒ } t | d ƒ } x- t D]% } | j |  | ƒ } | j | ƒ qK W| j ƒ  | j ƒ  t  j d ƒ | t  j d ƒ | t  j | ƒ d S(   sk   Replace in file in pathname name with val

    eg. replace("x1", "10", "constant/polymesh","blockMeshDict")s
   touch tempRM   s   rw+s	   rm -rf %ss
   mv temp %sN(	   R   R   R   R   R(   R|   Ro   R)   R*   (   t   namet   valR   R1   RJ   R{   R5   (    (    sA   /home/modyk/Desktop/Simulations/13Variations/ChamferX2/utility.pyRo     s    

(   t   __doc__R   R   t   scipyR    R	   R   R   R   R   R   R   R#   R&   R7   R=   R?   RL   R}   R„   R…   R†   R–   Ro   (    (    (    sA   /home/modyk/Desktop/Simulations/13Variations/ChamferX2/utility.pyt   <module>   s,   											&				W				/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2.3.0                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant/polyMesh";
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//Parametrized Orifice Geometry
//Run using:
//m4 -P blockMeshDict.m4 > blockMeshDict
//m4 definitions:
m4_changecom(//)m4_changequote([,])
m4_define(calc, [m4_esyscmd(perl -e 'use Math::Trig; printf ($1)')])
m4_define(VCOUNT, 0)
m4_define(vlabel, [[// ]Vertex $1 = VCOUNT m4_define($1, VCOUNT)m4_define([VCOUNT],
m4_incr(VCOUNT))])
//Mathematical constants:
m4_define(pi, 3.14159)
m4_define(rad, [calc($1*pi/180.0)])

//m4 spec:
m4_define(InletHeight, 3)
m4_define(OutletHeight, 3)
m4_define(ThroatHeight, 1)

m4_define(Chamfer, 0.2)

m4_define(ThroatLength, 1)
m4_define(BackLength, 0.9)
m4_define(BackAngle, rad(60.0))
m4_define(OutletLength, 73.8)
m4_define(InletLength, calc(80-OutletLength))

m4_define(Boundary, 0.05)

m4_define(ChamferY, calc(ThroatHeight+Chamfer))
m4_define(ChamferX, calc(-Chamfer))

m4_define(BackY, calc(ThroatHeight+(BackLength*tan(BackAngle))))
m4_define(BackX, calc(ThroatLength+BackLength))


//Arc E Inlet
m4_define(xE, calc(-InletLength))

//Arc F
m4_define(xFu, calc(-Chamfer))
m4_define(xFd, calc(-7.5*Chamfer))

//Arc G ChamferStart
m4_define(xGu, calc(-Chamfer))
m4_define(xGd, calc(-2.75*Chamfer))

//Arc H ThroatStart
m4_define(xHu, 0)
m4_define(xHd, calc(-1.25*Chamfer))

//Arc A ThroatEnd
m4_define(xAu, calc(ThroatLength))
m4_define(xAd, calc(1.25*ThroatLength))

//Arc B BackEnd
m4_define(xBu, calc(BackX))
m4_define(xBd, calc(1.25*BackX))

//Arc C
m4_define(xCu, calc(BackX))
m4_define(xCd, calc(1.43*BackX))

//Arc D Outlet
m4_define(xD, calc(OutletLength))

//plane z:
m4_define(z, 0.05)


/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2.3.0                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant/polyMesh";
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

convertToMeters 0.001;

vertices
(
	// A
	(xAd 0 -z) vlabel(A0)
	(xAd 0 z) vlabel(A1)
	(xAu ThroatHeight -z) vlabel(A2)
	(xAu ThroatHeight z) vlabel(A3)

	// B
	(xBd 0 -z) vlabel(B0)
	(xBd 0 z) vlabel(B1)
	(xBu BackY -z) vlabel(B2)
	(xBu BackY z) vlabel(B3)

	// C
	(xCd 0 -z) vlabel(C0)
	(xCd 0 z) vlabel(C1)
	(xCu OutletHeight -z) vlabel(C2)
	(xCu OutletHeight z) vlabel(C3)

	// D
	(xD 0 -z) vlabel(D0)
	(xD 0 z) vlabel(D1)
	(xD OutletHeight -z) vlabel(D2)
	(xD OutletHeight z) vlabel(D3)

	// E
	(xE 0 -z) vlabel(E0)
	(xE 0 z) vlabel(E1)
	(xE InletHeight -z) vlabel(E2)
	(xE InletHeight z) vlabel(E3)

	// F
	(xFd 0 -z) vlabel(F0)
	(xFd 0 z) vlabel(F1)
	(xFu InletHeight -z) vlabel(F2)
	(xFu InletHeight z) vlabel(F3)

	// G
	(xGd 0 -z) vlabel(G0)
	(xGd 0 z) vlabel(G1)
	(xGu ChamferY -z) vlabel(G2)
	(xGu ChamferY z) vlabel(G3)

	// H
	(xHd 0 -z) vlabel(H0)
	(xHd 0 z) vlabel(H1)
	(xHu ThroatHeight -z) vlabel(H2)
	(xHu ThroatHeight z) vlabel(H3)
);

edges
(
	arc F1 F3 (-1.4 0.6 z)
	arc F0 F2 (-1.4 0.6 -z)
	arc G1 G3 (-0.27 1.1 z)
	arc G0 G2 (-0.27 1.1 -z)
	arc H1 H3 (-0.15 0.6 z)
	arc H0 H2 (-0.15 0.6 -z)
	arc A1 A3 (1.15 0.6 z)
	arc A0 A2 (1.15 0.6 -z)
	arc B1 B3 (2.3 1 z)
	arc B0 B2 (2.3 1 -z)
	arc C1 C3 (2.6 1 z)
	arc C0 C2 (2.6 1 -z)
);

blocks
(
	hex (A0 B0 B2 A2 A1 B1 B3 A3) (100 120 1)
	simpleGrading
	(
		1
		(
			(0.95 0.8 1)
			(0.05 0.2 0.25)

		)
		1
	)

	hex (B0 C0 C2 B2 B1 C1 C3 B3) (30 120 1)
	simpleGrading
	(
		1
		(
			(0.95 0.8 1)
			(0.05 0.2 0.25)

		)
		1
	)
	hex (C0 D0 D2 C2 C1 D1 D3 C3) (1000 120 1)
	simpleGrading
	(
		7
		(
			(0.95 0.8 1)
			(0.05 0.2 0.25)

		)
		1
	)
	hex (E0 F0 F2 E2 E1 F1 F3 E3) (160 120 1)
	simpleGrading
	(
		0.3
		(
			(0.95 0.8 1)
			(0.05 0.2 0.25)

		)
		1
	)
	hex (F0 G0 G2 F2 F1 G1 G3 F3) (100 120 1)
	simpleGrading
	(
		1
		(
			(0.95 0.8 1)
			(0.05 0.2 0.25)

		)
		1
	)
	hex (G0 H0 H2 G2 G1 H1 H3 G3) (40 120 1)
	simpleGrading
	(
		1
		(
			(0.95 0.8 1)
			(0.05 0.2 0.25)

		)
		1
	)
	hex (H0 A0 A2 H2 H1 A1 A3 H3) (120 120 1)
	simpleGrading
	(
		1
		(
			(0.95 0.8 1)
			(0.05 0.2 0.25)

		)
		1
	)
);

boundary
(
	inlet
	{
            type inlet;

            faces
            (
                (E1 E3 E2 E0)
            );
        }
        
	outlet
	{
            type outlet;

            faces
            (
                (D0 D2 D3 D1)
            );
        }
        
	fixedwall
	{
            type wall;
            
            faces
            (
                (E3 F3 F2 E2)
                (F3 G3 G2 F2)
                (G3 H3 H2 G2)
                (H3 A3 A2 H2)
                (A3 B3 B2 A2)
                (B3 C3 C2 B2)
                (C3 D3 D2 C2)
            );
        }
	wedge1
	{
            type cyclic;
            
            faces
            (
                (E0 E2 F2 F0)
                (F0 F2 G2 G0)
                (G0 G2 H2 H0)
                (H0 H2 A2 A0)
                (A0 A2 B2 B0)
                (B0 B2 C2 C0)
                (C0 C2 D2 D0)
            );
            
            neighbourPatch  wedge2;
        }
        
	wedge2
	{
            type cyclic;
            
            faces
            (
                (F1 F3 E3 E1)
                (G1 G3 F3 F1)
                (H1 H3 G3 G1)
                (A1 A3 H3 H1)
                (B1 B3 A3 A1)
                (C1 C3 B3 B1)
                (D1 D3 C3 C1)
            );
            neighbourPatch  wedge1;
        }
	axis
	{
            type symmetry;
            
            faces
            (
                (E1 E0 F0 F1)
                (F1 F0 G0 G1)
                (G1 G0 H0 H1)
                (H1 H0 A0 A1)
                (A1 A0 B0 B1)
                (B1 B0 C0 C1)
                (C1 C0 D0 D1)
            );
        }    
);

mergePatchPairs
(
);
// ************************************************************************* //
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2.3.0                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant/polyMesh";
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//Parametrized Orifice Geometry
//Run using:
//m4 -P blockMeshDict.m4 > blockMeshDict
//m4 definitions:
m4_changecom(//)m4_changequote([,])
m4_define(calc, [m4_esyscmd(perl -e 'use Math::Trig; printf ($1)')])
m4_define(VCOUNT, 0)
m4_define(vlabel, [[// ]Vertex $1 = VCOUNT m4_define($1, VCOUNT)m4_define([VCOUNT],
m4_incr(VCOUNT))])
//Mathematical constants:
m4_define(pi, 3.14159)
m4_define(rad, [calc($1*pi/180.0)])

//m4 spec:
m4_define(InletHeight, 3)
m4_define(OutletHeight, 3)
m4_define(ThroatHeight, 1)

m4_define(Chamfer, 0.2)

m4_define(ThroatLength, 1)
m4_define(BackLength, 0.9)
m4_define(BackAngle, rad(60.0))
m4_define(OutletLength, 73.8)
m4_define(InletLength, calc(80-OutletLength))

m4_define(Boundary, 0.05)

m4_define(ChamferY, calc(ThroatHeight+Chamfer))
m4_define(ChamferX, calc(-Chamfer))

m4_define(BackY, calc(ThroatHeight+(BackLength*tan(BackAngle))))
m4_define(BackX, calc(ThroatLength+BackLength))


//Arc E Inlet
m4_define(xE, calc(-InletLength))

//Arc F
m4_define(xFu, calc(-Chamfer))
m4_define(xFd, calc(-7.5*Chamfer))

//Arc G ChamferStart
m4_define(xGu, calc(-Chamfer))
m4_define(xGd, calc(-2.75*Chamfer))

//Arc H ThroatStart
m4_define(xHu, 0)
m4_define(xHd, calc(-1.25*Chamfer))

//Arc A ThroatEnd
m4_define(xAu, calc(ThroatLength))
m4_define(xAd, calc(1.25*ThroatLength))

//Arc B BackEnd
m4_define(xBu, calc(BackX))
m4_define(xBd, calc(1.25*BackX))

//Arc C
m4_define(xCu, calc(BackX))
m4_define(xCd, calc(1.43*BackX))

//Arc D Outlet
m4_define(xD, calc(OutletLength))

//plane z:
m4_define(z, 0.2)


/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2.3.0                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant/polyMesh";
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

convertToMeters 0.001;

vertices
(
	// A
	(xAd 0 -z) vlabel(A0)
	(xAd 0 z) vlabel(A1)
	(xAu ThroatHeight -z) vlabel(A2)
	(xAu ThroatHeight z) vlabel(A3)

	// B
	(xBd 0 -z) vlabel(B0)
	(xBd 0 z) vlabel(B1)
	(xBu BackY -z) vlabel(B2)
	(xBu BackY z) vlabel(B3)

	// C
	(xCd 0 -z) vlabel(C0)
	(xCd 0 z) vlabel(C1)
	(xCu OutletHeight -z) vlabel(C2)
	(xCu OutletHeight z) vlabel(C3)

	// D
	(xD 0 -z) vlabel(D0)
	(xD 0 z) vlabel(D1)
	(xD OutletHeight -z) vlabel(D2)
	(xD OutletHeight z) vlabel(D3)

	// E
	(xE 0 -z) vlabel(E0)
	(xE 0 z) vlabel(E1)
	(xE InletHeight -z) vlabel(E2)
	(xE InletHeight z) vlabel(E3)

	// F
	(xFd 0 -z) vlabel(F0)
	(xFd 0 z) vlabel(F1)
	(xFu InletHeight -z) vlabel(F2)
	(xFu InletHeight z) vlabel(F3)

	// G
	(xGd 0 -z) vlabel(G0)
	(xGd 0 z) vlabel(G1)
	(xGu ChamferY -z) vlabel(G2)
	(xGu ChamferY z) vlabel(G3)

	// H
	(xHd 0 -z) vlabel(H0)
	(xHd 0 z) vlabel(H1)
	(xHu ThroatHeight -z) vlabel(H2)
	(xHu ThroatHeight z) vlabel(H3)
);

edges
(
	arc F1 F3 (-1.4 0.6 z)
	arc F0 F2 (-1.4 0.6 -z)
	arc G1 G3 (-0.27 1.1 z)
	arc G0 G2 (-0.27 1.1 -z)
	arc H1 H3 (-0.15 0.6 z)
	arc H0 H2 (-0.15 0.6 -z)
	arc A1 A3 (1.15 0.6 z)
	arc A0 A2 (1.15 0.6 -z)
	arc B1 B3 (2.3 1 z)
	arc B0 B2 (2.3 1 -z)
	arc C1 C3 (2.6 1 z)
	arc C0 C2 (2.6 1 -z)
);

blocks
(
	hex (A0 B0 B2 A2 A1 B1 B3 A3) (100 100 10)
	simpleGrading
	(
		1
		(
			(0.95 0.8 1)
			(0.05 0.2 0.25)

		)
		1
	)

	hex (B0 C0 C2 B2 B1 C1 C3 B3) (30 100 10)
	simpleGrading
	(
		1
		(
			(0.95 0.8 1)
			(0.05 0.2 0.25)

		)
		1
	)
	hex (C0 D0 D2 C2 C1 D1 D3 C3) (900 100 10)
	simpleGrading
	(
		7
		(
			(0.95 0.8 1)
			(0.05 0.2 0.25)

		)
		1
	)
	hex (E0 F0 F2 E2 E1 F1 F3 E3) (100 100 10)
	simpleGrading
	(
		0.3
		(
			(0.95 0.8 1)
			(0.05 0.2 0.25)

		)
		1
	)
	hex (F0 G0 G2 F2 F1 G1 G3 F3) (100 100 10)
	simpleGrading
	(
		1
		(
			(0.95 0.8 1)
			(0.05 0.2 0.25)

		)
		1
	)
	hex (G0 H0 H2 G2 G1 H1 H3 G3) (40 100 10)
	simpleGrading
	(
		1
		(
			(0.95 0.8 1)
			(0.05 0.2 0.25)

		)
		1
	)
	hex (H0 A0 A2 H2 H1 A1 A3 H3) (70 100 10)
	simpleGrading
	(
		1
		(
			(0.95 0.8 1)
			(0.05 0.2 0.25)

		)
		1
	)
);

boundary
(
	inlet
	{
            type inlet;

            faces
            (
                (E1 E3 E2 E0)
            );
        }
        
	outlet
	{
            type outlet;

            faces
            (
                (D0 D2 D3 D1)
            );
        }
        
	fixedwall
	{
            type wall;
            
            faces
            (
                (E3 F3 F2 E2)
                (F3 G3 G2 F2)
                (G3 H3 H2 G2)
                (H3 A3 A2 H2)
                (A3 B3 B2 A2)
                (B3 C3 C2 B2)
                (C3 D3 D2 C2)
            );
        }
	wedge1
	{
            type cyclic;
            
            faces
            (
                (E0 E2 F2 F0)
                (F0 F2 G2 G0)
                (G0 G2 H2 H0)
                (H0 H2 A2 A0)
                (A0 A2 B2 B0)
                (B0 B2 C2 C0)
                (C0 C2 D2 D0)
            );
            
            neighbourPatch  wedge2;
        }
        
	wedge2
	{
            type cyclic;
            
            faces
            (
                (F1 F3 E3 E1)
                (G1 G3 F3 F1)
                (H1 H3 G3 G1)
                (A1 A3 H3 H1)
                (B1 B3 A3 A1)
                (C1 C3 B3 B1)
                (D1 D3 C3 C1)
            );
            neighbourPatch  wedge1;
        }
	axis
	{
            type symmetry;
            
            faces
            (
                (E1 E0 F0 F1)
                (F1 F0 G0 G1)
                (G1 G0 H0 H1)
                (H1 H0 A0 A1)
                (A1 A0 B0 B1)
                (B1 B0 C0 C1)
                (C1 C0 D0 D1)
            );
        }    
);

mergePatchPairs
(
);
// ************************************************************************* //
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2.3.0                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant/polyMesh";
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//Parametrized Orifice Geometry
//Run using:
//m4 -P blockMeshDict.m4 > blockMeshDict
//m4 definitions:
m4_changecom(//)m4_changequote([,])
m4_define(calc, [m4_esyscmd(perl -e 'use Math::Trig; printf ($1)')])
m4_define(VCOUNT, 0)
m4_define(vlabel, [[// ]Vertex $1 = VCOUNT m4_define($1, VCOUNT)m4_define([VCOUNT],
m4_incr(VCOUNT))])
//Mathematical constants:
m4_define(pi, 3.14159)
m4_define(rad, [calc($1*pi/180.0)])

//m4 spec:
m4_define(InletHeight, 3)
m4_define(OutletHeight, 3)
m4_define(ThroatHeight, 1)

m4_define(Chamfer, 0.2)

m4_define(ThroatLength, 1)
m4_define(BackLength, 0.9)
m4_define(BackAngle, rad(60.0))
m4_define(OutletLength, 73.8)
m4_define(InletLength, calc(80-OutletLength))

m4_define(Boundary, 0.05)

m4_define(ChamferY, calc(ThroatHeight+Chamfer))
m4_define(ChamferX, calc(-Chamfer))

m4_define(BackY, calc(ThroatHeight+(BackLength*tan(BackAngle))))
m4_define(BackX, calc(ThroatLength+BackLength))


//Arc E Inlet
m4_define(xE, calc(-InletLength))

//Arc F
m4_define(xFu, calc(-Chamfer))
m4_define(xFd, calc(-7.5*Chamfer))

//Arc G ChamferStart
m4_define(xGu, calc(-Chamfer))
m4_define(xGd, calc(-2.75*Chamfer))

//Arc H ThroatStart
m4_define(xHu, 0)
m4_define(xHd, calc(-1.25*Chamfer))

//Arc A ThroatEnd
m4_define(xAu, calc(ThroatLength))
m4_define(xAd, calc(1.25*ThroatLength))

//Arc B BackEnd
m4_define(xBu, calc(BackX))
m4_define(xBd, calc(1.25*BackX))

//Arc C
m4_define(xCu, calc(BackX))
m4_define(xCd, calc(1.43*BackX))

//Arc D Outlet
m4_define(xD, calc(OutletLength))

//plane z:
m4_define(z, calc(InletHeight*tan(rad(2.5))))
m4_define(zHA, calc(ThroatHeight*tan(rad(2.5))))
m4_define(zG, calc(ChamferY*tan(rad(2.5))))
m4_define(zB, calc(BackY*tan(rad(2.5))))


/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2.3.0                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant/polyMesh";
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

convertToMeters 0.001;

vertices
(
	// A
	(xAd 0 0) vlabel(A0)
	(xAu ThroatHeight -zHA) vlabel(A2)
	(xAu ThroatHeight zHA) vlabel(A3)

	// B
	(xBd 0 0) vlabel(B0)
	(xBu BackY -zB) vlabel(B2)
	(xBu BackY zB) vlabel(B3)

	// C
	(xCd 0 0) vlabel(C0)
	(xCu OutletHeight -z) vlabel(C2)
	(xCu OutletHeight z) vlabel(C3)

	// D
	(xD 0 0) vlabel(D0)
	(xD OutletHeight -z) vlabel(D2)
	(xD OutletHeight z) vlabel(D3)

	// E
	(xE 0 0) vlabel(E0)
	(xE InletHeight -z) vlabel(E2)
	(xE InletHeight z) vlabel(E3)

	// F
	(xFd 0 0) vlabel(F0)
	(xFu InletHeight -z) vlabel(F2)
	(xFu InletHeight z) vlabel(F3)

	// G
	(xGd 0 0) vlabel(G0)
	(xGu ChamferY -zG) vlabel(G2)
	(xGu ChamferY zG) vlabel(G3)

	// H
	(xHd 0 0) vlabel(H0)
	(xHu ThroatHeight -zHA) vlabel(H2)
	(xHu ThroatHeight zHA) vlabel(H3)
);

edges
(
	arc F0 F3 (-1.4 1.0 0.0436609059828416)
	arc F0 F2 (-1.4 1.0 -0.0436609059828416)
	arc G0 G3 (-0.53 0.133 0.00580690540683)
	arc G0 G2 (-0.53 0.133 -0.00580690540683)
	arc H0 H3 (-0.05 0.9 0.0392948486177)
	arc H0 H2 (-0.05 0.9 -0.0392948486177)
	arc A0 A3 (1.05 0.9 0.0392948486177)
	arc A0 A2 (1.05 0.9 -0.0392948486177)
	arc B0 B3 (2.3 1.0 0.0436609059828416)
	arc B0 B2 (2.3 1.0 -0.0436609059828416)
	arc C0 C3 (2.6 1.0 0.0436609059828416)
	arc C0 C2 (2.6 1.0 -0.0436609059828416)
);

blocks
(
	hex (A0 B0 B2 A2 A0 B0 B3 A3) (100 120 1)
	simpleGrading
	(
		1
		(
			(0.95 0.8 1)
			(0.05 0.2 0.25)

		)
		1
	)

	hex (B0 C0 C2 B2 B0 C0 C3 B3) (30 120 1)
	simpleGrading
	(
		1
		(
			(0.95 0.8 1)
			(0.05 0.2 0.25)

		)
		1
	)
	hex (C0 D0 D2 C2 C0 D0 D3 C3) (1000 120 1)
	simpleGrading
	(
		7
		(
			(0.95 0.8 1)
			(0.05 0.2 0.25)

		)
		1
	)
	hex (E0 F0 F2 E2 E0 F0 F3 E3) (160 120 1)
	simpleGrading
	(
		0.3
		(
			(0.95 0.8 1)
			(0.05 0.2 0.25)

		)
		1
	)
	hex (F0 G0 G2 F2 F0 G0 G3 F3) (100 120 1)
	simpleGrading
	(
		1
		(
			(0.95 0.8 1)
			(0.05 0.2 0.25)

		)
		1
	)
	hex (G0 H0 H2 G2 G0 H0 H3 G3) (40 120 1)
	simpleGrading
	(
		1
		(
			(0.95 0.8 1)
			(0.05 0.2 0.25)

		)
		1
	)
	hex (H0 A0 A2 H2 H0 A0 A3 H3) (120 120 1)
	simpleGrading
	(
		1
		(
			(0.95 0.8 1)
			(0.05 0.2 0.25)

		)
		1
	)
);

boundary
(
	inlet
	{
            type inlet;

            faces
            (
                (E3 E2 E0 E0)
            );
        }
        
	outlet
	{
            type outlet;

            faces
            (
                (D0 D2 D3 D0)
            );
        }
        
	fixedwall
	{
            type wall;
            
            faces
            (
                (E3 F3 F2 E2)
                (F3 G3 G2 F2)
                (G3 H3 H2 G2)
                (H3 A3 A2 H2)
                (A3 B3 B2 A2)
                (B3 C3 C2 B2)
                (C3 D3 D2 C2)
            );
        }
	wedge1
	{
            type wedge;
            
            faces
            (
                (E0 E2 F2 F0)
                (F0 F2 G2 G0)
                (G0 G2 H2 H0)
                (H0 H2 A2 A0)
                (A0 A2 B2 B0)
                (B0 B2 C2 C0)
                (C0 C2 D2 D0)
            );
            
        }
        
	wedge2
	{
            type wedge;
            
            faces
            (
                (F0 F3 E3 E0)
                (G0 G3 F3 F0)
                (H0 H3 G3 G0)
                (A0 A3 H3 H0)
                (B0 B3 A3 A0)
                (C0 C3 B3 B0)
                (D0 D3 C3 C0)
            );
        }
	axis
	{
            type empty;
            
            faces
            (
                (E0 E0 F0 F0)
                (F0 F0 G0 G0)
                (G0 G0 H0 H0)
                (H0 H0 A0 A0)
                (A0 A0 B0 B0)
                (B0 B0 C0 C0)
                (C0 C0 D0 D0)
            );
        }   
);

mergePatchPairs
(
);
// ************************************************************************* //
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2.3.0                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant/polyMesh";
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//Parametrized Orifice Geometry
//Run using:
//m4 -P blockMeshDict.m4 > blockMeshDict
//m4 definitions:
m4_changecom(//)m4_changequote([,])
m4_define(calc, [m4_esyscmd(perl -e 'use Math::Trig; printf ($1)')])
m4_define(VCOUNT, 0)
m4_define(vlabel, [[// ]Vertex $1 = VCOUNT m4_define($1, VCOUNT)m4_define([VCOUNT],
m4_incr(VCOUNT))])
//Mathematical constants:
m4_define(pi, 3.14159)
m4_define(rad, [calc($1*pi/180.0)])

//m4 spec:
m4_define(InletHeight, 3)
m4_define(OutletHeight, 3)
m4_define(ThroatHeight, 1)

m4_define(ChamferX, 0.2)
m4_define(ChamferAngle, rad(45.0))

m4_define(ThroatLength, 0.5)
m4_define(BackLength, 0.9)
m4_define(BackAngle, rad(60.0))
m4_define(OutletLength, 73.8)
m4_define(InletLength, calc(80-OutletLength))

m4_define(Boundary, 0.05)

m4_define(ChamferYCoor, calc(ThroatHeight+(ChamferX*tan(ChamferAngle))))
m4_define(ChamferXCoor, calc(-ChamferX))

m4_define(BackY, calc(ThroatHeight+(BackLength*tan(BackAngle))))
m4_define(BackX, calc(ThroatLength+BackLength))


//Arc E Inlet
m4_define(xE, calc(-InletLength))

//Arc F
m4_define(xFu, calc(-ChamferX))
m4_define(xFd, calc(-7.5*ChamferX))

//Arc G ChamferStart
m4_define(xGu, calc(-ChamferX))
m4_define(xGd, calc(-3.9*ChamferX))

//Arc H ThroatStart
m4_define(xHu, 0)
m4_define(xHd, calc(-1.25*ChamferX))

//Arc A ThroatEnd
m4_define(xAu, calc(ThroatLength))
m4_define(xAd, calc(1.3*ThroatLength))

//Arc B BackEnd
m4_define(xBu, calc(BackX))
m4_define(xBd, calc(1.53*BackX))

//Arc C
m4_define(xCu, calc(BackX))
m4_define(xCd, calc(1.85*BackX))

//Arc D Outlet
m4_define(xD, calc(OutletLength))

//plane z:
m4_define(z, calc(InletHeight*tan(rad(2.5))))
m4_define(zHA, calc(ThroatHeight*tan(rad(2.5))))
m4_define(zG, calc(ChamferYCoor*tan(rad(2.5))))
m4_define(zB, calc(BackY*tan(rad(2.5))))

m4_define(ArcL, 0.1)

//Angles Arc
m4_define(Ang, rad(45.0))
m4_define(AngG, calc((rad(270)-ChamferAngle)/1.8))
m4_define(AngH, calc((rad(180)+ChamferAngle)/1.5))
m4_define(AngA, calc((rad(180)+BackAngle)/1.5))
m4_define(AngB, calc((rad(270)-BackAngle)/1.5))

//Arc Points X
m4_define(ArcFx, calc(ChamferXCoor-(ArcL*sin(Ang))))
m4_define(ArcGx, calc(ChamferXCoor-(ArcL*sin(AngG))))
m4_define(ArcHx, calc(-ArcL*sin(AngH)))
m4_define(ArcAx, calc(ThroatLength+(ArcL*sin(AngA))))
m4_define(ArcBx, calc(BackX+(ArcL*sin(AngB))))
m4_define(ArcCx, calc(BackX+(ArcL*sin(Ang))))

//Arc Points Y
m4_define(ArcFy, calc(InletHeight-(ArcL*cos(Ang))))
m4_define(ArcGy, calc(ChamferYCoor+(ArcL*cos(AngG))))
m4_define(ArcHy, calc(ThroatHeight+(ArcL*cos(AngH))))
m4_define(ArcAy, calc(ThroatHeight+(ArcL*cos(AngA))))
m4_define(ArcBy, calc(BackY+(ArcL*cos(AngB))))
m4_define(ArcCy, calc(OutletHeight-(ArcL*cos(Ang))))

//Arc Points Z
m4_define(zArcF, calc(ArcFy*tan(rad(2.5))))
m4_define(zArcG, calc(ArcGy*tan(rad(2.5))))
m4_define(zArcH, calc(ArcHy*tan(rad(2.5))))
m4_define(zArcA, calc(ArcAy*tan(rad(2.5))))
m4_define(zArcB, calc(ArcBy*tan(rad(2.5))))
m4_define(zArcC, calc(ArcCy*tan(rad(2.5))))

/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2.3.0                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant/polyMesh";
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

convertToMeters 0.001;

vertices
(
	// A
	(xAd 0 0) vlabel(A0)
	(xAu ThroatHeight -zHA) vlabel(A2)
	(xAu ThroatHeight zHA) vlabel(A3)

	// B
	(xBd 0 0) vlabel(B0)
	(xBu BackY -zB) vlabel(B2)
	(xBu BackY zB) vlabel(B3)

	// C
	(xCd 0 0) vlabel(C0)
	(xCu OutletHeight -z) vlabel(C2)
	(xCu OutletHeight z) vlabel(C3)

	// D
	(xD 0 0) vlabel(D0)
	(xD OutletHeight -z) vlabel(D2)
	(xD OutletHeight z) vlabel(D3)

	// E
	(xE 0 0) vlabel(E0)
	(xE InletHeight -z) vlabel(E2)
	(xE InletHeight z) vlabel(E3)

	// F
	(xFd 0 0) vlabel(F0)
	(xFu InletHeight -z) vlabel(F2)
	(xFu InletHeight z) vlabel(F3)

	// G
	(xGd 0 0) vlabel(G0)
	(xGu ChamferYCoor -zG) vlabel(G2)
	(xGu ChamferYCoor zG) vlabel(G3)

	// H
	(xHd 0 0) vlabel(H0)
	(xHu ThroatHeight -zHA) vlabel(H2)
	(xHu ThroatHeight zHA) vlabel(H3)
);

edges
(
	arc F0 F3 (ArcFx ArcFy zArcF)
	arc F0 F2 (ArcFx ArcFy -zArcF)
	arc G0 G3 (ArcGx ArcGy zArcG)
	arc G0 G2 (ArcGx ArcGy -zArcG)
	arc H0 H3 (ArcHx ArcHy zArcH)
	arc H0 H2 (ArcHx ArcHy -zArcH)
	arc A0 A3 (ArcAx ArcAy zArcA)
	arc A0 A2 (ArcAx ArcAy -zArcA)
	arc B0 B3 (ArcBx ArcBy zArcB)
	arc B0 B2 (ArcBx ArcBy -zArcB)
	arc C0 C3 (ArcCx ArcCy zArcC)
	arc C0 C2 (ArcCx ArcCy -zArcC)
);

blocks
(
	hex (A0 B0 B2 A2 A0 B0 B3 A3) (100 120 1)
	simpleGrading
	(
		1
		(
			(0.95 0.8 1)
			(0.05 0.2 0.25)

		)
		1
	)

	hex (B0 C0 C2 B2 B0 C0 C3 B3) (30 120 1)
	simpleGrading
	(
		1
		(
			(0.95 0.8 1)
			(0.05 0.2 0.25)

		)
		1
	)
	hex (C0 D0 D2 C2 C0 D0 D3 C3) (1000 120 1)
	simpleGrading
	(
		7
		(
			(0.95 0.8 1)
			(0.05 0.2 0.25)

		)
		1
	)
	hex (E0 F0 F2 E2 E0 F0 F3 E3) (160 120 1)
	simpleGrading
	(
		0.3
		(
			(0.95 0.8 1)
			(0.05 0.2 0.25)

		)
		1
	)
	hex (F0 G0 G2 F2 F0 G0 G3 F3) (100 120 1)
	simpleGrading
	(
		1
		(
			(0.95 0.8 1)
			(0.05 0.2 0.25)

		)
		1
	)
	hex (G0 H0 H2 G2 G0 H0 H3 G3) (40 120 1)
	simpleGrading
	(
		1
		(
			(0.95 0.8 1)
			(0.05 0.2 0.25)

		)
		1
	)
	hex (H0 A0 A2 H2 H0 A0 A3 H3) (120 120 1)
	simpleGrading
	(
		1
		(
			(0.95 0.8 1)
			(0.05 0.2 0.25)

		)
		1
	)
);

boundary
(
	inlet
	{
            type inlet;

            faces
            (
                (E3 E2 E0 E0)
            );
        }
        
	outlet
	{
            type outlet;

            faces
            (
                (D0 D2 D3 D0)
            );
        }
        
	fixedwall
	{
            type wall;
            
            faces
            (
                (E3 F3 F2 E2)
                (F3 G3 G2 F2)
                (G3 H3 H2 G2)
                (H3 A3 A2 H2)
                (A3 B3 B2 A2)
                (B3 C3 C2 B2)
                (C3 D3 D2 C2)
            );
        }
	wedge1
	{
            type wedge;
            
            faces
            (
                (E0 E2 F2 F0)
                (F0 F2 G2 G0)
                (G0 G2 H2 H0)
                (H0 H2 A2 A0)
                (A0 A2 B2 B0)
                (B0 B2 C2 C0)
                (C0 C2 D2 D0)
            );
            
        }
        
	wedge2
	{
            type wedge;
            
            faces
            (
                (F0 F3 E3 E0)
                (G0 G3 F3 F0)
                (H0 H3 G3 G0)
                (A0 A3 H3 H0)
                (B0 B3 A3 A0)
                (C0 C3 B3 B0)
                (D0 D3 C3 C0)
            );
        }
	axis
	{
            type empty;
            
            faces
            (
                (E0 E0 F0 F0)
                (F0 F0 G0 G0)
                (G0 G0 H0 H0)
                (H0 H0 A0 A0)
                (A0 A0 B0 B0)
                (B0 B0 C0 C0)
                (C0 C0 D0 D0)
            );
        }   
);

mergePatchPairs
(
);
// ************************************************************************* //
