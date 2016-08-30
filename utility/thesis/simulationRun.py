# -*- coding: utf-8 -*-
"""
Created on Tue Aug 18 15:48:22 2015

@author:modyk

Run Simulation Script

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
