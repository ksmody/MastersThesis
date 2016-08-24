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
