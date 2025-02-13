

"""
http://www.scipy.org/Cookbook/Least_Squares_Circle
"""

from numpy import *


import numpy as np
import matplotlib.pyplot as plt
import sys
from dump import dump
from scipy.optimize import curve_fit
import math
from math import sqrt, atan, pi



infile = "/home/devargya/Desktop/electro-wetting/50-emim-0.00/tmp.lammpstrj"


#nd = 1 # number of donor
nd = 2



donor = [1,0,1,2]  # pair of D and H  of water

#donor = [0,1]       #pair of D and H for cation [EMIM]



w = dump(infile)          # for water
w.map(1,"id",2,"type",3,"x",4,"y",5,"z")


w.aselect.test("$type == 12 or $type == 13")  #water in [emim][bf4]
#w.aselect.test("$type == 2 or $type == 5")  #cation in [emim][bf4]
w.delete()
w.sort()
w.write("temp/water.lammpstrj")





a = dump(infile)         # for acceptor of water
a.map(1,"id",2,"type",3,"x",4,"y",5,"z")
a.aselect.test("$type == 11")  #acceptor F in bf4
#a.aselect.test("$type == 1") #acceptor N in cation 
#a.aselect.test("$type == 13")  #acceptor O in H20 for bf4
#a.aselect.test("$type == 14")  #acceptor O in ntf2
#a.aselect.test("$type == 16")  #acceptor O in H2O for ntf2
#a.aselect.test("$type == 2")  #acceptor O in H20 for water
#a.aselect.test("$type == 15")  #acceptor O in H20 for bf4
#a.aselect.test("$type == 13")  #acceptor F in bf4[bmim]
#a.aselect.test("$type == 5")  #acceptor hcr
a.delete()
a.sort()
a.write("temp/acceptor.lammpstrj")



a = dump("temp/acceptor.lammpstrj",1000000)
a.map(1,"id",2,"type",3,"x",4,"y",5,"z")
w = dump("temp/water.lammpstrj",1000000)
w.map(1,"id",2,"type",3,"x",4,"y",5,"z")

rbin = 200   #number of radial bins
zbin = 100   #number of z bins
atom = np.zeros((zbin,rbin)) 
hbond = np.zeros((zbin,rbin))
xdata = np.zeros((zbin-1))

r = [] #np.zeros((zbin-1))
ydata = np.zeros((rbin-1))
zdata = np.zeros((rbin-1))

da = 95
nsnaps = 0
success = True
outfile = "temp/out.dat"
outfile2 = "temp/out2.dat"
cost = 0
#donor = [2,7,8,15,8,16]  # D-H

while 1:
	time1 = w.next()
	time2 = a.next()
	#d.scale(time)
	if (time1 == -1 or time2 == -1) : break
	#print time1,time2
	hb = 0
	mx = w.vecs(time1,"x")
	my = w.vecs(time1,"y")
	mz = w.vecs(time1,"z")
	nx = a.vecs(time2,"x")
	ny = a.vecs(time2,"y")
	nz = a.vecs(time2,"z")
	itype = w.vecs(time1,"type")
	matom = len(mx)
	natom = len(nx)
	minz,maxz = w.minmax("z")
	zl = minz
        zl = 0.0                     #coordinate of 1st layer atom of graphite sheet
	zu = minz + zbin
	xcom = sum(mx)/matom
	ycom = sum(my)/matom
	#dist = len(x)
	#dist = ((x-xcom)**2 + (y-ycom)**2)**(1/2)
	for kk in xrange(0,len(donor),2):
		#print 'k',k,donor[k],donor[k+1]
		#raw_input()

		for i in xrange(0,matom,3):
			x1 = mx[i+donor[kk]]    #  donor atom
        		y1 = my[i+donor[kk]]      
        		z1 = mz[i+donor[kk]]    
        		x2 = mx[i+donor[kk+1]]      #H1 hydrogen attached to donor
       	 		y2 = my[i+donor[kk+1]]
        		z2 = mz[i+donor[kk+1]]
			xmol = sum(mx[i:i+3])/3
			ymol = sum(my[i:i+3])/3		

			hb = 0
			success = True	
			#print matom, i, i+donor[kk+1]
                	#raw_input()
			for j in xrange(0,natom,1):
				x4 = nx[j]  #  acceptor atom
                		y4 = ny[j]  #
                		z4 = nz[j]  #
				roij = ((x4-x1)**2 + (y4-y1)**2 + (z4-z1)**2)**0.5  #distance between donor and acceptor
				if(roij <= 3.5) :
					m1=(((x2-x4)*(x2-x4))+((y2-y4)*(y2-y4))+((z2-z4)*(z2-z4)))**0.5   # H1-O2
                			k1=(((x1-x2)*(x1-x2))+((y1-y2)*(y1-y2))+((z1-z2)*(z1-z2)))**0.5   # O1-H1
					hb0 = 0
					if(m1 <= 2.5 and hb0 == 0) :
						cost=((x2-x1)*(x4-x1))+((y2-y1)*(y4-y1))+((z2-z1)*(z4-z1))   #angle between H1-O1 and O2-O1 vector
                				cost=cost/(k1*roij)

						if(cost > 0.86) :
							hb = hb + 1
							hb0 = hb0+1
		
			#raw_input()
				#xin = input('proceed??:')
			#print 'hb',hb, roij
			
			dist = (xmol-xcom)**2 + (ymol-ycom)**2
			dist = dist**0.5
			dist1 = (rbin*da/3.14)
			dist1 = dist1**0.5
			if(dist <= dist1):
        			k = 0
        			while success:
          				k = k + 1
          				dist2= ((k*da/3.14))**0.5
          				if(dist <= dist2): success = False     
          				#print 'dist and dist2', dist, dist2           
      				success = True
      				z7 = z1 - zl
      				z8 = int(z7/5) + 1
				#z7 = int(z7) + 1
      				#if(z1 > 200): print ' z and z1', i, z[i], z1
      				#print 'z8 j', z8, j
      				if(z8 < zbin and k < rbin): 
					atom[int(z8),k] = atom[int(z8),k] + 1
					hbond[int(z8),k] = hbond[int(z8),k] + hb
					#print 'hb z k atom', hb,z8,k, atom[int(z8),k]
					#raw_input()
	nsnaps += 1
	print time1,time2
hbond = hbond/nsnaps
atom = atom/nsnaps
#print 'atom', atom
#*****************************************printing radial and axial density*****************

fp = open(outfile,"w")
for i in xrange(zbin-1):
  for j in xrange(rbin-1):
    if(atom[i+1,j+1] > 0) : hbond[i+1,j+1]=hbond[i+1,j+1]/atom[i+1,j+1]	
    print >>fp, 5*i+2.5, (j*da/3.14)**0.5, atom[i+1,j+1],hbond[i+1,j+1]*nd
  print >>fp
fp.close()


