"""
Fixing an image using Ising's Model and Simulated Annealing

by
	Victor Hugo Vianna


Instructions:
Just put the file "img_disturbed.png" in the same folder as the script
"""

from PIL import Image, ImageFilter
#import gifmaker
import random as random
import math
import sys
import numpy as np

# constants
_BLACK = (0, 0, 0) 
_WHITE = (255, 255, 255)
_H = 0 # ising's model

_BETA = 100
_GAMMA=10


def theta(n):
	return math.exp(100)/n;

def spin(p):
	if(p==_BLACK):
		return 1;
	return -1;

def magnet(pix, sizeX, sizeY):
	resp=0
	for i in range (1, sizeX-1): ## only in the interior
		for j in range (1, sizeY-1):
			resp+=spin(pix[i,j])
	resp = 1.0*resp/ ((sizeX-2)*(sizeY-2))
	return resp;

def isValidPoint(x, y, sizeX, sizeY):
	return (x>=0 and x<sizeX and y>=0 and y<sizeY)
def isInteriorPoint(x, y, sizeX, sizeY):
	return isValidPoint(x,y,sizeX,sizeY) and x!=0 and x!=(sizeX-1) and y!=0 and y!=(sizeY-1) 


def metropolis(pixA, pixB, sizeX, sizeY, nIter):

	for k in range (1, nIter+1):

		# if (k%1000==0):  #if you want to see the evolution of the algorithm
		# 	img.save("img_corrected("+str(k/1000)+").png")

		x , y = random.randint(0,sizeX-1) , random.randint(0,sizeY-1)
		delta_ham =0; aux = 0;
	
		##calculating variation of the hamiltonian
		delta_ham = (_GAMMA**2)*0.5*(4*spin(pixA[x,y])*spin(pixB[x,y])) # 4*a*b == (-a-b)**2 -(a-b)**2

		p = x-1; q = y;
		coef= 2 if (isInteriorPoint(p,q, sizeX, sizeY) or isInteriorPoint(x,y, sizeX, sizeY)) else 1
		if isValidPoint(p,q,sizeX,sizeY):
			aux+=coef*spin(pixA[p,q])
		
		p = x+1; q = y;
		coef= 2 if (isInteriorPoint(p,q, sizeX, sizeY) or isInteriorPoint(x,y, sizeX, sizeY)) else 1
		if isValidPoint(p,q,sizeX,sizeY):
			aux+=coef*spin(pixA[p,q])

		p = x; q = y-1;
		coef= 2 if (isInteriorPoint(p,q, sizeX, sizeY) or isInteriorPoint(x,y, sizeX, sizeY)) else 1
		if isValidPoint(p,q,sizeX,sizeY):
			aux+=coef*spin(pixA[p,q])
		
		p = x; q = y+1;
		coef= 2 if (isInteriorPoint(p,q, sizeX, sizeY) or isInteriorPoint(x,y, sizeX, sizeY)) else 1
		if isValidPoint(p,q,sizeX,sizeY):
			aux+=coef*spin(pixA[p,q])

		delta_ham+= _BETA*spin(pixA[x,y])*aux


		if (np.log(random.uniform(0, 1)/theta(k)) <= -delta_ham ): 
			pixA[x, y] = _WHITE if (pixA[x, y]==_BLACK) else _BLACK # invert spin


'''
#code to disturb original image

img = Image.open("img.png")
sizeX = img.size[0]
sizeY = img.size[1]
pix = img.load()
for i in range(0, sizeX*sizeY/100): # change 10 percent of bits
	x , y = random.randint(0,sizeX-1) , random.randint(0,sizeY-1)
	if pix[x, y] == _BLACK:
		pix[x, y] = _WHITE
	else:
		pix[x, y] =_BLACK

img.save("img_disturbed.png")

'''


# correcting disturbed image

img = Image.open("img_disturbed.png") # initial candidate in the markov chain
disturbed = Image.open("img_disturbed.png")
pix = img.load() # pixels that will be gradually corrected
pix_dist = disturbed.load() #pixels of disturbed image
sizeX = img.size[0]
sizeY = img.size[1]

if len(sys.argv)==1:
	metropolis(pix, pix_dist, sizeX, sizeY, 1000000)
	img.save("B"+str(_BETA)+"G"+str(_GAMMA)+"N1mi.png")
else: # you can pass the number of iterations you want through the terminal
	metropolis(pix, pix_dist, sizeX, sizeY, int(sys.argv[1]))
	img.save("B"+str(_BETA)+"G"+str(_GAMMA)+"N"+sys.argv[1]+".png")

