import numpy as np
import matplotlib.pyplot as plt
import os

def do1(photodir,photostring):
	infile='%s/%s' % (photodir, photostring)

    # best way to extract header info?
	f=open(infile,'r')
	f.readline()
	f.readline()
	header_info = f.readline()
	f.close()
	nz, nm, nx, step, time, ra, pr, a, dt, count = header_info.split()
	nz=int(nz); nm=int(nm); nx=int(nx); step = int(step); time=float(time); ra=float(ra); pr=float(pr)
	a=float(a); dt=float(dt); count=int(count)

	k, n, temp, omega, stream = np.genfromtxt(infile,skip_header=5,unpack=True,usecols=(0,1,2,3,4))
       	    	
	x = a * np.arange(nx)/(nx-1)
	z = 1.0 * np.arange(nz)/(nz-1)
	
	# manually specify contours
	# levels = np.linspace(-2,2,num=30)

	print 'butts'
    	
   	# construct 2d solution array, currently temperature
	sol = np.zeros([nz,nx])
	for j in range(1,nz): ###################
		row = temp[np.where(k==j)[0]] # this is a block corresponding to one k value
		# row = omega[np.where(k==j)[0]] # this is a block corresponding to one k value
		for m in range(1,nm+1): # adds each fourier term
			sol[nz-j-1,:] = sol[nz-j-1,:] + row[m] * np.cos(m * 2*np.pi * x / a) # imshow wants k index flipped for some reason
					
	# make the actual plot
	plt.clf()
	# cp=plt.contourf(x,z,sol,100)
	# print sol
	
	cp=plt.imshow(sol,aspect=a,extent=(0,2*a,0,1),vmin=-1,vmax=1) # cmap=plt.get_cmap('hot')
	# cp=plt.contour(x,z,sol)
	
	plt.xlabel('x')
	plt.ylabel('z')
	cb=plt.colorbar(cp, shrink=0.9, extend='both')
	plt.title('step %i, time %e, photo no. %i' % (step, time, count))
	plt.draw()

def do_series(photodir='/Users/mank/Dropbox/mhd/photos'):
	full_list = os.listdir(photodir)
	photos = [x for x in full_list if not (('.' in x) or ('png' in x))]
	photos.sort()
	try:
		os.mkdir('%s/png' % photodir)
	except OSError:
		'directory %s/png already works' % photodir
	for photo in photos:
		if int(photo) % 2 == 0:
			continue
		do1(photodir,photo)
		plt.savefig('%s/png/%s.png' % (photodir,photo))
		print '%s: did %s out of %i' % (photodir, photo,len(photos))

def do_batch():
	for whichcase in ['tau_0.001','tau_0.0001','tau_0.000001']:
		do_series('/Users/mank/mhd/photos_%s' % whichcase)
		
