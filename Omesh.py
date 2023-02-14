################################################################################################################
#
#	Author: Zan Xu
#	This code generates a 2D Euler-type O-mesh for the NACA0012 airfoil.
# 	The grid is generated using Karman-Trefftz conformal transformation.
# 	
#	This code reproduces the result from the paper.
#   Vassberg, J. C., & Jameson, A. (2010). 
#	In pursuit of grid convergence for two-dimensional Euler solutions. Journal of Aircraft, 47(4), 1152-1166.
#
#	Equation reference is shown in the comments.
#
#################################################################################################################
#
#	Example Use:
#		python3 Omesh.py 
#			(This uses the default setting of generating 32 x 32 Omesh. The final result is shown as a plot)
#		python3 Omesh.py -NC 64 -p 0
#			(This generates a 64 x 64 Omesh, without showing the final plot)
#		python3 Omesh.py -NC 64 -t 1e-10 -s 1
#			(This changes the tolerance for the iterative procedure to 1e-10 and save the final result into
#				a plot3d mesh file (require the plot3d python module))
#
#################################################################################################################

import numpy as np
import matplotlib.pyplot as plt
import argparse

# Define some global constants for the NACA0012 airfoil
tau = 0.2818725
P = np.pi/(2*np.pi - tau)

z_1 = 1.0089304115
z_2 = 0.0079337

zeta_1 = 0.77043505
zeta_2 = 0.24642903
zeta_c = (zeta_1 + zeta_2) / 2

class CommandLineParser:
	def __init__(self):
		# Read in arguments from command line
		parser = argparse.ArgumentParser(description = "This script takes in few options for the script to generate 2D "
			" Euler-type O-meshes for the NACA0012 airfoil")
		parser.add_argument("-NC", "--NumCell", help = "Number of cell in each direction", required = False, default = "32")
		parser.add_argument("-t",  "--Tol",     help = "Tolerance for the iterative procedure", required = False, default = "1e-13")
		parser.add_argument("-p",  "--Plot",    help = "Whether to plot the final airfoil", required = False, default = 1)
		parser.add_argument("-s",  "--Save",    help = "Whether to save the grid file (require plot3d module)", required = False, default = 0)

		args = parser.parse_args()

		print("---------------------- Input arguments ----------------------\n")
		print(">>> Number of cells  :		{} (default: 32)".format(args.NumCell))
		print(">>> Tolerance        :		{} (default: 1e-13)".format(args.Tol))
		print(">>> Plotting (0 or 1):		{} (default: 1 (True))".format(int(args.Plot)))
		print(">>> Saving   (0 or 1):		{} (default: 0 (False))".format(int(args.Save)))
		print("\n-------------------------------------------------------------")

		#store the arguments
		self.args = args

		return

def ForwardTransformation(z):

	# Eq.(9), z -> zeta
	
	rhs = np.power((z - z_1) / (z - z_2), P)

	zeta = (zeta_1 - zeta_2*rhs) / (1 - rhs)
	
	return zeta

def BackwardTransformation(zeta):

	# Eq.(9), zeta -> z
	
	lhs = np.power((zeta - zeta_1) / (zeta - zeta_2), 1.0/P)

	z = (z_1 - z_2*lhs) / (1 - lhs)

	return z

def ComplexLinearInterpolate(zeta):

	# Linear interpolation for iterative procedure

	angle = np.angle(zeta - zeta_c)
	# The following two if statements are added to make sure 
	# __angle__ is decreasing for interpolation function later
	# since angle of real number can be pi or -pi
	if (np.absolute(angle[0]) - np.pi == 0):
		angle[0] = np.pi
	if (np.absolute(angle[-1]) - np.pi == 0):
		angle[-1] = -np.pi

	mag = np.absolute(zeta - zeta_c)

	angle_linear = np.linspace(-np.pi, np.pi, np.size(zeta))

	# __angle__ is inverted as np.interp require x-coordinate to be increasing
	mag_linear = np.interp(angle_linear, angle[::-1], mag)[::-1]

	zeta_interp = mag_linear * np.exp(1j*angle_linear) + zeta_c

	return zeta_interp

def GetInitialAirfoil(NC):

	# Generate a NACA0012 airfoil coordinate pair based on NC

	x_LE = 0.0
	x_TE = 1.0089304115
	x_u = np.linspace(x_LE, x_TE, NC//2+1)
	x_l = x_u[::-1][1::]
	y_u =  0.12/0.2*(0.2969*np.sqrt(x_u) - 0.1260*x_u - 0.3516*x_u*x_u + 0.2843*x_u*x_u*x_u - 0.1015*x_u*x_u*x_u*x_u)
	y_l = -0.12/0.2*(0.2969*np.sqrt(x_l) - 0.1260*x_l - 0.3516*x_l*x_l + 0.2843*x_l*x_l*x_l - 0.1015*x_l*x_l*x_l*x_l)

	# Concatenate the upper and lower surface into a loop such that
	# the starting point is the same as the end point
	x = np.concatenate((x_u, x_l))
	y = np.concatenate((y_u, y_l))

	return x, y

def GetAirfoilWithZ(z, NC):

	# Generate a NACA0012 airfoil coordinate pair based on input coordinate (projection)

	x_u = z.real[:NC//2+1]
	x_l = z.real[NC//2+1:]

	y_u =  0.12/0.2*(0.2969*np.sqrt(x_u, out=np.zeros_like(x_u), where=np.abs(x_u)>1e-14)  
		- 0.1260*x_u - 0.3516*x_u*x_u + 0.2843*x_u*x_u*x_u - 0.1015*x_u*x_u*x_u*x_u)
	y_l = -0.12/0.2*(0.2969*np.sqrt(x_l, out=np.zeros_like(x_l), where=np.abs(x_l)>1e-14)
		- 0.1260*x_l - 0.3516*x_l*x_l + 0.2843*x_l*x_l*x_l - 0.1015*x_l*x_l*x_l*x_l)

	# Concatenate the upper and lower surface into a loop
	x = np.concatenate((x_u, x_l))
	y = np.concatenate((y_u, y_l))

	return x, y


if __name__ == '__main__':
	
	CLParser = CommandLineParser()
	
	# Number of cells. Rename for brevity.
	NC = int(CLParser.args.NumCell)

	# Get an intial airfoil based on analytical definition of NACA airfoils
	x, y = GetInitialAirfoil(NC) # Eq.(1)

	# Iteratively find the airfoil coordinates with constant theta spacing
	tol = 1
	while (tol > float(CLParser.args.Tol)):

		zeta = ForwardTransformation(x+y*1j)  # Eq.(9)

		zeta_interp = ComplexLinearInterpolate(zeta)

		z = BackwardTransformation(zeta_interp) # Eq.(9)

		x_new, y_new = GetAirfoilWithZ(z, NC) # Eq.(1)

		tol = np.linalg.norm((x_new + y_new * 1j) - (x + y*1j), 2)
		
		x = x_new
		y = y_new

		# print("tolerance = ", tol)
		# Plot the old (z.real, z.imag) pair of coordinates vs the new (x_new, y_new) pair if necessary
		# plt.scatter(z.real, z.imag)
		# plt.scatter(x_new, y_new)
		# plt.show()

	# The resulting z, zeta after the iteraitve prcedure
	z = x + 1j * y
	zeta = zeta_interp

	# Find R_1 as the discrete arc length of the quasi circle
	R_1 = np.sum(np.absolute(zeta[1:] - zeta[:-1])) / (2*np.pi)

	idx_j = np.arange(NC+1)
	R_j = R_1*np.exp(idx_j * 2 * np.pi / NC) # Eq.(13)
	R_NC = R_j[-1]

	r = np.zeros((NC+1, NC+1))
	r[:,0] = np.absolute(zeta - zeta_c)
	for j in range(NC+1):
		r[:,j] = (r[:,0]*(R_NC - R_1) + R_NC*(R_j[j] - R_1))/(R_NC - R_1) # Eq.(14)

	# In the mapped plane, transform r, theta into zeta = xi + i*eta
	theta = np.linspace(np.pi, -np.pi, NC+1)
	zeta_ij = np.multiply(r, np.exp(1j*theta)[:,None]) + zeta_c

	# Transform back into the physical plane z = x + i*y
	z = BackwardTransformation(zeta_ij)

	# Enforce symmetry
	x = np.zeros((NC+1, NC+1))
	y = np.zeros((NC+1, NC+1))
	for i in range(NC//2+1):  # Eq.(15)
		ic = NC-i
		x[i,:] = 0.5*(z.real[i,:] + z.real[ic,:])
		y[i,:] = 0.5*(z.imag[i,:] - z.imag[ic,:])

	for ic in range(NC//2+1, NC+1):  # Eq.(16)
		i = NC - ic
		x[ic,:] = x[i,:]
		y[ic,:] = -y[i,:]

	# Post-processing 1: Plotting
	if int(CLParser.args.Plot) != 0:
		# Plot the final result in terms of (x, y) coordinates
		plt.plot(x, y, 'k-')
		plt.plot(x.T, y.T, 'k-')
		plt.xlim([-1.2, 2.2])
		plt.ylim([-1.5, 1.5])
		plt.show()

	# Post-processing 2: Saving the mesh
	if int(CLParser.args.Save) != 0:
		# Write to Plot3D meshfile. Require the module plot3d
		# https://nasa.github.io/Plot3D_utilities/build/html/
		from plot3d import write_plot3D, read_plot3D
		from plot3d.block import Block

		# Convert into 3D format with redundent 3rd axis.
		x_plot3d = x[:,:,np.newaxis]
		y_plot3d = y[:,:,np.newaxis]
		z_plot3d = np.zeros(np.shape(x_plot3d))

		blocks = [Block(x_plot3d, y_plot3d, z_plot3d)]

		write_plot3D('mesh_{:d}x{:d}.xyz'.format(NC, NC),blocks, binary=False)
