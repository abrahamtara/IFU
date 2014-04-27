#!/usr/bin/python

# plot3.py: work in progress
# this is the same as 'plot3_040914.py' the astro computers

# packages
from __future__ import division
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
from matplotlib import collections, transforms
import math
import pyfits
import numpy as np
import numpy.ma as ma
import datetime

# modules
from map import coords
from input_list import parse_input

# query: function that figures out how many files we want to plot
# and whether we will be plotting values from the continuum, 
# or values of certain emission lines from a data file
def query(line_list, datafiles):
	object_name = raw_input("Object name: > ")
	raw_num_plots = raw_input("Number of plots to display: > ")
	num_plots = int(raw_num_plots) # got to turn this into an int
	print "Type 'c' to plot the continuum, or give the name of the line you want to plot"
	print "(Your choices for the lines are: \n" + line_list
	data_type = raw_input("\n> ")
	# if we're going to be plotting the continuum, we don't need a separate data file
	if data_type == 'c':
		data_file = 'none'
	# but if we're going to be plotting a spectral line, then we need data files
	# the names of these files will be appended to the 'datafiles' list
	elif data_type in line_list:
		print "What are the names of the data files you want to plot?\n"
		print "Type them below in this order: Home, Off1, Off2"
		for i in range(num_plots):
			datafiles.append(raw_input("> "))
	else:
		print "You have not provided a viable option. Try again"
		exit
	return num_plots, data_type, datafiles, object_name

# parsedec: function that takes in a string 
# '+/DEG:ARCMIN:ARCSEC.ARCSEC' and turns it into arcseconds
# we will be using this to convert the central declination value given in the header
# into arcseconds
def parsedec(string): 
	array = string.split(":")
	deg = float(array[0])*3600
	arcminutes = float(array[1])*(60)
	arcseconds = float(array[2])
	if deg < 0:
		total_arcsecs = deg - arcminutes - arcseconds
	else:
		total_arcsecs = deg + arcminutes + arcseconds
	return total_arcsecs

# shift_coords: function to shift the default coordinates given 
# in the 'coords' dictionary according to the rotation angle
# and the offset distances for each offset position
def shift_coords(central_dec, ra_shift, dec_shift, angle):
	offset = math.radians(angle) # offset angle
	new_coords = {}
	for a in range(1, len(coords) + 1): # usually 1 to 83, excludes 83
		rotated_dec = coords[a][0]*math.sin(offset)+coords[a][1]*math.cos(offset)
		rotated_ra = coords[a][0]*math.cos(offset)-coords[a][1]*math.sin(offset)

		shifted_dec = rotated_dec + dec_shift
		shifted_ra = rotated_ra + ra_shift
		
		new_coords[a] = [shifted_ra, shifted_dec]
	return new_coords
 	
# load_files: function that uses the pyfits module load a fits file
# and return the file object and header
def load_files(num):
	if num == 1:
		position = "Offset 1"
	elif num == 2:
		position = "Offset 2"
	else:
		position = "    Home" # if num is 0
		
	position_string = "%s position fits file? > " % position
	file_name = str(raw_input(position_string))
	file_path = fitsdir + file_name
	this_fits = pyfits.open(file_path)
	this_obj = this_fits[0].data
	this_header = this_fits[0].header
	return this_obj, this_header, file_name

# load_values: depending on whether we're plotting the continuum or a spectral line
# this pulls values from the fits file, or uses the parse_input module to pull
# data from the desired spectral line. Returns an updated "colors" array 	
def load_values(plotnum, xlow, xhigh, object, header, colors, data_type, datafiles):
	if data_type == 'c':
		crval1 = header['CRVAL1']
		cd1_1 = header['CD1_1']
		xval = crval1 + cd1_1*np.arange(object[0,:].size)
	
		for ind in range(82):
			xinds = np.where((xval[:] > xlow ) & (xval[:] < xhigh))
			fluxval = np.median(np.median(object[ind,xinds])) 
			colors[plotnum][ind] = fluxval # here, ind goes from 0 to 81 inclusive

	else:
		colors[plotnum] = parse_input(datafiles[plotnum], data_type)

	return colors # now we've populated colors[i]

# get_header_data: A function to get the central declination and rotation angle
# from the fits header
def get_header_data(header):
	central_dec = parsedec(header['TARGDEC'])*(1/3600) # unsure if we should be using DEC or TARGDEC
	degrees_offset = parsedec(header['ROTOFF'])*(1/3600)
	return central_dec, degrees_offset

# average_colors: this takes the average of certain fibres, based on which fibres overlap when we fill the grid with 2 offsets
# note that fibre number j for position i is represented as colors[i][j-1] because the colors list
# goes from 0 to 81 inclusive, for all 82 fibres. this is how the data is organized in the 'fits' file. 
# To find the average, we first create a masked array that masks 'NaN' values of -100 from the mean. 
def average_colors(colors, i, num):
		
	if i == 2: # if we're on off2 
		colors[i][61] = ma.masked_values([colors[0][46], colors[1][43], colors[i][61]], -100).mean()
		colors[i][42] = ma.masked_values([colors[0][51], colors[1][46], colors[i][42]], -100).mean()
		colors[i][65] = ma.masked_values([colors[0][30], colors[1][51], colors[i][65]], -100).mean()
		colors[i][46] = ma.masked_values([colors[0][25], colors[1][40], colors[i][46]], -100).mean()
		colors[i][51] = ma.masked_values([colors[0][45], colors[1][25], colors[i][51]], -100).mean()
		colors[i][30] = ma.masked_values([colors[0][24], colors[1][45], colors[i][30]], -100).mean()
		colors[i][44] = ma.masked_values([colors[0][65], colors[1][42], colors[i][44]], -100).mean()
		colors[i][43] = ma.masked_values([colors[0][40], colors[i][43]], -100).mean()
		colors[i][50] = ma.masked_values([colors[0][47], colors[1][30], colors[i][50]], -100).mean()
		colors[i][47] = ma.masked_values([colors[1][24], colors[1][47]], -100).mean()
		colors[i][45] = ma.masked_values([colors[1][0], colors[i][45]], -100).mean()

	elif i == 1: # if we're mapping off1
		colors[i][41] = ma.masked_values([colors[0][43], colors[i][41]], -100).mean()
		colors[i][47] = ma.masked_values([colors[0][78], colors[i][47]], -100).mean()
		colors[i][61] = ma.masked_values([colors[0][42], colors[i][61]], -100).mean()
		colors[i][65] = ma.masked_values([colors[0][50], colors[i][65]], -100).mean()


		# and the rest stay the same! 

	##### FOR A376M150 THE SHIFTS ARE DIFFERENT: Uncomment below if we are mapping A376
#	if i == 2: 
#		colors[i][43] = np.average([colors[0][40], colors[1][46], colors[i][44]])
#		colors[i][61] = np.average([colors[0][46], colors[1][42], colors[i][61]])
#		colors[i][46] = np.average([colors[0][25], colors[1][51], colors[i][46]])
#		colors[i][42] = np.average([colors[0][51], colors[1][65], colors[i][42]])
#		colors[i][25] = np.average([colors[0][0], colors[1][45], colors[i][25]])
#		colors[i][51] = np.average([colors[0][45], colors[1][30], colors[i][51]])
#		colors[i][65] = np.average([colors[0][30], colors[1][50], colors[i][65]])
#		colors[i][30] = np.average([colors[0][24], colors[1][47], colors[i][30]])
#		colors[i][44] = np.average([colors[1][65], colors[i][44]])
#		colors[i][50] = np.average([colors[1][47], colors[i][50]])

#		colors[i][61] = np.average([colors[0][43], colors[i][61]])		
#		colors[i][44] = np.average([colors[0][42], colors[i][44]])	
#
	return colors

##
	if i == 2: # if we're on off2 
		colors[i][61] = ma.masked_values([colors[0][46], colors[1][43], colors[i][61]], -100).mean()
		colors[i][42] = ma.masked_values([colors[0][51], colors[1][46], colors[i][42]], -100).mean()
		colors[i][65] = ma.masked_values([colors[0][30], colors[1][51], colors[i][65]], -100).mean()
		colors[i][46] = ma.masked_values([colors[0][25], colors[1][40], colors[i][46]], -100).mean()
		colors[i][51] = ma.masked_values([colors[0][45], colors[1][25], colors[i][51]], -100).mean()
		colors[i][30] = ma.masked_values([colors[0][24], colors[1][45], colors[i][30]], -100).mean()
		colors[i][44] = ma.masked_values([colors[0][65], colors[1][42], colors[i][44]], -100).mean()
		colors[i][43] = ma.masked_values([colors[0][40], colors[i][43]], -100).mean()
		colors[i][50] = ma.masked_values([colors[0][47], colors[1][30], colors[i][50]], -100).mean()
		colors[i][47] = ma.masked_values([colors[1][24], colors[i][47]], -100).mean()
		colors[i][45] = ma.masked_values([colors[1][0], colors[i][45]], -100).mean()

	elif i == 1: # if we're mapping off1
		colors[i][41] = ma.masked_values([colors[0][43], colors[i][41]], -100).mean()
		colors[i][47] = ma.masked_values([colors[0][78], colors[i][47]], -100).mean()
		colors[i][61] = ma.masked_values([colors[0][42], colors[i][61]], -100).mean()
		colors[i][65] = ma.masked_values([colors[0][50], colors[i][65]], -100).mean()
##
	

###################################################################################### 
##                                                                                  ##
##                                                                                  ##
##                                     MAIN                                         ##
##                                                                                  ##
##                                                                                  ##
######################################################################################

# Initialization section: 
# defining the directory with the fits files, the low and high angstrom values to map
# the list of spectral lines we could plot, and the shift amounts for each offset in arcseconds. 

# should check for presence of necessary files and directories now! 
fitsdir = '/net/edwards/ta263/fits/' 
xlow = 5500.0  
xhigh = 6000.0  
# line_list = 'HdA HdF CN1 CN2 Ca4227 G4330 HgA HgF Fe4383 Ca4455 Fe4531\n C4668 Hb Fe5015 Mg1 Mg2 Mgb Fe5270 Fe5335 Fe5406\n Fe5709 Fe5782 NaD TiO1 TiO2'
line_list = 'Age  [Fe/H]  [Mg/Fe]  [C/Fe]  [N/Fe]  [Ca/Fe]  [O/Fe]  [Na/Fe]  [Si/Fe]  [Cr/Fe]  [Ti/Fe]  iso   x(imf)  Age_Hb  Age_Hg  Age_Hd'
actual_shifts = [(0,0),(-5.63,0),(-2.82, -4.93)] # incorporate this into the query! 407 shifts
# actual_shifts = [(0,0),(2.815, -4.8757), (-2.8645, -4.8985)] # incorporate this into the query! 376 shifts

datafiles = []

# Query: figure out what we're going to be plotting, and how many times we should loop the next section
num, data_type, datafiles, object_name = query(line_list, datafiles) 

# Initializing some lists
totalcoords = []
totalcolors = []
colors = [[0 for j in range(82)] for i in range(num)] # list comprehension!

# Loop
for i in range(0, num): # i goes from 0 to num-1
	object, header, file_name = load_files(i) # now we have the spectral data and the header data

	# get data from the file
	central_dec, degrees_offset = get_header_data(header)
	colors_raw = load_values(i, xlow, xhigh, object, header, colors, data_type, datafiles) # populates one of the colors lists
	colors_averaged = average_colors(colors_raw, i, num)
	actual_ra_shift = actual_shifts[i][0]
	actual_dec_shift = actual_shifts[i][1]
	shifted_coords = shift_coords(central_dec, actual_ra_shift, actual_dec_shift, degrees_offset)

	# populate 'totalcoords' list
	for j in range(1, len(coords) + 1): # usually 1 to 83, excludes 83
		totalcoords.append([(i, j), (shifted_coords[j][0], shifted_coords[j][1])]) # appending 2 tuples, one with info and one with coords

	# populate 'totalcolors' list
	for j in range(0, len(coords)): # 0 to 82, excludes 82
		this_color = colors_averaged[i][j]
		totalcolors.append(this_color)
	
	# now after this loop, we've added new coordinates to our totalcoords, and 
	# added new intensity values to our totalcolors

# split up total
fibre_info, just_coordinates = zip(*totalcoords)

# output data to a data file
data_text_name = object_name + ' data created at ' + str(datetime.datetime.now()) + '.txt'
data_text_file = open(data_text_name, 'w')
data_text_file.write('(Position, Fibre), (X-coord, Y-coord), Intensity\n')

for item in range(0, len(totalcolors)):
	#data_text_file.write(str(totalcoords[item]) + str(totalcolors[item]))
	data_text_file.write(str(fibre_info[item]).replace('[','').replace(']',''))
	if fibre_info[item] in o:
		data_text_file.write('*')
	data_text_file.write(", ")
	data_text_file.write(str(just_coordinates[item]).replace('[','').replace(']',''))
	data_text_file.write(", ")
	data_text_file.write(str(totalcolors[item]))
	data_text_file.write('\n')

data_text_file.close()

# make the figure!
fig = plt.figure()
ax = fig.add_subplot(111, xlim=(-80,80), ylim=(-80,80))

# point sizes for our circles? point size = 82, for 82 fibres
#pointsizes = [17.253]*82*num 
pointsizes = [82]*82*num 

# make a collection of circles
col = collections.CircleCollection(sizes=pointsizes, offsets=just_coordinates,
                                        transOffset=ax.transData)
trans = transforms.Affine2D().scale(fig.dpi/72.0) # one point is 1/72 inches
col.set_transform(trans)  # the points to pixels transform (from example)
ax.add_collection(col, autolim=True)
ax.autoscale_view()

# our colormap will be 'hot' in reverse, with outliers (-100) in green (but only if we're plotting line data)
cust_cm = cm.hot_r
cust_cm.set_bad('gray', -50) # -100 here is just the 'alpha' value for the colour, it isn't the outlier value

col.set_cmap(cust_cm)
col.set_norm(matplotlib.colors.Normalize())
col.set_edgecolors('black') # edge colours are an off-black hue 262626

# make a masked array from the data
if data_type != 'c':
	masked_totalcolors = ma.masked_equal(totalcolors, -100.0)
	col.set_array(masked_totalcolors)
else:
	col.set_array(np.array(totalcolors))

# make a colorbar
cbar = plt.colorbar(col)

if data_type == 'c':
	title_string = "%r continuum from %r to %r angstroms" % (object_name, xlow, xhigh)
else:
	title_string = "Plot of %r for object %r" % (data_type, object_name)

plt.title(title_string)

plt.show()

