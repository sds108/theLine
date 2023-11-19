import numpy as np
import matplotlib.pyplot as plt
import os, glob
from array import *
import imageio

borders_path = '/home/chad/Documents/projects/theLine/borders'
paths_path = '/home/chad/Documents/projects/theLine/paths'
outside_path = '/home/chad/Documents/projects/theLine/outside'
pPixels_path = '/home/chad/Documents/projects/theLine/pixelsPath'
bPixels_path = '/home/chad/Documents/projects/theLine/pixelsBorder'

borders = []
paths = []
outside = []
bPixels = []
pPixels = []

step = 0.001

steps = int(1 / step)

def bezier(t, a, b, c, d):
	return (pow(1 - t, 3) * a) + (3 * t * pow(1 - t, 2) * b) + (3 * t * t * (1 - t) * c) + (pow(t, 3) * d)

# Extract Borders
for filename in glob.glob(os.path.join(borders_path, '*.txt')):
  with open(filename, 'r') as f:
    borders.append(f.read())
	
# Extract Paths
for filename in glob.glob(os.path.join(paths_path, '*.txt')):
  with open(filename, 'r') as f:
    paths.append(f.read())
	
# Extract Outside
for filename in glob.glob(os.path.join(outside_path, '*.txt')):
  with open(filename, 'r') as f:
    outside.append(f.read())

# Extract pixels
for filename in glob.glob(os.path.join(pPixels_path, '*.txt')):
  with open(filename, 'r') as f:
    pPixels.append(f.read())
	
for filename in glob.glob(os.path.join(bPixels_path, '*.txt')):
  with open(filename, 'r') as f:
    bPixels.append(f.read())

# Extract border segments and draw them
for border in borders:
	segments = border.split('\n')
	bezierx = []
	beziery = []
	for params in segments:
		param = params.split(' ')
		for point in range(0, steps, 1):
			# plotting the points
			ax = float(param[0])
			ay = float(param[1])
			bx = float(param[2])
			by = float(param[3])
			cx = float(param[4])
			cy = float(param[5])
			dx = float(param[6])
			dy = float(param[7])
			#controlPointsX = [ax, bx, cx, dx]
			#controlPointsY = [ay, by, cy, dy]
			#plt.plot(controlPointsX, controlPointsY, 'b:', linewidth = 0.2)
			bezierx.append(bezier(point * step, ax, bx, cx, dx))
			beziery.append(bezier(point * step, ay, by, cy, dy))
			#plt.plot(bezier(point * step, ax, bx, cx, dx), bezier(point * step, ay, by, cy, dy), 'b.')
	plt.plot(bezierx, beziery, 'b-', linewidth = 1)
	
# Extract path segments and draw them
for path in paths:
	segments = path.split('\n')
	bezierx = []
	beziery = []
	for params in segments:
		param = params.split(' ')
		for point in range(0, steps, 1):
			# plotting the points
			ax = float(param[0])
			ay = float(param[1])
			bx = float(param[2])
			by = float(param[3])
			cx = float(param[4])
			cy = float(param[5])
			dx = float(param[6])
			dy = float(param[7])
			controlPointsX = [ax, bx, cx, dx]
			controlPointsY = [ay, by, cy, dy]
			#plt.plot(controlPointsX, controlPointsY, 'k:', linewidth = 0.2)
			plt.plot(ax, ay, 'b+')
			plt.plot(dx, dy, 'b+')
			bezierx.append(bezier(point * step, ax, bx, cx, dx))
			beziery.append(bezier(point * step, ay, by, cy, dy))
			#plt.plot(bezier(point * step, ax, bx, cx, dx), bezier(point * step, ay, by, cy, dy), 'r.')
	#plt.plot(bezierx, beziery, 'r.')
	plt.plot(bezierx, beziery, 'r-', linewidth = 1)

"""
# Extract outside points and draw them
for out in outside:
	point = out.split('\n')
	pointx = []
	pointy = []
	for params in point:
		param = params.split(' ')
		#pointx.append(float(param[0]))
		#pointy.append(float(param[1]))
		plt.plot(float(param[0]), float(param[1]), 'gx', linewidth = 0.5)
	#plt.plot(pointx, pointy, 'w.')
"""
	
"""
# Extract pixels and draw them
for out in bPixels:
	point = out.split('\n')
	pointx = []
	pointy = []
	for params in point:
		param = params.split(' ')
		#pointx.append(float(param[0]))
		#pointy.append(float(param[1]))
		plt.plot(int(param[0]), int(param[1]), 'b.', linewidth = 1)
		
for out in pPixels:
	point = out.split('\n')
	pointx = []
	pointy = []
	for params in point:
		param = params.split(' ')
		#pointx.append(float(param[0]))
		#pointy.append(float(param[1]))
		plt.plot(int(param[0]), int(param[1]), 'r.', linewidth = 1)
"""
			
# naming the x axis
plt.xlabel('x - axis')
# naming the y axis
plt.ylabel('y - axis')
# function to show the plot
plt.show()

