# theLine

To briefly summarize this repository, this is a WIP project I began in the summer of 2023.

The main idea behind the project was to try and make a reliable approximation method for measuring the best (fastest depending on many different conditions) possible path along a race track for a car to take in a race. 

There are of course a multitude of methods to attack this problem and achieve the same or similar results. However, somewhere along the planning of the project I became really interested in Bezier Curves. These curves are quite unique by the fact that rather than being a regular injective function, they are a set of time based equations, defining as many dimensions as one wants (2D, 3D, 4D...etc). On top of this, these dimensional functions can be easily defined by a set of points (that control the shape of the curve).

In theory, if one uses enough points to define a Bezier Curve, one can find the optimal placing of these points that leads to a curve (lap around the track) which provides the best overall race time.

Unfortunately, I am still working on a function that can accurately represent the distribution of weight (causing friction in the tyres) over the chassis of the vehicle based on the rate of turning, initial velocity, acceleration mid turn, drag, and the inertia of the chassis. This is a statically indeterminate problem that I have been trying to tackle on and off (currently quite busy with college semester). 

However, I have made a quite functional (few hiccups to fix) Bezier Curve adjusting algorithm that will adjust the positioning of the characteristic point coordinates to maximize or minimize a chosen parameter (track time is the goal, but curvature and distance work). 

So far, there is also a lot of redundant code in the main.cpp file and a lot of clean up is required. The computed coordinates of the characteristic points of the Bezier Curve are also saved to a text file, and the graphics are currently plotted using matplotlib in python (I will switch to SDL2 soon).
