# SimplexDIC
A very rudimentary global DIC implementation for Matlab using T3 or T6 triangles or T4 tetrahedra

## Overview
This code contains three global DIC implementations using Simplices (i.e. Trangles and Tetrahedra). The goal of this code is mainly teaching purposes.
As such, it is entirely in Matlab, and thus is very (very) slow. For the 2D case, this is less of an issue, but for the 3D (volume) case you need a lot of RAM to be able to handle reasonable size images (volumes). The 2D linear triangle (T3) DIC is nice, but the defining features of this library are that it contains a 2D quadratic triangle (T6) DIC version and a 3D tetrahedron (T4) DVC version. The T6 version, is limited to elements that have T3 geometries with their quadratic nodes exactly in the middle of the T3 edges. So it is not a "true" T6, but it does allow quadratic displacement fields and thus linear strain fields.

The T3/T6 version is the one used in this paper:
B Blaysat, J Neggers, M Grediac, F Sur, Towards criteria characterizing the metrological performance of full-field measurement techniques: application to the comparison between local and global versions of DIC, Experimental Mechanics 60, 393-407

## Help
The library contains example code for each DIC version (T3, T6, T4) in the main folder aptly named:
- run_me_T3.m
- run_me_T6.m
- run_me_T4.m
Have a look at these to see how to use the code.

### Library
The three DIC functions come with a few helper functions, for instance to compute the area of a triangle, or to flip triangles upside-up, etc. Have a look at the lib folder inside to find these library functions which each have some documentation in the header of the functions.

