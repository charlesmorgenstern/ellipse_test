# ellipse_test
Test gp only integration for the case of ellipses of the form r^2=x^2+2y^2

Functions in ellipse_test.jl
-------------------------------------------------------------------------------------------------------

plotsetup() plots the isosurfaces, gradient paths, and intersection points.

getintersections() calculates the intersections of the gradient paths and isosurfaces and stores results in an array. Columns 2i-1 and 2i contain the x and y coordinates for the ith ellipse counting from smallest to largest.

getparameters() calculates the parameter value for the intersections assuming the parameterization x(t)=acos(t) y(t)=bsin(t).

getarclengths() calculates the arc lengths of the segments of the ellipses between intersections using a path integral and parameterization of the ellipse. The ith column has the arc lengths of the ith ellipse counting from smallest to largest.

getgplengths() calculates the arc length of gradient paths between the intersections using a path integral and parameterization of the parabola. The ith row has the arc lengths of the ith gradient path counting from 0 to pi/2.

getcurvatures() calculates the curvature of the isosurfaces at intersection points. The ith column contains the curvatures for the ith ellipse.

getcurvderivative() calculates the first derivative of the curvature of isosurfaces with respect to the arc length of the isosurface as intersection points using a finite difference method. The ith column contains the dervatives of the ith ellipse.

getarclengthsgponly() calculates the arc lengths of the segments of the ellipses between intersections using the recurrence relation. The ith column has the arc lengths of the ith ellipse counting from smallest to largest.

getareatrap() uses the arc lengths from getarclengths() and trapezoidal rule to calculate the area of the largest ellipse.

checkareatrap() checks the relative error of getareatrap() with an exact area.

getareagponlytrap() uses the arc lengths from the recurrence relation and trapezoidal rule to calculate the area of the largest ellipse.

checkareagponlytrap() checks the relative area of checkareagponlytrap() with an exact area.

checkarclengths() checks that the arc lengths from getarclengths() add up to the total arc lengths of each ellipse and returns relative errors.

checkarclengthsgponly() hecks that the arc lengths from the recurrence relation add up to the total arc lengths of each ellipse and returns relative errors.

checkarclengthserr() computes the relative error of arc lengths computed with the recurrence relation using the arc lengths from getarclengths() as a reference solution since it is verified by the other functions. 

Functions in ellipse_test_varymesh.jl
-------------------------------------------------------------------------------------------------------
Throughout nrho is the number of ellipses, nc is the number of parabolas.

The following functions are the same as above but are modified to allow for arbirtrary grids:

plotsetup(nrho,nc), getintersections(nrho,nc), getparameters(nrho,nc), getarclengths(nrho,nc), getgplengths(nrho,nc), getcurvatures(nrho,nc), getcurvderivative(nrho,nc), getarclengthsgponly(nrho,nc), getareatrap(nrho,nc), checkareatrap(nrho,nc), getareagponlytrap(nrho,nc), checkareagponlytrap(nrho,nc), checkarclengths(nrho,nc), checkarclengthsgponly(nrho,nc), checkarclengthserr(nrho,nc)

The following functions are unique to ellipse_test_varymesh.jl:

getrhos(nrho) get rho parameter for nrho equally spaced ellipses

getc(nc) get c parameter for nc equally spaced parabolas by arc length on the largest ellipse.

maxarclengtherr(nrho,nc) returns the maximum relative error of arc length for the grid.

errortablearea() generates a table with relative error and estimated order of convergence as the grid is refined using the reference arc lengths.

errortableareagponly() generates a table with relative error and estimated order of convergence as the grid is refined using the recurrence relation.

errortablearclength() generates a table with maximum relative error of arc length and estimated order of convergence as the grid is refined.

animategrids() creates an animation of the grid refinement used for the tables.

