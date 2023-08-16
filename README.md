# ellipse_test
Test gp only integration for the case of ellipses of the form r^2=x^2+2y^2

plotsetup() plots the isosurfaces, gradient paths, and intersection points.

getintersections() calculates the intersections of the gradient paths and isosurfaces and stores results in an array. Columns 2i-1 and 2i contain the x and y coordinates for the ith ellipse counting from smallest to largest.

getparameters() calculates the parameter value for the intersections assuming the parameterization x(t)=acos(t) y(t)=bsin(t).

getarclengths() calculates the arc lengths of the segments of the ellipses between intersections using a path integral and parameterization of the ellipse. The ith column has the arc lengths of the ith ellipse counting from smallest to largest.

getgplengths() calculates the arc length of gradient paths between the intersections using a path integral and parameterization of the parabola. The ith row has the arc lengths of the ith gradient path counting from 0 to pi/2.

getcurvatures() calculates the curvature of the isosurfaces at intersection points. The ith column contains the curvatures for the ith ellipse.

getcurvderivative() calculates the first derivative of the curvature of isosurfaces with respect to the arc length of the isosurface as intersection points using a finite difference method. The ith column contains the dervatives of the ith ellipse.

getarclengthsgponly() calculates the arc lengths of the segments of the ellipses between intersections using the reccurence relation. The ith column has the arc lengths of the ith ellipse counting from smallest to largest.

getareatrap() uses the arc lengths from getarclengths() and trapezoidal rule to calculate the area of the largest ellipse.

checkareatrap() checks the relative error of getareatrap() with an exact area.

getareagponlytrap() uses the arc lengths from the recurrence relation and trapezoidal rule to calculate the area of the largest ellipse.

checkareagponlytrap() checks the relative area of checkareagponlytrap() with an exact area.

checkarclengths() checks that the arc lengths from getarclengths() add up to the total arc lengths of each ellipse and returns relative errors.

checkarclengthsgponly() hecks that the arc lengths from the reccurence relation add up to the total arc lengths of each ellipse and returns relative errors.

checkarclengthserr() computes the relative error of arc lengths computed with the recurrence relation using the arc lengths from getarclengths() as a reference solution since it is verified by the other functions. 

