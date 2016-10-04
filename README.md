# Calculate Integral of a given function

At input file, on the first line we have the limits in which we want
to calculate the integral value - a, and b,
and on the second line we have n, the number of iterations.
With bigger n comes a better precision of the integral value.
The function that we want to compute the integral is
represented by the f function in the .cpp file.

### Example:

##### integral.in


```
-1 1
100000
```
##### integral.out


```
--------- n = 100000 --------
```

###### Rectangle Method
* 0.0238566370

###### Trapeze Method
* 0.0238566370

###### Simpson's Method
* 0.0238566427

###### Boole's Method
* 0.0238566373

###### Newton-Cotes with 2 nodes method
* 0.0238566519

###### Newton-Cotes with 3 nodes method
* 0.0238566358

###### Gauss with 2 nodes method
* 0.0238566472

###### Gauss with 3 nodes method
* 0.0238566466
