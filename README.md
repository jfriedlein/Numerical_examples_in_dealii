# Numerical_examples_in_dealii
A collection of different numerical examples implemented in deal.ii

"A numerical example is worth a thousand equations." [Frans Paul VAN DER MEER, Dissertation-Propositions]

... and takes about as many lines of code.


## Changelog
From the 24.07.2020 we obey a new standard (numExS07).

Major changes:
* cleanup of boundary ids
* option to specify the loading direction and the loaded faces


## ToDo
* add the enumerator_list
* add an exemplary call to make_constraints, make_grid
* change the interfaces of the make_grid functions to obtain the same setup as the deal.ii geometries (all dimensions and parameters as explicit input arguments)
* add a note on the difference between the qplate and the infinite plate (maybe add analytical solution code from PA)
* add a few pictures of the meshes here
* add a documentation (input arguments, interface)
* update the QPlate grid with the new deal.ii function plate_with_a_hole
* add optional parameters/switches (e.g. apply z-sym on pos z-face) as global variables into a standardised framework

## Interface
* make_grid(*):

Creates the triangulation (2D or 3D), requires additional parameters (can also be hardcoded into the function) and boundary ids

* make_constraints(*):

Apply the boundary conditions onto the faces, e.g. symmetry BC.


### Rod
Parameters

Be aware that you can also choose the geometry of the notch (round, sharp), which might affect your results.

* Total length of the rod
* Radius of the cylindrical part
* length of the notch
* Radius of the leftover material in the notch (or specify the tool radius)
* Number of elements in the notched area in y-direction (number of repetitions)
* Number of elements in the cylindrical part in y-direction

<img src="https://github.com/jfriedlein/Numerical_examples_in_dealii/blob/master/images/Rod%20-%20geometry%20notch60.jpg" width="500">

### Quarter/Eight of a plate with hole

### Bar

### One element test (also distorted) & Distorted 8 elements
By setting the number of global refinements to zero, you obtain the one-element test.
For number of global refinements equal to 1 (in 3D) you obtain 8 Elements arranged as a cube, which can also be distorted internally (overall still a cube, but distorted but matching elements inside).

<img src="https://github.com/jfriedlein/Numerical_examples_in_dealii/blob/master/images/Dis8El.png" width="250">


### Tensile specimen SEP1230 (parameterised)
It's parameterised so you can also change the geometry.

<img src="https://github.com/jfriedlein/Numerical_examples_in_dealii/blob/master/images/tensileSpecimen_SEP1230.png" width="500">

### Three-point beam
Bending of a notched beam

!!! Not developed anymore, version frozen in March 2020 !!!

