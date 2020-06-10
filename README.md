# Numerical_examples_in_dealii
A collection of different numerical examples implemented in deal.ii

"A numerical example is worth a thousand equations." [Frans Paul VAN DER MEER, Dissertation-Propositions]

... and takes about as many lines of code.

## ToDo
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

* Total length of the rod
* Radius of the cylindrical part
* length of the notch
* Radius of the leftover material in the notch (or specify the tool radius)
* Number of elements in the notched area in y-direction (number of repetitions)
* Number of elements in the cylindrical part in y-direction

<img src="https://github.com/jfriedlein/Numerical_examples_in_dealii/blob/master/images/Rod%20-%20geometry%20notch60.jpg" width="500">
