# Numerical_examples_in_dealii
A collection of different numerical examples implemented in deal.ii

## ToDo
* add a few pictures of the meshes here
* add a documentation (input arguments, interface)
* update the QPlate grid with the new deal.ii function plate_with_a_hole
* add optional parameters/switches (e.g. apply z-sym on pos z-face) as global variables into a standardised framework

## Interface
make_grid(*):

Creates the triangulation, requires additional parameters (can also be hardcoded into the function) and boundary ids

### Rod

<img src="https://github.com/jfriedlein/Numerical_examples_in_dealii/blob/master/images/Rod%20-%20geometry%20notch60.jpg" width="500">
