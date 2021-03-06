# Numerical_examples_in_dealii
A collection of different numerical examples implemented in deal.ii

"A numerical example is worth a thousand equations." [Frans Paul VAN DER MEER, Dissertation-Propositions]

... and takes about as many lines of code.


## Changelog
Proposal for updated standard:
- group additional parameter
- create standalone make_grid, with all parameters are explicit input argument

Still being implemented consistently:
From the 04.01.2021 we obey a new standard (numExS11).

- contain member variable named body_dimensions (see Bar_model), which contains characteristic dimensions (length, width, thickness) and maybe also some for paths (e.g. the location of the evaluation path);
- incorporating contact and many more numerical examples (updated soon)

Major changes:
* clean-up
* contact
* member variable named body_dimensions, which contains characteristic dimensions (length, width, thickness) 


From the 24.07.2020 we obey a new standard (numExS07).

Major changes:
* clean-up of boundary ids
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
*** @todo *** document the interfaces

* make_grid(*):

```
    QuarterHyperCube_Merged::make_grid<dim> ( triangulation, parameter );
    loading_direction = QuarterHyperCube_Merged::loading_direction;
    boundary_id_collection[enums::id_primary_load] = QuarterHyperCube_Merged::id_boundary_load;
    boundary_id_collection[enums::id_secondary_load] = QuarterHyperCube_Merged::id_boundary_secondaryLoad;
```

Creates the triangulation (2D or 3D), requires additional parameters (can also be hardcoded into the function) and boundary ids (`std::vector<enums::enum_boundary_ids> boundary_id_collection`).

* make_constraints(*):

```
			QuarterHyperCube_Merged::make_constraints<dim> ( constraints, fe, n_components, dof_handler_ref, apply_dirichlet_bc,
															 load_increment, parameter);
```

Apply the boundary conditions onto the faces, e.g. symmetry BC.

## Available numerical examples
The names try to be as general as possible. Look closely, there are many options, so e.g. the HyperRectangle is not just a rectangle, but can be notched multiple times by defined round notches.

### HyperCube: one-element (also distorted) & Distorted 8 element patch test
By setting the number of global refinements to zero, you obtain the one-element test.
For number of global refinements equal to 1 (in 3D) you obtain 8 Elements arranged as a cube, which can also be distorted internally (overall still a cube, but distorted but matching elements inside). With the option (twitch) the Dis8El version can also be twitched, so it is no longer a cube

@todo add picture of OET and distored OET

<img src="https://github.com/jfriedlein/Numerical_examples_in_dealii/blob/master/images/Dis8El.png" width="250">

### HyperRectangle: 
For the bar or sheet strip (just a matter of the width) initial inhomogenous refinements are available. For finite plasticity, we recommend making the finer mesh part approximately a square to also capture the shear bands (develop for isotropy at 55°).




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


### Tensile specimen SEP1230 (parameterised)
It's parameterised so you can also change the geometry.

<img src="https://github.com/jfriedlein/Numerical_examples_in_dealii/blob/master/images/tensileSpecimen_SEP1230.png" width="500">

### SphereRigid-cube_contact: Pushing a rigid sphere into a cube in 2D
For this example to work you require the assemble routines for contact that are NOT yet available online.

### Three-point beam
Bending of a notched beam

!!! Not developed anymore, version frozen in March 2020 !!!



