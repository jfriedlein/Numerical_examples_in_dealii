// deal.II headers
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>

// C++ headers for some math operations
// @todo-assure Does this list contain everything that is needed and not more?
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>

// Numerical example helper function (required, can be downloaded from https://github.com/jfriedlein/Numerical_examples_in_dealii)
// also contains enumerators as part of "enums::"
#include "./numEx-helper_fnc.h"

using namespace dealii;

namespace Unconstrained_elastoplastic_test
/*
 * A single (or multiple) element(s) under shear load in x-direction, dimensions widthxdim
 * By selecting 1 global refinement and distortion we can create the distorted patch element test in 2D or 3D
 * @todo The refine special enums need to be local values for each numEx
 *
 * CERTIFIED TO STANDARD numExS11 (210104)
 */
{
	// @todo-optimize I guess the following variables will never be destroyed, so the memory remains filled?

	// Name of the numerical example
	// @todo-optimize Can we extract the name of the namespace as string? (however, that would limit the naming freedom)
	 std::string numEx_name = "Unconstrained_elastoplastic_test";
	 
	// The loading direction: \n
	// In which coordinate direction the load shall be applied, so x/y/z. We assume the positive axis direction, where
	// changes in the direction (+ -) are possible with positive and negative loads.
	 const unsigned int loading_direction = enums::x;

	// The loaded faces:
	 const enums::enum_boundary_ids id_boundary_load = enums::id_boundary_none;
	 const enums::enum_boundary_ids id_boundary_secondaryLoad = enums::id_boundary_none;

	// Characteristic body dimensions
	 std::vector<double> body_dimensions (5);
	 
	// Evaluation point
	// @note We cannot init the point yet, because we don't have the geometry dimensions and geometry
	 Point<3> eval_point;
	 
	// Evaluation path
	 Point<3> eval_path_start;
	 Point<3> eval_path_end;

	// Some internal parameters
	 struct parameterCollection
	 {
		// We declare the search tolerance as "static constexpr" so we don't need
		// to create an instance of this class just to get access to this variable.
		// If you append this struct, then you should probably create an instance and avoid the "static ..." blabla.
		static constexpr double search_tolerance = 1e-12;
	 };
	 
	// Evaluation points: \n
	// @todo We need \a dim here instead of 2, but dim is unkown at this place -> redesign
	 std::vector< numEx::EvalPointClass<3> > eval_points_list (2, numEx::EvalPointClass<3>() );


	// All additional parameters
	// @todo Group them somehow
	  
	/**
	 * Apply the boundary conditions (support and load) on the given AffineConstraints \a constraints. \n
	 * For the HyperCube that are three symmetry constraints on each plane (x=0, y=0, z=0) and the load on the \a id_boundary_load (for Dirichlet).
	 * Alternatively, we apply the rigid body motion for contact.
	 * @todo Change setup (see HyperR) where we define boundary for BC_xMinus and then use if ( BC_xMinus==... ), use the value of BC_xMinus directly
	 */
	template<int dim>
	void make_constraints ( AffineConstraints<double> &constraints, const FESystem<dim> &fe, DoFHandler<dim> &dof_handler_ref,
							const bool &apply_dirichlet_bc, double &current_load_increment, const Parameter::GeneralParameters &parameter,
							const unsigned int current_load_step )
	{
		// BC for the load ...
		 if ( parameter.driver == enums::Dirichlet )  // ... as Dirichlet only for Dirichlet as driver, alternatively  ...
		 {
			 const double load_increment = apply_dirichlet_bc ? current_load_increment : 0;
			 // First, add lines to the constraints matrix
			 // Then, we fill all desired non-zero entries, here we fill every entry even the entries that are zero
			  constraints.add_line(0);
			  constraints.set_inhomogeneity(0,0);

			  constraints.add_line(3);
			  constraints.set_inhomogeneity(3,0);

			  constraints.add_line(5);
			  constraints.set_inhomogeneity(5,0);

			  constraints.add_line(6);
			  constraints.set_inhomogeneity(6,load_increment);

			  constraints.add_line(7);
			  constraints.set_inhomogeneity(7,0);
		 }
	}

	// HyperCube grid: 2D and 3D
	template<int dim>
	void make_grid ( Triangulation<dim> &triangulation, const Parameter::GeneralParameters &parameter )
	{
		const double search_tolerance = parameterCollection::search_tolerance;

		const double width = parameter.width;
		
		// Assign the characteristic dimensions of the cube
		 body_dimensions[enums::x] = width;
		 body_dimensions[enums::y] = width;
		 body_dimensions[enums::z] = width;
		
		// Set the evaluation point
		 eval_point[enums::x] = 0;
		 eval_point[enums::y] = width*std::sqrt(2);
		 eval_point[enums::z] = 0;

		// Set the evaluation path points
		 eval_path_start = Point<3> (0,0,0);
		 eval_path_start = eval_point;

		// Create the triangulation
		 GridGenerator::hyper_cube(triangulation,0,width);

		// Clear all existing boundary ID's
		 numEx::clear_boundary_IDs( triangulation );
		 
		// Set the boundary IDs
		 for ( typename Triangulation<dim>::active_cell_iterator
			  cell = triangulation.begin_active();
			  cell != triangulation.end(); ++cell )
		 {
			for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
			  if (cell->face(face)->at_boundary())
			  {
				if (std::abs(cell->face(face)->center()[enums::x] - 0.0) < search_tolerance)
				{
					cell->face(face)->set_boundary_id(enums::id_boundary_xMinus);
				}
				else if (std::abs(cell->face(face)->center()[enums::x] - body_dimensions[enums::x]) < search_tolerance)
				{
					cell->face(face)->set_boundary_id(enums::id_boundary_xPlus);
				}
				else if (std::abs(cell->face(face)->center()[enums::y] - 0.0) < search_tolerance)
				{
					cell->face(face)->set_boundary_id(enums::id_boundary_yMinus);
				}
				else if (std::abs(cell->face(face)->center()[enums::y] - body_dimensions[enums::y]) < search_tolerance)
				{
					cell->face(face)->set_boundary_id(enums::id_boundary_yPlus);
				}
				else if (dim==3 && std::abs(cell->face(face)->center()[enums::z] - 0.0) < search_tolerance)
				{
					cell->face(face)->set_boundary_id(enums::id_boundary_zMinus);
				}
				else if (dim==3 && std::abs(cell->face(face)->center()[enums::z] - body_dimensions[enums::z]) < search_tolerance)
				{
					cell->face(face)->set_boundary_id(enums::id_boundary_zPlus);
				}
				else
				{
					// There are only 6 faces for a cube in 3D, so if we missed one, something went terribly wrong
					AssertThrow(false, ExcMessage( numEx_name+" - make_grid 3D<< Found an unidentified face at the boundary. "
												   "Maybe it slipt through the assignment or that face is simply not needed. "
												   "So either check the implementation or comment this line in the code") );
				}
			  }
		 }

		 // rotation by 45° = pi/4 in counter-clockwise direction
		  const double rotation_angle_in_radian = (std::atan(1)*4.) / 4.;
		  GridTools::rotate(rotation_angle_in_radian,enums::z,triangulation);

		// Refinement
		 triangulation.refine_global( parameter.nbr_global_refinements );

		// Evaluation points and the related list of them
		 numEx::EvalPointClass<3> eval_topLeftX ( eval_point, enums::x );
		 numEx::EvalPointClass<3> eval_topLeftY ( eval_point, enums::y );

		 eval_points_list = {eval_topLeftX,eval_topLeftY};

		// Output the triangulation as eps or inp
		 //numEx::output_triangulation( triangulation, enums::output_eps, numEx_name );
	}
}
