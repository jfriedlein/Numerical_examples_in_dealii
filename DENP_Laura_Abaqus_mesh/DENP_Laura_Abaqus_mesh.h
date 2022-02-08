// deal.II headers
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>

// C++ headers for some math operations
#include <iostream>
#include <fstream>
#include <cmath>

// Numerical example helper function (required, can be downloaded from https://github.com/jfriedlein/Numerical_examples_in_dealii)
// also contains enumerators as part of "enums::"
#include "../numEx-helper_fnc.h"

using namespace dealii;

namespace DENP_Laura_Abaqus_mesh
/*
 * xxx
 *
 * CERTIFIED TO STANDARD xxx
 */
{
	// Name of the numerical example
	 std::string numEx_name = "DENP_Laura_Abaqus_mesh";
 
	// The loading direction: \n
	// In which coordinate direction the load shall be applied, so x/y/z.
	 const unsigned int loading_direction = enums::y;

	// The loaded faces:
	 const enums::enum_boundary_ids id_boundary_load = enums::id_boundary_yPlus;
	 const enums::enum_boundary_ids id_boundary_secondaryLoad = enums::id_boundary_xPlus;

	// Characteristic body dimensions
	 std::vector<double> body_dimensions (5);

	// Some internal parameters
	 struct parameterCollection
	 {
		const double search_tolerance = 1e-12;

		const types::manifold_id manifold_id_right_radius = 11;
		const types::manifold_id manifold_id_left_radius = 10;
	 };

	  const enums::enum_BC BC_xMinus = enums::BC_none;
	  const enums::enum_BC BC_yPlus  = enums::BC_x0_z0;
	  const bool constrain_sideways_sliding_of_loaded_face = false;
	  const enums::enum_BC BC_yMinus = enums::BC_fix;
	  const enums::enum_BC BC_zMinus = enums::BC_sym;
	  const enums::enum_notch_type notch_type = enums::notch_round;
	  const bool notch_twice = true;
	  const bool DENP_Laura = true;
	  const bool DENP_Hagen = false;
	  const bool SheStrip = false;

	// Evaluation points: \n
	// We want points, one for the contraction of the center
	// and one for the contraction of the top face.
	// We don't know the coordinates yet, because the mesh has not yet been created.
	// So we fill the data in make_grid.
	// @todo We need \a dim here instead of 3, but dim is unkown at this place -> redesign
	 std::vector< numEx::EvalPointClass<3> > eval_points_list (2, numEx::EvalPointClass<3>() );
	 
	/**
	 * Apply the boundary conditions (support and load) on the given AffineConstraints \a constraints. \n
	 * For the HyperRectangle that are three symmetry constraints on each plane (x=0, y=0, z=0) and the load on the \a id_boundary_load (for Dirichlet).
	 */
	template<int dim>
	void make_constraints ( AffineConstraints<double> &constraints, const FESystem<dim> &fe, DoFHandler<dim> &dof_handler_ref,
							const bool &apply_dirichlet_bc, double &current_load_increment,
							const Parameter::GeneralParameters &parameter )
	{
		parameterCollection parameters_internal;

		// BC on x0 plane
		 if ( BC_xMinus==enums::BC_sym )	
			numEx::BC_apply( enums::id_boundary_xMinus, enums::x, 0, apply_dirichlet_bc, dof_handler_ref, fe, constraints );
		 else if ( BC_xMinus==enums::BC_fix )
			numEx::BC_apply_fix( enums::id_boundary_xMinus, dof_handler_ref, fe, constraints );
		 
		// BC on y0 plane
		 if ( BC_yMinus==enums::BC_fix )
		 {
			// For compression we fix/clamp the Y0 plane, so it does not run away
			 numEx::BC_apply_fix( enums::id_boundary_yMinus, dof_handler_ref, fe, constraints );
		 }
		 else if ( BC_yMinus==enums::BC_sym)
			numEx::BC_apply( enums::id_boundary_yMinus, enums::y, 0, apply_dirichlet_bc, dof_handler_ref, fe, constraints );
		 
		 if ( BC_yPlus==enums::BC_x0_z0 )
		 {
			 numEx::BC_apply( enums::id_boundary_yPlus, enums::x, 0, apply_dirichlet_bc, dof_handler_ref, fe, constraints );
			 numEx::BC_apply( enums::id_boundary_yPlus, enums::z, 0, apply_dirichlet_bc, dof_handler_ref, fe, constraints );
		 }
		 else if ( BC_yPlus==enums::BC_fix)
			numEx::BC_apply_fix( enums::id_boundary_yPlus, dof_handler_ref, fe, constraints );

		// BC on z0 plane ...
		 if ( dim==3 ) // ... only for 3D
		 {
			// For compression we don't fix anything in the third direction, because y0 was already clamped.
			// @todo However, what about the upper part?
			if ( BC_zMinus == enums::BC_sym )
				numEx::BC_apply( enums::id_boundary_zMinus, enums::z, 0, apply_dirichlet_bc, dof_handler_ref, fe, constraints );
		 }
		 
		// BC for the yPlus
		 if ( constrain_sideways_sliding_of_loaded_face && BC_yPlus==enums::BC_x0 )
			numEx::BC_apply( enums::id_boundary_yPlus, enums::x, 0, apply_dirichlet_bc, dof_handler_ref, fe, constraints );

		// BC for the load ...
		 if ( parameter.driver == enums::Dirichlet )  // ... as Dirichlet only for Dirichlet as driver
			numEx::BC_apply( id_boundary_load, loading_direction, current_load_increment, apply_dirichlet_bc, dof_handler_ref, fe, constraints );
	}

	// 2D grid
	template <int dim>
	void make_grid( Triangulation<2> &triangulation, const Parameter::GeneralParameters &parameter )
	{

	}
	
	// 3D grid
	template <int dim>
	void make_grid( Triangulation<3> &triangulation, const Parameter::GeneralParameters &parameter )
	{
		parameterCollection parameters_internal;
		const double search_tolerance = parameters_internal.search_tolerance;

		// Assign the dimensions of the hyper rectangle and store them as characteristic lengths
		 const double width = 18;
		 body_dimensions[enums::x] = width;
		 const double length = 60;
		 body_dimensions[enums::y] = length;
		 const double thickness = 1.5;
		 body_dimensions[enums::z] = thickness;

		  GridIn<dim> grid_in;
		  grid_in.attach_triangulation(triangulation);
		  std::string path_to_inp;
		  switch( parameter.refine_special )
		  {
		  	  case 1:
				  path_to_inp = "../Numerical_examples_in_dealii/DENP_Laura_Abaqus_mesh/DENP_Laura_Abaqus_mesh_m1_792el.inp";
				  break;
		  	  case 2:
				  path_to_inp = "../Numerical_examples_in_dealii/DENP_Laura_Abaqus_mesh/DENP_Laura_Abaqus_mesh_m2_2316el.inp";
				  break;
		  	  case 3:
				  path_to_inp = "../Numerical_examples_in_dealii/DENP_Laura_Abaqus_mesh/DENP_Laura_Abaqus_mesh_m3_3584el.inp";
				  break;
		  	  case 4:
				  path_to_inp = "../Numerical_examples_in_dealii/DENP_Laura_Abaqus_mesh/DENP_Laura_Abaqus_mesh_m4_3808el.inp";
				  break;
		  	  case 5:
				  path_to_inp = "../Numerical_examples_in_dealii/DENP_Laura_Abaqus_mesh/DENP_Laura_Abaqus_mesh_m5_6040el.inp";
				  break;
		  	  default:
		  		  AssertThrow(false,ExcMessage("DENP_Laura_Abaqus_mesh::make_grid<< Not available mesh chosen via refine_special."));
		  		  break;
		  }
		  std::ifstream input_file(path_to_inp);

		  grid_in.read_abaqus(input_file);

		// Clear all existing boundary ID's
		 numEx::clear_boundary_IDs( triangulation );

		// Set boundary IDs
		for (typename Triangulation<dim>::active_cell_iterator
			 cell = triangulation.begin_active();
			 cell != triangulation.end(); ++cell)
		{
			for (unsigned int face=0; face < GeometryInfo<dim>::faces_per_cell; ++face)
			{
			  if (cell->face(face)->at_boundary())
			  {
				// Set boundary IDs
				if (std::abs(cell->face(face)->center()[enums::x] - 0.0) < search_tolerance)
					cell->face(face)->set_boundary_id(enums::id_boundary_xMinus);
				else if ( std::abs(cell->face(face)->center()[enums::x] - body_dimensions[enums::x] ) < search_tolerance)
					cell->face(face)->set_boundary_id(enums::id_boundary_xPlus);
				else if (std::abs(cell->face(face)->center()[enums::y] - 0.0) < search_tolerance)
					cell->face(face)->set_boundary_id(enums::id_boundary_yMinus);
				else if (std::abs(cell->face(face)->center()[enums::y] - body_dimensions[enums::y]) < search_tolerance)
					cell->face(face)->set_boundary_id(enums::id_boundary_yPlus);
				else if ( std::abs(cell->face(face)->center()[enums::z] - 0.0) < search_tolerance)
					cell->face(face)->set_boundary_id(enums::id_boundary_zMinus);
				else if ( std::abs(cell->face(face)->center()[enums::z] - body_dimensions[enums::z]) < search_tolerance)
					cell->face(face)->set_boundary_id(enums::id_boundary_zPlus);
			  }
			}
		} // end for(cell)

		// Evaluation points and the related list of them
		 numEx::EvalPointClass<3> eval_center ( Point<3>(width-2.5,0,0), enums::x );
		 numEx::EvalPointClass<3> eval_top ( Point<3>(body_dimensions[enums::x],body_dimensions[enums::y],0), enums::x );

		 eval_points_list = {eval_center,eval_top};

		// Output the triangulation as eps or inp
//		 numEx::output_triangulation( triangulation, enums::output_eps_as_2D, numEx_name, true );
	}
}
