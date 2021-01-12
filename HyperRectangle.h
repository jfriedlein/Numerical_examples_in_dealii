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
#include "./numEx-helper_fnc.h"

using namespace dealii;

namespace HyperRectangle
/*
 * A hyper rectangle with three symmetry constraints, loaded in y-direction
 *
 * CERTIFIED TO STANDARD xxx
 */
{
	// Name of the numerical example
	 std::string numEx_name = "HyperRectangle";
 
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

	// All additional parameters
	const bool trigger_localisation_by_notching = true;
	const enums::enum_coord notched_face = enums::x;
	const bool use_fine_and_coarse_brick = true;
	const bool refine_local_isotropic = false;
	const bool notch_round = true;
	
	// Boundary conditions
	 const bool constrain_sideways_sliding_of_loaded_face = false;
	 const bool apply_sym_constraint_on_top_face = false; // to simulate plane strain for 3D
	 
	 const enums::enum_BC BC_xMinus = enums::BC_sym; // standard, tension
	 const enums::enum_BC BC_yMinus = enums::BC_sym; // standard, tension
	 const enums::enum_BC BC_zMinus = enums::BC_sym; // standard, tension

//	 const enums::enum_BC BC_xMinus = enums::BC_none; // compression, Seupel et al
//	 const enums::enum_BC BC_yMinus = enums::BC_fix; // compression, Seupel et al
//	 const enums::enum_BC BC_zMinus = enums::BC_none; // 3D compression, Seupel et al

	// Notching
	 const types::manifold_id manifold_id_notch_left = 10;
	 
	// Loading type and required modifications
	// @note "compression": \n
	// We use different boundary conditions and notch the body in the middle of its length and not a y=0.
	// @note "tension" or "standard": \n
	// We model 1/8 of the entire body and notch the body at y=0 (equals the middle of the entire body).
	 const enums::enum_loading_type loading_type = enums::Brick_Seupel_etal_a;

	 
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
		 
		// BC on y0 plane
		 if ( BC_yMinus==enums::BC_fix )
		 {
			// For compression we fix/clamp the Y0 plane, so it does not run away
			 numEx::BC_apply_fix( enums::id_boundary_yMinus, dof_handler_ref, fe, constraints );
		 }
		 else if ( BC_yMinus==enums::BC_sym)
			numEx::BC_apply( enums::id_boundary_yMinus, enums::y, 0, apply_dirichlet_bc, dof_handler_ref, fe, constraints );
		 
		// BC on z0 plane ...
		 if ( dim==3 ) // ... only for 3D
		 {
			// For compression we don't fix anything in the third direction, because y0 was already clamped.
			// @todo However, what about the upper part?
			if ( BC_zMinus == enums::BC_sym )
				numEx::BC_apply( enums::id_boundary_zMinus, enums::z, 0, apply_dirichlet_bc, dof_handler_ref, fe, constraints );

			if ( apply_sym_constraint_on_top_face )
				numEx::BC_apply( enums::id_boundary_zPlus, enums::z, 0, apply_dirichlet_bc, dof_handler_ref, fe, constraints );
		 }
		 
		// BC for the loaded face
		 if ( constrain_sideways_sliding_of_loaded_face )
			numEx::BC_apply( id_boundary_load, enums::x, 0, apply_dirichlet_bc, dof_handler_ref, fe, constraints );
		 
		// BC for the load ...
		 if ( parameter.driver == enums::Dirichlet )  // ... as Dirichlet only for Dirichlet as driver
			numEx::BC_apply( id_boundary_load, loading_direction, current_load_increment, apply_dirichlet_bc, dof_handler_ref, fe, constraints );
	}


	void make_grid_flat( Triangulation<2> &tria_flat,
						 const double &length, const double &width, const numEx::NotchClass<2> &notch,
						 const unsigned int n_elements_in_x_for_coarse_mesh, const unsigned int n_refine_global, const unsigned int n_refine_local )
	{
		parameterCollection parameters_internal;
		const double search_tolerance = parameters_internal.search_tolerance;

		// The ratio of the y and x lengths (typically greater than one)
		 const double edge_length_ratio = length / width;
		 
		// Using \a ceil(*), because elements are typically elongated in this direction (at least under tension)
		 const unsigned int n_elements_in_y_for_homogeneous_mesh = n_elements_in_x_for_coarse_mesh * std::ceil( edge_length_ratio );
		 const unsigned int n_elements_in_y_overhead = n_elements_in_x_for_coarse_mesh * std::ceil( std::ceil( edge_length_ratio ) - edge_length_ratio ) ;

		// Refine at least a square part (widthxwidth) if desired. In case the plate is wider than long, we limit the length to the 0.9 *length,
		// so we still leave a coarse part. If you want to limit the size of the refined fraction, just reduce the (1.*width).
		// @todo-optimize The 0.1 coarse part is rather silly, maybe switch to uniform refinement instead for such a case.
		 const double length_refined = std::min( 1.*width, 0.9 * length );
		 
		// Created the base mesh from a brick, either as ...
		 if ( use_fine_and_coarse_brick ) // ... a fine and a coarse part or ...
		 {
			// The bricks are spanned by three points (p1,p2,p3). The bar is created from two bricks, 
			// where the first will be meshed very fine (p1->p2) and the second remains coarse (p2->p3).
			 Point<2> p1 (0,0);
			 Point<2> p2 (width, length_refined); // extends in y-direction its length (loaded in y-direction as the other models)
			 Point<2> p3 (0, length);
			
			// Splitting the elements up into the coarse and fine parts
			 const unsigned int n_elements_in_y_for_fine_mesh = n_elements_in_x_for_coarse_mesh + n_elements_in_y_overhead;
			 const unsigned int n_elements_in_y_for_coarse_mesh = n_elements_in_y_for_homogeneous_mesh - n_elements_in_y_for_fine_mesh;
			 
			// Vector containing the number of elements in each dimension
			// The coarse segment consists of the set number of elements in the y-direction
			// @note The number of global refinements are introduced "softer", which means that we do not use powers of 2 as for the local refinements.
			 std::vector<unsigned int> repetitions_coarse (2);
			 repetitions_coarse[enums::x] = n_elements_in_x_for_coarse_mesh * ( n_refine_global + 1 );
			 repetitions_coarse[enums::y] = n_elements_in_y_for_coarse_mesh * ( n_refine_global + 1 );

			 // @todo Add option for automatic adjustment of the nbr of elements to the given notch geometry
			 
			// The fine segment consists of at least 2 elements plus possible refinements
			 std::vector<unsigned int> repetitions_fine (2);
			 repetitions_fine[enums::x] = repetitions_coarse[enums::x]; // equal to the coarse part, because we are not able to introduce hanging nodes like that
			 // In y-direction we refine the mesh 2^n times, only if we don't want the refinements to be isotropic (refine_local_isotropic==false)
			  repetitions_fine[enums::y] = n_elements_in_y_for_fine_mesh * std::pow(2., n_refine_local * !refine_local_isotropic) * (n_refine_global+1);

			Triangulation<2> triangulation_fine, triangulation_coarse;
			{
				// The fine brick
				 GridGenerator::subdivided_hyper_rectangle ( triangulation_fine,
															 repetitions_fine,
															 p1,
															 p2 );
	
				// The coarse brick
				 GridGenerator::subdivided_hyper_rectangle ( triangulation_coarse,
															 repetitions_coarse,
															 p2,
															 p3 );
			}

			// Merging fine and coarse brick
			// @note The interface between the two bricks needs to be meshed identically.
			// deal.II cannot detect hanging nodes there.
			 GridGenerator::merge_triangulations( triangulation_fine,
												  triangulation_coarse,
												  tria_flat,
												  1e-9 * length /*merge tolerance at interface*/ );
		 }
		 else // ... using a uniform brick with xy refinements
		 {
			 std::vector<unsigned int> repetitions (2);
			 repetitions[enums::x] = n_elements_in_x_for_coarse_mesh * (n_refine_global+1);
			 repetitions[enums::y] = n_elements_in_y_for_homogeneous_mesh * (n_refine_global+1);

			 Point<2> p1 (0,0);
			 Point<2> p2 (width, length); // extends in y-direction its length (loaded in y-direction as the other models)

			 GridGenerator::subdivided_hyper_rectangle ( tria_flat,
														 repetitions,
														 p1,
														 p2 );
		 }

		// Clear all existing boundary ID's
		 numEx::clear_boundary_IDs( tria_flat );

		// Set boundary IDs
		for (typename Triangulation<2>::active_cell_iterator
			 cell = tria_flat.begin_active();
			 cell != tria_flat.end(); ++cell)
		{
			for (unsigned int face=0; face < GeometryInfo<2>::faces_per_cell; ++face)
			  if (cell->face(face)->at_boundary())
			  {
				// Set boundary IDs
				if (std::abs(cell->face(face)->center()[enums::x] - 0.0) < search_tolerance)
				{
					cell->face(face)->set_boundary_id(enums::id_boundary_xMinus);
				}
				else if ( std::abs(cell->face(face)->center()[enums::x] - body_dimensions[enums::x] ) < search_tolerance)
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
				else
				{
					// There are only 6 faces for a cube in 3D, so if we missed one, something went terribly wrong
					AssertThrow(false, ExcMessage( numEx_name+" - make_grid 3D<< Found an unidentified face at the boundary. "
												   "Maybe it slipt through the assignment or that face is simply not needed. "
												   "So either check the implementation or comment this line in the code") );
				}
			  }
		} // end for(cell)

		
		// Notch the brick
		 if ( trigger_localisation_by_notching )
		 {
			 if ( notch_round )
			 {
				 // prepare the mesh
				  numEx::prepare_tria_for_notching( tria_flat, notch );
				 numEx::notch_body( tria_flat, notch );
				 
				 const Point<2> cyl_center_2D ( notch.cyl_center[0], notch.cyl_center[1] );
				 
				 // The manifolds are only usable for the 2D mesh not for 3D
				  SphericalManifold<2> spherical_manifold ( cyl_center_2D );
				  tria_flat.set_manifold( manifold_id_notch_left, spherical_manifold );
			 }
//			 else
//			 {
//				 Point<3> face_normal(-1.,0,0);
//				 Point<3> notch_reference_point ( 0, width*0.6, 0);
//				 double notch_depth = (1.-notch_reduction)*width;
//				 numEx::NotchClass<2> notch ( enums::notch_linear, notch_length, notch_depth, notch_reference_point, enums::id_boundary_xMinus,
//												face_normal, enums::y );
//				 
//				 numEx::prepare_tria_for_notching( triangulation, notch );
//				 numEx::notch_body( triangulation, notch );
//			 }
		 }
	}

	
	// 2D grid
	template <int dim>
	void make_grid( Triangulation<2> &triangulation, const Parameter::GeneralParameters &parameter )
	{
		parameterCollection parameters_internal;
		const double search_tolerance = parameters_internal.search_tolerance;

		// Assign the dimensions of the hyper rectangle and store them as characteristic lengths
		 const double width = parameter.width;
		 body_dimensions[enums::x] = width;
		 const double length = parameter.height;
		 body_dimensions[enums::y] = length;
		
		const unsigned int n_elements_in_x_for_coarse_mesh = parameter.grid_y_repetitions;
		const double notch_reduction = parameter.ratio_x;
		 Point<3> notch_reference_point ( width, width, 0);
		 Point<3> face_normal(1.,0,0);
		Point<3> axis_dir (0,0,1);

		 double notch_depth = (1.-notch_reduction)*width;
		 
		 // notching
		 numEx::NotchClass<2> notch ( enums::notch_round, parameter.notchWidth, notch_depth, notch_reference_point, enums::id_boundary_xPlus,
										face_normal, enums::y, manifold_id_notch_left );
		 
		make_grid_flat( triangulation, length, width, notch,
						n_elements_in_x_for_coarse_mesh, parameter.nbr_global_refinements, parameter.nbr_holeEdge_refinements );
			
	}
	

// 3d grid
	template<int dim>
	void make_grid( Triangulation<3> &triangulation, const Parameter::GeneralParameters &parameter )
	{
		parameterCollection parameters_internal;
		const double search_tolerance = parameters_internal.search_tolerance;

		// Assign the dimensions of the hyper rectangle and store them as characteristic lengths
		 const double width = parameter.width;
		 body_dimensions[enums::x] = width;
		 const double length = parameter.height;
		 body_dimensions[enums::y] = length;
		 const double thickness = parameter.thickness;
		 body_dimensions[enums::z] = thickness;
		
		const unsigned int n_elements_in_x_for_coarse_mesh = parameter.grid_y_repetitions;
		const double notch_reduction = parameter.ratio_x;
		 Point<3> notch_reference_point ( width, width, 0);
		 Point<3> face_normal(1.,0,0);
		Point<3> axis_dir (0,0,1);

		 double notch_depth = (1.-notch_reduction)*width;
		 
		 // notching
		 numEx::NotchClass<2> notch ( enums::notch_round, parameter.notchWidth, notch_depth, notch_reference_point, enums::id_boundary_xPlus,
										face_normal, enums::y, manifold_id_notch_left );
		 
		Triangulation<2> tria_flat;
		make_grid_flat( tria_flat, length, width, notch,
						n_elements_in_x_for_coarse_mesh, parameter.nbr_global_refinements, parameter.nbr_holeEdge_refinements );
		
		GridGenerator::extrude_triangulation( tria_flat, parameter.nbr_elementsInZ, parameter.thickness, triangulation, true );

		// Redo the manifold for 3D
		 CylindricalManifold<3> cylindrical_manifold (axis_dir, notch.cyl_center);
		 triangulation.set_manifold( manifold_id_notch_left, cylindrical_manifold );
		
		// Set boundary IDs
		for (typename Triangulation<3>::active_cell_iterator
			 cell = triangulation.begin_active();
			 cell != triangulation.end(); ++cell)
		{
			for (unsigned int face=0; face < GeometryInfo<2>::faces_per_cell; ++face)
			  if (cell->face(face)->at_boundary())
			  {
				if ( std::abs(cell->face(face)->center()[enums::z] - 0.0) < search_tolerance)
					cell->face(face)->set_boundary_id(enums::id_boundary_zMinus);
				else if ( std::abs(cell->face(face)->center()[enums::z] - body_dimensions[enums::z]) < search_tolerance)
					cell->face(face)->set_boundary_id(enums::id_boundary_zPlus);
			  }
		}
		 
		// In case we want local refinements to be isotropic (so refined in y and x and z).
		// We do this after the notching, so we can get away with slightly deeper notches.
//		 if ( refine_local_isotropic )
//		 {
//			for ( unsigned int nbr_local_ref=0; nbr_local_ref < parameter.nbr_holeEdge_refinements; nbr_local_ref++ )
//			{
//				for (typename Triangulation<dim>::active_cell_iterator
//							 cell = triangulation.begin_active();
//							 cell != triangulation.end(); ++cell)
//				{
//					for (unsigned int vertex=0; vertex < GeometryInfo<dim>::vertices_per_cell; ++vertex)
//					{
//						// Find all cells that lay in an exemplary damage band with size 1/4 from the y=0 face
//						if ( cell->vertex(vertex)[enums::y] < length_refined )
//						{
//							cell->set_refine_flag();
//							break; // break the for(vertex)
//						}
//					} // end for(vertex)
//				} // end for(cell)
//				triangulation.execute_coarsening_and_refinement();
//			}
//		 } // end if(refine_local_isotropic)
		
		// Output the triangulation as eps or inp
		 //numEx::output_triangulation( triangulation, enums::output_eps, numEx_name );
	}
}
