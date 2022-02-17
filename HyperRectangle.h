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

#include "numEx-baseClass.h"

/**
 * @brief 
 * CERTIFIED TO STANDARD S22
 * 
 * @tparam dim 
 */
template<int dim>
class numEx_HyperRectangle : public numEx_class<dim>
{
  public:
    std::string numEx_name() {
    	return "HyperRectangle";
    };
	
    unsigned int loading_direction() {
		return enums::y;
	};	
    
    std::vector< enums::enum_boundary_ids > id_boundary_loads()
	{
		std::vector< enums::enum_boundary_ids > id_boundary_loads_list (2);
		// The loaded faces:
		 id_boundary_loads_list[enums::id_primary_load] = enums::id_boundary_yPlus;
		 id_boundary_loads_list[enums::id_secondary_load] = enums::id_boundary_xPlus;

		return id_boundary_loads_list;
	};	

    std::vector< types::manifold_id > manifold_ids()
	{
		const types::manifold_id manifold_id_right_radius = 11;
		const types::manifold_id manifold_id_left_radius = 10;

		return { manifold_id_left_radius, manifold_id_right_radius };
	}

    void make_grid ( /*input-> */ const Parameter::GeneralParameters &parameter,
    				 /*output->*/ Triangulation<dim> &triangulation, 
    				 	 	 	  std::vector<double> &body_dimensions,
    							  std::vector< numEx::EvalPointClass<3> > &eval_points_list,
								  const std::string relativePath
    				);
    
    void make_constraints ( AffineConstraints<double> &constraints, const FESystem<dim> &fe, DoFHandler<dim> &dof_handler_ref,
    						const bool &apply_dirichlet_bc, double &current_load_increment, const Parameter::GeneralParameters &parameter );

	void make_grid_flat
	( 
		Triangulation<2> &tria_flat,
		const double &length, const double &width, 
		const std::vector< numEx::NotchClass<2> > &notch_list,
		const unsigned int n_elements_in_x_for_coarse_mesh, const unsigned int n_refine_global, 
		const unsigned int n_refine_local, std::vector<double> &body_dimensions, 
		const unsigned int aux_var
	);

};

	// All additional parameters
	const bool trigger_localisation_by_notching = true;
	const enums::enum_coord notched_face = enums::x;
	enums::enum_refine_special refine_special = enums::Mesh_HyperRectangle_coarse_and_fine_brick;
	const bool refine_local_isotropic = false;
	
	// Boundary conditions
	 const bool apply_sym_constraint_on_top_face = false; // to simulate plane strain for 3D, top face refers to zPlus
	 
	 // standard, tension:
	//   const enums::enum_BC BC_xMinus = enums::BC_sym;
	//   const enums::enum_BC BC_yPlus  = enums::BC_none;
	//   const bool constrain_sideways_sliding_of_loaded_face = false;
	//   const enums::enum_BC BC_yMinus = enums::BC_sym;
	//   const enums::enum_BC BC_zMinus = enums::BC_sym;
	//   const enums::enum_notch_type notch_type = enums::notch_linear;
	//   const bool notch_twice = false;
	//   const bool DENP_Laura = false;
	//   const bool DENP_Hagen = false;
	//   const bool SheStrip = true;

	  const bool Neto_planeStrain = false;
	  const bool refine_globally = false;

	 // DENP COMPLAS21
	  const enums::enum_BC BC_xMinus = enums::BC_none;
	  const enums::enum_BC BC_yPlus  = enums::BC_none;
	  const bool constrain_sideways_sliding_of_loaded_face = false;
	  const enums::enum_BC BC_yMinus = enums::BC_fix;
	  const enums::enum_BC BC_zMinus = enums::BC_sym;
	  const enums::enum_notch_type notch_type = enums::notch_round;
	  const bool notch_twice = true;
	  const bool DENP_Laura = false;
	  const bool DENP_Hagen = false;
	  const bool SheStrip = false;
	  const bool DENP_planeStrain = true;


	 // compression, Seupel et al:
//	 const enums::enum_BC BC_xMinus = enums::BC_none;
//	 const enums::enum_BC BC_yPlus  = enums::BC_none;//enums::BC_x0; //enums::BC_none; // guide top face
//	 const bool constrain_sideways_sliding_of_loaded_face = false;
//	 const enums::enum_BC BC_yMinus = enums::BC_fix;
//	 const enums::enum_BC BC_zMinus = enums::BC_none; // 3D compression, Seupel et al
//	 const enums::enum_notch_type notch_type = enums::notch_round;
//	 const bool notch_twice = true;
//	 const bool DENP_Laura = false;
//	 const bool DENP_Hagen = true;
//	  const bool SheStrip = false;

	 // 3D strip, left clamped, right pulled without contraction
//		const enums::enum_BC BC_yMinus = enums::BC_fix;
//		const enums::enum_BC BC_yPlus  = enums::BC_x0_z0;//enums::BC_x0; //enums::BC_none; // guide top face
//		const bool constrain_sideways_sliding_of_loaded_face = false;
//		const enums::enum_BC BC_xMinus = enums::BC_none;
//		const enums::enum_BC BC_zMinus = enums::BC_none;
//		const enums::enum_notch_type notch_type = enums::notch_round;
//		const bool notch_twice = false;
//		const bool DENP_Laura = false;
//		const bool DENP_Hagen = false;
//		const bool SheStrip = false;

	// Loading type and required modifications
	// @note "compression": \n
	// We use different boundary conditions and notch the body in the middle of its length and not a y=0.
	// @note "tension" or "standard": \n
	// We model 1/8 of the entire body and notch the body at y=0 (equals the middle of the entire body).
	 //const enums::enum_loading_type loading_type = enums::Brick_Seupel_etal_a;
	 //const enums::enum_loading_type loading_type = enums::compression;

	/**
	 * Apply the boundary conditions (support and load) on the given AffineConstraints \a constraints. \n
	 * For the HyperRectangle that are three symmetry constraints on each plane (x=0, y=0, z=0) and the load on the \a id_boundary_load (for Dirichlet).
	 */
template<int dim>
void numEx_HyperRectangle<dim>::make_constraints
	( 
		AffineConstraints<double> &constraints, const FESystem<dim> &fe, DoFHandler<dim> &dof_handler_ref,
		const bool &apply_dirichlet_bc, double &current_load_increment,
		const Parameter::GeneralParameters &parameter
	)
{
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
		
	// BC for the yPlus
	 if ( constrain_sideways_sliding_of_loaded_face && BC_yPlus==enums::BC_x0 )
		numEx::BC_apply( enums::id_boundary_yPlus, enums::x, 0, apply_dirichlet_bc, dof_handler_ref, fe, constraints );

	// BC for the load ...
	 if ( parameter.driver == enums::Dirichlet )  // ... as Dirichlet only for Dirichlet as driver
	 {
		const std::vector< enums::enum_boundary_ids > id_boundary_loads_list = id_boundary_loads(); 

		numEx::BC_apply( id_boundary_loads_list[enums::id_primary_load], loading_direction(), current_load_increment, apply_dirichlet_bc, dof_handler_ref, fe, constraints );
	 }
}

template<int dim>
void numEx_HyperRectangle<dim>::make_grid_flat
	( 
		Triangulation<2> &tria_flat,
		const double &length, const double &width, 
		const std::vector< numEx::NotchClass<2> > &notch_list,
		const unsigned int n_elements_in_x_for_coarse_mesh, const unsigned int n_refine_global, 
		const unsigned int n_refine_local, std::vector<double> &body_dimensions, 
		const unsigned int aux_var
	)
{
	// @todo Remove hardcoding, load from numEx_baseClass
	const double search_tolerance = 1e-12;

	// The ratio of the y and x lengths (typically greater than one)
	 const double edge_length_ratio = length / width;
		
	// Using \a ceil(*), because elements are typically elongated in this direction (at least under tension)
	unsigned int n_elements_in_y_for_homogeneous_mesh = n_elements_in_x_for_coarse_mesh * std::ceil( edge_length_ratio );

	// hardcoded
	if ( refine_special == enums::Mesh_refine_none )
		n_elements_in_y_for_homogeneous_mesh = aux_var; //15;

	const unsigned int n_elements_in_y_overhead = n_elements_in_x_for_coarse_mesh * std::ceil( std::ceil( edge_length_ratio ) - edge_length_ratio ) ;

	// Refine at least a square part (widthxwidth) if desired. In case the plate is wider than long, we limit the length to the 0.9 *length,
	// so we still leave a coarse part. If you want to limit the size of the refined fraction, just reduce the (1.*width).
	// @todo-optimize The 0.1 coarse part is rather silly, maybe switch to uniform refinement instead for such a case.
	double length_refined = std::min( 1.*width, 0.9 * length );

	if ( Neto_planeStrain )
		length_refined=length/6.;
		
	// Created the base mesh from a brick, either as ...
		if ( refine_special==enums::Mesh_HyperRectangle_coarse_and_fine_brick /*use_fine_and_coarse_brick*/ ) // ... a fine and a coarse part or ...
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

			 if ( Neto_planeStrain )
			 {
				 repetitions_coarse[enums::x] = 10;
				 repetitions_coarse[enums::y] = 10;
				 repetitions_fine[enums::x] = 10;
				 repetitions_fine[enums::y] = 10;
			 }

			 if ( SheStrip )
			 {
				 repetitions_coarse[enums::x] = 5;
				 repetitions_coarse[enums::y] = 14;
				 repetitions_fine[enums::x] = 5;
				 repetitions_fine[enums::y] = 6 * std::pow(2., n_refine_local-1);
			 }

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
			 if ( notch_twice && DENP_Laura )
			 {
				 repetitions[enums::x] = 5;
				 repetitions[enums::y] = 16;
			 }
			 else if ( notch_twice && DENP_Hagen )
			 {
				 repetitions[enums::x] = 3;
				 repetitions[enums::y] = 18;
			 }
			 else if ( notch_twice && DENP_planeStrain )
			 {
				 repetitions[enums::x] = 6;
				 repetitions[enums::y] = 24;
			 }
			 else
			 {
				 repetitions[enums::x] = n_elements_in_x_for_coarse_mesh * (n_refine_global+1);
				 repetitions[enums::y] = n_elements_in_y_for_homogeneous_mesh * (n_refine_global+1);
			 }

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
			{
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
//					AssertThrow(false, ExcMessage( numEx_name+" - make_grid 2D<< Found an unidentified face at the boundary. "
//												   "Maybe it slipt through the assignment or that face is simply not needed. "
//												   "So either check the implementation or comment this line in the code") );
				}
			  }
			}
		} // end for(cell)

		  if ( false /*pre-refine to improve notch mesh*/ )
		  {
			  	const double notch_offset = 10.;
				for (typename Triangulation<2>::active_cell_iterator
							 cell = tria_flat.begin_active();
							 cell != tria_flat.end(); ++cell)
				{
						// Find all cells that lay in an exemplary damage band
						 if ( std::abs( cell->center()[enums::y] - ( notch_offset/width * cell->center()[enums::x] + notch_list[1].cyl_center[enums::y] ) )
							  < 1.5*notch_list[1].length/2. )
							cell->set_refine_flag();
				} // end for(cell)
				tria_flat.execute_coarsening_and_refinement();
		  }
		
		// Notch the brick
		 if ( trigger_localisation_by_notching && notch_list[0].depth > 1e-20 )
		 {
			 //if ( notch_type == enums::notch_round )
			 {
				 // prepare the mesh
				  if ( notch_type == enums::notch_round )
					 numEx::prepare_tria_for_notching( tria_flat, notch_list[0] );

				 numEx::notch_body( tria_flat, notch_list[0] );
				 
				const std::vector< types::manifold_id > manifold_ids_list = manifold_ids();
				const types::manifold_id manifold_id_notch_left = manifold_ids_list[0];
				const types::manifold_id manifold_id_notch_right = manifold_ids_list[1];

				 if ( notch_type == enums::notch_round )
				 {
					 const Point<2> cyl_center_2D ( notch_list[0].cyl_center[0], notch_list[0].cyl_center[1] );
					 // The manifolds are only usable for the 2D mesh not for 3D
					  SphericalManifold<2> spherical_manifold ( cyl_center_2D );
					  tria_flat.set_manifold( manifold_id_notch_right, spherical_manifold );
				 }
				 
				  if ( notch_twice )
				  {
					 numEx::prepare_tria_for_notching( tria_flat, notch_list[1] );
					 numEx::notch_body( tria_flat, notch_list[1] );

					 const Point<2> cyl_center_2D ( notch_list[1].cyl_center[0], notch_list[1].cyl_center[1] );

					 // The manifolds are only usable for the 2D mesh not for 3D
					  SphericalManifold<2> spherical_manifold ( cyl_center_2D );
					  tria_flat.set_manifold( manifold_id_notch_left, spherical_manifold );
				  }
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

	// Output the triangulation as eps or inp
	 //numEx::output_triangulation( tria_flat, enums::output_eps, numEx_name );
}

	
// 2D grid
template <int dim>
void numEx_HyperRectangle<dim>::make_grid
	( 
		/*input-> */ const Parameter::GeneralParameters &parameter,
		/*output->*/ Triangulation<dim> &triangulation,
					 std::vector<double> &body_dimensions,
					 std::vector< numEx::EvalPointClass<3> > &eval_points_list,
					 const std::string relativePath
	) 
{
	AssertThrow(dim==2, ExcMessage(numEx_name()+" << not yet available for 3D"));

	refine_special = enums::enum_refine_special(parameter.refine_special);

	// Assign the dimensions of the hyper rectangle and store them as characteristic lengths
		const double width = parameter.width;
		body_dimensions[enums::x] = width;
		const double length = parameter.height;
		body_dimensions[enums::y] = length;
		const unsigned int aux_var = int(parameter.holeRadius);

		double notch_offset;
		if ( DENP_Laura )
			notch_offset = 10;
		else if ( DENP_Hagen )
			notch_offset = 20;
		else
			notch_offset = width;

		 // double notch for compression or bottom notch for tension
		  const double notch_y_right = notch_twice ? (length/2.+notch_offset/2.) : 0;
		 const double notch_y_left = length/2. - notch_offset/2.;

		const std::vector< types::manifold_id > manifold_ids_list = manifold_ids();
		const types::manifold_id manifold_id_notch_left = manifold_ids_list[0];
		const types::manifold_id manifold_id_notch_right = manifold_ids_list[1];

		// notching
		 // First notch on the right
		  const unsigned int n_elements_in_x_for_coarse_mesh = parameter.grid_y_repetitions;
		  const double notch_reduction = parameter.ratio_x;
		  Point<3> notch_reference_point1 ( width, notch_y_right, 0);
		  Point<3> face_normal1(1.,0,0);
		  double notch_depth = (1.-notch_reduction)*width;

		  numEx::NotchClass<2> notch1 ( notch_type, parameter.notchWidth, notch_depth, notch_reference_point1, enums::id_boundary_xPlus,
										face_normal1, enums::y, manifold_id_notch_right );

		 // Second notch on the left
		  Point<3> notch_reference_point2 ( 0, notch_y_left, 0);
		  Point<3> face_normal2(-1.,0,0);

		  numEx::NotchClass<2> notch2 ( notch_type, parameter.notchWidth, notch_depth, notch_reference_point2, enums::id_boundary_xMinus,
										face_normal2, enums::y, manifold_id_notch_left );

		// Create the 2D base mesh
		 if ( notch_twice )
			make_grid_flat( triangulation, length, width, {notch1,notch2},
							n_elements_in_x_for_coarse_mesh, parameter.nbr_global_refinements, parameter.nbr_holeEdge_refinements, body_dimensions, aux_var );
		 else if (notch_twice==false)
			make_grid_flat( triangulation, length, width, {notch1},
							n_elements_in_x_for_coarse_mesh, parameter.nbr_global_refinements, parameter.nbr_holeEdge_refinements, body_dimensions,  aux_var );
		 else // just to get rid of the unused variable warning
		 	std::cout << "relativePath " << relativePath << std::endl;
			 
		// Local refinements
		 if ( notch_twice )
		 {
			 for ( unsigned int nbr_local_ref=0; nbr_local_ref < parameter.nbr_holeEdge_refinements; nbr_local_ref++ )
			 {
				for (typename Triangulation<dim>::active_cell_iterator
							 cell = triangulation.begin_active();
							 cell != triangulation.end(); ++cell)
				{
						// Find all cells that lay in an exemplary damage band
						 if ( std::abs( cell->center()[enums::y] - ( notch_offset/width * cell->center()[enums::x] + notch_y_left ) )
						 	  < 1.75*parameter.notchWidth/2. )
							cell->set_refine_flag();
						 if ( refine_globally )
							 cell->set_refine_flag();
				} // end for(cell)
				triangulation.execute_coarsening_and_refinement();
			 }

			 // special case:
			 if ( DENP_Laura )
				 triangulation.refine_global(parameter.nbr_global_refinements);
		 }
		 else
		 {
			 for ( unsigned int nbr_local_ref=0; nbr_local_ref < parameter.nbr_holeEdge_refinements; nbr_local_ref++ )
			 {
				for (typename Triangulation<dim>::active_cell_iterator
							 cell = triangulation.begin_active();
							 cell != triangulation.end(); ++cell)
				{
						// Find all cells that lay in an exemplary damage band with size 1/4 from the y=0 face
						if ( cell->center()[enums::y] < width )
							cell->set_refine_flag();
				} // end for(cell)
				triangulation.execute_coarsening_and_refinement();
			 }
		}

	// Evaluation points and the related list of them
	//numEx::EvalPointClass<3> eval_center ( Point<3>(width-notch_depth,0,0), enums::x );
	const numEx::EvalPointClass<3> eval_top ( Point<3>(body_dimensions[enums::x],body_dimensions[enums::y],0), enums::y );

	eval_points_list[enums::eval_point_0] = eval_top;
}