// Numerical example helper function (required, can be downloaded from https://github.com/jfriedlein/Numerical_examples_in_dealii)
// also contains enumerators as part of "enums::"
#include "numEx-helper_fnc.h"

#include "numEx-baseClass.h"

/**
 * @brief 
 * CERTIFIED TO STANDARD S22
 * 
 * @tparam dim 
 */
template<int dim>
class numEx_HyperCube : public numEx_class<dim>
{
  public:
    std::string numEx_name() {
    	return "HyperCube";
    };
	
    unsigned int loading_direction() {
		return enums::y;
	};	
    
    std::vector< enums::enum_boundary_ids > id_boundary_loads()
	{
    	// The loaded faces:
    	const enums::enum_boundary_ids id_boundary_load = enums::id_boundary_yPlus;
    	const enums::enum_boundary_ids id_boundary_secondaryLoad = enums::id_boundary_xPlus;
    	
		std::vector< enums::enum_boundary_ids > id_boundary_loads_list (2);
		id_boundary_loads_list[enums::id_primary_load] = id_boundary_load;
		id_boundary_loads_list[enums::id_secondary_load] = id_boundary_secondaryLoad;

		return id_boundary_loads_list;
	};	

    void make_grid ( /*input-> */ const Parameter::GeneralParameters &parameter,
					 /*output->*/ Triangulation<dim> &triangulation, 
								  std::vector<double> &body_dimensions,
								  std::vector< numEx::EvalPointClass<3> > &eval_points_list,
								  const std::string relativePath
					);
    
    void make_constraints ( AffineConstraints<double> &constraints, const FESystem<dim> &fe, DoFHandler<dim> &dof_handler_ref,
    						const bool &apply_dirichlet_bc, double &current_load_increment, const Parameter::GeneralParameters &parameter );

};


// HyperCube grid: 2D and 3D
template<int dim>
void numEx_HyperCube<dim>::make_grid 
	(
		/*input-> */ const Parameter::GeneralParameters &parameter,
		/*output->*/ Triangulation<dim> &triangulation, 
					 std::vector<double> &body_dimensions,
					 std::vector< numEx::EvalPointClass<3> > &eval_points_list,
					 const std::string relativePath
	)
{
	const double search_tolerance = numEx_class<dim>::search_tolerance;
	const double width = parameter.width;

	// Assign the characteristic dimensions of the cube
	 body_dimensions[enums::x] = width;
	 body_dimensions[enums::y] = width;
	 body_dimensions[enums::z] = width;
	
	// Evaluation points and the related list of them
	 const numEx::EvalPointClass<3> eval_top ( Point<3>(0,body_dimensions[enums::y],0), enums::y );
	 eval_points_list[enums::eval_point_0] = eval_top;

	 const numEx::EvalPointClass<3> eval_path_start ( Point<3>(0,0,0), enums::coord_none );
	 const numEx::EvalPointClass<3> eval_path_end = eval_top;
	 eval_points_list[enums::eval_path0_start] = eval_path_start;
	 eval_points_list[enums::eval_path0_end  ] = eval_path_end;

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
				AssertThrow(false, ExcMessage( numEx_name()+" - make_grid<< Found an unidentified face at the boundary. "
											   "Maybe it slipt through the assignment or that face is simply not needed. "
											   "So either check the implementation or comment this line in the code") );
				// Just to get rid of the unused variable warning
				 std::cout << "Relative path " << relativePath << std::endl;
			}
		  }
	 }

	// Refinement
	 triangulation.refine_global( parameter.nbr_global_refinements );

	// Output the triangulation as eps or inp
	 //numEx::output_triangulation( triangulation, enums::output_eps, numEx_name );
}


/**
 * Apply the boundary conditions (support and load) on the given AffineConstraints \a constraints. \n
 * For the HyperCube that are three symmetry constraints on each plane (x=0, y=0, z=0) and the load on the \a id_boundary_load (for Dirichlet).
 * Alternatively, we apply the rigid body motion for contact.
 * @todo Change setup (see HyperR) where we define boundary for BC_xMinus and then use if ( BC_xMinus==... ), use the value of BC_xMinus directly
 */
template<int dim>
void numEx_HyperCube<dim>::make_constraints ( AffineConstraints<double> &constraints, const FESystem<dim> &fe, DoFHandler<dim> &dof_handler_ref,
											  const bool &apply_dirichlet_bc, double &current_load_increment, const Parameter::GeneralParameters &parameter )
{
	// BC on x0 plane
	 numEx::BC_apply( enums::id_boundary_xMinus, enums::x, 0, apply_dirichlet_bc, dof_handler_ref, fe, constraints );

	// BC on y0 plane
	 numEx::BC_apply( enums::id_boundary_yMinus, enums::y, 0, apply_dirichlet_bc, dof_handler_ref, fe, constraints );

	// Fix the bottom face
	 //numEx::BC_apply_fix( enums::id_boundary_yMinus, dof_handler_ref, fe, constraints );
	 
	// BC on z0 plane ...
	 if ( dim==3 ) // ... only for 3D
		numEx::BC_apply( enums::id_boundary_zMinus, enums::z, 0, apply_dirichlet_bc, dof_handler_ref, fe, constraints );

	// BC for the load ...
	 if ( parameter.driver == enums::Dirichlet )  // ... as Dirichlet only for Dirichlet as driver, alternatively  ...
	 {
		const std::vector< enums::enum_boundary_ids > id_boundary_loads_list = id_boundary_loads(); 
		const unsigned int loading_dir = loading_direction();

		numEx::BC_apply( id_boundary_loads_list[enums::id_primary_load], loading_dir, current_load_increment, apply_dirichlet_bc, dof_handler_ref, fe, constraints );
		//numEx::BC_apply( id_boundary_loads_list[1], enums::x, current_load_increment, apply_dirichlet_bc, dof_handler_ref, fe, constraints );
	 }
	 else if ( parameter.driver == enums::Contact ) // ... as contact
	 {
		//if (apply_dirichlet_bc == true )
		//	rigid_wall->move( current_load_increment );
 	 }
}
