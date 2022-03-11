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

// Append the \a enums namespace
namespace enums
{
	/**
	 * Clamping types
	 */
	 enum enum_Beam_clamping
	 {
		beam_clamped_free = 0,   //!< clamped_free: left face is fixed in all directions, right face is unconstrained
		beam_clamped_sliding = 1,//!< clamped_sliding: left face is fixed in all directions, right face can slide in y-direction
		beam_clamped_clamped = 2 //!< clamped_clamped: @todo Not yet implemented
	 };
}


/**
 * @brief 
 * CERTIFIED TO STANDARD S22
 * 
 * @tparam dim 
 */
template<int dim>
class numEx_Beam : public numEx_class<dim>
{
  public:
      std::string numEx_name() {
    	return "Beam";
    };

	// The loading direction: \n
	// In which coordinate direction the load shall be applied, so x/y/z.
    unsigned int loading_direction() {
		return enums::y; // enums::y, enums::x
	};	

    std::vector< enums::enum_boundary_ids > id_boundary_loads()
	{
    	// The loaded faces:
		std::vector< enums::enum_boundary_ids > id_boundary_loads_list (2);
		id_boundary_loads_list[enums::id_primary_load] = enums::id_boundary_xPlus;
		id_boundary_loads_list[enums::id_secondary_load] = enums::id_boundary_none;

		return id_boundary_loads_list;
	};	

	void make_grid ( /*input-> */ const Parameter::GeneralParameters &parameter,
    				 /*output->*/ Triangulation<2> &triangulation, 
    				 	 	 	  std::vector<double> &body_dimensions,
    							  std::vector< numEx::EvalPointClass<3> > &eval_points_list,
								  const std::string relativePath
    				);
	void make_grid ( /*input-> */ const Parameter::GeneralParameters &parameter,
    				 /*output->*/ Triangulation<3> &triangulation, 
    				 	 	 	  std::vector<double> &body_dimensions,
    							  std::vector< numEx::EvalPointClass<3> > &eval_points_list,
								  const std::string relativePath
    				);
    
    void make_constraints ( AffineConstraints<double> &constraints, const FESystem<dim> &fe, DoFHandler<dim> &dof_handler_ref,
    						const bool &apply_dirichlet_bc, double &current_load_increment, const Parameter::GeneralParameters &parameter );

};


	template<int dim>
	void numEx_Beam<dim>::make_constraints ( AffineConstraints<double> &constraints, const FESystem<dim> &fe, DoFHandler<dim> &dof_handler_ref,
							const bool &apply_dirichlet_bc, double &current_load_increment,
							const Parameter::GeneralParameters &parameter )
	{
		// USER PARAMETER
		const unsigned int beam_type = enums::beam_clamped_free;

		// clamping on X0 plane: set x, y and z displacements on x0 plane to zero
		if ( beam_type == enums::beam_clamped_free )
 			numEx::BC_apply_fix( enums::id_boundary_xMinus, dof_handler_ref, fe, constraints );
		else
			numEx::BC_apply( enums::id_boundary_xMinus, enums::x, 0, apply_dirichlet_bc, dof_handler_ref, fe, constraints );

		// BC on z0 plane ...
		//  if ( dim==3 ) // ... only for 3D
		// 	numEx::BC_apply( enums::id_boundary_zMinus, enums::z, 0, apply_dirichlet_bc, dof_handler_ref, fe, constraints );

//		// Point-constraint
//		 // First, add lines to the constraints matrix
//		 const unsigned int constrained_dof = ( (parameter.degree==1) ? 19 : 43 );
//		 constraints.add_line(constrained_dof);
//		 // Then, we fill all desired non-zero entries, here we fill every entry even the entries that are zero
//		 // symmetric shear
//		  constraints.set_inhomogeneity(constrained_dof,0);

		// BC for the load ...
		 if ( parameter.driver == enums::Neumann )
		 {
			numEx::BC_apply_fix( enums::id_boundary_xMinus, dof_handler_ref, fe, constraints );
		 }
		 else if ( parameter.driver == enums::Dirichlet )  // ... as Dirichlet only for Dirichlet as driver, alternatively  ...
		 {
			 	 const std::vector< enums::enum_boundary_ids > id_boundary_loads_list = id_boundary_loads(); 

			//numEx::BC_apply_fix( enums::id_boundary_xMinus, dof_handler_ref, fe, constraints );
			numEx::BC_apply( id_boundary_loads_list[enums::id_primary_load], loading_direction(), current_load_increment, apply_dirichlet_bc, dof_handler_ref, fe, constraints );

			if ( beam_type==enums::beam_clamped_sliding )
				numEx::BC_apply( id_boundary_loads_list[enums::id_primary_load], enums::x, 0, apply_dirichlet_bc, dof_handler_ref, fe, constraints );

			// classical
//			 numEx::BC_apply( id_boundary_load, loading_direction, current_load_increment, apply_dirichlet_bc, dof_handler_ref, fe, constraints );

			// shear-free
//			{
//				const double lambda_n1 = lambda_n + current_load_increment;
//				const double current_load_increment_x = body_dimensions[enums::x] * ( std::sin(lambda_n1)/lambda_n1 - 1. ) - body_dimensions[enums::x] * ( std::sin(lambda_n+1e-20)/(lambda_n+1e-20) - 1. );
//				const double current_load_increment_y = body_dimensions[enums::x] * ( std::cos(lambda_n1)/lambda_n1 - 1./lambda_n1 ) - body_dimensions[enums::x] * ( std::cos(lambda_n+1e-10)/(lambda_n+1e-10) - 1./(lambda_n+1e-10) );
//
//				std::cout << "current_load_increment_x " << current_load_increment_x << std::endl;
//				std::cout << "current_load_increment_y " << current_load_increment_y << std::endl;
//
//				numEx::BC_apply( id_boundary_load, enums::x, current_load_increment_x, apply_dirichlet_bc, dof_handler_ref, fe, constraints );
//				numEx::BC_apply( id_boundary_load, enums::y, current_load_increment_y, apply_dirichlet_bc, dof_handler_ref, fe, constraints );
//			}

//			// prescribed pseudo pure bending load
//			numEx::BeamEnd<dim> beamEnd (body_dimensions[enums::x], body_dimensions[enums::y], lambda_n, current_load_increment, n_components);
//
//			if (apply_dirichlet_bc == true )
//			{
//				// Apply the given load
//				 VectorTools::interpolate_boundary_values(
//						 dof_handler_ref,
//															id_boundary_load,
//															beamEnd,
//															constraints
//														 );
//			}
//			else
//			{
//				VectorTools::interpolate_boundary_values(
//						dof_handler_ref,
//															id_boundary_load,
//															ZeroFunction<dim> ( n_components ),
//															constraints
//														);
//			}
		 }
		 else if ( parameter.driver == enums::Contact ) // ... as contact
		 {
	 	 }
	}


	// 2D grid
	template <int dim>
	void numEx_Beam<dim>::make_grid( /*input-> */ const Parameter::GeneralParameters &parameter,
									/*output->*/ Triangulation<2> &triangulation, 
												std::vector<double> &body_dimensions,
												std::vector< numEx::EvalPointClass<3> > &eval_points_list,
												const std::string relativePath )
	{
		const double search_tolerance = numEx_class<dim>::search_tolerance;

		// ToDo: use the values from the parameter file
		const double width = parameter.width; // use thickness=width for square bottom area
		const double thickness = parameter.thickness;
		const double length = parameter.height;

		// USER PARAMETERS
		 double refined_fraction=length/4.;

		body_dimensions[enums::x] = length;
		body_dimensions[enums::y] = width;
		body_dimensions[enums::z] = thickness;

		// The bar is created from two bricks, where the first will be meshed very fine
		// and the second remains coarse. The bricks are spanned by three points.
		 Point<dim> p1 (0,0);
		 Point<dim> p2 (length * refined_fraction, width); // extends in y-direction its height (loaded in y-direction as the othe models)
		 Point<dim> p3 (length, 0); // extends in y-direction its height (loaded in y-direction as the othe models)
		 Point<dim> p4 (length, width);


		 if ( parameter.refine_special == enums::Mesh_refine_beam_fineCoarseBrick )
		 {
			// Vector containing the number of elements in each dimension
			// The coarse segment consists of the set number of elements in the y-direction
			 std::vector<unsigned int> repetitions (dim);
			 repetitions[enums::x] = 6 * (parameter.nbr_global_refinements+1);
			 repetitions[enums::y] = 3; // y

			// The fine segment consists of at least 2 elements plus possible refinements
			 std::vector<unsigned int> repetitions_fine (dim);
			 repetitions_fine[enums::x] = 6 * (parameter.nbr_global_refinements+1);
			 repetitions_fine[enums::y] = 3; // y

			Triangulation<dim> triangulation_fine, triangulation_coarse;
			// The fine brick
			 GridGenerator::subdivided_hyper_rectangle ( triangulation_fine,
														 repetitions_fine,
														 p1,
														 p2 );

			// The coarse brick
			 GridGenerator::subdivided_hyper_rectangle ( triangulation_coarse,
														 repetitions,
														 p2,
														 p3 );

			// Merging fine and coarse brick
			// @note The interface between the two bricks needs to be meshed identically.
			// deal.II cannot detect hanging nodes there.
			GridGenerator::merge_triangulations( triangulation_fine,
												 triangulation_coarse,
												 triangulation,
												 1e-9 * length );

			// Local refinement
			if ( parameter.nbr_holeEdge_refinements >= 0 && parameter.nbr_global_refinements==0 )
			{
				for (typename Triangulation<dim>::active_cell_iterator
							 cell = triangulation.begin_active();
							 cell != triangulation.end(); ++cell)
				{
					for ( unsigned int face=0; face < GeometryInfo<dim>::faces_per_cell; face++ )
						if ( cell->center()[loading_direction()] < length * refined_fraction )
						{
							cell->set_refine_flag();
							break;
						}
				}
				triangulation.execute_coarsening_and_refinement();
			}
		 }
		else if ( parameter.refine_special == enums::Mesh_refine_uniform ) // use uniform brick with xy refinements
		{
			 std::vector<unsigned int> repetitions (dim);
			 repetitions[enums::x] = std::pow(2.,parameter.nbr_holeEdge_refinements);
			 repetitions[enums::y] = std::pow(2.,parameter.nbr_holeEdge_refinements); // y

			 GridGenerator::subdivided_hyper_rectangle ( triangulation,
														 repetitions,
														 p1,
														 p4 );
		}
		 else if ( parameter.refine_special == enums::Mesh_refine_beam_nx1 )
		 {
			 std::vector<unsigned int> repetitions (dim);
			 repetitions[enums::x] = 4;
			 repetitions[enums::y] = 1; // y

			 GridGenerator::subdivided_hyper_rectangle ( triangulation,
														 repetitions,
														 p1,
														 p4 );
		 }
		 else if ( parameter.refine_special == enums::Mesh_refine_beam_nxm )
		 {
			 std::vector<unsigned int> repetitions (dim);
			 repetitions[enums::x] = parameter.grid_y_repetitions;
			 repetitions[enums::y] = parameter.nbr_elementsInZ;

			 GridGenerator::subdivided_hyper_rectangle ( triangulation,
														 repetitions,
														 p1,
														 p4 );
		 }
		 else
			 AssertThrow(false,ExcMessage("Beam - make_grid<< Undefined refine_special case."))

		// Clear all existing boundary ID's
		 numEx::clear_boundary_IDs( triangulation );

		// Set boundary IDs and manifolds
		for (typename Triangulation<dim>::active_cell_iterator
			 cell = triangulation.begin_active();
			 cell != triangulation.end(); ++cell)
		{
			for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
			  if (cell->face(face)->at_boundary())
			  {
				//Set boundary IDs
				if (std::abs(cell->face(face)->center()[0] - 0.0) < search_tolerance)
				{
					cell->face(face)->set_boundary_id(enums::id_boundary_xMinus);
//					for ( unsigned int vertex=0; vertex < GeometryInfo<dim>::vertices_per_face; ++vertex)
//					{
//						if ( cell->face(face)->vertex(vertex).distance(Point<dim> (0,0)) < 1e-10 )
//						{
//							cell->face(face)->set_boundary_id(enums::id_boundary_xMinus2);
//						}
//					}
				}
				else if (std::abs(cell->face(face)->center()[enums::x] - body_dimensions[enums::x]) < search_tolerance)
				{
					cell->face(face)->set_boundary_id(enums::id_boundary_xPlus);
				}
				else if (std::abs(cell->face(face)->center()[1] - 0.0) < search_tolerance)
				{
					cell->face(face)->set_boundary_id(enums::id_boundary_yMinus);
//						cell->set_material_id( enums::tracked_QP );
				}
				else if (std::abs(cell->face(face)->center()[enums::y] - body_dimensions[enums::y]) < search_tolerance)
				{
					cell->face(face)->set_boundary_id(enums::id_boundary_yPlus);
				}
//					else
//					{
//						// There are just eight faces, so if we missed one, something went clearly terribly wrong
//						 AssertThrow(false, ExcMessage("Beam - make_grid 2D<< Found an unidentified face at the boundary. Maybe it slipt through the assignment or that face is simply not needed. So either check the implementation or comment this line in the code"));
//					}
			  }
		}
//
//			triangulation.refine_global(parameter.nbr_global_refinements);	// ... Parameter.prm file
//
//			// refine in a special manner only cells around the origin
//			{
//				for ( unsigned int nbr_local_ref=0; nbr_local_ref<parameter.nbr_holeEdge_refinements; nbr_local_ref++ )
//				{
//					for (typename Triangulation<dim>::active_cell_iterator
//								 cell = triangulation.begin_active();
//								 cell != triangulation.end(); ++cell)
//					{
//						for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
//						  {
//							// Find all cells that lay in an exemplary damage band with size 1.5 mm from the y=0 face
//							if (std::abs(cell->face(face)->center()[enums::y] ) < length/8. )
//							{
//								cell->set_refine_flag();
//								break;
//							}
//						  }
//					}
//					triangulation.execute_coarsening_and_refinement();
//				}
//			}
//
//			// Mark the innermost cell(s) for softening to trigger the damage development
//			if ( triangulation.n_active_cells()>1)
//			{
//				bool found_cell=false;
//				for (typename Triangulation<dim>::active_cell_iterator
//							 cell = triangulation.begin_active();
//							 cell != triangulation.end(); ++cell)
//				{
//					for (unsigned int vertex=0; vertex < GeometryInfo<dim>::vertices_per_cell; ++vertex)
//					 // Find the cell that has one vertex at the origin
//	//				  if ( (cell->vertex(vertex)).distance(centre)<1e-12 )
//					 // Find cells that lay on the xz-plane (y=0)
//					  if ( std::abs(cell->vertex(vertex)[1]) < 1 ) // vs 1e-12 used previously
//					  {
//						  // Avoid overwriting the tracked cell
//						  if ( cell->material_id() != enums::tracked_QP )
//							  cell->set_material_id(1);
//						  found_cell = true;
//						  break;
//					  }
//
//					if ( /*only weaken one cell:*/ false  && found_cell )	// ToDo: check: does a break leave only the inner for-loop?
//						break;
//				}
//
//				AssertThrow(found_cell, ExcMessage("Beam<< Was not able to identify the cell at the origin(0,0,0). Please recheck the triangulation or adapt the code."));
//			}
//
//
//			// Notch the specimen by moving some nodes inwards to form a notch
//			if ( damage_trigger_by_notching )
//			{
////				Point<2> edge_point (1,0);
////				Point<2> notching (-0.1,0);
////				// A quick assurance variable to assure that at least a single vertex has been found,
////				// so our search criterion where to look for the vertices is not completely off
////				 bool found_vertex=false;
////				// Looping over all cells to notch the to-be-notched cells
////				 for ( typename Triangulation<dim>::active_cell_iterator
////							 cell = triangulation.begin_active();
////							 cell != triangulation.end(); ++cell )
////				 {
////					for (unsigned int vertex=0; vertex < GeometryInfo<dim>::vertices_per_cell; ++vertex)
////					 // Find vertices that are in the first 1/16 of the entire length
////					  if (cell->vertex(vertex).distance(edge_point) < search_tolerance )
////					  {
////						  // The found vertex is moved by the \a notching vector
////						  // The notching shall be linear, hence a vertex in the notch is fully notched and the farther you
////						  // move away from the notch the lower the notching gets (linearly).
////						   cell->vertex(vertex) += notching;
////						  found_vertex = true;
////					  }
////				 }
//
//				// Declare the shift vector for the notching
//				 Point<2> notching; // initially zero
//				// Depending on the desired notching direction (notched_face),
//				// we set the according shift component to the overall reduction
//				 notching[notched_face] = - body_dimensions[notched_face] * ( 1.-notch_reduction );
//
//				// A quick assurance variable to assure that at least a single vertex has been found,
//				// so our search criterion where to look for the vertices is not completely off
//				 bool found_vertex=false;
//				// Looping over all cells to notch the to-be-notched cells
//				 for ( typename Triangulation<dim>::active_cell_iterator
//							 cell = triangulation.begin_active();
//							 cell != triangulation.end(); ++cell )
//				 {
//					for (unsigned int vertex=0; vertex < GeometryInfo<dim>::vertices_per_cell; ++vertex)
//					 // Find vertices that are in the first 1/16 of the entire length
//					  if ( std::abs(cell->vertex(vertex)[loading_direction]) <  notch_length )
//						  if ( std::abs( cell->vertex(vertex)[notched_face] - body_dimensions[notched_face]) < search_tolerance )
//						  {
//							  // The found vertex is moved by the \a notching vector
//							  // The notching shall be linear, hence a vertex in the notch is fully notched and the farther you
//							  // move away from the notch the lower the notching gets (linearly).
//							   cell->vertex(vertex) += notching * ( notch_length - cell->vertex(vertex)[loading_direction] ) / notch_length;
//							  found_vertex = true;
//						  }
//				 }
//
//				AssertThrow(found_vertex, ExcMessage("Beam<< We weren't able to find at least a single vertex to be notched."));
//			}
//
		// Output the triangulation as eps or inp
		 //numEx::output_triangulation( triangulation, enums::output_eps, numEx_name );
		}



// 3d grid
	template <int dim>
	void numEx_Beam<dim>::make_grid( /*input-> */ const Parameter::GeneralParameters &parameter,
									/*output->*/ Triangulation<3> &triangulation, 
												std::vector<double> &body_dimensions,
												std::vector< numEx::EvalPointClass<3> > &eval_points_list,
												const std::string relativePath )
	{
		const double search_tolerance = numEx_class<dim>::search_tolerance;

		// ToDo: use the values from the parameter file
		const double width = parameter.width; // use thickness=width for square bottom area
		const double thickness = parameter.thickness;
		const double length = parameter.height;

		// USER PARAMETERS
		 double refined_fraction=length/4.;

		body_dimensions[enums::x] = length;
		body_dimensions[enums::y] = width;
		body_dimensions[enums::z] = thickness;

		// The bar is created from two bricks, where the first will be meshed very fine
		// and the second remains coarse. The bricks are spanned by three points.
		 Point<dim> p1 (0,0,0);
		 Point<dim> p2 (length * refined_fraction, width, thickness); // extends in y-direction its height (loaded in y-direction as the othe models)
		 Point<dim> p3 (length, 0, 0); // extends in y-direction its height (loaded in y-direction as the othe models)
		 Point<dim> p4 (length, width, thickness);

		// Create the evaluation point
		{
			const numEx::EvalPointClass<3> eval_point_0 ( p4, loading_direction() );
			eval_points_list[enums::eval_point_0] = eval_point_0;
		}

		// Use fine and coarse brick
		if ( parameter.refine_special==enums::Mesh_refine_special_standard )
		{
			// Vector containing the number of elements in each dimension
			// The coarse segment consists of the set number of elements in the y-direction
			 std::vector<unsigned int> repetitions (3);
			 repetitions[enums::x] = 6 * (parameter.nbr_global_refinements+1);
			 repetitions[enums::y] = 3; // y
			 repetitions[enums::z] = parameter.nbr_elementsInZ;

			// The fine segment consists of at least 2 elements plus possible refinements
			 std::vector<unsigned int> repetitions_fine (3);
			 repetitions_fine[enums::x] = 6 * (parameter.nbr_global_refinements+1);
			 repetitions_fine[enums::y] = 3; // y
			 repetitions_fine[enums::z] = parameter.nbr_elementsInZ;

			Triangulation<3> triangulation_fine, triangulation_coarse;
			// The fine brick
			 GridGenerator::subdivided_hyper_rectangle ( triangulation_fine,
														 repetitions_fine,
														 p1,
														 p2 );

			// The coarse brick
			 GridGenerator::subdivided_hyper_rectangle ( triangulation_coarse,
														 repetitions,
														 p2,
														 p3 );

			// Merging fine and coarse brick
			// @note The interface between the two bricks needs to be meshed identically.
			// deal.II cannot detect hanging nodes there.
			GridGenerator::merge_triangulations( triangulation_fine,
												 triangulation_coarse,
												 triangulation,
												 1e-9 * length );

			// Local refinement
			if ( parameter.nbr_holeEdge_refinements >= 0 && parameter.nbr_global_refinements==0 )
			{
				for (typename Triangulation<dim>::active_cell_iterator
							 cell = triangulation.begin_active();
							 cell != triangulation.end(); ++cell)
				{
					for ( unsigned int face=0; face < GeometryInfo<dim>::faces_per_cell; face++ )
						if ( cell->center()[loading_direction()] < length * refined_fraction )
						{
							cell->set_refine_flag();
							break;
						}
				}
				triangulation.execute_coarsening_and_refinement();
			}
		}
		// use uniform brick with xy refinements
		else if ( parameter.refine_special==enums::Rod_refine_special_uniform )
		{
			 std::vector<unsigned int> repetitions (3);
			 repetitions[enums::x] = std::pow(2.,parameter.nbr_holeEdge_refinements);
			 repetitions[enums::y] = std::pow(2.,parameter.nbr_holeEdge_refinements); // y
			 repetitions[enums::z] = parameter.nbr_elementsInZ;

			 GridGenerator::subdivided_hyper_rectangle ( triangulation,
														 repetitions,
														 p1,
														 p4 );
		}
		else if ( parameter.refine_special==enums::Mesh_refine_special_innermost )
		{
			 std::vector<unsigned int> repetitions (3);
			 repetitions[enums::x] = parameter.grid_y_repetitions*2.;
			 repetitions[enums::y] = parameter.grid_y_repetitions; // y
			 repetitions[enums::z] = parameter.nbr_elementsInZ;

			 GridGenerator::subdivided_hyper_rectangle ( triangulation,
														 repetitions,
														 p1,
														 p4 );
		}

		// Clear boundary ID's
		for (typename Triangulation<dim>::active_cell_iterator
			 cell = triangulation.begin_active();
			 cell != triangulation.end(); ++cell)
		{
			for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
			  if (cell->face(face)->at_boundary())
			  {
				  cell->face(face)->set_all_boundary_ids(0);
			  }
		}

		// Set boundary IDs and and manifolds
		const Point<dim> direction (0,0,1);
		const Point<dim> centre (0,0,0);
		for (typename Triangulation<dim>::active_cell_iterator
			 cell = triangulation.begin_active();
			 cell != triangulation.end(); ++cell)
		{
			for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
			  if (cell->face(face)->at_boundary())
			  {
				// Set boundary IDs
				if (std::abs(cell->face(face)->center()[0] - 0.0) < search_tolerance)
				{
					cell->face(face)->set_boundary_id(enums::id_boundary_xMinus);
				}
				else if ( std::abs(cell->face(face)->center()[enums::x] - body_dimensions[enums::x] ) < search_tolerance)
				{
					cell->face(face)->set_boundary_id(enums::id_boundary_xPlus);
				}
				else if (std::abs(cell->face(face)->center()[1] - 0.0) < search_tolerance)
				{
					cell->face(face)->set_boundary_id(enums::id_boundary_yMinus);
				}
				else if (std::abs(cell->face(face)->center()[enums::y] - body_dimensions[enums::y]) < search_tolerance)
				{
					cell->face(face)->set_boundary_id(enums::id_boundary_yPlus);
				}
				else if (std::abs(cell->face(face)->center()[2] - 0.0) < search_tolerance)
				{
					cell->face(face)->set_boundary_id(enums::id_boundary_zMinus);
				}
				else if (std::abs(cell->face(face)->center()[enums::z] - body_dimensions[enums::z]) < search_tolerance)
				{
					cell->face(face)->set_boundary_id(enums::id_boundary_zPlus);
				}
				else
				{
					// There are just eight faces, so if we missed one, something went clearly terribly wrong
					 AssertThrow(false, ExcMessage("Beam - make_grid 3D<< Found an unidentified face at the boundary. Maybe it slipt through the assignment or that face is simply not needed. So either check the implementation or comment this line in the code"));
				}
			  }
		}

		triangulation.refine_global(parameter.nbr_global_refinements);	// ... Parameter.prm file

		// Refine in a special manner only cells around the origin
		{
			for ( unsigned int nbr_local_ref=0; nbr_local_ref < parameter.nbr_holeEdge_refinements; nbr_local_ref++ )
			{
				for (typename Triangulation<dim>::active_cell_iterator
							 cell = triangulation.begin_active();
							 cell != triangulation.end(); ++cell)
				{
					for (unsigned int vertex=0; vertex < GeometryInfo<dim>::vertices_per_cell; ++vertex)
					{
						// Find all cells that lay in an exemplary damage band with size 1/4 from the y=0 face
						if ( cell->vertex(vertex)[enums::x] < length/4. )
						{
							cell->set_refine_flag();
							break;
						}
					}
				}
				triangulation.execute_coarsening_and_refinement();
			}
		}



		// include the following two scopes to see directly how the variation of the input parameters changes the geometry of the grid
//		{
//			std::ofstream out ("grid-3d_quarter_plate_merged.eps");
//			GridOut grid_out;
//			grid_out.write_eps (triangulation, out);
//			std::cout << "Grid written to grid-3d_quarter_plate_merged.eps" << std::endl;
//			AssertThrow(false,ExcMessage("done"));
//		}
//		{
//			std::ofstream out_ucd("Grid-3d_quarter_plate_merged.inp");
//			GridOut grid_out;
//			GridOutFlags::Ucd ucd_flags(true,true,true);
//			grid_out.set_flags(ucd_flags);
//			grid_out.write_ucd(triangulation, out_ucd);
//			std::cout<<"Mesh written to Grid-3d_quarter_plate_merged.inp "<<std::endl;
//		}
	}