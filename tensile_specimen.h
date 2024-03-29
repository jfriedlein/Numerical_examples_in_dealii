#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>

#include <iostream>
#include <fstream>
#include <cmath>


// Numerical example helper function (required, can be downloaded from https://github.com/jfriedlein/Numerical_examples_in_dealii)
// also contains enumerators as part of "enums::"
#include "./numEx-helper_fnc.h"

using namespace dealii;

namespace tensileSpecimen
{
	// Name of the numerical example
	 std::string numEx_name = "tensileSpecimen";

	// The loading direction: \n
	// In which coordinate direction the load shall be applied, so x/y/z.
	 const unsigned int loading_direction = enums::x;

	// The loaded faces:
	 const enums::enum_boundary_ids id_boundary_load = enums::id_boundary_xPlus;
	 //const enums::enum_boundary_ids id_boundary_secondaryLoad = enums::id_boundary_yPlus;

	// Desired overall length
	// @todo Currently hardcoded, should be a parameter
	 const double desired_length = 70.; // SEP1230: 70, A50: 84

	// Some internal parameters
	 struct parameterCollection
	 {
		const types::boundary_id boundary_id_radius_lower = 10;
		const types::manifold_id manifold_id_radius_lower = 10;
		
		const types::boundary_id boundary_id_radius_upper = 11;
		const types::manifold_id manifold_id_radius_upper = 11;	

		const double search_tolerance = 1e-12;
	 };


	template<int dim>
	void make_constraints ( AffineConstraints<double> &constraints, const FESystem<dim> &fe, DoFHandler<dim> &dof_handler_ref,
							const bool &apply_dirichlet_bc, double &current_load_increment,
							const Parameter::GeneralParameters &parameter )
	{
		/* inputs:
		 * dof_handler_ref,
		 * fe
		 * apply_dirichlet_bc
		 * constraints
		 * current_load_increment
		 */

		// Symmetry constraints:
		// Update and apply new constraints
		//		on x0_plane for symmetry (displacement_in_x = 0)
		//		on y0_plane for symmetry (displacement_in_y = 0)
		//		on z0_plane for symmetry (displacement_in_z = 0)

		parameterCollection parameters_internal;

		const FEValuesExtractors::Vector displacement(0);
		const FEValuesExtractors::Scalar x_displacement(0);
		const FEValuesExtractors::Scalar y_displacement(1);

		// on X0 plane
		// like a symmetry constraint, because the left part is longer (continuous on)

		// clamp the left face
		 numEx::BC_apply_fix( enums::id_boundary_xMinus, dof_handler_ref, fe, constraints );

		// on Y0 edge
		// symmetry constraint (only 1/4 of the model is considered)
		 numEx::BC_apply( enums::id_boundary_yMinus, enums::y, 0, apply_dirichlet_bc, dof_handler_ref, fe, constraints );

		// on Z0 plane
		// symmetry constraint in thickness direction
		if ( dim==3 )
		{
			 numEx::BC_apply( enums::id_boundary_zMinus, enums::z, 0, apply_dirichlet_bc, dof_handler_ref, fe, constraints );
		}

		// BC for the load ...
		 if ( parameter.driver == enums::Dirichlet )  // ... as Dirichlet only for Dirichlet as driver
			numEx::BC_apply( id_boundary_load, loading_direction, current_load_increment, apply_dirichlet_bc, dof_handler_ref, fe, constraints );
	}


	// grid generation for 2D and 3D
	template <int dim>
	void make_2D_grid( Triangulation<dim> &triangulation_2D, const Parameter::GeneralParameters &parameter )
	{
		parameterCollection parameters_internal;

		/*
		 * parameters:
		 * * transition_radius
		 * * hwidth_b
		 * * hwidth_B
		 * * hthickness
		 * * length_parallel
		 * * extension_lower
		 * * extension_upper
		 * * n_grid_repetitions_central_part
		 * * n_global_refinements
		 */

		const double search_tolerance = parameters_internal.search_tolerance;

		// radius of the transition between width b and width B
		 const double transition_radius = parameter.holeRadius;

		// half widths for the quarter model (split in y and z direction)
		 const double hwidth_b = parameter.notchWidth/2.;
		 const double hwidth_B = parameter.width/2.;

		// half thickness of the specimen
		 //const double hthickness = parameter.thickness/2.;

		// length of the parallel thinner sectin with width b
		 const double length_parallel = parameter.height;

		// possiblity to extend the model to the left and right by these additional lengths
		 const double extension_lower = 0.;
		 const double extension_upper = 0.;

		AssertThrow( hwidth_B > hwidth_b + ( 1 - std::sin(22.5 * std::atan(1)*4/180)) * transition_radius,
					 ExcMessage("tensileSpecimen<< We cannot yet model meshes with a very soft radius. "
							    "This would require us to delete some of the uppermost cells and shift "
							    "significantly more vertices."));

	  // @note
	  // The mesh is first created in 2D and then extruded in 3D before some additional modifcations are applied.


	  // ************************************************************************************************************
	  // Central brick (currently still 2D rectangle)
	  // The central thinner parallel section of the specimen with width b and length \a length_parallel
		// The points spanning the brick
		 Point<2> p1 (-length_parallel/2.,0);
		 Point<2> p2 (length_parallel/2., hwidth_b);

		// Vector containing the number of elements in each dimension
		 std::vector<unsigned int> repetitions (2);
		 // @todo Rename \a grid_y_repetitions
		 // @todo Maybe we want to align the specimen in y-direction, so we would have to change the enums::x .. accessors accordingly
		 repetitions[enums::x] = parameter.grid_y_repetitions;
		 repetitions[enums::y] = 2; // minimum for compatibility with prerefined radial parts

		Triangulation<2> central_rectangular_part;
		GridGenerator::subdivided_hyper_rectangle
				( 	central_rectangular_part,
					repetitions,
					p1,
					p2 );
	  // ************************************************************************************************************

	  // ************************************************************************************************************
	  // radial part lower
	  // To get the radius, we cut a plate with hole into pieces
		   // Center point of the lower radius
			Point<2> lower_radius_center;
		   // create the quarter plate with hole
			Triangulation<2> tria_quarter_plate_hole_lower;
			{
			   // Set the center point
				lower_radius_center[enums::x] = - length_parallel/2.;
				lower_radius_center[enums::y] = hwidth_b + transition_radius; // this makes it tangential to the rectangular part in the middle

			   // Create the full plate with hole
				Triangulation<2> tria_plate_hole;
				GridGenerator::plate_with_a_hole
					(
						tria_plate_hole,
						transition_radius,
						transition_radius + hwidth_b,
						0,
						0,
						extension_lower, // possible extension for clamping
						0,
						lower_radius_center
					);

				// Adapt the x-length of the radius to get the defined overall length
				// We do this right here before we refine the mesh to avoid the distortion of inner cells,
				// which don't exist before the following refine_global.
				// The following is very simple, but
				// @todo This only works for minor changes in the length, for larger we need the pad for \a plate_with_a_hole, or
				// abort if the desired length gets to small
				 {
					double initial_pos = - transition_radius - hwidth_b - length_parallel/2.;
					double new_pos = - desired_length/2.;

					numEx::shift_vertex_layer( tria_plate_hole, initial_pos, new_pos, enums::x );
				 }

			   // Prerefine the plate
			   // This is necessary at this very moment. When we later apply our own ids and manifolds and there is
			   // only a single element between the outer straight edge of the plate and the hole edge, the cylindrical
			   // manifold is applied on both edges. So in the end we would have a plate with inner and outer cylindrical,
			   // which is not what I want.
			   // @todo-optimize: Nevertheless, all of that makes no sense. There must be an easier way.
				tria_plate_hole.refine_global(1);

			   // Remove 3/4 of the plate, leaving only the desired quadrant
				std::set<typename Triangulation<2>::active_cell_iterator > cells_to_remove;
				for (typename Triangulation<2>::active_cell_iterator
					 cell = tria_plate_hole.begin_active();
					 cell != tria_plate_hole.end(); ++cell)
				{
					// Remove all cells that are not in the third quadrant
					if (cell->center()[enums::x] > lower_radius_center[enums::x] || cell->center()[enums::y] > lower_radius_center[enums::y] )
						cells_to_remove.insert(cell);
				}

				Assert(cells_to_remove.size() > 0, ExcInternalError());
				Assert(cells_to_remove.size() != tria_plate_hole.n_active_cells(), ExcInternalError());
				GridGenerator::create_triangulation_with_removed_cells(tria_plate_hole,cells_to_remove,tria_quarter_plate_hole_lower);
			}

			// Adapt the width B (y-dimension) of the lower radial part
			// By moving all the vertices that are too high (\a initial_pos) downwards to the new position \a new_pos that
			// was defined as the width B.
			 {
				double initial_pos = hwidth_b + transition_radius;
				double new_pos = hwidth_B;
				const unsigned int direction = enums::y;
				{
					bool shifted_node = false;
					for (typename Triangulation<2>::active_cell_iterator
					   cell = tria_quarter_plate_hole_lower.begin_active();
					   cell != tria_quarter_plate_hole_lower.end(); ++cell)
					{
					  for (unsigned int vertex=0; vertex<GeometryInfo<2>::vertices_per_cell; ++vertex )
					  {
						  if ( std::abs( cell->vertex(vertex)[direction] - initial_pos) < 1e3 * search_tolerance /*requires this softer tolerance*/ )
						  {
							  Point<2> shift_vector;
							  // If this point is on the radius, we have to move it tangential (so down and to the RIGHT)
							   if ( std::abs( cell->vertex(vertex).distance(lower_radius_center) - transition_radius) < 1e3 * search_tolerance )
							   {
								  shift_vector[enums::y] = - (hwidth_b + transition_radius - hwidth_B);
								  shift_vector[enums::x] = transition_radius
														   - std::sqrt( transition_radius*transition_radius - shift_vector[enums::y]*shift_vector[enums::y] );
							   }
							  // If not, then we simply move it down by the y-difference
							   else
							   {
								  shift_vector[direction] = (new_pos-initial_pos);
							   }
							  cell->vertex(vertex) += shift_vector;
							  shifted_node = true; // -> We have shifted at least a single node
						  }
					  }
					}
					// Ensure that we shifted at least a single node
					 AssertThrow( shifted_node==true,
								  ExcMessage("tensileSpecimen<< You haven't moved a single node. Please check the selection criterion initial_pos."));
				}
			 }
		// ************************************************************************************************************

		// ************************************************************************************************************
		// radial part upper
		// See above for the similar lower part. Still, be aware of the minor but critical differences.
		 Point<2> upper_radius_center;
		 Triangulation<2> tria_quarter_plate_hole_upper;
		 {
			upper_radius_center[enums::x] = + length_parallel/2.; // center point is now on the RIGHT-hand side of the origin
			upper_radius_center[enums::y] = hwidth_b + transition_radius;

			Triangulation<2> tria_plate_hole;
			GridGenerator::plate_with_a_hole
				(
					tria_plate_hole,
					transition_radius,
					transition_radius + hwidth_b,
					0,
					0,
					0,
					extension_upper, // possible extension for clamping to the RIGHT
					upper_radius_center
				);

			// Adapt the x-length of the radius to get the defined overall length
			// We do this right here before we refine the mesh to avoid the distortion of inner cells,
			// which don't exist before the following refine_global.
			 {
				double initial_pos = transition_radius + hwidth_b + length_parallel/2.;
				double new_pos = desired_length/2.;

				numEx::shift_vertex_layer( tria_plate_hole, initial_pos, new_pos, enums::x );
			 }

			tria_plate_hole.refine_global(1);

			std::set<typename Triangulation<2>::active_cell_iterator > cells_to_remove;
			for (typename Triangulation<2>::active_cell_iterator
				 cell = tria_plate_hole.begin_active();
				 cell != tria_plate_hole.end(); ++cell)
			{
				// Remove all cells that are not in the FOURTH quadrant
				if (cell->center()[enums::x] < upper_radius_center[enums::x] || cell->center()[enums::y] > upper_radius_center[enums::y] )
					cells_to_remove.insert(cell);
			}

			Assert(cells_to_remove.size() > 0, ExcInternalError());
			Assert(cells_to_remove.size() != tria_plate_hole.n_active_cells(), ExcInternalError());
			GridGenerator::create_triangulation_with_removed_cells(tria_plate_hole,cells_to_remove,tria_quarter_plate_hole_upper);
		 }

		 // Adapt the width B of the upper radial part
		  {
			double initial_pos = hwidth_b + transition_radius;
			double new_pos = hwidth_B;
			const unsigned int direction = enums::y;
			{
				bool shifted_node = false;
				for (typename Triangulation<2>::active_cell_iterator
				   cell = tria_quarter_plate_hole_upper.begin_active();
				   cell != tria_quarter_plate_hole_upper.end(); ++cell)
				{
				  for (unsigned int vertex=0; vertex<GeometryInfo<2>::vertices_per_cell; ++vertex )
				  {
					  if ( std::abs( cell->vertex(vertex)[direction] - initial_pos) < (1e3 * search_tolerance) /*requires this softer tolerance*/ )
					  {
						  Point<2> shift_vector;
						  if ( std::abs( cell->vertex(vertex).distance(upper_radius_center) - transition_radius) < (1e3 * search_tolerance) )
						  {
							  // Move it tho the bottom LEFT
							  shift_vector[enums::y] = - (hwidth_b + transition_radius - hwidth_B);
							  shift_vector[enums::x] = transition_radius
													   - std::sqrt( transition_radius*transition_radius - shift_vector[enums::y]*shift_vector[enums::y] );
							  shift_vector[enums::x] *= -1; // for the upper radius we need to move the shifted node to the left when moving it downwards
						  }
						  else
						  {
							  shift_vector[direction] = (new_pos-initial_pos);
						  }
						  cell->vertex(vertex) += shift_vector;
						  shifted_node = true; // -> We have shifted at least a single node
					  }
				  }
				}
				// Ensure that we shifted at least a single node
				 AssertThrow( shifted_node==true, ExcMessage("tensileSpecimen<< You haven't moved a single node. Please check the selection criterion initial_pos."));
			}
		  }
		// ************************************************************************************************************


		// ************************************************************************************************************
		// Merge the three mesh parts
		// Here we also use a softer tolerance (1e-8), because else it will not be properly merged
			GridGenerator::merge_triangulations( tria_quarter_plate_hole_lower,
												 central_rectangular_part,
												 triangulation_2D, 1e3 * search_tolerance );
			GridGenerator::merge_triangulations( triangulation_2D,
												 tria_quarter_plate_hole_upper,
												 triangulation_2D, 1e3 * search_tolerance );

 			// Clear boundary ID's
 			 for (typename Triangulation<2>::active_cell_iterator
 				 cell = triangulation_2D.begin_active();
 				 cell != triangulation_2D.end(); ++cell)
 			 {
 				for (unsigned int face=0; face < GeometryInfo<2>::faces_per_cell; ++face)
 				  if (cell->face(face)->at_boundary())
 				  {
 					  cell->face(face)->set_all_boundary_ids(0);
 				  }
 			 }
		// ************************************************************************************************************
	}


	// 2D grid
	template <int dim>
	void make_grid( Triangulation<2> &triangulation, const Parameter::GeneralParameters &parameter )
	{
		make_2D_grid( triangulation, parameter );

		// @todo-optimize Somehow avoid this doubling of parameters
		parameterCollection parameters_internal;

		const double search_tolerance = parameters_internal.search_tolerance;

		// radius of the transition between width b and width B
		 const double transition_radius = parameter.holeRadius;

		// half widths for the quarter model (split in y and z direction)
		 const double hwidth_b = parameter.notchWidth/2.;
		 //const double hwidth_B = parameter.width/2.;

		// half thickness of the specimen
		 //const double hthickness = parameter.thickness/2.;

		// length of the parallel thinner sectin with width b
		 const double length_parallel = parameter.height;

		// possiblity to extend the model to the left and right by these additional lengths
		 //const double extension_lower = 0.;
		 //const double extension_upper = 0.;

		// Set the lower center point
		 Point<2> lower_radius_center;
		 lower_radius_center[enums::x] = - length_parallel/2.;
		 lower_radius_center[enums::y] = hwidth_b + transition_radius; // this makes it tangential to the rectangular part in the middle

		// Set the upper center point
		 Point<2> upper_radius_center;
		 upper_radius_center[enums::x] = - lower_radius_center[enums::x];
		 upper_radius_center[enums::y] = lower_radius_center[enums::y]; // this makes it tangential to the rectangular part in the middle

		// Compute the x-coordinate of the left- and rightmost vertices
		 //double lower_end = - ( length_parallel/2. + transition_radius + hwidth_b + extension_lower );
		 double lower_end = - desired_length/2.;
		 //double upper_end = + ( length_parallel/2. + transition_radius + hwidth_b + extension_upper );
		 double upper_end = desired_length/2.;

		// Set boundary IDs and and manifolds
		for (typename Triangulation<dim>::active_cell_iterator
			 cell = triangulation.begin_active();
			 cell != triangulation.end(); ++cell)
		{
			for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
			  if (cell->face(face)->at_boundary())
			  {
				// Set boundary IDs
				// lower end
				if (std::abs(cell->face(face)->center()[enums::x] - lower_end) < search_tolerance)
				{
					cell->face(face)->set_boundary_id(enums::id_boundary_xMinus);
				}
				// upper end
				else if (std::abs(cell->face(face)->center()[enums::x] - upper_end) < search_tolerance)
				{
					cell->face(face)->set_boundary_id(enums::id_boundary_xPlus);
				}
				// y0 (symmetry)
				else if (std::abs(cell->face(face)->center()[enums::y] - 0.0) < search_tolerance)
				{
					cell->face(face)->set_boundary_id(enums::id_boundary_yMinus);
				}
				else
				{
					for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_face; ++vertex)
					{
					 // Project the cell vertex to the XY plane and test the distance from the cylinder axis
						Point<dim> current_vertex = cell->face(face)->vertex(vertex);
						// If the point lies above the inner parallel area ...
						if ( current_vertex[enums::x] > ( length_parallel/2. + 5.*search_tolerance ) )
						{
							// ... and the radius fits to transition radius
							if ( std::abs(current_vertex.distance(upper_radius_center) - transition_radius) < search_tolerance )
							{
								cell->face(face)->set_boundary_id(parameters_internal.boundary_id_radius_upper);
								cell->face(face)->set_all_manifold_ids(parameters_internal.manifold_id_radius_upper);
								break;
							}
						}
						// If the point lies below the inner parallel area ...
						else if ( current_vertex[enums::x] < - ( length_parallel/2. + 5.*search_tolerance ) )
						{
							// ... and the radius fits to transition radius
							if ( std::abs(current_vertex.distance(lower_radius_center) - transition_radius) < search_tolerance )
							{
								cell->face(face)->set_boundary_id(parameters_internal.boundary_id_radius_lower);
								cell->face(face)->set_all_manifold_ids(parameters_internal.manifold_id_radius_lower);
								break;
							}
						}
					}
				}
			  }
		}

		// Apply cylindrical manifolds to both radii
		 SphericalManifold<dim> spherical_manifold_upper ( upper_radius_center);
		 triangulation.set_manifold(parameters_internal.manifold_id_radius_upper, spherical_manifold_upper);

		 SphericalManifold<dim> spherical_manifold_lower ( lower_radius_center);
		 triangulation.set_manifold(parameters_internal.manifold_id_radius_lower, spherical_manifold_lower);

		// Notch the parallel area in the middle
		 // Find the nodes (plural because of thickness) at x=0 and shift them down by 0.5% of the hwidth_b
		 if ( true/*notch the tensile specimen*/ )
		 {
			for (typename Triangulation<dim>::active_cell_iterator
				 cell = triangulation.begin_active();
				 cell != triangulation.end(); ++cell)
			{
					for (unsigned int vertex=0; vertex < GeometryInfo<dim>::vertices_per_cell; ++vertex)
						if ( ( std::abs(cell->vertex(vertex)[enums::x]) < 5.*search_tolerance ) )
							if (( std::abs(cell->vertex(vertex)[enums::y] - hwidth_b) < 5.*search_tolerance ) )
							  cell->vertex(vertex)[enums::y] -= 0.03 * hwidth_b;
			}
		 }

		// Global refinement
		 triangulation.refine_global(parameter.nbr_global_refinements);

		// Refine the cells in the parallel part
		 for ( unsigned int nbr_local_ref=0; nbr_local_ref < parameter.nbr_holeEdge_refinements; nbr_local_ref++ )
		 {
			for (typename Triangulation<dim>::active_cell_iterator
						 cell = triangulation.begin_active();
						 cell != triangulation.end(); ++cell)
			{
				for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
					// Find all cells that lay in an exemplary damage band with size 1.5 mm from the y=0 face
					if ( std::abs( cell->face(face)->center()[enums::x] ) <= ( length_parallel/(4.+3.*double(nbr_local_ref) ) ) )
					{
						if ( nbr_local_ref==1 || nbr_local_ref==3 || nbr_local_ref==5 || nbr_local_ref==7 ) // even
							cell->set_refine_flag();
						else
							cell->set_refine_flag(RefinementCase<dim>::cut_x); // refine only in the x-direction
						break;
					}
			}
			triangulation.execute_coarsening_and_refinement();
		 }

		// include the following two scopes to see directly how the variation of the input parameters changes the geometry of the grid
		/*
		{
			std::ofstream out ("grid-2d_quarter_plate_merged.eps");
			GridOut grid_out;
			GridOutFlags::Eps<2> eps_flags;
			eps_flags.line_width = 0.1;
			grid_out.set_flags (eps_flags);
			grid_out.write_eps (triangulation, out);
			std::cout << "Grid written to grid-2d_quarter_plate_merged.eps" << std::endl;
			std::cout << "nElem: " << triangulation.n_active_cells() << std::endl;
			AssertThrow(false,ExcMessage("ddd"));
		}

		{
			std::ofstream out_ucd("Grid-2d_quarter_plate_merged.inp");
			GridOut grid_out;
			GridOutFlags::Ucd ucd_flags(true,true,true);
			grid_out.set_flags(ucd_flags);
			grid_out.write_ucd(triangulation, out_ucd);
			std::cout<<"Mesh written to Grid-2d_quarter_plate_merged.inp "<<std::endl;
		}
		*/
	}



// 3d grid
	template <int dim>
	void make_grid( Triangulation<3> &triangulation, const Parameter::GeneralParameters &parameter )
	{
		parameterCollection parameters_internal;
		
		/*
		 * parameters:
		 * * transition_radius
		 * * hwidth_b
		 * * hwidth_B
		 * * hthickness
		 * * length_parallel
		 * * extension_lower
		 * * extension_upper
		 * * n_grid_repetitions_central_part
		 * * n_global_refinements
		 */
		
		const double search_tolerance = parameters_internal.search_tolerance;

		// radius of the transition between width b and width B
		 const double transition_radius = parameter.holeRadius;
		
		// half widths for the quarter model (split in y and z direction)
		 const double hwidth_b = parameter.notchWidth/2.;
		 const double hwidth_B = parameter.width/2.;
		
		// half thickness of the specimen
		 const double hthickness = parameter.thickness/2.;
		
		// length of the parallel thinner sectin with width b
		 const double length_parallel = parameter.height;
		 
		// possiblity to extend the model to the left and right by these additional lengths
		 //const double extension_lower = 0.;
		 //const double extension_upper = 0.;
		
		AssertThrow( hwidth_B > hwidth_b + ( 1 - std::sin(22.5 * std::atan(1)*4/180)) * transition_radius,
					 ExcMessage("tensileSpecimen<< We cannot yet model meshes with a very soft radius. "
							    "This would require us to delete some of the uppermost cells and shift "
							    "significantly more vertices."));

		Triangulation<2> triangulation_2D;
		make_2D_grid(triangulation_2D, parameter);

		// Set the lower center point
		 Point<2> lower_radius_center;
		 lower_radius_center[enums::x] = - length_parallel/2.;
		 lower_radius_center[enums::y] = hwidth_b + transition_radius; // this makes it tangential to the rectangular part in the middle

		// Set the upper center point
		 Point<2> upper_radius_center;
		 upper_radius_center[enums::x] = - lower_radius_center[enums::x];
		 upper_radius_center[enums::y] = lower_radius_center[enums::y]; // this makes it tangential to the rectangular part in the middle

		// ************************************************************************************************************		
		// Extrude 2D grid to 3D
		 GridGenerator::extrude_triangulation( triangulation_2D,
										   	   parameter.nbr_elementsInZ + 1,
											   hthickness,
											   triangulation );
		
		// ************************************************************************************************************	
	    // From now on 3D
		// Compute the x-coordinate of the left- and rightmost vertices
		 //double lower_end = - ( length_parallel/2. + transition_radius + hwidth_b + extension_lower );
		 double lower_end = - desired_length/2.;
		 //double upper_end = + ( length_parallel/2. + transition_radius + hwidth_b + extension_upper );
		 double upper_end = desired_length/2.;

		// Expand the center points to 3D
		 Point<3> upper_radius_center_3D;
		 Point<3> lower_radius_center_3D;
		 for ( unsigned int i=0; i<2; i++ )
		 {
			upper_radius_center_3D[i] = upper_radius_center[i];
			lower_radius_center_3D[i] = lower_radius_center[i];
		 }
		
		// Set boundary IDs and and manifolds
		for (typename Triangulation<dim>::active_cell_iterator
			 cell = triangulation.begin_active();
			 cell != triangulation.end(); ++cell)
		{
			for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
			  if (cell->face(face)->at_boundary())
			  {
				// Set boundary IDs
				// lower end
				if (std::abs(cell->face(face)->center()[enums::x] - lower_end) < search_tolerance)
				{
					cell->face(face)->set_boundary_id(enums::id_boundary_xMinus);
				}
				// upper end
				else if (std::abs(cell->face(face)->center()[enums::x] - upper_end) < search_tolerance)
				{
					cell->face(face)->set_boundary_id(enums::id_boundary_xPlus);
				}
				// y0 (symmetry)
				else if (std::abs(cell->face(face)->center()[enums::y] - 0.0) < search_tolerance)
				{
					cell->face(face)->set_boundary_id(enums::id_boundary_yMinus);
				}
				else if (std::abs(cell->face(face)->center()[2] - 0.0) < search_tolerance)
				{
					cell->face(face)->set_boundary_id(enums::id_boundary_zMinus);
				}
				else if (std::abs(cell->face(face)->center()[2] - hthickness) < search_tolerance)
				{
					cell->face(face)->set_boundary_id(enums::id_boundary_zPlus);
				}
				else
				{
					for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_face; ++vertex)
					{
					 // Project the cell vertex to the XY plane and test the distance from the cylinder axis
						Point<dim> vertex_proj = cell->face(face)->vertex(vertex);
						vertex_proj[enums::z] = 0.0;
						// If the point lies above the inner parallel area ...
						if ( vertex_proj[enums::x] > ( length_parallel/2. + 5.*search_tolerance ) )
						{
							// ... and the radius fits to transition radius
							if ( std::abs(vertex_proj.distance(upper_radius_center_3D) - transition_radius) < search_tolerance )
							{
								cell->face(face)->set_boundary_id(parameters_internal.boundary_id_radius_upper);
								cell->face(face)->set_all_manifold_ids(parameters_internal.manifold_id_radius_upper);
								break;
							}
						}
						// If the point lies below the inner parallel area ...
						else if ( vertex_proj[enums::x] < - ( length_parallel/2. + 5.*search_tolerance ) )
						{
							// ... and the radius fits to transition radius
							if ( std::abs(vertex_proj.distance(lower_radius_center_3D) - transition_radius) < search_tolerance )
							{
								cell->face(face)->set_boundary_id(parameters_internal.boundary_id_radius_lower);
								cell->face(face)->set_all_manifold_ids(parameters_internal.manifold_id_radius_lower);
								break;
							}
						}
					}
				}
			  }
		}

		// Apply cylindrical manifolds to both radii
		 Tensor<1,dim> cylinder_axis;
		 cylinder_axis[enums::z] = 1;

		 CylindricalManifold<dim> cylindrical_manifold_upper (cylinder_axis, upper_radius_center_3D);
		 triangulation.set_manifold(parameters_internal.manifold_id_radius_upper, cylindrical_manifold_upper);

		 CylindricalManifold<dim> cylindrical_manifold_lower (cylinder_axis, lower_radius_center_3D);
		 triangulation.set_manifold(parameters_internal.manifold_id_radius_lower, cylindrical_manifold_lower);

		 // "Global" refinements are done solely in the xy-plane
		  for ( unsigned int nbr_gl_ref=0; nbr_gl_ref < parameter.nbr_global_refinements; nbr_gl_ref++ )
		  {
			for (typename Triangulation<dim>::active_cell_iterator
						 cell = triangulation.begin_active();
						 cell != triangulation.end(); ++cell)
			{
				cell->set_refine_flag(RefinementCase<dim>::cut_xy); // refine in x and y-direction
			}
			triangulation.execute_coarsening_and_refinement();
		  }
		 //triangulation.refine_global(parameter.nbr_global_refinements);

		 // @todo check the use of only anisotropic xy refinements to keep the thickness direction
		// Refine the cells in the parallel part
		// Once refine by cut_x
		 for (typename Triangulation<dim>::active_cell_iterator
					 cell = triangulation.begin_active();
					 cell != triangulation.end(); ++cell)
		 {
			for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
				// Find all cells that lay in an exemplary damage band with size 1.5 mm from the y=0 face
				if ( std::abs( cell->face(face)->center()[enums::x] ) <= parameter.referenceLength/2. )//length_parallel/(4.+2.*double(nbr_local_ref)) )
				{
					// @todo Multiple local anisotropic refinements cause DII to fail, Why?
//						if ( nbr_local_ref==1 || nbr_local_ref==3 || nbr_local_ref==5 || nbr_local_ref==7 ) // even
						cell->set_refine_flag(RefinementCase<dim>::cut_x); // refine in x and y-direction
//						else
//							cell->set_refine_flag(RefinementCase<dim>::cut_x); // refine only in the x-direction
					break;
				}
		 }
		 triangulation.execute_coarsening_and_refinement();

		// Refine innermost part by cut_xy
		 for ( unsigned int nbr_local_ref=0; nbr_local_ref < parameter.nbr_holeEdge_refinements; nbr_local_ref++ )
		 {
			for (typename Triangulation<dim>::active_cell_iterator
						 cell = triangulation.begin_active();
						 cell != triangulation.end(); ++cell)
			{
				for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
					// Find all cells that lay in an exemplary damage band with size 1.5 mm from the y=0 face
					if ( std::abs( cell->face(face)->center()[enums::x] ) <= 2. )//length_parallel/(4.+2.*double(nbr_local_ref)) )
					{
						// @todo Multiple local anisotropic refinements cause DII to fail, Why?
//						if ( nbr_local_ref==1 || nbr_local_ref==3 || nbr_local_ref==5 || nbr_local_ref==7 ) // even
							cell->set_refine_flag(RefinementCase<dim>::cut_x); // refine in x and y-direction
//						else
//							cell->set_refine_flag(RefinementCase<dim>::cut_x); // refine only in the x-direction
						break;
					}
			}
			triangulation.execute_coarsening_and_refinement();
		 }

		// Notch the parallel area in the middle
		 // Find the nodes (plural because of thickness) at x=0 and shift them down by 0.5% of the hwidth_b
		 if ( false/*notch the tensile specimen*/ )
		 {
			for (typename Triangulation<dim>::active_cell_iterator
				 cell = triangulation.begin_active();
				 cell != triangulation.end(); ++cell)
			{
					for (unsigned int vertex=0; vertex < GeometryInfo<dim>::vertices_per_cell; ++vertex)
						if ( ( std::abs(cell->vertex(vertex)[enums::x]) < 5.*search_tolerance ) )
							if (( std::abs(cell->vertex(vertex)[enums::y] - hwidth_b) < 5.*search_tolerance ) )
							  cell->vertex(vertex)[enums::y] -= 0.03 * hwidth_b;
			}
		 }

		// Output the triangulation as eps or inp
		 //numEx::output_triangulation( triangulation, enums::output_eps, numEx_name );
	}
}

