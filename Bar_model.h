#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>

#include <iostream>
#include <fstream>
#include <cmath>

#include "../MA-Code/enumerator_list.h"

using namespace dealii;

namespace BarModel
/*
 * A bar with three symmetry constraints, loaded in y-direction
 *
 * CERTIFIED TO STANDARD numExS07 (200724)
 */
{
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

	// USER PARAMETERS
	 const bool constrain_sideways_sliding_of_loaded_face = false;

	// Loading type and required modifications
	// @note "compression": \n
	// We use different boundary conditions and notch the body in the middle of its length and not a y=0.
	// @note "tension" or "standard": \n
	// We model 1/8 of the entire body and notch the body at y=0 (equals the middle of the entire body).
	 const enums::enum_loading_type loading_type = enums::compression;

	template<int dim>
	void make_constraints ( AffineConstraints<double> &constraints, const FESystem<dim> &fe, unsigned int &n_components, DoFHandler<dim> &dof_handler_ref,
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
		if ( loading_type==enums::compression )
		{
			// For compression we cannot apply a symmetry constraint on the X0 plane
		}
		else
		{
			if (apply_dirichlet_bc == true )
			{
				VectorTools::interpolate_boundary_values(
															dof_handler_ref,
															enums::id_boundary_xMinus,
															ZeroFunction<dim> (n_components),
															constraints,
															fe.component_mask(x_displacement)
														);
			}
			else	// in the exact same manner
			{
				VectorTools::interpolate_boundary_values(
															dof_handler_ref,
															enums::id_boundary_xMinus,
															ZeroFunction<dim> (n_components),
															constraints,
															fe.component_mask(x_displacement)
														);
			}
			}

		// on Y0 edge
		if ( loading_type==enums::compression )
		{
			// For compression we fix/clamp the Y0 plane, so it does not run away
			FEValuesExtractors::Vector displacements(0);
			ComponentMask disp_mask = fe.component_mask (displacements);
			if (apply_dirichlet_bc == true )
			{
				VectorTools::interpolate_boundary_values(
															dof_handler_ref,
															enums::id_boundary_yMinus,
															ZeroFunction<dim> (n_components),
															constraints,
															disp_mask // all disp components (not the global damage)
														);
			}
			else	// in the exact same manner
			{
				VectorTools::interpolate_boundary_values(
															dof_handler_ref,
															enums::id_boundary_yMinus,
															ZeroFunction<dim> (n_components),
															constraints,
															disp_mask // all disp components (not the global damage)
														);
			}
		}
		else
		{
			if (apply_dirichlet_bc == true )
			{
				VectorTools::interpolate_boundary_values(
															dof_handler_ref,
															enums::id_boundary_yMinus,
															ZeroFunction<dim> (n_components),
															constraints,
															fe.component_mask(y_displacement)
														);
			}
			else	// in the exact same manner
			{
				VectorTools::interpolate_boundary_values(
															dof_handler_ref,
															enums::id_boundary_yMinus,
															ZeroFunction<dim> (n_components),
															constraints,
															fe.component_mask(y_displacement)
														);
			}
		}

		// on Z0 plane
		if ( loading_type==enums::compression )
		{
			// for compression we don't fix anything in the third direction, because y0 was already clamped.
			// @todo However, what about the upper part?
			if ( dim==3 )
				AssertThrow(false, ExcMessage("Bar_model<< Compression of bar yet untested for 3D. Take care of the z-boundary conditions."));
		}
		else
		{
			if ( dim==3 )
			{
				const FEValuesExtractors::Scalar z_displacement(2);

				if (apply_dirichlet_bc == true )
				{
					VectorTools::interpolate_boundary_values(
																dof_handler_ref,
																enums::id_boundary_zMinus,
																ZeroFunction<dim> (n_components),
																constraints,
																fe.component_mask(z_displacement)
															);
				}
				else	// in the exact same manner
				{
					VectorTools::interpolate_boundary_values(
																dof_handler_ref,
																enums::id_boundary_zMinus,
																ZeroFunction<dim> (n_components),
																constraints,
																fe.component_mask(z_displacement)
															);
				}

				if ( false/*apply sym BC on positive z-face also*/ )
				{
					if (apply_dirichlet_bc == true )
					{
						VectorTools::interpolate_boundary_values(
																	dof_handler_ref,
																	enums::id_boundary_zPlus,
																	ZeroFunction<dim> (n_components),
																	constraints,
																	fe.component_mask(z_displacement)
																);
					}
					else	// in the exact same manner
					{
						VectorTools::interpolate_boundary_values(
																	dof_handler_ref,
																	enums::id_boundary_zPlus,
																	ZeroFunction<dim> (n_components),
																	constraints,
																	fe.component_mask(z_displacement)
																);
					}
				}
			}
		}

		if ( parameter.driver == 2/*Dirichlet*/ ) // ToDo-optimize: use string in parameterfile denoting "Dirichlet" so the enumerator is not undermined
		{
			// on top edge
			if (apply_dirichlet_bc == true )
			{
				VectorTools::interpolate_boundary_values(
															dof_handler_ref,
															id_boundary_load,
															ConstantFunction<dim> (current_load_increment/*add only the increment*/, n_components),
															constraints,
															fe.component_mask(y_displacement)
														);
				if ( constrain_sideways_sliding_of_loaded_face )
					VectorTools::interpolate_boundary_values(
																dof_handler_ref,
																id_boundary_load,
																ZeroFunction<dim> (n_components),
																constraints,
																fe.component_mask(x_displacement)
															);
			}
			else
			{
				VectorTools::interpolate_boundary_values(
															dof_handler_ref,
															id_boundary_load,
															ZeroFunction<dim> (n_components),
															constraints,
															fe.component_mask(y_displacement)
														);
				if ( constrain_sideways_sliding_of_loaded_face )
					VectorTools::interpolate_boundary_values(
																dof_handler_ref,
																id_boundary_load,
																ZeroFunction<dim> (n_components),
																constraints,
																fe.component_mask(x_displacement)
															);
			}
		}
	}


	// 2D grid
		template <int dim>
		void make_grid( Triangulation<2> &triangulation, const Parameter::GeneralParameters &parameter )
		{
			// ToDo-assure: the use the values from the parameter file
			const double width = parameter.width; // use thickness=width for square bottom area
			const double length = parameter.height;

			const bool damage_trigger_by_notching = true;
			const enums::enum_coord notched_face = enums::x;
			const double notch_reduction = parameter.ratio_x;
			// The notch length is set (for consistency) s.t. its mesh discretisation is exact for 1 local refinement (start)
			 double notch_length = parameter.notchWidth;//length/10.;//(2.*parameter.grid_y_repetitions);

			 const bool notch_rounded = true;

			body_dimensions[enums::x] = width;
			body_dimensions[enums::y] = length;

			parameterCollection parameters_internal;

			const double search_tolerance = parameters_internal.search_tolerance;


			// The points that span the brick
			 Point<dim> p1 (0,0);
			 Point<dim> p2 (width, length); // extends in y-direction its length (loaded in y-direction as the other models)

			// vector containing the number of elements in each dimension
			 std::vector<unsigned int> repetitions (3);
			 repetitions[0]=1; // x
			 repetitions[1]=parameter.grid_y_repetitions; // y

			GridGenerator::subdivided_hyper_rectangle 	( 	triangulation,
															repetitions,
															p1,
															p2
														);

			// ToDo-optimize: The following is similar for 2D and 3D, maybe merge it
			//Clear boundary ID's
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

			//Set boundary IDs and and manifolds
			for (typename Triangulation<dim>::active_cell_iterator
				 cell = triangulation.begin_active();
				 cell != triangulation.end(); ++cell)
			{
				for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
				  if (cell->face(face)->at_boundary())
				  {
					//Set boundary IDs
					if (std::abs(cell->face(face)->center()[0] - 0.0) < search_tolerance )//&& (false || cell->face(face)->center()[enums::y] < 0.5) ) // option: not symmetric in x
					{
						cell->face(face)->set_boundary_id(enums::id_boundary_xMinus);
					}
					else if (std::abs(cell->face(face)->center()[0] - width) < search_tolerance)
					{
						cell->face(face)->set_boundary_id(enums::id_boundary_xPlus);
					}
					else if (std::abs(cell->face(face)->center()[1] - 0.0) < search_tolerance)
					{
						cell->face(face)->set_boundary_id(enums::id_boundary_yMinus);
//						cell->set_material_id( enums::tracked_QP );
					}
					else if (std::abs(cell->face(face)->center()[1] - length) < search_tolerance)
					{
						cell->face(face)->set_boundary_id(enums::id_boundary_yPlus);
					}
					else
					{
						// There are just eight faces, so if we missed one, something went clearly terribly wrong
						 AssertThrow(false, ExcMessage("BarModel - make_grid 2D<< Found an unidentified face at the boundary. Maybe it slipt through the assignment or that face is simply not needed. So either check the implementation or comment this line in the code"));
					}
				  }
			}

			if ( notch_rounded/*use round notch*/ && loading_type == enums::compression )
				triangulation.refine_global(2);	// ... Parameter.prm file
			else
				triangulation.refine_global(parameter.nbr_global_refinements);	// ... Parameter.prm file

			// @warning When implementing the compression example we moved the notching in front of the local refinements.
			// @todo Does that do something?
			// Notch the specimen by moving some nodes inwards to form a notch
			if ( damage_trigger_by_notching )
			{
				// Declare the shift vector for the notching
				 Point<2> notching; // initially zero
				// Depending on the desired notching direction (notched_face),
				// we set the according shift component to the overall reduction
				 notching[notched_face] = - body_dimensions[notched_face] * ( 1.-notch_reduction );

				// A quick assurance variable to assure that at least a single vertex has been found,
				// so our search criterion where to look for the vertices is not completely off
				 bool found_vertex=false;
				 double notch_location = ( (loading_type==enums::compression) ? ((length-width)/2.+width) : 0. );

				// Looping over all cells to notch the to-be-notched cells
				 for ( typename Triangulation<dim>::active_cell_iterator
							 cell = triangulation.begin_active();
							 cell != triangulation.end(); ++cell )
				 {
					for (unsigned int vertex=0; vertex < GeometryInfo<dim>::vertices_per_cell; ++vertex)
					 // Find vertices that are in the first 1/16 of the entire length
					  if ( std::abs(cell->vertex(vertex)[loading_direction]-notch_location) <=  notch_length )
						  if ( std::abs( cell->vertex(vertex)[notched_face] - body_dimensions[notched_face]) < search_tolerance )
						  {
							  // The found vertex is moved by the \a notching vector
							  // The notching shall be linear, hence a vertex in the notch is fully notched and the farther you
							  // move away from the notch the lower the notching gets (linearly).
							   cell->vertex(vertex) += notching * ( notch_length - std::abs(cell->vertex(vertex)[loading_direction]-notch_location) ) / notch_length;
							  found_vertex = true;
						  }
				 }

				AssertThrow(found_vertex, ExcMessage("BarModel<< We weren't able to find at least a single vertex to be notched."));

				// For compression we also notch the left edge
				if ( loading_type == enums::compression )
				{
					notch_location = (length-width)/2.;

					// A quick assurance variable to assure that at least a single vertex has been found,
					// so our search criterion where to look for the vertices is not completely off
					 found_vertex=false;

					// Looping over all cells to notch the to-be-notched cells
					 for ( typename Triangulation<dim>::active_cell_iterator
								 cell = triangulation.begin_active();
								 cell != triangulation.end(); ++cell )
					 {
						for (unsigned int vertex=0; vertex < GeometryInfo<dim>::vertices_per_cell; ++vertex)
						 // Find vertices that are in the first 1/16 of the entire length
						  if ( std::abs(cell->vertex(vertex)[loading_direction]-notch_location) <=  notch_length )
							  if ( std::abs( cell->vertex(vertex)[notched_face] - 0.) < search_tolerance )
							  {
								  // The found vertex is moved by the \a notching vector
								  // The notching shall be linear, hence a vertex in the notch is fully notched and the farther you
								  // move away from the notch the lower the notching gets (linearly).
								   cell->vertex(vertex) -= notching * ( notch_length - std::abs(cell->vertex(vertex)[loading_direction]-notch_location) ) / notch_length;
								  found_vertex = true;
							  }
					 }

					AssertThrow(found_vertex, ExcMessage("BarModel<< We weren't able to find at least a single vertex to be notched."));
				}
			}

			// Apply the manifolds
			if ( notch_rounded/*use round notch*/ && loading_type == enums::compression )
			{
				// Find the inclined faces
				for (typename Triangulation<dim>::active_cell_iterator
							 cell = triangulation.begin_active();
							 cell != triangulation.end(); ++cell)
				{
					for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
				    {
						Point<dim> face_center = cell->face(face)->center();
						if ( face_center[enums::x] > width/100. && face_center[enums::x] < (width/8.-width/100.) )
							cell->face(face)->set_all_manifold_ids(parameters_internal.manifold_id_left_radius);
						else if ( face_center[enums::x] < (width*99./100.) && face_center[enums::x] > (width*7./8.+width/100.) )
							cell->face(face)->set_all_manifold_ids(parameters_internal.manifold_id_right_radius);
				    }
				}

				const double actual_notch_hlength = 1.25;
				const double indent = body_dimensions[notched_face] * ( 1.-notch_reduction );
				const double offset_of_center_from_face = (actual_notch_hlength*actual_notch_hlength - indent*indent)/(2.*indent);

				// For the left radius
				 Point<dim> centre_left_radius (-offset_of_center_from_face, (length-width)/2.);
				 static SphericalManifold<dim> spherical_manifold_left (centre_left_radius);
				 triangulation.set_manifold(parameters_internal.manifold_id_left_radius,spherical_manifold_left);

				// For the right radius
				 Point<dim> centre_right_radius (width+offset_of_center_from_face, (length-width)/2.+width );
				 static SphericalManifold<dim> spherical_manifold_right (centre_right_radius);
				 triangulation.set_manifold(parameters_internal.manifold_id_right_radius,spherical_manifold_right);

				triangulation.refine_global(parameter.nbr_global_refinements);	// ... Parameter.prm file
			}

			// Refine in a special manner locally, if the number of local refinements is >0
			if ( parameter.nbr_holeEdge_refinements > 0 )
			{
				if ( loading_type == enums::compression )
				{
					for ( unsigned int nbr_local_ref=0; nbr_local_ref<parameter.nbr_holeEdge_refinements; nbr_local_ref++ )
					{
						for (typename Triangulation<dim>::active_cell_iterator
									 cell = triangulation.begin_active();
									 cell != triangulation.end(); ++cell)
						{
							for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
							  {
								Point<dim> face_center = cell->face(face)->center();
								// Find all cells that lay in an exemplary damage band with size 2xnotch_length along the diagonal
								if (    face_center[enums::y] < ((length-width)/2. + notch_length + face_center[enums::x])
									 && face_center[enums::y] > ((length-width)/2. - notch_length + face_center[enums::x]) )
								{
									cell->set_refine_flag();
									break;
								}
							  }
						}
						triangulation.execute_coarsening_and_refinement();
					}
				}
				else
				{
					for ( unsigned int nbr_local_ref=0; nbr_local_ref < parameter.nbr_holeEdge_refinements; nbr_local_ref++ )
					{
						for (typename Triangulation<dim>::active_cell_iterator
									 cell = triangulation.begin_active();
									 cell != triangulation.end(); ++cell)
						{
							for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
							  {
								// Find all cells that lay in an exemplary damage band with size 1.5 mm from the y=0 face
								if (std::abs(cell->face(face)->center()[enums::y] ) < length/5. )
								{
									cell->set_refine_flag();
									break;
								}
							  }
						}
						triangulation.execute_coarsening_and_refinement();
					}
				}
			}

			// Mark the innermost cell(s) for softening to trigger the damage development
			if ( triangulation.n_active_cells()>1)
			{
				bool found_cell=false;
				for (typename Triangulation<dim>::active_cell_iterator
							 cell = triangulation.begin_active();
							 cell != triangulation.end(); ++cell)
				{
					for (unsigned int vertex=0; vertex < GeometryInfo<dim>::vertices_per_cell; ++vertex)
					 // Find the cell that has one vertex at the origin
	//				  if ( (cell->vertex(vertex)).distance(centre)<1e-12 )
					 // Find cells that lay on the xz-plane (y=0)
					  if ( std::abs(cell->vertex(vertex)[1]) < 1 ) // vs 1e-12 used previously
					  {
						  // Avoid overwriting the tracked cell
						  if ( cell->material_id() != enums::tracked_QP )
							  cell->set_material_id(1);
						  found_cell = true;
						  break;
					  }

					if ( /*only weaken one cell:*/ false  && found_cell )	// ToDo: check: does a break leave only the inner for-loop?
						break;
				}

				AssertThrow(found_cell, ExcMessage("BarModel<< Was not able to identify the cell at the origin(0,0,0). Please recheck the triangulation or adapt the code."));
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

		const double search_tolerance = parameters_internal.search_tolerance;

		// USER PARAMETERS
		//const bool refinement_only_globally = false;
		//const bool damage_trigger_by_materialParameters = false;
		const bool damage_trigger_by_notching = true;
		const enums::enum_coord notched_face = enums::x;
		// The refined fraction consists of \a nbr_coarse_in_fine_section coarse elements that make up the fine section
		 const unsigned int nbr_coarse_in_fine_section = 2;
		const double refined_fraction = double(nbr_coarse_in_fine_section)/parameter.grid_y_repetitions;
		const bool use_fine_and_coarse_brick = true;
		const bool hardcoded_repetitions = false;

		// ToDo: use the values from the parameter file
		const double width = parameter.width; // use thickness=width for square bottom area
		const double thickness = parameter.thickness;
		const double length = parameter.height;
		const double notch_reduction = parameter.ratio_x;
		// The notch length is set (for consistency) s.t. its mesh discretisation is exact for 1 local refinement (start)
		 double notch_length = length/10.;//(2.*parameter.grid_y_repetitions);

		body_dimensions[enums::x] = width;
		body_dimensions[enums::y] = length;
		body_dimensions[enums::z] = thickness;

		// The bar is created from two bricks, where the first will be meshed very fine
		// and the second remains coarse. The bricks are spanned by three points.
		 Point<dim> p1 (0,0,0);
		 Point<dim> p2 (width, length * refined_fraction, thickness); // extends in y-direction its height (loaded in y-direction as the othe models)
		 Point<dim> p3 (0, length, 0); // extends in y-direction its height (loaded in y-direction as the othe models)
		 Point<dim> p4 (width, length, thickness);

		if ( use_fine_and_coarse_brick )
		{
			// Vector containing the number of elements in each dimension
			// The coarse segment consists of the set number of elements in the y-direction
			 std::vector<unsigned int> repetitions (3);
			 repetitions[enums::x] = 6 * (parameter.nbr_global_refinements+1);
			 repetitions[enums::y] = parameter.grid_y_repetitions * (1 - refined_fraction) * (parameter.nbr_global_refinements+1); // y
			 repetitions[enums::z] = parameter.nbr_elementsInZ;

			// The fine segment consists of at least 2 elements plus possible refinements
			 std::vector<unsigned int> repetitions_fine (3);
			 repetitions_fine[enums::x] = 6 * (parameter.nbr_global_refinements+1);
			 repetitions_fine[enums::y] = (parameter.grid_y_repetitions * refined_fraction) * std::pow(2.,parameter.nbr_holeEdge_refinements) * (parameter.nbr_global_refinements+1); // y
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
						if ( cell->center()[loading_direction] < length * refined_fraction )
						{
							cell->set_refine_flag();
							break;
						}
				}
				triangulation.execute_coarsening_and_refinement();
			}
		}
		else // use uniform brick with xy refinements
		{
			 std::vector<unsigned int> repetitions (3);

			 if ( hardcoded_repetitions )
			 {
				 repetitions[enums::x] = 15;
				 repetitions[enums::y] = 10;
				 repetitions[enums::z] = 4;
			 }
			 else
			 {
				 repetitions[enums::x] = std::pow(2.,parameter.nbr_holeEdge_refinements);
				 repetitions[enums::y] = std::pow(2.,parameter.nbr_holeEdge_refinements); // y
				 repetitions[enums::z] = parameter.nbr_elementsInZ;
			 }

			 GridGenerator::subdivided_hyper_rectangle ( triangulation,
														 repetitions,
														 p1,
														 p4 );

			 notch_length = length/8.;
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
				else if (std::abs(cell->face(face)->center()[1] - length) < search_tolerance)
				{
					cell->face(face)->set_boundary_id(enums::id_boundary_yPlus);
				}
				else if (std::abs(cell->face(face)->center()[2] - 0.0) < search_tolerance)
				{
					cell->face(face)->set_boundary_id(enums::id_boundary_zMinus);
				}
				else if (std::abs(cell->face(face)->center()[2] - thickness) < search_tolerance)
				{
					cell->face(face)->set_boundary_id(enums::id_boundary_zPlus);
				}
				else
				{
					// There are just eight faces, so if we missed one, something went clearly terribly wrong
					 AssertThrow(false, ExcMessage("BarModel - make_grid 3D<< Found an unidentified face at the boundary. Maybe it slipt through the assignment or that face is simply not needed. So either check the implementation or comment this line in the code"));
				}
			  }
		}

//		if ( refinement_only_globally )
//			triangulation.refine_global(parameter.nbr_global_refinements);	// ... Parameter.prm file
//		else // Refine in a special manner only cells around the origin
//		{
//			for ( unsigned int nbr_local_ref=0; nbr_local_ref < parameter.nbr_holeEdge_refinements; nbr_local_ref++ )
//			{
//				for (typename Triangulation<dim>::active_cell_iterator
//							 cell = triangulation.begin_active();
//							 cell != triangulation.end(); ++cell)
//				{
//					for (unsigned int vertex=0; vertex < GeometryInfo<dim>::vertices_per_cell; ++vertex)
//					{
//						// Find all cells that lay in an exemplary damage band with size 1/4 from the y=0 face
//						if ( cell->vertex(vertex)[enums::y] < length/4. )
//						{
//							cell->set_refine_flag();
//							break;
//						}
//					}
//				}
//				triangulation.execute_coarsening_and_refinement();
//			}
//		}
//
//		// Mark the innermost cell(s) for softening to trigger the damage development
//		if ( triangulation.n_active_cells()>1 && damage_trigger_by_materialParameters )
//		{
//			bool found_cell=false;
//			for (typename Triangulation<dim>::active_cell_iterator
//						 cell = triangulation.begin_active();
//						 cell != triangulation.end(); ++cell)
//			{
//				for (unsigned int vertex=0; vertex < GeometryInfo<dim>::vertices_per_cell; ++vertex)
//				 // Find the cell that has one vertex at the origin
////				  if ( (cell->vertex(vertex)).distance(centre)<1e-12 )
//				 // Find cells that lay on the xz-plane (y=0)
//				  if ( std::abs(cell->vertex(vertex)[enums::y]) < 1 ) // mark cells in the first millimeter as "to be softened"
//				  {
//					  cell->set_material_id(1);
//					  found_cell = true;
//					  break;
//				  }
//
//				if ( /*only weaken one cell:*/ false  && found_cell )	// ToDo: check: does a break leave only the inner for-loop?
//					break;
//			}
//
//			AssertThrow(found_cell, ExcMessage("BarModel<< Was not able to identify the cell at the origin(0,0,0). Please recheck the triangulation or adapt the code."));
//		}

		// Notch the specimen by moving some nodes inwards to form a notch
		if ( triangulation.n_active_cells() > 1 && damage_trigger_by_notching )
		{
			// Declare the shift vector for the notching
			 Point<3> notching; // initially zero
			// Depending on the desired notching direction (notched_face),
			// we set the according shift component to the overall reduction
			 notching[notched_face] = - body_dimensions[notched_face] * ( 1.-notch_reduction );

			// A quick assurance variable to assure that at least a single vertex has been found,
			// so our search criterion where to look for the vertices is not completely off
			 bool found_vertex=false;
			// Looping over all cells to notch the to-be-notched cells
			 for ( typename Triangulation<dim>::active_cell_iterator
						 cell = triangulation.begin_active();
						 cell != triangulation.end(); ++cell )
			 {
				for (unsigned int vertex=0; vertex < GeometryInfo<dim>::vertices_per_cell; ++vertex)
				 // Find vertices that are in the first 1/16 of the entire length
				  if ( std::abs(cell->vertex(vertex)[loading_direction]) <  notch_length )
					  if ( std::abs( cell->vertex(vertex)[notched_face] - body_dimensions[notched_face]) < search_tolerance )
					  {
						  // The found vertex is moved by the \a notching vector
						  // The notching shall be linear, hence a vertex in the notch is fully notched and the farther you
						  // move away from the notch the lower the notching gets (linearly).
						   cell->vertex(vertex) += notching * ( notch_length - cell->vertex(vertex)[loading_direction] ) / notch_length;
						  found_vertex = true;
					  }
			 }

			AssertThrow(found_vertex, ExcMessage("BarModel<< We weren't able to find at least a single vertex to be notched."));
		}


		// include the following two scopes to see directly how the variation of the input parameters changes the geometry of the grid
		/*
		{
			std::ofstream out ("grid-3d_quarter_plate_merged.eps");
			GridOut grid_out;
			grid_out.write_eps (triangulation, out);
			std::cout << "Grid written to grid-3d_quarter_plate_merged.eps" << std::endl;
		}
		{
			std::ofstream out_ucd("Grid-3d_quarter_plate_merged.inp");
			GridOut grid_out;
			GridOutFlags::Ucd ucd_flags(true,true,true);
			grid_out.set_flags(ucd_flags);
			grid_out.write_ucd(triangulation, out_ucd);
			std::cout<<"Mesh written to Grid-3d_quarter_plate_merged.inp "<<std::endl;
		}
		*/
	}
}
