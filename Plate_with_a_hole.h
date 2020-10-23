#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>

#include "../MA-Code/enumerator_list.h"

#include <iostream>
#include <fstream>
#include <cmath>

using namespace dealii;

namespace PlateWithAHole
/*
 * A quarter of a plate with hole in 2D or 1/8 in 3D
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

	// Some internal parameters
	 struct parameterCollection
	 {
		const types::boundary_id boundary_id_hole = 10;
		const types::manifold_id manifold_id_hole = 10;

		const double search_tolerance = 1e-12;
	 };

	template<int dim>
	void make_constraints ( AffineConstraints<double>  &constraints, const FESystem<dim> &fe, unsigned int &n_components, DoFHandler<dim> &dof_handler_ref,
							const bool &apply_dirichlet_bc, const double &current_load_increment,
							const Parameter::GeneralParameters &parameter )
	{
		/* inputs:
		 * dof_handler_ref,
		 * fe
		 * apply_dirichlet_bc
		 * constraints
		 * current_load_increment
		 */

		const FEValuesExtractors::Vector displacement(0);
		const FEValuesExtractors::Scalar x_displacement(0);
		const FEValuesExtractors::Scalar y_displacement(1);

		// Update and apply new constraints
		//		on y0_plane (bottom) for counter bearing (displacement_in_y = 0)
		//		to stop the plate from sliding to the left and right, we fix a single node in the center in x=0

		// on bottom edge
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

		// fixed node in x-direction
		// @todo Currently we fix an entire cell face at x=0 to x=0
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

		// on thickness symmetry plane
		if ( dim==3 )
		{
			const FEValuesExtractors::Scalar z_displacement(2);

			if ( false/*free 3D deformation*/ )
			{

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

				// For the exact 2D plane strain condition (infinite thickness)
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
			else
			{
				if (apply_dirichlet_bc == true )
				{
					VectorTools::interpolate_boundary_values(
																dof_handler_ref,
																enums::id_boundary_xMinus,
																ZeroFunction<dim> (n_components),
																constraints,
																fe.component_mask(z_displacement)
															);
				}
				else	// in the exact same manner
				{
					VectorTools::interpolate_boundary_values(
																dof_handler_ref,
																enums::id_boundary_xMinus,
																ZeroFunction<dim> (n_components),
																constraints,
																fe.component_mask(z_displacement)
															);
				}
			}
		}

		if ( parameter.driver == 2/*Dirichlet*/ )
		{
			// on top/loaded edge
			if (apply_dirichlet_bc == true )
			{
				VectorTools::interpolate_boundary_values(
															dof_handler_ref,
															id_boundary_load,
															// ToDo: adapt this to also work for the load_history
															ConstantFunction<dim> (current_load_increment/*add only the increment*/, n_components),
															//ConstantFunction<dim> (parameter.pressure_load / nbr_loadsteps/*add only the increment*/, n_components),
															constraints,
															fe.component_mask(y_displacement)
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
			}
		}
	}

	// ToDo-optimize: use existing DII command	void GridGenerator::plate_with_a_hole

	// to see the effects of the inputs (lengths, refinements, etc) consider using the output (.eps, etc) below
	void make_2d_plate_with_hole( Triangulation<2> &tria_2d,
										  const double half_length,
										  //const double half_width,
										  const double hole_radius,
										  //const double hole_division_fraction,
										  const Parameter::GeneralParameters &parameter )
	{
		//const double width =  2.0*half_width;
		const double hole_diameter = 2.0*hole_radius;
		//const double internal_width = hole_diameter + hole_division_fraction*(width - hole_diameter);

		Point<2> centre_2d(0,0);
		const types::manifold_id  	polar_manifold_id = 0;
		const types::manifold_id  	tfi_manifold_id = 1;

		const double height2Width_ratio = 1.;

		GridGenerator::plate_with_a_hole 	( 	tria_2d,
												hole_radius,
												half_length/height2Width_ratio, // Width of the plate is \a length, the height is 2*length
												half_length*(1.-1./height2Width_ratio),
												half_length*(1.-1./height2Width_ratio),
												0.,
												0.,
												centre_2d,
												polar_manifold_id,
												tfi_manifold_id
											);

		// Attach a manifold to the curved boundary and refine
		// @note We can only guarantee that the vertices sit on the curve, so we must test with their position instead of the cell centre.
		for (typename Triangulation<2>::active_cell_iterator
		   cell = tria_2d.begin_active();
		   cell != tria_2d.end(); ++cell)
		{
		  for (unsigned int face=0; face<GeometryInfo<2>::faces_per_cell; ++face)
			if (cell->face(face)->at_boundary())
			  for (unsigned int vertex=0; vertex<GeometryInfo<2>::vertices_per_face; ++vertex)
			  {
				if (std::abs(cell->face(face)->vertex(vertex).distance(centre_2d) - hole_diameter/2.0) < 1e-12)
				 {
					cell->face(face)->set_manifold_id(10);
					break;
				 }
			  }
		}

		SphericalManifold<2> spherical_manifold_2d (centre_2d);
		tria_2d.set_manifold(10,spherical_manifold_2d);

		if ( parameter.stepwise_global_refinement==true ) // for the special case of step by step global refinement ...
			tria_2d.refine_global( 1 );	// ... only refine the initial grid once
		else	// for the standard case of AMR refine the grid as specified in the ...
			tria_2d.refine_global(parameter.nbr_global_refinements);	// ...Parameter.prm file; has to be refined before the manifolds are deleted again

//			AssertThrow(parameter.nbr_holeEdge_refinements==0, ExcMessage("QuarterPlate mesh creation: Sorry, right now you cannot use hole edge refinements."));

		tria_2d.reset_manifold(10); // Clear manifold

		//GridGenerator::flatten_triangulation(tria_2d,tria_2d);

		// include the following scope to see directly how the variation of the input parameters changes the geometry of the grid

//		{
//			std::ofstream out ("grid-tria_2d.eps");
//			GridOut grid_out;
//			grid_out.write_eps (tria_2d, out);
//			std::cout << "Grid written to grid-tria_2d.eps" << std::endl;
//		}
	}

// 2D grid
	template <int dim>
	void make_grid( Triangulation<2> &triangulation, const Parameter::GeneralParameters &parameter )
	{
		parameterCollection parameters_internal;

		// size of the plate divided by the size of the hole
		  double ratio_width_To_holeRadius = parameter.width;
		  double hwidth = parameter.width;
		  double holeRadius = parameter.holeRadius;

		  // size of the inner mesh (hypercube with hole) relative to size of the whole plate
		  //double ratio_x = parameter.ratio_x;

		const double search_tolerance = parameters_internal.search_tolerance;

		make_2d_plate_with_hole(
											triangulation,
											ratio_width_To_holeRadius,			// length
											//ratio_width_To_holeRadius * 1.0,	// width: *1.0 => square
											holeRadius,	// hole radius = diameter/2
										    //ratio_x,
											parameter
										);

		//Clear boundary ID's
		for ( typename Triangulation<dim>::active_cell_iterator
				cell = triangulation.begin_active();
					cell != triangulation.end(); ++cell )
		{
			for (unsigned int face=0; face < GeometryInfo<dim>::faces_per_cell; ++face)
			  if (cell->face(face)->at_boundary())
			  {
				  cell->face(face)->set_all_boundary_ids(0);
				  cell->face(face)->set_manifold_id(0);
			  }
		}

		//Set boundary IDs and manifolds
		const Point<dim> centre (0,0);
		for ( typename Triangulation<dim>::active_cell_iterator
				cell = triangulation.begin_active();
					cell != triangulation.end(); ++cell )
		{
			for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
			{
			// at boundary
			  if (cell->face(face)->at_boundary())
			  {
				//Set boundary IDs
				 // The bottom is now all the way down to -hwidth
				if (std::abs(cell->face(face)->center()[enums::y] + hwidth) < search_tolerance)
				{
					cell->face(face)->set_boundary_id(enums::id_boundary_yMinus);	// the bottom edge

					for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_face; ++vertex)
					 if (std::abs(cell->vertex(vertex)[enums::y] - 0) < search_tolerance)
					  if (std::abs(cell->vertex(vertex)[enums::x] - parameter.holeRadius) < search_tolerance)
					  {
						  // We found the cell that lies at the bottom edge next to the hole (bottom left corner)
						  cell->set_material_id( enums::tracked_QP );
						  break;
					  }
				}
				else if (std::abs(cell->face(face)->center()[1] - ratio_width_To_holeRadius) < search_tolerance)
				{
					cell->face(face)->set_boundary_id(enums::id_boundary_yPlus); // the top edge
				}
				else
				{
					// Be aware that we have to access the vertex as cell->face->vertex
					for (unsigned int vertex=0; vertex < GeometryInfo<dim>::vertices_per_face; ++vertex)
					  if (std::abs(cell->face(face)->vertex(vertex).distance(centre) - parameter.holeRadius) < search_tolerance)
					  {
						  cell->face(face)->set_boundary_id(parameters_internal.boundary_id_hole);	// the hole edge
						  break;
					  }
				}

				// Set manifold IDs
				for (unsigned int vertex=0; vertex < GeometryInfo<dim>::vertices_per_face; ++vertex)
				  if (std::abs(cell->vertex(vertex).distance(centre) - parameter.holeRadius) < search_tolerance)
				  {
					  cell->face(face)->set_manifold_id(parameters_internal.manifold_id_hole);
					  break;
				  }
			  }
			}
		}

		static SphericalManifold<dim> spherical_manifold (centre);
		triangulation.set_manifold(parameters_internal.manifold_id_hole,spherical_manifold);

		// The following does not work?
//		// pre-refinement of the damaged area (around y=0)
//		for (unsigned int refine_counter=0; refine_counter<parameter.nbr_holeEdge_refinements; refine_counter++)
//		{
//			for (typename Triangulation<dim>::active_cell_iterator
//						 cell = triangulation.begin_active();
//						 cell != triangulation.end(); ++cell)
//			{
//				if ( std::abs( cell->center()[enums::y]) < hwidth/4. )
//					cell->set_refine_flag();
//			}
//			triangulation.execute_coarsening_and_refinement();
//		}

		// refine hole boundary
		for (unsigned int refine_counter=0; refine_counter<parameter.nbr_holeEdge_refinements; refine_counter++)
		{
			for (typename Triangulation<dim>::active_cell_iterator
						 cell = triangulation.begin_active();
						 cell != triangulation.end(); ++cell)
			{
				for (unsigned int face=0; face < GeometryInfo<dim>::faces_per_cell; ++face)
					if ( cell->face(face)->boundary_id() == parameters_internal.boundary_id_hole )
					{
						cell->set_refine_flag();
						break;
					}
			}
			triangulation.execute_coarsening_and_refinement();
		}

		// The x-fixed faces are set now, to get the smallest possible face
		for ( typename Triangulation<dim>::active_cell_iterator
				cell = triangulation.begin_active();
					cell != triangulation.end(); ++cell )
		{
			for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
			{
			  if (cell->face(face)->at_boundary() )
			  {
				  if ( ( cell->face(face)->center()[enums::y] + hwidth ) < search_tolerance ) // lower face
				  {
					for (unsigned int vertex=0; vertex < GeometryInfo<dim>::vertices_per_face; ++vertex)
					 {
						 //Project the cell vertex to the XY plane and test the distance from the cylinder axis
							Point<2> vertex_v = cell->face(face)->vertex(vertex);
							Point<2> desired_mid_point (0,-hwidth);

							if (std::abs(vertex_v.distance(desired_mid_point)) < search_tolerance)
							{
								std::cout << "found l at " << vertex_v << std::endl;
								cell->face(face)->set_boundary_id(enums::id_boundary_xMinus);
								break;
							}
					 }
				  }
					// We also fix some face at the top of the plate to x=0, so it does not drift/shear sidewides
					// @todo-ensure Do we create stress peaks by doing so?
					else if ( ( cell->face(face)->center()[enums::y] - hwidth ) < search_tolerance ) // upper face
					{
						for (unsigned int vertex=0; vertex < GeometryInfo<dim>::vertices_per_face; ++vertex)
						 {
						 //Project the cell vertex to the XY plane and test the distance from the cylinder axis
							Point<2> vertex_v = cell->face(face)->vertex(vertex);
							Point<2> desired_mid_point (0,+hwidth);

							if (std::abs(vertex_v.distance(desired_mid_point)) < search_tolerance)
							{
								std::cout << "found u at " << vertex_v << std::endl;
								cell->face(face)->set_boundary_id(enums::id_boundary_xMinus);
								break;
							}
						}
					}
				}
			  }
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


	/**
	 * Make 3D mesh of plate with a hole (only positive z-part)
	 */
	template <int dim>
	void make_grid( Triangulation<3> &triangulation, const Parameter::GeneralParameters &parameter )
	{
		parameterCollection parameters_internal;

		// size of the plate divided by the size of the hole
		 double ratio_width_To_holeRadius = parameter.width;
		 double hwidth = parameter.width;
		 double holeRadius = parameter.holeRadius;

		  // size of the inner mesh (hypercube with hole) relative to size of the whole plate
		  //double ratio_x = parameter.ratio_x;

		const double search_tolerance = parameters_internal.search_tolerance;

		Triangulation<2> tria_2d;
				make_2d_plate_with_hole (
											tria_2d,
											ratio_width_To_holeRadius,			// length
											//ratio_width_To_holeRadius * 1.0,	// width: *1.0 => square
											holeRadius,	// hole radius = diameter/2
											//ratio_x,
											parameter
										);

		// only relevant for 3d grid:
		  const unsigned int n_repetitions_z = parameter.nbr_elementsInZ;			// nbr of Unterteilungen in z-direction for 3d meshing; 1=one element in z; 2=two el.s in z; ...

		GridGenerator::extrude_triangulation(tria_2d,
										   n_repetitions_z+1,
										   parameter.thickness/2.0,
										   triangulation);

		// Clear boundary ID's
		for (typename Triangulation<dim>::active_cell_iterator
			 cell = triangulation.begin_active();
			 cell != triangulation.end(); ++cell )
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
			{
				// Set the INTERNAL (not at the boundary) face to the xminus id
				// Find the cell face at x=0, containing the node (0,-hwidth,0) pointing in x-direction
				if ( (std::abs(cell->face(face)->center()[enums::x] - 0.0) < search_tolerance) )
				{
					// @todo-optimize: We only look for single cell-face, so we could stop checking this, if we found something (e.g. bool found_xMinus)
					if ( cell->face(face)->center()[enums::y] < 0 )
					{
						for (unsigned int vertex=0; vertex < GeometryInfo<dim>::vertices_per_face; ++vertex)
						 {
						 //Project the cell vertex to the XY plane and test the distance from the cylinder axis
							Point<3> vertex_v = cell->face(face)->vertex(vertex);
							Point<3> desired_mid_point (0,-hwidth,0);

							if (std::abs(vertex_v.distance(desired_mid_point)) < search_tolerance)
							{
								cell->face(face)->set_boundary_id(enums::id_boundary_xMinus);
								break;
							}
						}
					}
					// We also fix some face at the top of the plate to x=0, so it does not drift/shear sidewides
					// @todo-ensure Do we create stress peaks by doing so?
					else if ( cell->face(face)->center()[enums::y] > 0 )
					{
						for (unsigned int vertex=0; vertex < GeometryInfo<dim>::vertices_per_face; ++vertex)
						 {
						 //Project the cell vertex to the XY plane and test the distance from the cylinder axis
							Point<3> vertex_v = cell->face(face)->vertex(vertex);
							Point<3> desired_mid_point (0,+hwidth,0);

							if (std::abs(vertex_v.distance(desired_mid_point)) < search_tolerance)
							{
								cell->face(face)->set_boundary_id(enums::id_boundary_xMinus);
								break;
							}
						}
					}
				}
			  // Set boundary IDs
			  if (cell->face(face)->at_boundary())
			  {
				// The bottom face (yMinus) is now at the negative hwidth (plate reaches from y = -hwidth to +hwidth)
				if (std::abs(cell->face(face)->center()[enums::y] + hwidth) < search_tolerance)
				{
					cell->face(face)->set_boundary_id(enums::id_boundary_yMinus);
					// Tracked QP
					for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_face; ++vertex)
					{
					 if (std::abs(cell->vertex(vertex)[enums::y] - 0) < search_tolerance)
						if (std::abs(cell->vertex(vertex)[enums::z] - 0) < search_tolerance)
						  if (std::abs(cell->vertex(vertex)[enums::x] - parameter.holeRadius) < search_tolerance)
						  {
							  // We found the cell that lies at the bottom edge next to the hole (bottom left corner)
							  cell->set_material_id( enums::tracked_QP );
							  break;
						  }
					}
				}
				else if (std::abs(cell->face(face)->center()[1] - ratio_width_To_holeRadius) < search_tolerance)
				{
					cell->face(face)->set_boundary_id(enums::id_boundary_yPlus);
				}
				else if (std::abs(cell->face(face)->center()[2] - 0.0) < search_tolerance)
				{
					cell->face(face)->set_boundary_id(enums::id_boundary_zMinus);
				}
				else if (std::abs(cell->face(face)->center()[2] - parameter.thickness/2.0) < search_tolerance)
				{
					cell->face(face)->set_boundary_id(enums::id_boundary_zPlus);
				}
				else
				{
					for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_face; ++vertex)
					 {
					 //Project the cell vertex to the XY plane and test the distance from the cylinder axis
						Point<dim> vertex_proj = cell->vertex(vertex);
						vertex_proj[2] = 0.0;
						if (std::abs(vertex_proj.distance(centre) - parameter.holeRadius) < search_tolerance)
						{
							cell->face(face)->set_boundary_id(parameters_internal.boundary_id_hole);
							break;
						}
					}
				}

				//Set manifold IDs
				for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_face; ++vertex)
				{
					//Project the cell vertex to the XY plane and test the distance from the cylinder axis
					Point<dim> vertex_proj = cell->vertex(vertex);
					vertex_proj[2] = 0.0;
					if (std::abs(vertex_proj.distance(centre) - parameter.holeRadius) < search_tolerance)
					{
						//Set manifold ID on face and edges
						cell->face(face)->set_all_manifold_ids(parameters_internal.manifold_id_hole);
						break;
					  }
				  }
			  }
			}
		}

		// pre-refinement of the possibly damaged area (around y=0)
//		for (unsigned int refine_counter=0; refine_counter<parameter.nbr_holeEdge_refinements; refine_counter++)
//		{
//			for (typename Triangulation<dim>::active_cell_iterator
//						 cell = triangulation.begin_active();
//						 cell != triangulation.end(); ++cell)
//			{
//				double distance2D = std::sqrt( cell->center()[0]*cell->center()[0] + cell->center()[1]*cell->center()[1] );
//
//				for ( unsigned int face=0; face < GeometryInfo<dim>::faces_per_cell; face++ )
//				{
//					// We use the absolute y-coordinate here to refine cells above and below the horizontal centre line
//					if ( std::abs( cell->center()[loading_direction] ) < 30 )
//					{
//						cell->set_refine_flag(RefinementCase<dim>::cut_xy);
//						break;
//					}
//				}
//			}
//			triangulation.execute_coarsening_and_refinement();
//		}

		// refine hole boundary
		for (unsigned int refine_counter=0; refine_counter<parameter.nbr_holeEdge_refinements; refine_counter++)
		{
			for (typename Triangulation<dim>::active_cell_iterator
						 cell = triangulation.begin_active();
						 cell != triangulation.end(); ++cell)
			{
				for (unsigned int face=0; face < GeometryInfo<dim>::faces_per_cell; ++face)
					if ( cell->face(face)->boundary_id() == parameters_internal.boundary_id_hole )
					{
						cell->set_refine_flag();
						break;
					}
			}
			triangulation.execute_coarsening_and_refinement();
		}

		// The x-fixed faces are set now, to get the smallest possible face
		for ( typename Triangulation<dim>::active_cell_iterator
				cell = triangulation.begin_active();
					cell != triangulation.end(); ++cell )
		{
			for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
			{
				// Set the INTERNAL (not at the boundary) face to the xminus id
				// Find the cell face at x=0, containing the node (0,-hwidth,0) pointing in x-direction
				if ( (std::abs(cell->face(face)->center()[enums::x] - 0.0) < search_tolerance) )
				{
					// @todo-optimize: We only look for single cell-face, so we could stop checking this, if we found something (e.g. bool found_xMinus)
					if ( cell->face(face)->center()[enums::y] < 0 )
					{
						for (unsigned int vertex=0; vertex < GeometryInfo<dim>::vertices_per_face; ++vertex)
						 {
						 //Project the cell vertex to the XY plane and test the distance from the cylinder axis
							Point<3> vertex_v = cell->face(face)->vertex(vertex);
							Point<3> desired_mid_point (0,-hwidth);

							if (std::abs(vertex_v.distance(desired_mid_point)) < search_tolerance)
							{
								cell->face(face)->set_boundary_id(enums::id_boundary_xMinus);
								break;
							}
						}
					}
					// We also fix some face at the top of the plate to x=0, so it does not drift/shear sidewides
					// @todo-ensure Do we create stress peaks by doing so?
					else if ( cell->face(face)->center()[enums::y] > 0 )
					{
						for (unsigned int vertex=0; vertex < GeometryInfo<dim>::vertices_per_face; ++vertex)
						 {
						 //Project the cell vertex to the XY plane and test the distance from the cylinder axis
							Point<3> vertex_v = cell->face(face)->vertex(vertex);
							Point<3> desired_mid_point (0,+hwidth);

							if (std::abs(vertex_v.distance(desired_mid_point)) < search_tolerance)
							{
								cell->face(face)->set_boundary_id(enums::id_boundary_xMinus);
								break;
							}
						}
					}
				}
			}
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
