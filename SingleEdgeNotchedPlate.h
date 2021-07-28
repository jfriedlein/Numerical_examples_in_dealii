#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>

#include <iostream>
#include <fstream>
#include <cmath>

using namespace dealii;

namespace SingleEdgeNotchedPlate
/*
 * A single element with three symmetry constraints, loaded in y-direction, dimensions 1x1x1
 * By selecting 1 global refinement we can create the distorted 8 element test
 *
 * STILL UNCERTIFIED ( since introduction of contact, requires update of standard to incorporate contact numEx)
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
	 };


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

		const FEValuesExtractors::Vector displacement(0);
		const FEValuesExtractors::Scalar x_displacement(0);
		const FEValuesExtractors::Scalar y_displacement(1);

		// on Y0 edge
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

		if (apply_dirichlet_bc == true )
		{
			VectorTools::interpolate_boundary_values(
														dof_handler_ref,
														enums::id_boundary_fix,
														ZeroFunction<dim> (n_components),
														constraints,
														fe.component_mask(x_displacement)
													);
			VectorTools::interpolate_boundary_values(
														dof_handler_ref,
														enums::id_boundary_fix,
														ZeroFunction<dim> (n_components),
														constraints,
														fe.component_mask(y_displacement)
													);
		}
		else	// in the exact same manner
		{
			VectorTools::interpolate_boundary_values(
														dof_handler_ref,
														enums::id_boundary_fix,
														ZeroFunction<dim> (n_components),
														constraints,
														fe.component_mask(x_displacement)
													);
			VectorTools::interpolate_boundary_values(
														dof_handler_ref,
														enums::id_boundary_fix,
														ZeroFunction<dim> (n_components),
														constraints,
														fe.component_mask(y_displacement)
													);
		}


		// on Z0 plane
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
		}

		if ( parameter.driver == enums::Dirichlet ) // ToDo-optimize: use string in parameterfile denoting "Dirichlet" so the enumerator is not undermined
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


	template<int dim>
	void shift_vertex_by_vector ( Triangulation<dim> &tria, const std::vector< Point<dim> > &points, const std::vector< Point<dim> > &shift)
	{
		unsigned int shifted_node = 0;
		const unsigned int n_points = points.size();
		for (typename Triangulation<dim>::active_cell_iterator
		   cell = tria.begin_active();
		   cell != tria.end(); ++cell)
		{
		  for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_cell; ++vertex )
			  for ( unsigned int i=0; i<n_points; i++)
				  if ( cell->vertex(vertex).distance(points[i]) < 1e-12/*search_tolerance*/ )
				  {
					  cell->vertex(vertex) += shift[i];
					  shifted_node += 1; // -> We have shifted at least a single node
				  }
		}
		AssertThrow( shifted_node == n_points, ExcMessage("OneElementTest<< Distortion, we only shifted "+std::to_string(shifted_node)+
														  " instead of "+std::to_string(n_points)+" vertices."));
	}


// 3d grid
	template <int dim>
	void make_grid( Triangulation<3> &triangulation, const Parameter::GeneralParameters &parameter )
	{
		parameterCollection parameters_internal;

		const double search_tolerance = parameters_internal.search_tolerance;

		const double width = 1; // unit cube

        // USER PARAMETERS:
        const bool twitch_Dis8El_cube = true;

		GridGenerator::hyper_cube(triangulation);

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
		const Point<dim> direction (0,0,1);
		const Point<dim> centre (0,0,0);
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
				}
				else if (std::abs(cell->face(face)->center()[0] - width) < search_tolerance)
				{
					cell->face(face)->set_boundary_id(enums::id_boundary_xPlus);
				}
				else if (std::abs(cell->face(face)->center()[1] - 0.0) < search_tolerance)
				{
					cell->face(face)->set_boundary_id(enums::id_boundary_yMinus);
				}
				else if (std::abs(cell->face(face)->center()[1] - width) < search_tolerance)
				{
					cell->face(face)->set_boundary_id(enums::id_boundary_yPlus);
				}
				else if (std::abs(cell->face(face)->center()[2] - 0.0) < search_tolerance)
				{
					cell->face(face)->set_boundary_id(enums::id_boundary_zMinus);
				}
				else if (std::abs(cell->face(face)->center()[2] - width) < search_tolerance)
				{
					cell->face(face)->set_boundary_id(enums::id_boundary_zPlus);
				}
				else
				{
					AssertThrow(false, ExcMessage("OneElement - make_grid 3D<< Found an unidentified face at the boundary. Maybe it slipt through the assignment or that face is simply not needed. So either check the implementation or comment this line in the code"));
				}
			  }
		}

		triangulation.refine_global(parameter.nbr_global_refinements);	// ... Parameter.prm file

		// Distortion
		if ( /*distortion*/ false )
		{
			if ( triangulation.n_active_cells()==1 )
			{
				std::vector< Point<3> > points_xyz (2);
				std::vector< Point<3> > shift_dxdydz (2);

				Point<3> x1y1z1 (1,1,1);
				Point<3> shift_of_x1y1z1 (-0.5,0,-0.5);
				Point<3> x0y1z1 (0,1,1);
				Point<3> shift_of_x0y1z1 (0,0,-0.25);

				points_xyz[0] = x1y1z1;
				points_xyz[1] = x0y1z1;
				shift_dxdydz[0] = shift_of_x1y1z1;
				shift_dxdydz[1] = shift_of_x0y1z1;

				shift_vertex_by_vector( triangulation, points_xyz, shift_dxdydz );
			}
			// Distorted 8 elements (Dis8El)
			else if ( triangulation.n_active_cells()==8 )
			{
				std::vector< Point<3> > points_xyz (16);
				std::vector< Point<3> > shift_dxdydz (16);

				Point<3> x0y05z0 (0,0.5,0);
				Point<3> shift_of_x0y05z0(0,0.1,0);
				points_xyz[0] = x0y05z0;
				shift_dxdydz[0] = shift_of_x0y05z0;

				Point<3> x05y0z0 (0.5,0,0);
				Point<3> shift_of_x05y0z0(0,-0.125,0);
				points_xyz[1] = x05y0z0;
				shift_dxdydz[1] = shift_of_x05y0z0;

				Point<3> x05y05z0 (0.5,0.5,0);
				Point<3> shift_of_x05y05z0 (0,-0.3,0);
				points_xyz[2] = x05y05z0;
				shift_dxdydz[2] = shift_of_x05y05z0;

				Point<3> x1y05z0 (1,0.5,0);
				Point<3> shift_of_x1y05z0(0,-0.2,0);
				points_xyz[3] = x1y05z0;
				shift_dxdydz[3] = shift_of_x1y05z0;

				Point<3> x05y1z0 (0.5,1,0);
				Point<3> shift_of_x05y1z0(0.15,0,0);
				points_xyz[4] = x05y1z0;
				shift_dxdydz[4] = shift_of_x05y1z0;

				Point<3> x0y05z05 (0,0.5,0.5);
				Point<3> shift_of_x0y05z05(0,0.2,0);
				points_xyz[5] = x0y05z05;
				shift_dxdydz[5] = shift_of_x0y05z05;

				Point<3> x05y05z05 (0.5,0.5,0.5);
				Point<3> shift_of_x05y05z05(0.1,-0.1,-0.15);
				points_xyz[6] = x05y05z05;
				shift_dxdydz[6] = shift_of_x05y05z05;

				Point<3> x1y05z05 (1,0.5,0.5);
				Point<3> shift_of_x1y05z05(0,0.1,-0.1);
				points_xyz[7] = x1y05z05;
				shift_dxdydz[7] = shift_of_x1y05z05;

				Point<3> x0y1z05 (0,1,0.5);
				Point<3> shift_of_x0y1z05(0,0,0.1);
				points_xyz[8] = x0y1z05;
				shift_dxdydz[8] = shift_of_x0y1z05;

				Point<3> x05y1z05 (0.5,1,0.5);
				Point<3> shift_of_x05y1z05(0.1,0,-0.05);
				points_xyz[9] = x05y1z05;
				shift_dxdydz[9] = shift_of_x05y1z05;

				Point<3> x1y1z05 (1,1,0.5);
				Point<3> shift_of_x1y1z05(0,0,0.2);
				points_xyz[10] = x1y1z05;
				shift_dxdydz[10] = shift_of_x1y1z05;

				Point<3> x05y0z1 (0.5,0,1);
				Point<3> shift_of_x05y0z1(-0.1,0,0);
				points_xyz[11] = x05y0z1;
				shift_dxdydz[11] = shift_of_x05y0z1;

				Point<3> x0y05z1 (0,0.5,1);
				Point<3> shift_of_x0y05z1(0,-0.125,0);
				points_xyz[12] = x0y05z1;
				shift_dxdydz[12] = shift_of_x0y05z1;

				Point<3> x05y05z1 (0.5,0.5,1);
				Point<3> shift_of_x05y05z1(0,-0.2,0);
				points_xyz[13] = x05y05z1;
				shift_dxdydz[13] = shift_of_x05y05z1;

				Point<3> x1y05z1 (1,0.5,1);
				Point<3> shift_of_x1y05z1(0,0.25,0);
				points_xyz[14] = x1y05z1;
				shift_dxdydz[14] = shift_of_x1y05z1;

				Point<3> x05y1z1 (0.5,1,1);
				Point<3> shift_of_x05y1z1(0.1,0,0);
				points_xyz[15] = x05y1z1;
				shift_dxdydz[15] = shift_of_x05y1z1;

				shift_vertex_by_vector( triangulation, points_xyz, shift_dxdydz );

				if ( twitch_Dis8El_cube /*to be twitched*/)
				{
					std::vector< Point<3> > point_xyz (1);
					std::vector< Point<3> > shift_point (1);

					Point<3> x1y1z1 (1,1,1);
					Point<3> shift_of_x1y1z1(0.01,0,0);
					point_xyz[0] = x1y1z1;
					shift_point[0] = shift_of_x1y1z1;

					shift_vertex_by_vector( triangulation, point_xyz, shift_point );
				}
			}
		}

		if ( triangulation.n_active_cells()>1)
		{
			bool found_cell=false;
			for (typename Triangulation<dim>::active_cell_iterator
						 cell = triangulation.begin_active();
						 cell != triangulation.end(); ++cell)
			{
				for (unsigned int vertex=0; vertex < GeometryInfo<dim>::vertices_per_cell; ++vertex)
				  if ( (cell->vertex(vertex)).distance(centre)<1e-12 )
				  {
					  cell->set_material_id(1);
					  found_cell = true;
					  break;
				  }

				if ( found_cell )
					break;
			}

			AssertThrow(found_cell, ExcMessage("OneElement<< Was not able to identify the cell at the origin(0,0,0). Please recheck the triangulation or adapt the code."));
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


// 2d grid
	template <int dim>
	void make_grid( Triangulation<2> &triangulation, const Parameter::GeneralParameters &parameter )
	{
		parameterCollection parameters_internal;

		const double search_tolerance = parameters_internal.search_tolerance;

		const double width = parameter.width;
		const double height = parameter.width/2.;

		body_dimensions[enums::x] = width;
		body_dimensions[enums::y] = height;

		// The points that span the brick
		 Point<dim> p1 (0,0);
		 Point<dim> p2 (width, height); // extends in y-direction its length (loaded in y-direction as the other models)

		 // @note We always need an even number of elements in x-direction for the half notch and the half support on yMinus
		// vector containing the number of elements in each dimension
		 std::vector<unsigned int> repetitions (2);
		 repetitions[0]=2; // x
		 repetitions[1]=1; // y

		GridGenerator::subdivided_hyper_rectangle 	( 	triangulation,
														repetitions,
														p1,
														p2
													);
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
		const Point<dim> centre (0,0);

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
					//cell->face(face)->set_boundary_id(enums::id_boundary_xMinus);
				}
				else if (std::abs(cell->face(face)->center()[0] - width) < search_tolerance)
				{
					//cell->face(face)->set_boundary_id(enums::id_boundary_xPlus);
				}
				else if (std::abs(cell->face(face)->center()[1] - 0.0) < search_tolerance)
				{
					if ( cell->face(face)->center()[enums::x] > width/2. )
						cell->face(face)->set_boundary_id(enums::id_boundary_yMinus);
				}
				else if (std::abs(cell->face(face)->center()[1] - height) < search_tolerance)
				{
					cell->face(face)->set_boundary_id(enums::id_boundary_yPlus);
				}
//				else
//				{
//					AssertThrow(false, ExcMessage("SENP - make_grid 2D<< Found an unidentified face at the boundary. Maybe it slipt through the assignment or that face is simply not needed. So either check the implementation or comment this line in the code"));
//				}
			  }
		}

		// Distort
//		for (typename Triangulation<dim>::active_cell_iterator
//					 cell = triangulation.begin_active();
//					 cell != triangulation.end(); ++cell)
//		{
//			//cell->vertex(3)[enums::y] -= 0.15;
//			Point<dim> corner_leftBottom (width,0);
//			for (unsigned int vertex=0; vertex < GeometryInfo<dim>::vertices_per_cell; ++vertex)
//			  if ( (cell->vertex(vertex)).distance(corner_leftBottom)<1e-12 )
//			  {
//				  cell->vertex(vertex)[enums::x] -= 0.1;
//				  break;
//			  }
//		}

		triangulation.refine_global(parameter.nbr_global_refinements);	// ... Parameter.prm file

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
					if (    face_center[enums::x] > width/2.
						 && face_center[enums::y] < height/4. )
					{
						cell->set_refine_flag();
						break;
					}
				  }
			}
			triangulation.execute_coarsening_and_refinement();
		}

		// After the refinements the fixed cell is small enough so we can mark it indiviually
		const Point<dim> point_notch (width/2.,0);

		for (typename Triangulation<dim>::active_cell_iterator
			 cell = triangulation.begin_active();
			 cell != triangulation.end(); ++cell)
		{
			for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
			  if (cell->face(face)->boundary_id() == enums::id_boundary_yMinus)
			  {
					// If there is a point right at the \a point_notch we mark this cell as fixed.
					for (unsigned int vertex=0; vertex < GeometryInfo<dim>::vertices_per_cell; ++vertex)
					  if ( cell->vertex(vertex).distance(point_notch) < search_tolerance )
					  {
						  // Found the correct cell, now we select the face that is actually at the boundary and mark it
							for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
							  if (cell->face(face)->at_boundary())
							  {
								  cell->face(face)->set_boundary_id(enums::id_boundary_fix);
								  break;
							  }
					  }
				}
		}

		// Output the triangulation as eps or inp
		 //numEx::output_triangulation( triangulation, enums::output_eps, numEx_name );
	}
}
