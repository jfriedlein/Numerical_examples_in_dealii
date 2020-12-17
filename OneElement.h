#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>

#include <iostream>
#include <fstream>
#include <cmath>

#include "../contact-rigidBody-dealii/contact-bodies.cc"
#include "../contact-rigidBody-dealii/contact-rigid.cc"

using namespace dealii;

namespace OneElement
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

	// Some internal parameters
	 struct parameterCollection
	 {
		const double search_tolerance = 1e-12;
	 };

	// Wall
//	 Point<2> wall_point_on_plane = Point<2>(0.0,1.0);
//	 const Point<2> wall_normal_unit_vector = Point<2>(0,-1.0);
//	 std::shared_ptr<WallRigid<2>> rigid_wall = std::shared_ptr<WallRigid<2>>(new WallRigid<2>( {wall_point_on_plane,wall_normal_unit_vector,wall_normal_unit_vector} , {} ));

//	 Point<2> wall_point_on_plane = Point<2>(0.5,1.0);
//	 const Point<2> wall_normal_unit_vector = Point<2>(0,-1.0);
//	 std::shared_ptr<HalfWallRigid<2>> rigid_wall = std::shared_ptr<HalfWallRigid<2>>(new HalfWallRigid<2>( {wall_point_on_plane,wall_normal_unit_vector,wall_normal_unit_vector} , {-1} ));


	// Punch
	 Point<2> punch_center = Point<2>(0.0,3.0);
	 const Point<2> punch_loading_vector = Point<2>(0.,-1.);
	 std::shared_ptr<SphereRigid<2>> rigid_wall = std::shared_ptr<SphereRigid<2>>(new SphereRigid<2>( {punch_center,punch_loading_vector,punch_loading_vector}, {1.0,0,0} ));


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

		// on X0 plane
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
		else if ( parameter.driver == enums::Contact )
		{
			if (apply_dirichlet_bc == true )
				rigid_wall->move(current_load_increment);
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

		GridGenerator::hyper_cube(triangulation,0,width);

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
				else
				{
					AssertThrow(false, ExcMessage("OneElement - make_grid 2D<< Found an unidentified face at the boundary. Maybe it slipt through the assignment or that face is simply not needed. So either check the implementation or comment this line in the code"));
				}
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

			AssertThrow(found_cell, ExcMessage("OneElement<< Was not able to identify the cell at the origin(0,0). Please recheck the triangulation or adapt the code."));
		}


		// Include the following two scopes to see directly how the variation of the input parameters changes the geometry of the grid
//		{
//			std::ofstream out ("grid-3d_quarter_plate_merged.eps");
//			GridOut grid_out;
//			grid_out.write_eps (triangulation, out);
//			std::cout << "Grid written to grid-3d_quarter_plate_merged.eps" << std::endl;
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

	template <int dim>
	void assemble_contact
	(
			const typename DoFHandler<dim>::active_cell_iterator &cell,
			const double &penalty_stiffness,
			const FESystem<dim> &fe,
			FEFaceValues<dim> &fe_face_values_ref,
			const FEValuesExtractors::Vector u_fe,
			const unsigned int n_q_points_f,
			const Vector<double> &current_solution,
			std::vector< std::shared_ptr< PointHistory<dim> > > lqph,
			const std::vector<types::global_dof_index> local_dof_indices,
			FullMatrix<double> &cell_matrix,
			Vector<double> &cell_rhs
	)
	{
		// Assemble the contact pair
		assemble_contact(
							cell,
							rigid_wall,
							enums::id_boundary_yPlus,
							penalty_stiffness,
							fe,
							fe_face_values_ref,
							u_fe,
							n_q_points_f,
							current_solution,
							lqph,
							local_dof_indices,
							cell_matrix,
							cell_rhs
						 );
	}
}
