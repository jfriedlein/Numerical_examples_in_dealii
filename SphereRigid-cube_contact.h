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

#include "../contact-rigidBody-dealii/contact-bodies.cc"
#include "../contact-rigidBody-dealii/contact-rigid.cc"


using namespace dealii;

// Append the \a enums namespace
namespace enums
{
	/**
	 * Clamping types
	 */
	 enum enum_SRC_clamping
	 {
		SRC_clamped_free = 0,   //!< clamped_free: left face is fixed in all directions, right face is unconstrained
		SRC_clamped_sliding = 1,//!< clamped_sliding: left face is fixed in all directions, right face can slide in y-direction
		SRC_clamped_clamped = 2 //!< clamped_clamped: @todo Not yet implemented
	 };

	 enum enum_loading
	 {
		 loading_prescribed = 0,
		 loading_by_contact = 1
	 };
}


namespace SphereRigid_Cube
/*
 * A beam along x with one symmetry constraint in z, loaded in y-direction
 *
 * STILL UNCERTIFIED (requires update of standard to incorporate contact numEx)
 */
{
	// The loading direction: \n
	// In which coordinate direction the load shall be applied, so x/y/z.
	 const unsigned int loading_direction = enums::y;

	// The loaded faces:
	 const enums::enum_boundary_ids id_boundary_load = enums::id_boundary_yPlus;
	 const enums::enum_boundary_ids id_boundary_secondaryLoad = enums::id_boundary_none;
	 const enums::enum_boundary_ids id_contact = enums::id_boundary_yPlus;

	// Characteristic body dimensions
	 std::vector<double> body_dimensions (5);

	// Some internal parameters
	 struct parameterCollection
	 {
		const double search_tolerance = 1e-12;
	 };

	// USER PARAMETER
	 const unsigned int beam_type = enums::SRC_clamped_sliding;
	 const unsigned int loading_type = enums::loading_by_contact;

	// @note We cannot add contact faces that point away from the punch and are on the other side of the body
	 // @todo add "dim" instead of "2"
	// Punch
	 Point<2> punch_center = Point<2>(0.0,2.04);
	 const Point<2> punch_loading_vector = Point<2>(0.,-1.);
	 const double punch_radius = 1.;
	 std::shared_ptr<SphereRigid<2>> rigid_punch = std::shared_ptr<SphereRigid<2>>(new SphereRigid<2>( {punch_center,punch_loading_vector,punch_loading_vector} , {punch_radius} ));
//	 Point<2> punch_center = Point<2>(0.0,1.01);
//	 const Point<2> punch_loading_vector = Point<2>(0.,-1.);
//	 const double punch_radius = 0.75;
//	 const double angle_cone = 20; // 20° cone angle to stabilise the computation
//	 std::shared_ptr<RoundedPunch<2>> rigid_punch = std::shared_ptr<RoundedPunch<2>>(new RoundedPunch<2>( {punch_center,punch_loading_vector,punch_loading_vector} , {punch_radius,angle_cone} ));

	// Wall on the bottom
	 Point<2> wall_point_on_plane = Point<2>(0.0,-0.5);
	 const Point<2> wall_normal_unit_vector = Point<2>(0,1.0);
	 std::shared_ptr<WallRigid<2>> rigid_wall = std::shared_ptr<WallRigid<2>>(new WallRigid<2>( {wall_point_on_plane,wall_normal_unit_vector,wall_normal_unit_vector} , {} ));

	// Die on the left
	 Point<2> die_point_on_plane = Point<2>(1.5,-0.1);
	 const Point<2> die_normal_unit_vector = Point<2>(-0.1,.5);
	 std::shared_ptr<HalfWallRigid<2>> rigid_die = std::shared_ptr<HalfWallRigid<2>>(new HalfWallRigid<2>( {die_point_on_plane,die_normal_unit_vector,die_normal_unit_vector} , {} ));

	template<int dim>
	void make_constraints ( AffineConstraints<double> &constraints, const FESystem<dim> &fe, unsigned int &n_components, DoFHandler<dim> &dof_handler_ref,
							const bool &apply_dirichlet_bc, double &current_load_increment )
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

		// set x, y and z displacements on x0 plane to zero
		// Only set the displacement dofs to zero not all dofs (hence e.g. not the damage dofs)
		if (apply_dirichlet_bc == true )
		{
			VectorTools::interpolate_boundary_values(
														dof_handler_ref,
														enums::id_boundary_xPlus,
														ZeroFunction<dim> (n_components),
														constraints,
														fe.component_mask(displacement)
													);
		}
		else	// in the exact same manner
		{
			VectorTools::interpolate_boundary_values(
														dof_handler_ref,
														enums::id_boundary_xPlus,
														ZeroFunction<dim> (n_components),
														constraints,
														fe.component_mask(displacement)
													);
		}

		if (apply_dirichlet_bc == true )
		{
			if ( loading_type==enums::loading_prescribed )
				VectorTools::interpolate_boundary_values(
															dof_handler_ref,
															id_boundary_load,
															ConstantFunction<dim> (current_load_increment/*add only the increment*/, n_components),
															constraints,
															fe.component_mask(y_displacement)
														);
			else if ( loading_type==enums::loading_by_contact )
				rigid_punch->move(current_load_increment);

			if ( beam_type==enums::SRC_clamped_sliding )
				VectorTools::interpolate_boundary_values(
															dof_handler_ref,
															enums::id_boundary_xMinus,
															ZeroFunction<dim> (n_components),
															constraints,
															fe.component_mask(x_displacement)
														);
		}
		else
		{
			if ( loading_type==enums::loading_prescribed)
				VectorTools::interpolate_boundary_values(
															dof_handler_ref,
															id_boundary_load,
															ZeroFunction<dim> (n_components),
															constraints,
															fe.component_mask(y_displacement)
														);
			if ( beam_type==enums::SRC_clamped_sliding )
				VectorTools::interpolate_boundary_values(
															dof_handler_ref,
															enums::id_boundary_xMinus,
															ZeroFunction<dim> (n_components),
															constraints,
															fe.component_mask(x_displacement)
														);
		}
	}


	// 2D grid
		template <int dim>
		void make_grid( Triangulation<2> &triangulation )
		{
			//create a vector of begin and end positions of the blocks
			std::vector<double> x_position{0.,2.};
			std::vector<double> y_position{0.,1.};

			body_dimensions[enums::x] = x_position[1];
			body_dimensions[enums::x] = y_position[1];

			// Only one block
		//	GridGenerator::hyper_cube (triangulation,x_position[0],x_position[1]);
			Point<dim> p1 (x_position[0],y_position[0]);
			Point<dim> p2 (x_position[1],y_position[1]);

			//GridGenerator::subdivided_hyper_rectangle( triangulation, {1*parameter.nbr_global_refinements,4*parameter.nbr_global_refinements}, p1, p2 );
			GridGenerator::subdivided_hyper_rectangle( triangulation, {12*2,12}, p1, p2 );

			// Add the punch as a dummy body
			// That is a bad way of doing it, because we add many dofs, we loop over in the assembly (expensive, useless)
//			{
//				Triangulation<dim> tria_punch;
//
//				// rectangle
//		//		{
//		//			Point<dim> left (0.4,1.01);
//		//			Point<dim> right (0.6,1.1);
//		//			GridGenerator::hyper_rectangle (tria_punch,left,right);
//		//		}
//
//				// sphere
//				{
//					GridGenerator::hyper_ball( tria_punch , punch_center, punch_radius );
//				}
//
//				for (auto cell: tria_punch.active_cell_iterators())
//				{
//					cell->set_material_id(99);
//				}
//
//				GridGenerator::merge_triangulations(triangulation,tria_punch,triangulation);
//
//				for (auto cell: triangulation.active_cell_iterators())
//				{
//					if ( cell->material_id()==99 )
//					{
//						for(unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
//							if ( cell->face(face)->at_boundary() )
//								cell->face(face)->set_all_manifold_ids(1);
//					}
//				}
//				 static SphericalManifold<dim> spherical_manifold_upper (punch_center);
//				 triangulation.set_manifold(1,spherical_manifold_upper);
//
//				for ( unsigned int nbr_local_ref=0; nbr_local_ref < 4; nbr_local_ref++ )
//				{
//					for (auto cell: triangulation.active_cell_iterators())
//						if ( cell->material_id()==99 )
//							cell->set_refine_flag();
//					triangulation.execute_coarsening_and_refinement();
//				}
//			}

			// global refinement
			 //triangulation.refine_global(parameter.nbr_global_refinements);

			//set boundary ids
			for (auto cell: triangulation.active_cell_iterators())
			{
				for(unsigned int j=0; j<GeometryInfo<dim>::faces_per_cell; ++j)
				{
					//contact surface: rigid punch side - block
					if ( cell->face(j)->at_boundary()
							&& ( std::abs(cell->face(j)->center()[1] - y_position[1]) < 1e-10 )) // top
						cell->face(j)->set_boundary_id(enums::id_boundary_yPlus);
					else if ( cell->face(j)->at_boundary()
							&& ( std::abs(cell->face(j)->center()[0] - x_position[0]) < 1e-10 )) // left
						cell->face(j)->set_boundary_id(enums::id_boundary_xMinus);
					else if ( cell->face(j)->at_boundary()
							&& ( std::abs(cell->face(j)->center()[0] - x_position[1]) < 1e-10 )) // right
						cell->face(j)->set_boundary_id(enums::id_boundary_xPlus);
					//fixed surface
					else if (cell->face(j)->at_boundary() && std::abs(cell->face(j)->center()[1] - y_position[0]) < 1e-10)
						cell->face(j)->set_boundary_id(enums::id_boundary_yMinus);
				}
				if ( cell->material_id()==99 )
					for(unsigned int j=0; j<GeometryInfo<dim>::faces_per_cell; ++j)
						cell->face(j)->set_boundary_id(enums::id_body_dummy);
			}

//
//			// include the following two scopes to see directly how the variation of the input parameters changes the geometry of the grid
//			/*
//			{
//				std::ofstream out ("grid-2d_quarter_plate_merged.eps");
//				GridOut grid_out;
//				GridOutFlags::Eps<2> eps_flags;
//				eps_flags.line_width = 0.1;
//				grid_out.set_flags (eps_flags);
//				grid_out.write_eps (triangulation, out);
//				std::cout << "Grid written to grid-2d_quarter_plate_merged.eps" << std::endl;
//				std::cout << "nElem: " << triangulation.n_active_cells() << std::endl;
//				AssertThrow(false,ExcMessage("ddd"));
//			}
//
//			{
//				std::ofstream out_ucd("Grid-2d_quarter_plate_merged.inp");
//				GridOut grid_out;
//				GridOutFlags::Ucd ucd_flags(true,true,true);
//				grid_out.set_flags(ucd_flags);
//				grid_out.write_ucd(triangulation, out_ucd);
//				std::cout<<"Mesh written to Grid-2d_quarter_plate_merged.inp "<<std::endl;
//			}
//			*/
		}



// 3d grid
	template <int dim>
	void make_grid( Triangulation<3> &triangulation, const Parameter::GeneralParameters &parameter )
	{

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
			const std::vector<types::global_dof_index> local_dof_indices,
			FullMatrix<double> &cell_matrix,
			Vector<double> &cell_rhs
	)
	{
		// Assemble the bottom wall
		assemble_contact(
							cell,
							rigid_wall,
							enums::id_boundary_yMinus,
							penalty_stiffness,
							fe,
							fe_face_values_ref,
							u_fe,
							n_q_points_f,
							current_solution,
							local_dof_indices,
							cell_matrix,
							cell_rhs
						 );

		// Assemble the left die
//		assemble_contact(
//							cell,
//							rigid_die,
//							enums::id_boundary_yMinus,
//							penalty_stiffness,
//							fe,
//							fe_face_values_ref,
//							u_fe,
//							n_q_points_f,
//							current_solution,
//							local_dof_indices,
//							cell_matrix,
//							cell_rhs
//						 );

		// Assemble the punch
		assemble_contact(
							cell,
							rigid_punch,
							SphereRigid_Cube::id_contact,
							penalty_stiffness,
							fe,
							fe_face_values_ref,
							u_fe,
							n_q_points_f,
							current_solution,
							local_dof_indices,
							cell_matrix,
							cell_rhs
						 );
	}
}
