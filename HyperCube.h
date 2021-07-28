// deal.II headers
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>

// C++ headers for some math operations
// @todo-assure Does this list contain everything that is needed and not more?
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>

// Numerical example helper function (required, can be downloaded from https://github.com/jfriedlein/Numerical_examples_in_dealii)
// also contains enumerators as part of "enums::"
#include "./numEx-helper_fnc.h"

// Contact (if unused, you can remove it together with the \a assemble_contact(*) function at the end)
#include "../contact-rigidBody-dealii/contact-bodies.cc"
#include "../contact-rigidBody-dealii/contact-rigid.cc"

using namespace dealii;

namespace HyperCube
/*
 * A single (or multiple) element(s) with three symmetry constraints, loaded in y-direction, dimensions widthxdim
 * By selecting 1 global refinement and distortion we can create the distorted patch element test in 2D or 3D
 * @todo The refine special enums need to be local values for each numEx
 *
 * CERTIFIED TO STANDARD numExS11 (210104)
 */
{
	// @todo-optimize I guess the following variables will never be destroyed, so the memory remains filled?

	// Name of the numerical example
	// @todo-optimize Can we extract the name of the namespace as string? (however, that would limit the naming freedom)
	 std::string numEx_name = "HyperCube";
	 
	// The loading direction: \n
	// In which coordinate direction the load shall be applied, so x/y/z. We assume the positive axis direction, where
	// changes in the direction (+ -) are possible with positive and negative loads.
	 const unsigned int loading_direction = enums::y; // standard

	// The loaded faces:
	 const enums::enum_boundary_ids id_boundary_load = enums::id_boundary_yPlus;
	 const enums::enum_boundary_ids id_boundary_secondaryLoad = enums::id_boundary_xPlus;

	// Characteristic body dimensions
	 std::vector<double> body_dimensions (5);
	 
	// Evaluation point
	// @note We cannot init the point yet, because we don't have the geometry dimensions and geometry
	 Point<3> eval_point;
	 
	// Evaluation path
	 Point<3> eval_path_start;
	 Point<3> eval_path_end;

	// Some internal parameters
	 struct parameterCollection
	 {
		// We declare the search tolerance as "static constexpr" so we don't need
		// to create an instance of this class just to get access to this variable.
		// If you append this struct, then you should probably create an instance and avoid the "static ..." blabla.
		static constexpr double search_tolerance = 1e-12;
	 };
	 
	// All additional parameters
	// @todo Group them somehow
	 const bool element_distortion = false;
     const bool twitch_Dis8El_cube = false;
	// Sphere size
	// Sphere position
	// ...

	// Contact bodies and geometries
	 // Wall: Pushing down on the cube
//	  Point<2> wall_point_on_plane = Point<2>(0.0,1.0);
//	  const Point<2> wall_normal_unit_vector = Point<2>(0,-1.0);
//	  std::shared_ptr<WallRigid<2>> rigid_wall = std::shared_ptr<WallRigid<2>>(new WallRigid<2>( {wall_point_on_plane,wall_normal_unit_vector,wall_normal_unit_vector} , {} ));

	 // HalfWall: Pushing on the right half of the top face
//	  Point<2> wall_point_on_plane = Point<2>(0.5,1.0);
//	  const Point<2> wall_normal_unit_vector = Point<2>(0,-1.0);
//	  std::shared_ptr<HalfWallRigid<2>> rigid_wall = std::shared_ptr<HalfWallRigid<2>>(new HalfWallRigid<2>( {wall_point_on_plane,wall_normal_unit_vector,wall_normal_unit_vector} , {-1} ));

	 // Sphere
	  Point<2> punch_center = Point<2>(0.0,2.0);
	  const Point<2> punch_loading_vector = Point<2>(0.,-1.);
	  std::shared_ptr<SphereRigid<2>> rigid_wall = std::shared_ptr<SphereRigid<2>>(new SphereRigid<2>( {punch_center,punch_loading_vector,punch_loading_vector}, {1.0,0,0} ));

	  
	/**
	 * Apply the boundary conditions (support and load) on the given AffineConstraints \a constraints. \n
	 * For the HyperCube that are three symmetry constraints on each plane (x=0, y=0, z=0) and the load on the \a id_boundary_load (for Dirichlet).
	 * Alternatively, we apply the rigid body motion for contact.
	 * @todo Change setup (see HyperR) where we define boundary for BC_xMinus and then use if ( BC_xMinus==... ), use the value of BC_xMinus directly
	 */
	template<int dim>
	void make_constraints ( AffineConstraints<double> &constraints, const FESystem<dim> &fe, DoFHandler<dim> &dof_handler_ref,
							const bool &apply_dirichlet_bc, double &current_load_increment, const Parameter::GeneralParameters &parameter )
	{
		// BC on x0 plane
		 numEx::BC_apply( enums::id_boundary_xMinus, enums::x, 0, apply_dirichlet_bc, dof_handler_ref, fe, constraints );

		// BC on y0 plane
		 numEx::BC_apply( enums::id_boundary_yMinus, enums::y, 0, apply_dirichlet_bc, dof_handler_ref, fe, constraints );

		// BC on z0 plane ...
		 if ( dim==3 ) // ... only for 3D
			numEx::BC_apply( enums::id_boundary_zMinus, enums::z, 0, apply_dirichlet_bc, dof_handler_ref, fe, constraints );

//		// Fix the bottom face
//		 numEx::BC_apply_fix( enums::id_boundary_yMinus, dof_handler_ref, fe, constraints );

		// BC for the load ...
		 if ( parameter.driver == enums::Dirichlet )  // ... as Dirichlet only for Dirichlet as driver, alternatively  ...
		 {
			numEx::BC_apply( id_boundary_load, loading_direction, current_load_increment, apply_dirichlet_bc, dof_handler_ref, fe, constraints );
			//numEx::BC_apply( id_boundary_secondaryLoad, enums::x, current_load_increment, apply_dirichlet_bc, dof_handler_ref, fe, constraints );
		 }
		 else if ( parameter.driver == enums::Contact ) // ... as contact
		 {
			if (apply_dirichlet_bc == true )
				rigid_wall->move( current_load_increment );
	 	 }
	}

	// HyperCube grid: 2D and 3D
	template<int dim>
	void make_grid ( Triangulation<dim> &triangulation, const Parameter::GeneralParameters &parameter )
	{
		const double search_tolerance = parameterCollection::search_tolerance;

		const double width = parameter.width;
		
		// Assign the characteristic dimensions of the cube
		 body_dimensions[enums::x] = width;
		 body_dimensions[enums::y] = width;
		 body_dimensions[enums::z] = width;
		
		// Set the evaluation point
		 eval_point[enums::x] = 0;
		 eval_point[enums::y] = width;
		 eval_point[enums::z] = 0;

		// Set the evaluation path points
		 eval_path_start = Point<3> (0,0,0);
		 eval_path_start = eval_point;

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
					AssertThrow(false, ExcMessage( numEx_name+" - make_grid 3D<< Found an unidentified face at the boundary. "
												   "Maybe it slipt through the assignment or that face is simply not needed. "
												   "So either check the implementation or comment this line in the code") );
				}
			  }
		 }

		// Refinement
		 triangulation.refine_global( parameter.nbr_global_refinements );

		// Mark the cell at the centre \a centre by a material id of 1 ...
		 if ( triangulation.n_active_cells()>1 ) // ... only if there is more than one cell
		 {
			bool found_cell=false;
			Point<dim> centre ( 0, 0, 0 );
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
					break; // break out of the entire for(cell) loop
			} // end for(cell)

			AssertThrow(found_cell, ExcMessage(numEx_name+"<< Was not able to identify the cell at the origin(0,0,0). Please recheck the triangulation or adapt the code."));
		 }

		// Distortion
		if ( element_distortion )
		{
			if ( dim==3 )
			{
				// Single distorted element
				if ( triangulation.n_active_cells()==1 )
				{
					std::vector< Point<dim> > points_xyz (2);
					std::vector< Point<dim> > shift_dxdydz (2);
	
					Point<dim> x1y1z1 (1,1,1);
					Point<dim> shift_of_x1y1z1 (-0.5,0,-0.5);
					Point<dim> x0y1z1 (0,1,1);
					Point<dim> shift_of_x0y1z1 (0,0,-0.25);
	
					points_xyz[0] = x1y1z1;
					points_xyz[1] = x0y1z1;
					shift_dxdydz[0] = shift_of_x1y1z1;
					shift_dxdydz[1] = shift_of_x0y1z1;
	
					numEx::shift_vertex_by_vector( triangulation, points_xyz, shift_dxdydz, numEx_name );
				}
				// Distorted 8 elements (Dis8El)
				else if ( triangulation.n_active_cells()==8 )
				{
					std::vector< Point<dim> > points_xyz (16);
					std::vector< Point<dim> > shift_dxdydz (16);
	
					// @todo-optimize The following is super ugly and should be done in a 
					// summarised fashion.
					Point<dim> x0y05z0 (0,0.5,0);
					Point<dim> shift_of_x0y05z0(0,0.1,0);
					points_xyz[0] = x0y05z0;
					shift_dxdydz[0] = shift_of_x0y05z0;
	
					Point<dim> x05y0z0 (0.5,0,0);
					Point<dim> shift_of_x05y0z0(0,-0.125,0);
					points_xyz[1] = x05y0z0;
					shift_dxdydz[1] = shift_of_x05y0z0;
	
					Point<dim> x05y05z0 (0.5,0.5,0);
					Point<dim> shift_of_x05y05z0 (0,-0.3,0);
					points_xyz[2] = x05y05z0;
					shift_dxdydz[2] = shift_of_x05y05z0;
	
					Point<dim> x1y05z0 (1,0.5,0);
					Point<dim> shift_of_x1y05z0(0,-0.2,0);
					points_xyz[3] = x1y05z0;
					shift_dxdydz[3] = shift_of_x1y05z0;
	
					Point<dim> x05y1z0 (0.5,1,0);
					Point<dim> shift_of_x05y1z0(0.15,0,0);
					points_xyz[4] = x05y1z0;
					shift_dxdydz[4] = shift_of_x05y1z0;
	
					Point<dim> x0y05z05 (0,0.5,0.5);
					Point<dim> shift_of_x0y05z05(0,0.2,0);
					points_xyz[5] = x0y05z05;
					shift_dxdydz[5] = shift_of_x0y05z05;
	
					Point<dim> x05y05z05 (0.5,0.5,0.5);
					Point<dim> shift_of_x05y05z05(0.1,-0.1,-0.15);
					points_xyz[6] = x05y05z05;
					shift_dxdydz[6] = shift_of_x05y05z05;
	
					Point<dim> x1y05z05 (1,0.5,0.5);
					Point<dim> shift_of_x1y05z05(0,0.1,-0.1);
					points_xyz[7] = x1y05z05;
					shift_dxdydz[7] = shift_of_x1y05z05;
	
					Point<dim> x0y1z05 (0,1,0.5);
					Point<dim> shift_of_x0y1z05(0,0,0.1);
					points_xyz[8] = x0y1z05;
					shift_dxdydz[8] = shift_of_x0y1z05;
	
					Point<dim> x05y1z05 (0.5,1,0.5);
					Point<dim> shift_of_x05y1z05(0.1,0,-0.05);
					points_xyz[9] = x05y1z05;
					shift_dxdydz[9] = shift_of_x05y1z05;
	
					Point<dim> x1y1z05 (1,1,0.5);
					Point<dim> shift_of_x1y1z05(0,0,0.2);
					points_xyz[10] = x1y1z05;
					shift_dxdydz[10] = shift_of_x1y1z05;
	
					Point<dim> x05y0z1 (0.5,0,1);
					Point<dim> shift_of_x05y0z1(-0.1,0,0);
					points_xyz[11] = x05y0z1;
					shift_dxdydz[11] = shift_of_x05y0z1;
	
					Point<dim> x0y05z1 (0,0.5,1);
					Point<dim> shift_of_x0y05z1(0,-0.125,0);
					points_xyz[12] = x0y05z1;
					shift_dxdydz[12] = shift_of_x0y05z1;
	
					Point<dim> x05y05z1 (0.5,0.5,1);
					Point<dim> shift_of_x05y05z1(0,-0.2,0);
					points_xyz[13] = x05y05z1;
					shift_dxdydz[13] = shift_of_x05y05z1;
	
					Point<dim> x1y05z1 (1,0.5,1);
					Point<dim> shift_of_x1y05z1(0,0.25,0);
					points_xyz[14] = x1y05z1;
					shift_dxdydz[14] = shift_of_x1y05z1;
	
					Point<dim> x05y1z1 (0.5,1,1);
					Point<dim> shift_of_x05y1z1(0.1,0,0);
					points_xyz[15] = x05y1z1;
					shift_dxdydz[15] = shift_of_x05y1z1;
	
					numEx::shift_vertex_by_vector( triangulation, points_xyz, shift_dxdydz, numEx_name );
	
					// Up to now the geometry is still a cube, but the elements are
					// internally distorted. If we also want to twitch the entire cube
					// we can do this in the following
					if ( twitch_Dis8El_cube )
					{
						std::vector< Point<dim> > point_xyz (1);
						std::vector< Point<dim> > shift_point (1);
	
						Point<dim> x1y1z1 (1,1,1);
						Point<dim> shift_of_x1y1z1(0.01,0,0);
						point_xyz[0] = x1y1z1;
						shift_point[0] = shift_of_x1y1z1;
	
						numEx::shift_vertex_by_vector( triangulation, point_xyz, shift_point, numEx_name );
					}
				} // end if(Dis8El)
			} // end if(dim==3)
			else if ( dim==2 )
			{
				for (typename Triangulation<dim>::active_cell_iterator
							 cell = triangulation.begin_active();
							 cell != triangulation.end(); ++cell)
				{
					//cell->vertex(3)[enums::y] -= 0.15;
					Point<dim> corner_rightBottom (width,0);
					for (unsigned int vertex=0; vertex < GeometryInfo<dim>::vertices_per_cell; ++vertex)
					  if ( (cell->vertex(vertex)).distance(corner_rightBottom)<1e-12 )
					  {
						  cell->vertex(vertex)[enums::x] += 0.5;
						  break;
					  }
				}
			} // end if(dim==2)
		} // end if(element_distortion)
		
		// Output the triangulation as eps or inp
		 //numEx::output_triangulation( triangulation, enums::output_eps, numEx_name );
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
							{enums::id_boundary_yPlus},
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
