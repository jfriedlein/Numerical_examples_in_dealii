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

namespace Rod
{
	class parameterCollection
	{
	public:
		parameterCollection( std::vector<unsigned int> Vec_boundary_id_collection /*[5,3,6,4,1]*/)
		:
			boundary_id_minus_X(Vec_boundary_id_collection[0]),
			boundary_id_minus_Y(Vec_boundary_id_collection[1]),
			boundary_id_plus_X (Vec_boundary_id_collection[2]),
			boundary_id_plus_Y (Vec_boundary_id_collection[3]),
			boundary_id_minus_Z(Vec_boundary_id_collection[4])
		{
		}

		const types::boundary_id boundary_id_minus_X;// = 5;
		const types::boundary_id boundary_id_minus_Y;// = 3;
		const types::boundary_id boundary_id_plus_X; // = 6;
		const types::boundary_id boundary_id_plus_Y; // = 4;

		const types::boundary_id boundary_id_minus_Z;// = 1;
		const types::boundary_id boundary_id_plus_Z =  2;

		const double search_tolerance = 1e-12;
	};


	template<int dim>
	void make_constraints ( AffineConstraints<double> &constraints, const FESystem<dim> &fe, unsigned int &n_components, DoFHandler<dim> &dof_handler_ref,
							const bool &apply_dirichlet_bc, double &current_load_increment,
							const Parameter::GeneralParameters &parameter, std::vector<unsigned int> Vec_boundary_id_collection )
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

		parameterCollection parameters_internal ( Vec_boundary_id_collection );

		const FEValuesExtractors::Vector displacement(0);
		const FEValuesExtractors::Scalar x_displacement(0);
		const FEValuesExtractors::Scalar y_displacement(1);

		// on X0 plane
		const int boundary_id_X0 = parameters_internal.boundary_id_minus_X;

		if (apply_dirichlet_bc == true )
		{
			VectorTools::interpolate_boundary_values(
														dof_handler_ref,
														boundary_id_X0,
														ZeroFunction<dim> (n_components),
														constraints,
														fe.component_mask(x_displacement)
													);
		}
		else	// in the exact same manner
		{
			VectorTools::interpolate_boundary_values(
														dof_handler_ref,
														boundary_id_X0,
														ZeroFunction<dim> (n_components),
														constraints,
														fe.component_mask(x_displacement)
													);
		}

		// on Y0 edge
		const int boundary_id_Y0 = parameters_internal.boundary_id_minus_Y;

		if (apply_dirichlet_bc == true )
		{
			VectorTools::interpolate_boundary_values(
														dof_handler_ref,
														boundary_id_Y0,
														ZeroFunction<dim> (n_components),
														constraints,
														fe.component_mask(y_displacement)
													);
		}
		else	// in the exact same manner
		{
			VectorTools::interpolate_boundary_values(
														dof_handler_ref,
														boundary_id_Y0,
														ZeroFunction<dim> (n_components),
														constraints,
														fe.component_mask(y_displacement)
													);
		}

		// on Z0 plane
		if ( dim==3 )
		{
			const FEValuesExtractors::Scalar z_displacement(2);
			const int boundary_id_Z0 = parameters_internal.boundary_id_minus_Z;

			if (apply_dirichlet_bc == true )
			{
				VectorTools::interpolate_boundary_values(
															dof_handler_ref,
															boundary_id_Z0,
															ZeroFunction<dim> (n_components),
															constraints,
															fe.component_mask(z_displacement)
														);
			}
			else	// in the exact same manner
			{
				VectorTools::interpolate_boundary_values(
															dof_handler_ref,
															boundary_id_Z0,
															ZeroFunction<dim> (n_components),
															constraints,
															fe.component_mask(z_displacement)
														);
			}
		}

		if ( parameter.driver == 2/*Dirichlet*/ ) // ToDo-optimize: use string in parameterfile denoting "Dirichlet" so the enumerator is not undermined
		{
			const int boundary_id_top = parameters_internal.boundary_id_plus_Y;

			// on top edge
			if (apply_dirichlet_bc == true )
			{
				VectorTools::interpolate_boundary_values(
															dof_handler_ref,
															boundary_id_top,
															ConstantFunction<dim> (current_load_increment/*add only the increment*/, n_components),
															constraints,
															fe.component_mask(y_displacement)
														);
			}
			else
			{
				VectorTools::interpolate_boundary_values(
															dof_handler_ref,
															boundary_id_top,
															ZeroFunction<dim> (n_components),
															constraints,
															fe.component_mask(y_displacement)
														);
			}
		}
	}




// 3d grid
	template <int dim>
	void make_grid( Triangulation<3> &triangulation, const Parameter::GeneralParameters &parameter, std::vector<unsigned int> Vec_boundary_id_collection )
	{
		parameterCollection parameters_internal ( Vec_boundary_id_collection );

		const double search_tolerance = parameters_internal.search_tolerance;

		const double half_length = 53.34/2.;
		const double radius = 6.4135;

		{
			Triangulation<dim> tria_full_cylinder;
			GridGenerator::cylinder(tria_full_cylinder, radius, half_length);

			// Let's first refine the "cylinder" ones, because the initial mesh is a brick
			tria_full_cylinder.refine_global( 1 );

			GridTools::rotate( std::atan(1)*2, 2, tria_full_cylinder);

			std::set<typename Triangulation<dim>::active_cell_iterator > cells_to_remove;
			for (typename Triangulation<dim>::active_cell_iterator
				 cell = tria_full_cylinder.begin_active();
				 cell != tria_full_cylinder.end(); ++cell)
			{
				// Remove all cells that are not in the first quadrant.
				// The 1/8 shall reside in the positive x,y,z quadrant
				if (cell->center()[0] < 0.0 || cell->center()[1] < 0.0 || cell->center()[2] < 0.0 )
					cells_to_remove.insert(cell);
			}

			Assert(cells_to_remove.size() > 0, ExcInternalError());
			Assert(cells_to_remove.size() != tria_full_cylinder.n_active_cells(), ExcInternalError());
			GridGenerator::create_triangulation_with_removed_cells(tria_full_cylinder,cells_to_remove,triangulation);
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
					cell->face(face)->set_boundary_id(parameters_internal.boundary_id_minus_X);
				}
				else if (std::abs(cell->face(face)->center()[1] - 0.0) < search_tolerance)
				{
					cell->face(face)->set_boundary_id(parameters_internal.boundary_id_minus_Y);
				}
				else if (std::abs(cell->face(face)->center()[1] - half_length) < search_tolerance)
				{
					cell->face(face)->set_boundary_id(parameters_internal.boundary_id_plus_Y);
				}
				else if (std::abs(cell->face(face)->center()[2] - 0.0) < search_tolerance)
				{
					cell->face(face)->set_boundary_id(parameters_internal.boundary_id_minus_Z);
				}
			  }
		}

		// Attach a manifold to the curved boundary and refine
		// @note We can only guarantee that the vertices sit on the curve, so we must test with their position instead of the cell centre.
		for (typename Triangulation<dim>::active_cell_iterator
		   cell = triangulation.begin_active();
		   cell != triangulation.end(); ++cell)
		{
		  for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
			if (cell->face(face)->at_boundary())
			  for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_face; ++vertex)
			  {
				 double distance_2d_xz = std::sqrt( cell->vertex(vertex)[0]*cell->vertex(vertex)[0] + cell->vertex(vertex)[2]*cell->vertex(vertex)[2] );
				 if (std::abs(distance_2d_xz - radius) < 1e-12)
				 {
					cell->face(face)->set_all_manifold_ids(10);
					break;
				 }
			  }
		}

		CylindricalManifold<dim> cylindrical_manifold_3d (1); // y-axis
		triangulation.set_manifold(10,cylindrical_manifold_3d);

		triangulation.refine_global(parameter.nbr_global_refinements);	// ... Parameter.prm file

		// Mark the cells at the center for softening (similar to reduction in cross sectional area
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
			}

			AssertThrow(found_cell, ExcMessage("Rod: Was not able to identify the cell at the origin(0,0,0). Please recheck the triangulation or adapt the code."));
		}

		// add some local refinements
		for (unsigned int refine_counter=0; refine_counter<parameter.nbr_holeEdge_refinements; refine_counter++)
		{
			//  You probably cannot use this local refinements because the new nodes on the boundary are not compatible to the old ones. Maybe you have to deactivate the cylindrical manifold.
//			if ( parameter.nbr_holeEdge_refinements>0 )
//				triangulation.reset_manifold(10); // Clear manifold

			for (typename Triangulation<dim>::active_cell_iterator
			   cell = triangulation.begin_active();
			   cell != triangulation.end(); ++cell)
			{
				for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
					if (cell->face(face)->at_boundary())
						if ( cell->face(face)->center()[1] <= half_length/(5.*(refine_counter+1)) )
						{
							cell->set_refine_flag();
							break;
						}
			}
			triangulation.execute_coarsening_and_refinement();
		}

		// include the following two scopes to see directly how the variation of the input parameters changes the geometry of the grid
		/*
		{
			std::ofstream out ("grid-3d_quarter_plate_merged.eps");
			GridOut grid_out;
			grid_out.write_eps (triangulation, out);
			std::cout << "Grid written to grid-3d_quarter_plate_merged.eps" << std::endl;
		}
		*/
		/*
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
