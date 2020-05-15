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

		const double hole_diameter = 2; // so the hole_radius (a) = 1

		const types::boundary_id boundary_id_minus_X;// = 5;
		const types::boundary_id boundary_id_minus_Y;// = 3;
		const types::boundary_id boundary_id_plus_X; // = 6;
		const types::boundary_id boundary_id_plus_Y; // = 4;

		const types::boundary_id boundary_id_minus_Z;// = 1;
		const types::boundary_id boundary_id_plus_Z =  2;

		const double search_tolerance = 1e-12;

		// only relevant for 3d grid:
		  const unsigned int n_repetitions_z = 2;			// nbr of Unterteilungen in z-direction for 3d meshing
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

			if ( false/*apply sym BC on positive z-face also*/ )
			{
				if (apply_dirichlet_bc == true )
				{
					VectorTools::interpolate_boundary_values(
																dof_handler_ref,
																parameters_internal.boundary_id_plus_Z,
																ZeroFunction<dim> (n_components),
																constraints,
																fe.component_mask(z_displacement)
															);
				}
				else	// in the exact same manner
				{
					VectorTools::interpolate_boundary_values(
																dof_handler_ref,
																parameters_internal.boundary_id_plus_Z,
																ZeroFunction<dim> (n_components),
																constraints,
																fe.component_mask(z_displacement)
															);
				}
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



	// 2D grid
		template <int dim>
		void make_grid( Triangulation<2> &triangulation, const Parameter::GeneralParameters &parameter, std::vector<unsigned int> Vec_boundary_id_collection )
		{
			parameterCollection parameters_internal ( Vec_boundary_id_collection );

			const double search_tolerance = parameters_internal.search_tolerance;

			// ToDo: use the values from the parameter file
			const double width = 1; // equals depth, square bottom area
			const double height = 8;

			// The points that span the brick
			 Point<dim> p1 (0,0);
			 Point<dim> p2 (width, height); // extends in y-direction its height (loaded in y-direction as the othe models)

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
					if (std::abs(cell->face(face)->center()[0] - 0.0) < search_tolerance)
					{
						cell->face(face)->set_boundary_id(parameters_internal.boundary_id_minus_X);
					}
					else if (std::abs(cell->face(face)->center()[0] - width) < search_tolerance)
					{
						cell->face(face)->set_boundary_id(parameters_internal.boundary_id_plus_X);
					}
					else if (std::abs(cell->face(face)->center()[1] - 0.0) < search_tolerance)
					{
						cell->face(face)->set_boundary_id(parameters_internal.boundary_id_minus_Y);
//						cell->set_material_id( enums::tracked_QP );
					}
					else if (std::abs(cell->face(face)->center()[1] - height) < search_tolerance)
					{
						cell->face(face)->set_boundary_id(parameters_internal.boundary_id_plus_Y);
					}
					else
					{
						// There are just eight faces, so if we missed one, something went clearly terribly wrong
						 AssertThrow(false, ExcMessage("BarModel - make_grid 2D: Found an unidentified face at the boundary. Maybe it slipt through the assignment or that face is simply not needed. So either check the implementation or comment this line in the code"));
					}
				  }
			}

			if ( true/*only refine globally*/ )
				triangulation.refine_global(parameter.nbr_global_refinements);	// ... Parameter.prm file
			else // refine in a special manner only cells around the origin
			{
				for ( unsigned int nbr_local_ref=0; nbr_local_ref<parameter.nbr_holeEdge_refinements; nbr_local_ref++ )
				{
					for (typename Triangulation<dim>::active_cell_iterator
								 cell = triangulation.begin_active();
								 cell != triangulation.end(); ++cell)
					{
						for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
						  if (cell->face(face)->at_boundary())
						  {
							// Find all cells that lay in an exemplary damage band with size 1.5 mm from the y=0 face
							if (std::abs(cell->face(face)->center()[1] ) < 1./(nbr_local_ref+1.))
							{
								cell->set_refine_flag();
								break;
							}
						  }
					}
					triangulation.execute_coarsening_and_refinement();
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

				AssertThrow(found_cell, ExcMessage("BarModel: Was not able to identify the cell at the origin(0,0,0). Please recheck the triangulation or adapt the code."));
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
	void make_grid( Triangulation<3> &triangulation, const Parameter::GeneralParameters &parameter, std::vector<unsigned int> Vec_boundary_id_collection )
	{
		parameterCollection parameters_internal ( Vec_boundary_id_collection );

		const double search_tolerance = parameters_internal.search_tolerance;

		// ToDo: use the values from the parameter file
		const double width = 1; // equals depth, square bottom area
		const double height = 8;

		// The points that span the brick
		 Point<dim> p1 (0,0,0);
		 Point<dim> p2 (width, height, width); // extends in y-direction its height (loaded in y-direction as the othe models)

		// vector containing the number of elements in each dimension
		 std::vector<unsigned int> repetitions (3);
		 repetitions[0]=1; // x
		 repetitions[1]=parameter.grid_y_repetitions; // y
		 repetitions[2]=1; // z

		GridGenerator::subdivided_hyper_rectangle
				( 	triangulation,
					repetitions,
					p1,
					p2 );

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
					cell->face(face)->set_boundary_id(parameters_internal.boundary_id_minus_X);
				}
				else if (std::abs(cell->face(face)->center()[0] - width) < search_tolerance)
				{
					cell->face(face)->set_boundary_id(parameters_internal.boundary_id_plus_X);
				}
				else if (std::abs(cell->face(face)->center()[1] - 0.0) < search_tolerance)
				{
					cell->face(face)->set_boundary_id(parameters_internal.boundary_id_minus_Y);
				}
				else if (std::abs(cell->face(face)->center()[1] - height) < search_tolerance)
				{
					cell->face(face)->set_boundary_id(parameters_internal.boundary_id_plus_Y);
				}
				else if (std::abs(cell->face(face)->center()[2] - 0.0) < search_tolerance)
				{
					cell->face(face)->set_boundary_id(parameters_internal.boundary_id_minus_Z);
				}
				else if (std::abs(cell->face(face)->center()[2] - width) < search_tolerance)
				{
					cell->face(face)->set_boundary_id(parameters_internal.boundary_id_plus_Z);
				}
				else
				{
					// There are just eight faces, so if we missed one, something went clearly terribly wrong
					 AssertThrow(false, ExcMessage("BarModel - make_grid 3D: Found an unidentified face at the boundary. Maybe it slipt through the assignment or that face is simply not needed. So either check the implementation or comment this line in the code"));
				}
			  }
		}

		if ( true/*only refine globally*/ )
			triangulation.refine_global(parameter.nbr_global_refinements);	// ... Parameter.prm file
		else // refine in a special manner only cells around the origin
		{
			for ( unsigned int nbr_local_ref=0; nbr_local_ref<parameter.nbr_holeEdge_refinements; nbr_local_ref++ )
			{
				for (typename Triangulation<dim>::active_cell_iterator
							 cell = triangulation.begin_active();
							 cell != triangulation.end(); ++cell)
				{
					for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
					  if (cell->face(face)->at_boundary())
					  {
						// Find all cells that lay in an exemplary damage band with size 1.5 mm from the y=0 face
						if (std::abs(cell->face(face)->center()[1] ) < 1./(nbr_local_ref+1.))
						{
							cell->set_refine_flag();
							break;
						}
					  }
				}
				triangulation.execute_coarsening_and_refinement();
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
					  cell->set_material_id(1);
					  found_cell = true;
					  break;
				  }

				if ( /*only weaken one cell:*/ false  && found_cell )	// ToDo: check: does a break leave only the inner for-loop?
					break;
			}

			AssertThrow(found_cell, ExcMessage("BarModel: Was not able to identify the cell at the origin(0,0,0). Please recheck the triangulation or adapt the code."));
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
