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

namespace Miyauchi_shear
/*
 * A Miyauchi specimen for shearing
 *
 * CERTIFIED TO STANDARD numExS07 (200724)
 */
{
	// Name of the numerical example
	 std::string numEx_name = "Miyauchi_shear";

	// The loading direction: \n
	// In which coordinate direction the load shall be applied, so x/y/z.
	 const unsigned int loading_direction = enums::x;

	// Evaluation point
	 Point<3> eval_point;

	// The loaded faces:
	 const enums::enum_boundary_ids id_boundary_load = enums::id_boundary_xPlus2;
	 const enums::enum_boundary_ids id_boundary_secondaryLoad = enums::id_boundary_none;

	// Characteristic body dimensions
	 std::vector<double> body_dimensions (5);

	// Some internal parameters
	 struct parameterCollection
	 {
		const types::manifold_id manifold_id_left_radius = 10;
		const types::manifold_id manifold_id_right_radius = 11;

		const double search_tolerance = 1e-12;
	 };

	template<int dim>
	void make_constraints ( AffineConstraints<double> &constraints, const FESystem<dim> &fe, DoFHandler<dim> &dof_handler_ref,
							const bool &apply_dirichlet_bc, double &current_load_increment,
							const Parameter::GeneralParameters &parameter )
	{
		// on yPlus plane constrain y dofs (symmetry BC)
		 numEx::BC_apply( enums::id_boundary_yPlus, enums::y, 0, apply_dirichlet_bc, dof_handler_ref, fe, constraints );

		// on xPlus1 plane constrain x dofs
		 numEx::BC_apply( enums::id_boundary_xPlus1, enums::x, 0, apply_dirichlet_bc, dof_handler_ref, fe, constraints );

		// on xMinus1 plane constrain x dofs
		 if ( false ) //guide_left_end
			 numEx::BC_apply( enums::id_boundary_xMinus1, enums::x, 0, apply_dirichlet_bc, dof_handler_ref, fe, constraints );

		// BC for the load ...
		 if ( parameter.driver == enums::Dirichlet )  // ... as Dirichlet only for Dirichlet as driver
			numEx::BC_apply( id_boundary_load, loading_direction, current_load_increment, apply_dirichlet_bc, dof_handler_ref, fe, constraints );
	}


	template <int dim>
	void make_half_grid( Triangulation<dim> &triangulation )
	{
		const double widthX = 81.25;
		const double heightY_half = 65./2.;
		const double widthX_notch = 15.;
		const double radius_notch = 3./2.;
		const double width_innerPlate = (widthX_notch/2. + radius_notch);

		// The points that span the left lower brick
		 Point<dim> p1 (-widthX/2.,-heightY_half/2.);
		 Point<dim> p2 (-width_innerPlate, -radius_notch*2.);
		 Triangulation<2> tria_brick_leftlow;
		 GridGenerator::subdivided_hyper_rectangle 	( 	tria_brick_leftlow,
														{8,4},
														p1,
														p2
													);

		// The points that span the left lower match brick
		 Point<dim> p11 (p1[enums::x], p2[enums::y]);
		 Point<dim> p21 (p2[enums::x], -radius_notch);
		 Triangulation<2> tria_brickMatch_leftlow;
		 GridGenerator::subdivided_hyper_rectangle 	( 	tria_brickMatch_leftlow,
														{8,1},
														p11,
														p21
													);

		// The points that span the left upper brick
		 Point<dim> p3 (p1[enums::x],-p1[enums::y]);
		 Point<dim> p4 (p2[enums::x],-p2[enums::y]);
		 Triangulation<2> tria_brick_leftup;
		 GridGenerator::subdivided_hyper_rectangle 	( 	tria_brick_leftup,
														{8,4},
														p3,
														p4
													);

		// The points that span the left lower match brick
		 Point<dim> p31 (p11[enums::x], -p11[enums::y]);
		 Point<dim> p41 (p21[enums::x], -p21[enums::y]);
		 Triangulation<2> tria_brickMatch_leftup;
		 GridGenerator::subdivided_hyper_rectangle 	( 	tria_brickMatch_leftup,
														{8,1},
														p31,
														p41
													);

		// Left half plate with hole
		 Point<dim> centre_left ( -width_innerPlate, 0 );
		 const types::manifold_id  	polar_manifold_id = 0;
		 const types::manifold_id  	tfi_manifold_id = 1;
		 Triangulation<2> tria_plateWithHole_left;
		 // Create the full plate with hole
		  GridGenerator::plate_with_a_hole 	( 	tria_plateWithHole_left,
												/*inner radius*/radius_notch,
												/*outer radius*/2*radius_notch, //width_innerPlate,
												/*pad bottom  */(heightY_half-4.*radius_notch)/2.,
												/*pad top     */(heightY_half-4.*radius_notch)/2.,
												/*pad left    */0.,
												/*pad right   */width_innerPlate-2.*radius_notch,
												centre_left,
												polar_manifold_id,
												tfi_manifold_id
											);

		  // Remove the left part of the plate for the left half plate
			Triangulation<2> tria_HalfPlateWithHole_left;
			std::set<typename Triangulation<dim>::active_cell_iterator > cells_to_remove;
			for (typename Triangulation<2>::active_cell_iterator
				 cell = tria_plateWithHole_left.begin_active();
				 cell != tria_plateWithHole_left.end(); ++cell)
			{
				// Remove all cells that are not in the first quadrant
				if (cell->center()[enums::x] < centre_left[enums::x] )
					cells_to_remove.insert(cell);
			}

			Assert(cells_to_remove.size() > 0, ExcInternalError());
			Assert(cells_to_remove.size() != tria_plateWithHole_left.n_active_cells(), ExcInternalError());
			GridGenerator::create_triangulation_with_removed_cells(tria_plateWithHole_left,cells_to_remove,tria_HalfPlateWithHole_left);


		 // Merge the left parts
		  GridGenerator::merge_triangulations(	{&tria_brick_leftup, &tria_brickMatch_leftup, &tria_HalfPlateWithHole_left, &tria_brickMatch_leftlow, &tria_brick_leftlow},
												 triangulation,
												 1.0e-10,
												 true );
	}


	// 2D grid
	template <int dim>
	void make_grid( Triangulation<2> &triangulation, const Parameter::GeneralParameters &parameter )
	{
		// ToDo-assure: the use the values from the parameter file
		const double widthX = 81.25;
		const double heightY_half = 65./2.;
		const double widthX_notch = 15.;
		const double radius_notch = 3./2.;
		const double width_innerPlate = (widthX_notch/2. + radius_notch);

		body_dimensions[enums::x] = widthX;
		body_dimensions[enums::y] = heightY_half;

//		// Set the evaluation point
//		if ( loading_direction == enums::y )
//		{
//			 eval_point[enums::x] = body_dimensions[enums::x];
//			 eval_point[enums::y] = body_dimensions[enums::y];
//		}
//		else if ( loading_direction == enums::x )
//		{
//			 eval_point[enums::x] = body_dimensions[enums::x];
//			 eval_point[enums::y] = body_dimensions[enums::y]/2.;
//		}

		parameterCollection parameters_internal;

		const double search_tolerance = parameters_internal.search_tolerance;

		Triangulation<dim> tria_left, tria_right;
		make_half_grid( tria_left );
		make_half_grid( tria_right );

		// @todo-optimize Better use function transform instead of rotate
		GridTools::rotate(std::atan(1) * 4./*180 degrees*/, tria_right);

		 // Merge the left and right part
		  GridGenerator::merge_triangulations(	{&tria_left, &tria_right},
												 triangulation,
												 1.0e-10,
												 true );

		if ( true ) // remove the left inner part
		{
			std::set<typename Triangulation<dim>::active_cell_iterator > cells_to_remove;
			for (typename Triangulation<2>::active_cell_iterator
				 cell = triangulation.begin_active();
				 cell != triangulation.end(); ++cell)
			{
				// Remove all cells that are not in the first quadrant
				if ( cell->center()[enums::x] < -width_innerPlate && cell->center()[enums::y] > 0 )
					cells_to_remove.insert(cell);
			}

			Assert(cells_to_remove.size() > 0, ExcInternalError());
			Assert(cells_to_remove.size() != triangulation.n_active_cells(), ExcInternalError());
			GridGenerator::create_triangulation_with_removed_cells(triangulation,cells_to_remove,triangulation);

		}

		// Clear all existing boundary ID's
		 numEx::clear_boundary_IDs( triangulation );

		//Set boundary IDs and and manifolds
		for (typename Triangulation<dim>::active_cell_iterator
			 cell = triangulation.begin_active();
			 cell != triangulation.end(); ++cell)
		{
			for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
			  if (cell->face(face)->at_boundary())
			  {
				// Set boundary IDs
				if ( std::abs(cell->face(face)->center()[enums::x] - body_dimensions[enums::x]/2.) < search_tolerance )
				{
					if ( cell->face(face)->center()[enums::y] > 0 )
					{
						cell->face(face)->set_boundary_id(enums::id_boundary_xPlus2);
					}
					else
					{
						cell->face(face)->set_boundary_id(enums::id_boundary_xPlus1);
					}
				}
				else if (std::abs(cell->face(face)->center()[enums::y] - body_dimensions[enums::y]/2.) < search_tolerance )
				{
					cell->face(face)->set_boundary_id(enums::id_boundary_yPlus);
				}
				else if ( std::abs(cell->face(face)->center()[enums::x] + body_dimensions[enums::x]/2.) < search_tolerance )
				{
					if ( cell->face(face)->center()[enums::y] > 0 )
					{
						cell->face(face)->set_boundary_id(enums::id_boundary_xMinus2);
					}
					else
					{
						cell->face(face)->set_boundary_id(enums::id_boundary_xMinus1);
					}
				}
			  }
		}

		// Attach the notch radius manifolds
		{
			 Point<dim> centre_left  ( -width_innerPlate, 0 );
			 Point<dim> centre_right (  width_innerPlate, 0 );

			for (typename Triangulation<dim>::active_cell_iterator
				 cell = triangulation.begin_active();
				 cell != triangulation.end(); ++cell)
			{
				for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
				  if ( cell->face(face)->at_boundary() && std::abs( cell->face(face)->center()[enums::x] ) < width_innerPlate )
					for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_face; ++vertex)
					{
						if ( std::abs( (cell->face(face)->vertex(vertex)).distance(centre_left) - radius_notch ) < search_tolerance )
						{
							cell->face(face)->set_all_manifold_ids(parameters_internal.manifold_id_left_radius);
						}
						else if ( std::abs( (cell->face(face)->vertex(vertex)).distance(centre_right) - radius_notch ) < search_tolerance )
						{
							cell->face(face)->set_all_manifold_ids(parameters_internal.manifold_id_right_radius);
						}
					}
			}
			// For the upper radius
			 static SphericalManifold<dim> spherical_manifold_left (centre_left);
			 triangulation.set_manifold(parameters_internal.manifold_id_left_radius,spherical_manifold_left);

			// For the lower radius
			 static SphericalManifold<dim> spherical_manifold_right (centre_right);
			 triangulation.set_manifold(parameters_internal.manifold_id_right_radius,spherical_manifold_right);
		}

		triangulation.refine_global(parameter.nbr_global_refinements);	// ... Parameter.prm file

		// local refinements
		 {
			for ( unsigned int nbr_local_ref=0; nbr_local_ref<parameter.nbr_holeEdge_refinements; nbr_local_ref++ )
			{
				for (typename Triangulation<dim>::active_cell_iterator
							 cell = triangulation.begin_active();
							 cell != triangulation.end(); ++cell)
				{
					// Find all cells that lay in an exemplary damage band with size 1.5 mm from the y=0 face
					if ( std::abs( cell->center()[enums::x] ) < width_innerPlate  && std::abs( cell->center()[enums::y] ) < 2*radius_notch )
						cell->set_refine_flag();
				}
				triangulation.execute_coarsening_and_refinement();
			}
		 }

//		 numEx::output_triangulation(triangulation,enums::output_eps,numEx_name);
	}



// 3d grid
	template <int dim>
	void make_grid( Triangulation<3> &triangulation, const Parameter::GeneralParameters &parameter )
	{
		std::cout << "number of active cells: " << triangulation.n_active_cells() << std::endl;
		std::cout << "Chosen numerical example: " << parameter.numExample << std::endl;

		AssertThrow(false, ExcMessage("Miyauchi_shear<< not yet implemented for 3D."));
	}
}
