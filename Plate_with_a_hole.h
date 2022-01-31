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
	// Name of the numerical example
	 std::string numEx_name = "PlateWithAHole";

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
	void make_constraints ( AffineConstraints<double>  &constraints, const FESystem<dim> &fe, DoFHandler<dim> &dof_handler_ref,
							const bool &apply_dirichlet_bc, const double &current_load_increment,
							const Parameter::GeneralParameters &parameter )
	{
		// clamping on Y0 plane: set x, y and z displacements on x0 plane to zero
		 numEx::BC_apply_fix( enums::id_boundary_yMinus, dof_handler_ref, fe, constraints );

//		// fixed node in x-direction
//		// @todo Currently we fix an entire cell face at x=0 to x=0
//		if (apply_dirichlet_bc == true )
//		{
//			VectorTools::interpolate_boundary_values(
//														dof_handler_ref,
//														enums::id_boundary_xMinus,
//														ZeroFunction<dim> (n_components),
//														constraints,
//														fe.component_mask(x_displacement)
//													);
//		}
//		else	// in the exact same manner
//		{
//			VectorTools::interpolate_boundary_values(
//														dof_handler_ref,
//														enums::id_boundary_xMinus,
//														ZeroFunction<dim> (n_components),
//														constraints,
//														fe.component_mask(x_displacement)
//													);
//		}

		if ( parameter.driver == 2/*Dirichlet*/ )
			numEx::BC_apply( id_boundary_load, loading_direction, current_load_increment, apply_dirichlet_bc, dof_handler_ref, fe, constraints );
	}

	// ToDo-optimize: use existing DII command	void GridGenerator::plate_with_a_hole

	// to see the effects of the inputs (lengths, refinements, etc) consider using the output (.eps, etc) below
	void make_2d_plate_with_hole( Triangulation<2> &tria_2d_out,
										  const double half_length,
										  //const double half_width,
										  const double hole_radius,
										  //const double hole_division_fraction,
										  const Parameter::GeneralParameters &parameter )
	{
		Triangulation<2> tria_2d;
		//const double width =  2.0*half_width;
		const double hole_diameter = 2.0*hole_radius;
		//const double internal_width = hole_diameter + hole_division_fraction*(width - hole_diameter);

		Point<2> centre_2d(0,0);
		const types::manifold_id  	polar_manifold_id = 0;
		const types::manifold_id  	tfi_manifold_id = 1;

		const double height2Width_ratio = 3.;

		GridGenerator::plate_with_a_hole 	( 	tria_2d,
												/*inner radius*/hole_radius,
												/*outer radius*/half_length/height2Width_ratio, // Width of the plate is \a length, the height is 2*length
												/*pad bottom  */half_length*(1.-1./height2Width_ratio),
												/*pad top     */half_length*(1.-1./height2Width_ratio),
												/*pad left    */0.,
												/*pad right   */0.,
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

		// For some reason the flatten_triangulation needs to be done, else we got problems for some types of local refinements
		 GridGenerator::flatten_triangulation(tria_2d,tria_2d_out);

		// Output the triangulation as eps or inp
		 numEx::output_triangulation( tria_2d_out, enums::output_eps, numEx_name );
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

		// Clear all existing boundary ID's
		 numEx::clear_boundary_IDs( triangulation );

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

//					for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_face; ++vertex)
//					 if (std::abs(cell->vertex(vertex)[enums::y] - 0) < search_tolerance)
//					  if (std::abs(cell->vertex(vertex)[enums::x] - parameter.holeRadius) < search_tolerance)
//					  {
//						  // We found the cell that lies at the bottom edge next to the hole (bottom left corner)
//						  cell->set_material_id( enums::tracked_QP );
//						  break;
//					  }
				}
				else if (std::abs(cell->face(face)->center()[1] - ratio_width_To_holeRadius) < search_tolerance)
				{
					cell->face(face)->set_boundary_id(enums::id_boundary_yPlus); // the top edge
				}
				else
				{
					// Be aware that we have to access the vertex as cell->face->vertex
					for (unsigned int vertex=0; vertex < GeometryInfo<dim>::vertices_per_face; ++vertex)
					  if (std::abs(cell->face(face)->vertex(vertex).distance(centre) - holeRadius) < search_tolerance)
					  {
						  cell->face(face)->set_boundary_id(parameters_internal.boundary_id_hole);	// the hole edge
						  break;
					  }
				}

				// Set manifold IDs
				for (unsigned int vertex=0; vertex < GeometryInfo<dim>::vertices_per_face; ++vertex)
				  if (std::abs(cell->vertex(vertex).distance(centre) - holeRadius) < search_tolerance)
				  {
					  cell->face(face)->set_manifold_id(parameters_internal.manifold_id_hole);
					  break;
				  }
			  }
			}
		}

		static SphericalManifold<dim> spherical_manifold (centre);
		triangulation.set_manifold(parameters_internal.manifold_id_hole,spherical_manifold);

		// refine hole boundary
//		for (unsigned int refine_counter=0; refine_counter < 2; refine_counter++)
//		{
//			for (typename Triangulation<dim>::active_cell_iterator
//						 cell = triangulation.begin_active();
//						 cell != triangulation.end(); ++cell)
//			{
//				for (unsigned int face=0; face < GeometryInfo<dim>::faces_per_cell; ++face)
//					if ( cell->face(face)->boundary_id() == parameters_internal.boundary_id_hole )
//					{
//						cell->set_refine_flag();
//						break;
//					}
//			}
//			triangulation.execute_coarsening_and_refinement();
//		}

		// The following does not work?
		// pre-refinement of the damaged area (around y=0)
		for (unsigned int refine_counter=0; refine_counter < parameter.nbr_holeEdge_refinements; refine_counter++)
		{
			for (typename Triangulation<dim>::active_cell_iterator
						 cell = triangulation.begin_active();
						 cell != triangulation.end(); ++cell)
			{
				if ( std::abs( cell->center()[enums::y]) < holeRadius*0.9 )
					cell->set_refine_flag();
			}
			triangulation.execute_coarsening_and_refinement();
		}

//		// The x-fixed faces are set now, to get the smallest possible face
//		for ( typename Triangulation<dim>::active_cell_iterator
//				cell = triangulation.begin_active();
//					cell != triangulation.end(); ++cell )
//		{
//			for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
//			{
//			  if (cell->face(face)->at_boundary() )
//			  {
//				  if ( ( cell->face(face)->center()[enums::y] + hwidth ) < search_tolerance ) // lower face
//				  {
//					for (unsigned int vertex=0; vertex < GeometryInfo<dim>::vertices_per_face; ++vertex)
//					 {
//						 //Project the cell vertex to the XY plane and test the distance from the cylinder axis
//							Point<2> vertex_v = cell->face(face)->vertex(vertex);
//							Point<2> desired_mid_point (0,-hwidth);
//
//							if (std::abs(vertex_v.distance(desired_mid_point)) < search_tolerance)
//							{
//								std::cout << "found l at " << vertex_v << std::endl;
//								cell->face(face)->set_boundary_id(enums::id_boundary_xMinus);
//								break;
//							}
//					 }
//				  }
//					// We also fix some face at the top of the plate to x=0, so it does not drift/shear sidewides
//					// @todo-ensure Do we create stress peaks by doing so?
//					else if ( ( cell->face(face)->center()[enums::y] - hwidth ) < search_tolerance ) // upper face
//					{
//						for (unsigned int vertex=0; vertex < GeometryInfo<dim>::vertices_per_face; ++vertex)
//						 {
//						 //Project the cell vertex to the XY plane and test the distance from the cylinder axis
//							Point<2> vertex_v = cell->face(face)->vertex(vertex);
//							Point<2> desired_mid_point (0,+hwidth);
//
//							if (std::abs(vertex_v.distance(desired_mid_point)) < search_tolerance)
//							{
//								std::cout << "found u at " << vertex_v << std::endl;
//								cell->face(face)->set_boundary_id(enums::id_boundary_xMinus);
//								break;
//							}
//						}
//					}
//				}
//			  }
//		}
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
										   parameter.thickness,
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
//				if ( (std::abs(cell->face(face)->center()[enums::x] - 0.0) < search_tolerance) )
//				{
//					// @todo-optimize: We only look for single cell-face, so we could stop checking this, if we found something (e.g. bool found_xMinus)
//					if ( cell->face(face)->center()[enums::y] < 0 )
//					{
//						for (unsigned int vertex=0; vertex < GeometryInfo<dim>::vertices_per_face; ++vertex)
//						 {
//						 //Project the cell vertex to the XY plane and test the distance from the cylinder axis
//							Point<3> vertex_v = cell->face(face)->vertex(vertex);
//							Point<3> desired_mid_point (0,-hwidth,0);
//
//							if (std::abs(vertex_v.distance(desired_mid_point)) < search_tolerance)
//							{
//								cell->face(face)->set_boundary_id(enums::id_boundary_xMinus);
//								break;
//							}
//						}
//					}
//					// We also fix some face at the top of the plate to x=0, so it does not drift/shear sidewides
//					// @todo-ensure Do we create stress peaks by doing so?
//					else if ( cell->face(face)->center()[enums::y] > 0 )
//					{
//						for (unsigned int vertex=0; vertex < GeometryInfo<dim>::vertices_per_face; ++vertex)
//						 {
//						 //Project the cell vertex to the XY plane and test the distance from the cylinder axis
//							Point<3> vertex_v = cell->face(face)->vertex(vertex);
//							Point<3> desired_mid_point (0,+hwidth,0);
//
//							if (std::abs(vertex_v.distance(desired_mid_point)) < search_tolerance)
//							{
//								cell->face(face)->set_boundary_id(enums::id_boundary_xMinus);
//								break;
//							}
//						}
//					}
//				}
			  // Set boundary IDs
			  if ( cell->face(face)->at_boundary() )
			  {
				// The bottom face (yMinus) is now at the negative hwidth (plate reaches from y = -hwidth to +hwidth)
				if (std::abs(cell->face(face)->center()[enums::y] + hwidth) < search_tolerance)
				{
					cell->face(face)->set_boundary_id(enums::id_boundary_yMinus);
//					// Tracked QP
//					for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_face; ++vertex)
//					{
//					 if (std::abs(cell->vertex(vertex)[enums::y] - 0) < search_tolerance)
//						if (std::abs(cell->vertex(vertex)[enums::z] - 0) < search_tolerance)
//						  if (std::abs(cell->vertex(vertex)[enums::x] - parameter.holeRadius) < search_tolerance)
//						  {
//							  // We found the cell that lies at the bottom edge next to the hole (bottom left corner)
//							  cell->set_material_id( enums::tracked_QP );
//							  break;
//						  }
//					}
				}
				else if (std::abs(cell->face(face)->center()[1] - ratio_width_To_holeRadius) < search_tolerance)
				{
					cell->face(face)->set_boundary_id(enums::id_boundary_yPlus);
				}
				else if (std::abs(cell->face(face)->center()[2] - 0.0) < search_tolerance)
				{
					cell->face(face)->set_boundary_id(enums::id_boundary_zMinus);
				}
				else if (std::abs(cell->face(face)->center()[2] - parameter.thickness) < search_tolerance)
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

		// pre-refinement of the damaged area (around y=0)
		for (unsigned int refine_counter=0; refine_counter < parameter.nbr_holeEdge_refinements; refine_counter++)
		{
			for (typename Triangulation<dim>::active_cell_iterator
						 cell = triangulation.begin_active();
						 cell != triangulation.end(); ++cell)
			{
				if ( std::abs( cell->center()[enums::y]) < holeRadius*0.5 )
					cell->set_refine_flag();
			}
			triangulation.execute_coarsening_and_refinement();
		}

//		// refine hole boundary
//		for (unsigned int refine_counter=0; refine_counter<parameter.nbr_holeEdge_refinements; refine_counter++)
//		{
//			for (typename Triangulation<dim>::active_cell_iterator
//						 cell = triangulation.begin_active();
//						 cell != triangulation.end(); ++cell)
//			{
//				for (unsigned int face=0; face < GeometryInfo<dim>::faces_per_cell; ++face)
//					if ( cell->face(face)->boundary_id() == parameters_internal.boundary_id_hole )
//					{
//						cell->set_refine_flag();
//						break;
//					}
//			}
//			triangulation.execute_coarsening_and_refinement();
//		}

		// The x-fixed faces are set now, to get the smallest possible face
//		for ( typename Triangulation<dim>::active_cell_iterator
//				cell = triangulation.begin_active();
//					cell != triangulation.end(); ++cell )
//		{
//			for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
//			{
//				// Set the INTERNAL (not at the boundary) face to the xminus id
//				// Find the cell face at x=0, containing the node (0,-hwidth,0) pointing in x-direction
//				if ( (std::abs(cell->face(face)->center()[enums::x] - 0.0) < search_tolerance) )
//				{
//					// @todo-optimize: We only look for single cell-face, so we could stop checking this, if we found something (e.g. bool found_xMinus)
//					if ( cell->face(face)->center()[enums::y] < 0 )
//					{
//						for (unsigned int vertex=0; vertex < GeometryInfo<dim>::vertices_per_face; ++vertex)
//						 {
//						 //Project the cell vertex to the XY plane and test the distance from the cylinder axis
//							Point<3> vertex_v = cell->face(face)->vertex(vertex);
//							Point<3> desired_mid_point (0,-hwidth);
//
//							if (std::abs(vertex_v.distance(desired_mid_point)) < search_tolerance)
//							{
//								cell->face(face)->set_boundary_id(enums::id_boundary_xMinus);
//								break;
//							}
//						}
//					}
//					// We also fix some face at the top of the plate to x=0, so it does not drift/shear sidewides
//					// @todo-ensure Do we create stress peaks by doing so?
//					else if ( cell->face(face)->center()[enums::y] > 0 )
//					{
//						for (unsigned int vertex=0; vertex < GeometryInfo<dim>::vertices_per_face; ++vertex)
//						 {
//						 //Project the cell vertex to the XY plane and test the distance from the cylinder axis
//							Point<3> vertex_v = cell->face(face)->vertex(vertex);
//							Point<3> desired_mid_point (0,+hwidth);
//
//							if (std::abs(vertex_v.distance(desired_mid_point)) < search_tolerance)
//							{
//								cell->face(face)->set_boundary_id(enums::id_boundary_xMinus);
//								break;
//							}
//						}
//					}
//				}
//			}
//		}
	}
}
