#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/base/table.h>
#include <deal.II/base/function_lib.h>

#include <iostream>
#include <fstream>
#include <cmath>

using namespace dealii;

namespace ThreePointBeam
{
	class parameterCollection
	{
	public:
		const types::boundary_id boundary_id_support1;// = 5;
		const types::boundary_id boundary_id_support2;// = 3;
		const types::boundary_id boundary_id_plus_X; // = 6;
		const types::boundary_id boundary_id_load_surface; // = 4;

		const types::boundary_id boundary_id_minus_Z;// = 1;
		const types::boundary_id boundary_id_plus_Z =  2;
		const types::boundary_id boundary_id_hole = 10;
		const types::manifold_id manifold_id_hole = 10;


		parameterCollection( std::vector<unsigned int> Vec_boundary_id_collection /*[5,3,6,4,1]*/)
		:
			boundary_id_support1(Vec_boundary_id_collection[0]), // has been minus_X
			boundary_id_support2(Vec_boundary_id_collection[1]), // has been minus_Y
			boundary_id_plus_X (Vec_boundary_id_collection[2]), // unused for this model, could be used for an additional symmetry constraint
			boundary_id_load_surface (Vec_boundary_id_collection[3]), // has been Plus_Y
			boundary_id_minus_Z(Vec_boundary_id_collection[4]) // for the symmetry
		{
		}

		const double search_tolerance = 1e-12;

		// only relevant for 3d grid:
		  const unsigned int n_repetitions_z = 2;			// nbr of Unterteilungen in z-direction for 3d meshing
	};


// 2D grid
	template <int dim>
	void make_grid( Triangulation<2> &triangulation, const Parameter::GeneralParameters &parameter, std::vector<unsigned int> Vec_boundary_id_collection )
	{
		AssertThrow( false, ExcMessage("The 3 point beam mesh has not yet been implemented for 2D. Use either 3D or simply implement it yourself."));

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


//	template <int dim>
//	class poly_BC : public Functions::Polynomial<dim>
//	{
//	public:
//		poly_BC ( Table<2,double> exponents, std::vector< double > coefficients, const unsigned int n_components )
//	: Functions::Polynomial<dim> (exponents,coefficients), exponents(exponents), coefficients(coefficients), n_components(n_components) {}
//		virtual ~poly_BC() {};
//
//		Table<2,double> exponents;
//		std::vector< double > coefficients;
//		const unsigned int n_components;
//
//		virtual double value (const Point<dim>   &p, const unsigned int component ) const override final
//		{
//			return this->value(p,0);
//		}
//	};


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

		// on support1 plane: constraint in all directions
		const int boundary_id_X0 = parameters_internal.boundary_id_support1;

		if (apply_dirichlet_bc == true )
		{
			VectorTools::interpolate_boundary_values(
														dof_handler_ref,
														boundary_id_X0,
														ZeroFunction<dim> (n_components),
														constraints
													);
		}
		else	// in the exact same manner
		{
			VectorTools::interpolate_boundary_values(
														dof_handler_ref,
														boundary_id_X0,
														ZeroFunction<dim> (n_components),
														constraints
													);
		}


		// on support 2: Loslager, only constraint in y-direction, rest is left free
		const int boundary_id_Y0 = parameters_internal.boundary_id_support2;

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

		// on Z0 plane; symmetry constraint, so no displacement in z-direction
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
			const int boundary_id_top = parameters_internal.boundary_id_load_surface;

			// on top/loaded edge
			if (apply_dirichlet_bc == true )
			{
//				if ( parameter.numExample!=3 )
				{
					VectorTools::interpolate_boundary_values(
																dof_handler_ref,
																boundary_id_top,
																ConstantFunction<dim> (current_load_increment/*add only the increment*/, n_components),
																constraints,
																fe.component_mask(y_displacement)
															);
				}
//				else
//				{
//					Table<2,double> exponents (3,3);
//					// x
//					exponents[0][0]=0; // constant term x^0
//					exponents[1][0]=0; // linear term
//					exponents[2][0]=2; // quadratic term
//					// y
//					exponents[0][1]=0;
//					exponents[1][1]=0;
//					exponents[2][1]=0;
//					// z
//					exponents[0][2]=0;
//					exponents[1][2]=0;
//					exponents[2][2]=0;
//
//					// TESTING
////					current_load_increment = 1.;
////
//					std::vector< double > coefficients (3);
//					coefficients[0]=-current_load_increment;
//					coefficients[1]=0;
//					coefficients[2]=current_load_increment;
////
////					Functions::Polynomial<dim> poly_BC (exponents,coefficients);
////
////					Point<dim> p(0.5,0.5,0.25);
////					double value_at_p = poly_BC.value(p,0);
////
////					std::cout << "value at p=" << value_at_p << std::endl;
//
//					VectorTools::interpolate_boundary_values(
//																dof_handler_ref,
//																boundary_id_top,
//																poly_BC<dim> (exponents,coefficients,n_components),
//																constraints,
//																fe.component_mask(y_displacement)
//															);



//					// Polynomial: i * x² - i -> zero boundary condition (-1,1) and prescribed negative displacement i @x=0
//					Table<2,double> exponents (2,3);
//					// x
//					exponents[0][0]=0; // constant term x^0
//					exponents[1][0]=2; // quadratic term
//					// rest is left as zero
//
//					// TESTING
////					current_load_increment = 1.;
////
//					std::vector< double > coefficients (2);
//					coefficients[0]=-current_load_increment;
//					coefficients[1]=current_load_increment;
////
//					Functions::Polynomial<dim> poly_BC (exponents,coefficients);
//
//					// TESTING
////					Point<dim> p(0.5,0.5,0.25);
////					double value_at_p = poly_BC.value(p,0);
////					std::cout << "value at p=" << value_at_p << std::endl;
//
//					VectorTools::interpolate_boundary_values(
//																dof_handler_ref,
//																boundary_id_top,
//																poly_BC,
//																constraints,
//																fe.component_mask(y_displacement)
//															);
//				}
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




// 3d grid:
	template <int dim>
	void make_grid( Triangulation<3> &triangulation, const Parameter::GeneralParameters &parameter, std::vector<unsigned int> Vec_boundary_id_collection )
	{
		parameterCollection parameters_internal ( Vec_boundary_id_collection );

		const double search_tolerance = parameters_internal.search_tolerance;

		// set the dimensions of the three point beam
		const double length = parameter.width;
		const double height = parameter.height;
		const double thickness = parameter.thickness/2.;
		const double notchWidth = parameter.notchWidth; // whereas the height is set to 1/4 of the height due to the mesh refinement and deletion of certain cells to create the notch

		const double left_point = -length/3.;
		const double right_point = length + left_point;

		// Define the points spanning the three bodies
		Point<dim> pL_left( left_point, 0, 0);
		Point<dim> pR_right(  right_point, height, thickness);
		Point<dim> pN_right( 4*notchWidth, 0,  0);	// the size of the notch body is twice a wide as the notch because only the innermost cells after the refinement will be deleted to create the notch
		Point<dim> pN_left( -4*notchWidth, height, thickness);
		
		std::vector<unsigned int> repetitions_notch (3);
		repetitions_notch[0]=16;
		repetitions_notch[1]=4;
		repetitions_notch[2]=2;

		std::vector<unsigned int> repetitions_rest (3);
		repetitions_rest[0]=4;
		repetitions_rest[1]=4;
		repetitions_rest[2]=2;

		Triangulation<dim> body_N, body_R, body_L;
		// Notch body N
		GridGenerator::subdivided_hyper_rectangle( body_N, repetitions_notch, pN_right, pN_left );
		// right hex
		GridGenerator::subdivided_hyper_rectangle( body_R, repetitions_rest, pN_right, pR_right );
		// left hex
		GridGenerator::subdivided_hyper_rectangle( body_L, repetitions_rest, pN_left, pL_left );
		
		// merge the three bodies
		GridGenerator::merge_triangulations( body_N, body_R, triangulation );
		GridGenerator::merge_triangulations( triangulation, body_L, triangulation );

		// remove the cell(s) where the notches shall be
		std::set<typename Triangulation<dim>::active_cell_iterator > cells_to_remove;

		for (typename Triangulation<dim>::active_cell_iterator
			 cell = triangulation.begin_active();
			 cell != triangulation.end(); ++cell)
		{
			for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
			  if (cell->face(face)->at_boundary())
			  {
				// on the bottom face for positive x
				if (std::abs(cell->face(face)->center()[1] - 0.0) < search_tolerance)
				{
					for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_face; ++vertex)
					  if ( std::abs(cell->vertex(vertex)[0] - 2*notchWidth) < search_tolerance) // found a cell which lays next to the crack plane
					  {
						  cells_to_remove.insert(cell);
						  break;
					  }

					if ( std::abs(cell->face(face)->center()[0] - 2*notchWidth) < search_tolerance ) // found a cell which lays right on the crack plane
						cells_to_remove.insert(cell);
				}
				// second notch on the top face for negative x
				else if (std::abs(cell->face(face)->center()[1] - height) < search_tolerance)
				{
					for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_face; ++vertex)
					  if ( std::abs(cell->vertex(vertex)[0] + 2*notchWidth) < search_tolerance) // found a cell which lays next to the crack plane
					  {
						  cells_to_remove.insert(cell);
						  break;
					  }

					if ( std::abs(cell->face(face)->center()[0] + 2*notchWidth) < search_tolerance ) // found a cell which lays right on the crack plane
						cells_to_remove.insert(cell);
				}
			  }
		}

		Assert(cells_to_remove.size() > 0, ExcInternalError());
		Assert(cells_to_remove.size() != triangulation.n_active_cells(), ExcInternalError());
		GridGenerator::create_triangulation_with_removed_cells(triangulation,cells_to_remove,triangulation);

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
				if (std::abs(cell->face(face)->center()[1] - 0.0) < search_tolerance)
				{
					// -> on the bottom face
					for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_face; ++vertex)
					  if (std::abs(cell->vertex(vertex)[0] - left_point ) < search_tolerance) // found a cell where the support1 is
					  {
						  cell->face(face)->set_boundary_id(parameters_internal.boundary_id_support1);
						  break;
					  }
					  else if (std::abs(cell->vertex(vertex)[0] - right_point ) < search_tolerance) // found a cell where the support2 is
					  {
						  cell->face(face)->set_boundary_id(parameters_internal.boundary_id_support2);
						  break;
					  }
				}
				else if (std::abs(cell->face(face)->center()[1] - height) < search_tolerance)
				{
					const double sym_width_of_load_surface = notchWidth/2;
					// -> on the top face
					for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_face; ++vertex)
					  if ( (cell->vertex(vertex)[0] >= - sym_width_of_load_surface)  &&
						   (cell->vertex(vertex)[0] <= sym_width_of_load_surface ) )
						  // found a cell which lays next to the crack plane or in the search area
					  {
						  cell->face(face)->set_boundary_id(parameters_internal.boundary_id_load_surface);
						  break;
					  }

//					if ( ( cell->face(face)->center()[0] >= - sym_width_of_load_surface) &&
//						  (cell->face(face)->center()[0] <= sym_width_of_load_surface ) )
//						// found a cell which lays right on the crack plane or in the search area
//						cell->face(face)->set_boundary_id(parameters_internal.boundary_id_load_surface);
				}
				else if (std::abs(cell->face(face)->center()[2] - 0.0) < search_tolerance)
				{
					cell->face(face)->set_boundary_id(parameters_internal.boundary_id_minus_Z);
				}
				else
				{
					//AssertThrow(false, ExcMessage("3PointBeam - make_grid 3D: Found an unidentified face at the boundary. Maybe it slipt through the assignment or that face is simply not needed. So either check the implementation or comment this line in the code"));
				}

				// if the cell is at the right or left edge, we mark it with the material_id 3 to be able to switch the damage off for these later on
				if ( std::abs(cell->face(face)->center()[0] - right_point) < ( length/4. + search_tolerance)
						||
						std::abs(cell->face(face)->center()[0] - left_point) < ( length/4. + search_tolerance) ) // found a cell where the left edge is
				{
					  cell->set_material_id(3);
				}
			  }
		}

		// Putting a cylindrical manifold on the notch faces does not seem to do anything. More importantly we would need a half cylinder to only smooth the inner edges
//		const Point<dim> centre_3d (0, height-0.5*notchWidth, 0);
//		for (typename Triangulation<dim>::active_cell_iterator
//		   cell = triangulation.begin_active();
//		   cell != triangulation.end(); ++cell)
//		{
//		  for (unsigned int face=0; face<GeometryInfo<2>::faces_per_cell; ++face)
//			if (cell->face(face)->at_boundary())
//			  for (unsigned int vertex=0; vertex<GeometryInfo<2>::vertices_per_face; ++vertex)
//				if (  std::abs( (cell->vertex(vertex))[1] - height/4. ) < 1e-12 && std::abs( (cell->vertex(vertex))[0] ) <= notchWidth/2.)
////						std::a bs(cell->vertex(vertex).distance(centre_2d) - hole_diameter/2.0) < 1e-12)
//				 {
//					std::cout << cell->vertex(vertex) << std::endl;
//					cell->face(face)->set_manifold_id(10);
//					break;
//				 }
//		}
//
//		Tensor<1,dim> axis_dir;
//		axis_dir[2]=1;  // axis point in the z-direction
//		CylindricalManifold<dim> cylindrical_manifold_3d (axis_dir, centre_3d);
//		triangulation.set_manifold(10,cylindrical_manifold_3d);

		triangulation.refine_global(parameter.nbr_global_refinements);	// ... Parameter.prm file


		// Refine the cells around the center (x=0) of the beam
		if ( parameter.nbr_holeEdge_refinements>0 )
		{
			for ( unsigned int refinements=0; refinements < parameter.nbr_holeEdge_refinements; ++refinements )
			{
				for (typename Triangulation<dim>::active_cell_iterator
				   cell = triangulation.begin_active();
				   cell != triangulation.end(); ++cell)
				{
				  for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
					  for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_face; ++vertex)
						if ( std::abs( cell->vertex(vertex)[0] ) <= 2.*notchWidth )
						 {
							cell->set_refine_flag();
							break;
						 }
				}

				triangulation.execute_coarsening_and_refinement();
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



	// 3d grid:
	template <int dim>
	void make_grid( bool &deactivated_true, Triangulation<3> &triangulation, const Parameter::GeneralParameters &parameter, std::vector<unsigned int> Vec_boundary_id_collection )
	{
		parameterCollection parameters_internal ( Vec_boundary_id_collection );

		const double search_tolerance = parameters_internal.search_tolerance;

		// set the dimensions of the three point beam
		const double length = parameter.width;
		const double height = parameter.height;
		const double thickness = parameter.thickness;
		const double notchWidth = parameter.notchWidth; // whereas the height is set to 1/4 of the height due to the mesh refinement and deletion of certain cells to create the notch
		const Point<dim> centre (0,0,thickness/4.);

		const double left_point = length/3. ; //length/2.;
		const double right_point = length*2./3.; //length/2.;

		Triangulation<3> tria_half_plate_hole;
		{
		 	const types::manifold_id  	polar_manifold_id = 0;
		 	const types::manifold_id  	tfi_manifold_id = 1;
			Triangulation<3> tria_plate_hole;
//			GridGenerator::plate_with_a_hole 	( 	tria_plate_hole ,
//													notchWidth/2.,
//													notchWidth,
//													1.,
//													height - notchWidth,
//													left_point - notchWidth,
//													right_point - notchWidth,
//													centre,
//													polar_manifold_id,
//													tfi_manifold_id,
//													thickness/2.,
//													parameter.nbr_elementsInZ
//												);

			std::set<typename Triangulation<3>::active_cell_iterator > cells_to_remove;
			for (typename Triangulation<3>::active_cell_iterator
				 cell = tria_plate_hole.begin_active();
				 cell != tria_plate_hole.end(); ++cell)
			{
				//Remove all cells that are not in the first or second quadrant (corresponds to cells that have a negative y-coord)
				if ( cell->center()[1] < 0.0)
					cells_to_remove.insert(cell);
			}

			Assert(cells_to_remove.size() > 0, ExcInternalError());
			Assert(cells_to_remove.size() != tria_plate_hole.n_active_cells(), ExcInternalError());
			GridGenerator::create_triangulation_with_removed_cells(tria_plate_hole,cells_to_remove,triangulation);//tria_half_plate_hole);
		}


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
				if (std::abs(cell->face(face)->center()[1] - 0.0) < search_tolerance)
				{
					// -> on the bottom face
					for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_face; ++vertex)
					  if (std::abs(cell->vertex(vertex)[0] + left_point) < search_tolerance) // found a cell where the support1 is
					  {
						  cell->face(face)->set_boundary_id(parameters_internal.boundary_id_support1);
						  break;
					  }
					  else if (std::abs(cell->vertex(vertex)[0] - right_point) < search_tolerance) // found a cell where the support2 is
					  {
						  cell->face(face)->set_boundary_id(parameters_internal.boundary_id_support2);
						  break;
					  }
				}
				else if (std::abs(cell->face(face)->center()[1] - height) < search_tolerance)
				{
					// -> on the top face
					for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_face; ++vertex)
					  if (std::abs(cell->vertex(vertex)[0] - notchWidth/2) < search_tolerance ||
						  std::abs(cell->vertex(vertex)[0] + notchWidth/2) < search_tolerance)
						  // found a cell which lays next to the crack plane or in the search area
					  {
						  cell->face(face)->set_boundary_id(parameters_internal.boundary_id_load_surface);
						  break;
					  }

					if ( std::abs(cell->face(face)->center()[0] - notchWidth/2) < search_tolerance ||
						 std::abs(cell->face(face)->center()[0] + notchWidth/2) < search_tolerance)
						// found a cell which lays right on the crack plane or in the search area
						cell->face(face)->set_boundary_id(parameters_internal.boundary_id_load_surface);
				}
				else if (std::abs(cell->face(face)->center()[2] - 0.0) < search_tolerance)
				{
					cell->face(face)->set_boundary_id(parameters_internal.boundary_id_minus_Z);
				}
				else
				{
					for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_face; ++vertex)
					 {
					 //Project the cell vertex to the XY plane and test the distance from the cylinder axis
						Point<dim> vertex_proj = cell->vertex(vertex);
						vertex_proj[2] = 0.0;
						Point<dim> centre_proj (0,0,0);
						if (std::abs(vertex_proj.distance(centre_proj) - notchWidth/2.0) < search_tolerance)
						{
							cell->face(face)->set_boundary_id(parameters_internal.boundary_id_hole);
							cell->face(face)->set_manifold_id(parameters_internal.manifold_id_hole);
							std::cout << "found" << std::endl;
							break;
						}
					}
					//AssertThrow(false, ExcMessage("3PointBeam - make_grid 3D: Found an unidentified face at the boundary. Maybe it slipt through the assignment or that face is simply not needed. So either check the implementation or comment this line in the code"));
				}
			  }
		}


		Tensor<1,dim> axis_dir;
		axis_dir[2]=1;  // axis point in the z-direction
		CylindricalManifold<dim> cylindrical_manifold_3d (axis_dir, centre);
		triangulation.set_manifold(parameters_internal.manifold_id_hole,cylindrical_manifold_3d);


		// Refine the cells around the center (x=0) of the beam
		if ( parameter.nbr_holeEdge_refinements>0 )
		{
			for (typename Triangulation<dim>::active_cell_iterator
			   cell = triangulation.begin_active();
			   cell != triangulation.end(); ++cell)
			{
			  for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
				  for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_face; ++vertex)
					if ( std::abs( cell->vertex(vertex)[0] ) <= notchWidth )
					 {
						cell->set_refine_flag();
						break;
					 }
			}

			triangulation.execute_coarsening_and_refinement();
		}

		triangulation.refine_global(parameter.nbr_global_refinements);	// ... Parameter.prm file

//		std::ofstream out ("grid-threepointbeam.eps");
//		GridOut grid_out;
//		GridOutFlags::Eps<2> eps_flags;
//		eps_flags.line_width = 0.1;
//		grid_out.set_flags (eps_flags);
//		grid_out.write_eps (triangulation, out);
//		std::cout << "Grid written to grid-threepointbeam.eps" << std::endl;
//		std::cout << "nElem: " << triangulation.n_active_cells() << std::endl;
//		AssertThrow(false,ExcMessage("ddd"));
	}


	// 2d grid:
//	template <int dim>
//	void make_grid( Triangulation<2> &triangulation, const Parameter::GeneralParameters &parameter, std::vector<unsigned int> Vec_boundary_id_collection )
//	{
//		parameterCollection parameters_internal ( Vec_boundary_id_collection );
//
//		const double search_tolerance = parameters_internal.search_tolerance;
//
//		// set the dimensions of the three point beam
//		const double length = parameter.width;
//		const double height = parameter.height;
//		const double notchWidth = parameter.notchWidth; // whereas the height is set to 1/4 of the height due to the mesh refinement and deletion of certain cells to create the notch
//
//		// Define the points spanning the three bodies
//		Point<dim> pL_left( -length/2., 0);
//		Point<dim> pR_right(  length/2., height);
//		Point<dim> pN_right( notchWidth, 0);	// the size of the notch body is twice a wide as the notch because only the innermost cells after the refinement will be deleted to create the notch
//		Point<dim> pN_left( -notchWidth, height);
//
//		Triangulation<dim> body_N, body_R, body_L;
//		// Notch body N
//		GridGenerator::hyper_rectangle( body_N, pN_right, pN_left );
//		// right hex
//		GridGenerator::hyper_rectangle( body_R, pN_right, pR_right );
//		// left hex
//		GridGenerator::hyper_rectangle( body_L, pN_left, pL_left );
//
//		// merge the three bodies
//		GridGenerator::merge_triangulations( body_N, body_R, triangulation );
//		GridGenerator::merge_triangulations( triangulation, body_L, triangulation );
//
//		// Requires some refinements:
//		triangulation.refine_global( 2 );	// to get a usable initial mesh with a decent notch
//
//		// remove the cell(s) where the notch shall be
//		std::set<typename Triangulation<dim>::active_cell_iterator > cells_to_remove;
//
//		for (typename Triangulation<dim>::active_cell_iterator
//			 cell = triangulation.begin_active();
//			 cell != triangulation.end(); ++cell)
//		{
//			for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
//			  if (cell->face(face)->at_boundary())
//			  {
//				// on the bottom face
//				if (std::abs(cell->face(face)->center()[1] - 0.0) < search_tolerance)
//				{
//					for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_face; ++vertex)
//					  if (std::abs(cell->vertex(vertex)[0] - 0.0) < search_tolerance) // found a cell which lays next to the crack plane
//					  {
//						  cells_to_remove.insert(cell);
//						  break;
//					  }
//
//					if ( std::abs(cell->face(face)->center()[0] - 0.0) < search_tolerance ) // found a cell which lays right on the crack plane
//						cells_to_remove.insert(cell);
//				}
//			  }
//		}
//
//		Assert(cells_to_remove.size() > 0, ExcInternalError());
//		Assert(cells_to_remove.size() != triangulation.n_active_cells(), ExcInternalError());
//		GridGenerator::create_triangulation_with_removed_cells(triangulation,cells_to_remove,triangulation);
//
//		//Clear boundary ID's
//		for (typename Triangulation<dim>::active_cell_iterator
//			 cell = triangulation.begin_active();
//			 cell != triangulation.end(); ++cell)
//		{
//			for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
//			  if (cell->face(face)->at_boundary())
//			  {
//				  cell->face(face)->set_all_boundary_ids(0);
//			  }
//		}
//
//		//Set boundary IDs and and manifolds
//		const Point<dim> centre (0,0);
//		for (typename Triangulation<dim>::active_cell_iterator
//			 cell = triangulation.begin_active();
//			 cell != triangulation.end(); ++cell)
//		{
//			for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
//			  if (cell->face(face)->at_boundary())
//			  {
//				//Set boundary IDs
//				if (std::abs(cell->face(face)->center()[1] - 0.0) < search_tolerance)
//				{
//					// -> on the bottom face
//					for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_face; ++vertex)
//					  if (std::abs(cell->vertex(vertex)[0] + length/2.) < search_tolerance) // found a cell where the support1 is
//					  {
//						  cell->face(face)->set_boundary_id(parameters_internal.boundary_id_support1);
//						  break;
//					  }
//					  else if (std::abs(cell->vertex(vertex)[0] - length/2.) < search_tolerance) // found a cell where the support2 is
//					  {
//						  cell->face(face)->set_boundary_id(parameters_internal.boundary_id_support2);
//						  break;
//					  }
//				}
//				else if (std::abs(cell->face(face)->center()[1] - height) < search_tolerance)
//				{
//					// -> on the top face
//					for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_face; ++vertex)
//					  if (std::abs(cell->vertex(vertex)[0] - notchWidth/2) < search_tolerance ||
//						  std::abs(cell->vertex(vertex)[0] + notchWidth/2) < search_tolerance)
//						  // found a cell which lays next to the crack plane or in the search area
//					  {
//						  cell->face(face)->set_boundary_id(parameters_internal.boundary_id_load_surface);
//						  break;
//					  }
//
//					if ( std::abs(cell->face(face)->center()[0] - notchWidth/2) < search_tolerance ||
//						 std::abs(cell->face(face)->center()[0] + notchWidth/2) < search_tolerance)
//						// found a cell which lays right on the crack plane or in the search area
//						cell->face(face)->set_boundary_id(parameters_internal.boundary_id_load_surface);
//				}
//				else
//				{
//					//AssertThrow(false, ExcMessage("3PointBeam - make_grid 3D: Found an unidentified face at the boundary. Maybe it slipt through the assignment or that face is simply not needed. So either check the implementation or comment this line in the code"));
//				}
//			  }
//		}
//
//		const Point<dim> centre_2d (0, height-0.5*notchWidth);
//		for (typename Triangulation<dim>::active_cell_iterator
//		   cell = triangulation.begin_active();
//		   cell != triangulation.end(); ++cell)
//		{
//		  for (unsigned int face=0; face<GeometryInfo<2>::faces_per_cell; ++face)
//			if (cell->face(face)->at_boundary())
//			  for (unsigned int vertex=0; vertex<GeometryInfo<2>::vertices_per_face; ++vertex)
//				if (  std::abs( (cell->vertex(vertex))[1] - height/4. ) < 1e-12 && std::abs( (cell->vertex(vertex))[0] - notchWidth/2. ) < 1e-12 )
////						std::a bs(cell->vertex(vertex).distance(centre_2d) - hole_diameter/2.0) < 1e-12)
//				 {
//					cell->face(face)->set_manifold_id(10);
//					break;
//				 }
//		}
//
//		SphericalManifold<2> spherical_manifold_2d (centre_2d);
//		triangulation.set_manifold(10,spherical_manifold_2d);
//
//		triangulation.refine_global(parameter.nbr_global_refinements);	// ... Parameter.prm file
//
//		std::ofstream out ("grid-threepointbeam.eps");
//		GridOut grid_out;
//		GridOutFlags::Eps<2> eps_flags;
//		eps_flags.line_width = 0.1;
//		grid_out.set_flags (eps_flags);
//		grid_out.write_eps (triangulation, out);
//		std::cout << "Grid written to grid-threepointbeam.eps" << std::endl;
//		std::cout << "nElem: " << triangulation.n_active_cells() << std::endl;
//		AssertThrow(false,ExcMessage("ddd"));
//	}
}

