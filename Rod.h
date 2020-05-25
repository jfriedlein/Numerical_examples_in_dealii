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

		const types::manifold_id manifold_id_surf = 10;

		const double search_tolerance = 1e-12;
	};


	/*
	 * Shift a layer of vertices of the triangulation at the coord position \a initial_pos to the position \a new_pos
	 * @param direction Gives the shift direction 0(x), 1(y), 2(z)
	 */
	template <int dim>
	void shift_vertex_layer( Triangulation<dim> &triangulation, double &initial_pos, double &new_pos, unsigned int direction )
	{
		bool shifted_node = false;
		for (typename Triangulation<dim>::active_cell_iterator
		   cell = triangulation.begin_active();
		   cell != triangulation.end(); ++cell)
		{
		  for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_cell; ++vertex )
		  {
			  if ( std::abs( cell->vertex(vertex)[direction] - initial_pos) < 1e-12/*search_tolerance*/ )
			  {
				  Point<dim> shift_vector;
				  shift_vector[direction] = (new_pos-initial_pos);
				  cell->vertex(vertex) += shift_vector;
				  shifted_node = true; // -> We have shifted at least a single node
			  }
		  }
		}
		// Ensure that we shifted at least a single node
		 AssertThrow( shifted_node==true, ExcMessage("shift_vertex_layer<< You haven't moved a single node. Please check the selection criterion initial_pos."));
	}


	double get_current_notch_radius( double &y_coord, const double &half_notch_length, const double &radius, const double &notch_radius, const double &R  )
	{
		// Here you can choose between a radial notch (smooth dent) and a sharp triangular notch (viewed in the cross section)
		const bool radial_notch = true;

		if ( radial_notch )
			return R + notch_radius - std::sqrt( R*R - y_coord*y_coord );
		else  /*linear notch*/
			return y_coord/half_notch_length * (radius-notch_radius) + notch_radius;
	}


	/*
	 * @param triangulation
	 * @param half_notch_length Half the length of the notch in y-direction. We only model 1/8 of the entire bar, hence only 1/2 of the notch length
	 * @param notch_radius The radius of the rod that is left at y=0 in the notch
	 * @param R The radius of the notch corresponds to the tool radius that could be used on a lathe to create the notch.
	 */
	template <int dim>
	void notch_body( Triangulation<dim> &triangulation, const double &half_notch_length, const double &radius, const double &notch_radius, const double &R  )
	{
		enum enum_coord_directions
		{
			x = 0, y = 1, z = 2
		};
		const double search_tolerance = 1e-12;

		if ( /*Deep notch: also adapt inner nodes in the notched area*/true )
		{
			// ToDo-optimize: Do we need the \a index_list anymore? (compare shallow notch below)
			std::vector<unsigned int> index_list;

			// Generate the notch
			 for (typename Triangulation<dim>::active_cell_iterator
			   cell = triangulation.begin_active();
			   cell != triangulation.end(); ++cell)
			 {
				  for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_cell; ++vertex)
				  {
					  // We look for all the points that lie in the notched area (y-coord in half length of notch)
					   double y_coord = cell->vertex(vertex)[y];
					   unsigned int index_vertex = cell->vertex_index(vertex);
					   if ( y_coord < half_notch_length && (std::find(index_list.begin(), index_list.end(), index_vertex) == index_list.end()) )
					   {
						  double x_coord = cell->vertex(vertex)[x];
						  double vertex_radius, z_coord;
						  if ( dim==2 )
							  vertex_radius=x_coord;
						  else if ( dim==3 )
						  {
							  z_coord = cell->vertex(vertex)[z];
							  vertex_radius = std::sqrt(x_coord*x_coord + z_coord*z_coord);
						  }
						  // The radius of the leftover notched material describes an arc along the y-coord.
						  // Hence, the radius of the notch changes with the y-coord
						   double current_notch_radius = get_current_notch_radius(y_coord, half_notch_length, radius, notch_radius, R );
						  // Set the shift vector that moves the vertex inwards (along its radius)
						   Point<dim> shift_vector;
						   shift_vector[x] = (current_notch_radius - radius) * std::sqrt(vertex_radius/radius) * x_coord/radius;
						   if ( dim==3 )
							   shift_vector[z] = (current_notch_radius - radius) * std::sqrt(vertex_radius/radius) * z_coord/radius;
						  // Apply the shift vector to the vertex
						   cell->vertex(vertex) += shift_vector;
						  // Add the current vertex to the list of already shifted vertices
						   index_list.push_back(index_vertex);
					   }
				  }
			 }
		}
		else /*only move the outer nodes of the notch inwards, this is limited to shallow notches and does not distort the inner cells*/
		{
			Assert( (radius - notch_radius) < radius/6.,
					ExcMessage("Rod<< You choose the shallow notching, but your notch seems to be very deep. Consider using the deep notch option above."));
			 for (typename Triangulation<dim>::active_cell_iterator
			   cell = triangulation.begin_active();
			   cell != triangulation.end(); ++cell)
			 {
			  for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
				if (cell->face(face)->at_boundary())
				  for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_face; ++vertex)
				  {
					  // We look for all the points that lie in the notched area (y-coord in half length of notch)
					   double y_coord = cell->face(face)->vertex(vertex)[y];
					   if ( y_coord < half_notch_length )
					   {
						  double x_coord = cell->face(face)->vertex(vertex)[x];
						  double z_coord = cell->face(face)->vertex(vertex)[z];
						  // Look for point that lies on the outer surface
						   if ( std::abs( std::sqrt(x_coord*x_coord + z_coord*z_coord) - radius ) < search_tolerance )
						   {
							  // The radius of the leftover notched material describes an arc along the y-coord.
							  // Hence, the radius of the notch changes with the y-coord
							   double current_notch_radius = R + notch_radius - std::sqrt( R*R - y_coord*y_coord );
							  // Set the shift vector that moves the vertex inwards
							   Point<dim> shift_vector;
							   shift_vector[0] = (current_notch_radius - radius) / radius * x_coord ;
							   shift_vector[2] = (current_notch_radius - radius) / radius * z_coord;
							  // Apply the shift vector to the vertex
							   cell->face(face)->vertex(vertex) += shift_vector;
						   }
					   }
				  }
			 }
		}
	}


// 3d grid
	/*
	 * @param triangulation
	 * @param length_of_the_entireRod The length of the entire rod
	 * @param radius_of_the_entireRod The outer radius of the rod
	 * @param length_of_the_entireNotchedArea The length of the notch in y-direction. We only model 1/8 of the entire bar, hence only 1/2 of the notch length
	 *
	 * @todo Add the remaining parameters with description
	 */
	template <int dim>
	void make_grid (
						Triangulation<3> &triangulation, std::vector<unsigned int> Vec_boundary_id_collection,
						const double &length_of_the_entireRod,
						const double &radius_of_the_entireRod,
						const double &length_of_the_entireNotchedArea,
						const double &radius_reductionFactor_in_notchedArea,
						const unsigned int n_additonal_refinements_in_y=1,
						const unsigned int n_global_refinements=0,
						const int n_max_of_elements_in_the_coarse_area = 6
				   )

	{
		// parameterCollection that contains the boundary ids
		 parameterCollection parameters_internal ( Vec_boundary_id_collection );

		const double search_tolerance = parameters_internal.search_tolerance;

		const double half_length = length_of_the_entireRod/2.;
		const double radius = radius_of_the_entireRod;
		const double half_notch_length = length_of_the_entireNotchedArea/2.;
		const double notch_radius = radius_reductionFactor_in_notchedArea * radius;
		const unsigned int n_additional_refinements = n_additonal_refinements_in_y;

		// The radius of the notch, e.g. the tool radius that was used to create the notch from the outside
		 const double R = ( half_notch_length*half_notch_length + (radius - notch_radius)*(radius - notch_radius) )
					      / ( 2.*(radius - notch_radius) );

		enum enum_coord_directions
		{
			x = 0, y = 1, z = 2
		};

		Assert(n_additional_refinements>0, ExcMessage("Rod<< Mesh not implemented for only 4 elements in total. Please increase the nbr_holeEdge_refinements to at least 1."));

		// Create in a first step the triangulation representing 1/8 of a cylinder
		 {
			// First we create a cylinder
			 Triangulation<dim> tria_full_cylinder;
			 GridGenerator::cylinder(tria_full_cylinder, radius, half_length);

			// Let's first refine the "cylinder" ones, because the initial mesh is a brick
			 tria_full_cylinder.refine_global( 1 );

			// We rotate the cylinder (oriented along x-axis by default) by 90° (=std::atan(1)*2 rad) around the z-axis
			 GridTools::rotate( std::atan(1)*2, z, tria_full_cylinder);

			// We only model 1/8 of the entire rod, hence we remove everything else
			 std::set<typename Triangulation<dim>::active_cell_iterator > cells_to_remove;
			 for (typename Triangulation<dim>::active_cell_iterator
				 cell = tria_full_cylinder.begin_active();
				 cell != tria_full_cylinder.end(); ++cell)
			 {
				// Remove all cells that are not in the first quadrant.
				// The 1/8 shall reside in the positive x,y,z quadrant
				if (cell->center()[x] < 0.0 || cell->center()[y] < 0.0 || cell->center()[z] < 0.0 )
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
				  cell->face(face)->set_all_boundary_ids(0);
		 }

		// Set boundary IDs and and manifolds
		 for (typename Triangulation<dim>::active_cell_iterator
			 cell = triangulation.begin_active();
			 cell != triangulation.end(); ++cell)
		 {
			for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
			 // Cells that describe the boundary can only describe the boundary when they possess a face that lies at the boundary:
			  if (cell->face(face)->at_boundary())
			  {
				// Cell at the x0-plane
				 if (std::abs(cell->face(face)->center()[x] - 0.0) < search_tolerance)
					cell->face(face)->set_boundary_id(parameters_internal.boundary_id_minus_X);
				// Cell at the y0-plane
				 else if (std::abs(cell->face(face)->center()[y] - 0.0) < search_tolerance)
					cell->face(face)->set_boundary_id(parameters_internal.boundary_id_minus_Y);
				// Cell at the z0-plane
				 else if (std::abs(cell->face(face)->center()[z] - 0.0) < search_tolerance)
					cell->face(face)->set_boundary_id(parameters_internal.boundary_id_minus_Z);
				// Cell at the other end of the rod
				 else if (std::abs(cell->face(face)->center()[y] - half_length) < search_tolerance)
					cell->face(face)->set_boundary_id(parameters_internal.boundary_id_plus_Y);
			  }
		 }

		// Attach a manifold to the curved boundary
		// @note We can only guarantee that the vertices sit on the curve, so we must test with their position instead of the cell centre.
		for (typename Triangulation<dim>::active_cell_iterator
		   cell = triangulation.begin_active();
		   cell != triangulation.end(); ++cell)
		{
		  for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
			if (cell->face(face)->at_boundary())
			  for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_face; ++vertex)
			  {
				 // Compute the projected radius in the xz-plane, so the distance between the vertex and the y-axis
				  double distance_2d_xz = std::sqrt( cell->vertex(vertex)[x]*cell->vertex(vertex)[x] + cell->vertex(vertex)[z]*cell->vertex(vertex)[z] );
				 if ( std::abs(distance_2d_xz - radius) < search_tolerance )
				 {
					// This vertex lies on the outer surface, hence the face belongs to the manifold
					// @note For some reason it is essential to use \a set_all_manifold_ids() instead of just set_manifold_id
					cell->face(face)->set_all_manifold_ids(parameters_internal.manifold_id_surf);
					break;
				 }
			  }
		}

		// Create a cylindrical manifold to be put on the outer cylindrical surface
		 CylindricalManifold<dim> cylindrical_manifold_3d (y); // y-axis
		 triangulation.set_manifold( parameters_internal.manifold_id_surf, cylindrical_manifold_3d );

		// Global refinement of the mesh to get a better approximation of the contour:\n
		// Previous: 2 elements for quarter arc; After global refinement: 4 elements
		 triangulation.refine_global( 1 );

		// Add some local refinements:
		// Cells are cut in y-direction, so we simple get some more cells that will be rearranged subsequently
		// @note For some reason I cannot cut_y two cells that lie next to each other. Hence, I only refine the cell at the y0-plane
		 for (unsigned int refine_counter=0; refine_counter < n_additional_refinements; refine_counter++)
		 {
			for (typename Triangulation<dim>::active_cell_iterator
			   cell = triangulation.begin_active();
			   cell != triangulation.end(); ++cell)
			{
				for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
					if (cell->face(face)->at_boundary())
						if ( std::abs(cell->face(face)->center()[y])<search_tolerance )
						{
							cell->set_refine_flag(RefinementCase<dim>::cut_y); // refine only in the y-direction
							break;
						}
			}
			triangulation.execute_coarsening_and_refinement();
		 }

	  // Shift the refinement layers in y-direction:
	  // This is a bit tricky and can best be comprehended on paper for specific example values.
		double initial_pos, new_pos;

		const unsigned int nbr_of_y_cells = 4 + n_additional_refinements;
		unsigned int nbr_of_coarse_y_cells = std::min(int(std::ceil(nbr_of_y_cells/2.)),n_max_of_elements_in_the_coarse_area);
		const unsigned int nbr_of_fine_y_cells = nbr_of_y_cells - nbr_of_coarse_y_cells;

		// Shift the coarsest cells such that the coarser outer area is uniformly discretised
		 for ( unsigned int i=1; i<=3; i++ )
		 {
			initial_pos = half_length * (4-i)/4.;
			new_pos = (nbr_of_coarse_y_cells - i)/double(nbr_of_coarse_y_cells) * (half_length - half_notch_length) + half_notch_length;
			shift_vertex_layer( triangulation, initial_pos, new_pos, y );
		 }

		// We have to grab a few more cells from the local refinements in case we want more than 9 cells in y-direction
		 if ( nbr_of_coarse_y_cells>4 )
			for ( unsigned int i=3; i<=(nbr_of_coarse_y_cells-2); i++ )
			{
				initial_pos = half_length * 1./(std::pow(2,i));
				new_pos = (nbr_of_coarse_y_cells-1 - i)/double(nbr_of_coarse_y_cells) * (half_length - half_notch_length) + half_notch_length;
				shift_vertex_layer( triangulation, initial_pos, new_pos, y );
			}

		// A small trick to get this general framework to operate even for the two lowest refinements 1 and 2
		 if ( n_additional_refinements <= 2 )
			 nbr_of_coarse_y_cells = 4;

		// Now we are down to the notch length
		 for ( unsigned int i=(nbr_of_coarse_y_cells-1); i<=(n_additional_refinements+2); i++ )
		 {
			initial_pos = half_length * 1./(std::pow(2,i));
			new_pos = (nbr_of_y_cells-1 - i)/double(nbr_of_fine_y_cells)  * half_notch_length;
			shift_vertex_layer( triangulation, initial_pos, new_pos, y );
		 }

	    // Generate the notch
		 notch_body( triangulation, half_notch_length, radius, notch_radius, R );

		// Possibly some additional global isotropic refinements
		 triangulation.refine_global(n_global_refinements);	// ... Parameter.prm file

//		// include the following two scopes to see directly how the variation of the input parameters changes the geometry of the grid
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
		void make_grid( Triangulation<3> &triangulation, const Parameter::GeneralParameters &parameter, std::vector<unsigned int> Vec_boundary_id_collection )

		{
			// parameterCollection that contains the boundary ids
			 parameterCollection parameters_internal ( Vec_boundary_id_collection );

			const double search_tolerance = parameters_internal.search_tolerance;

			const double half_length = parameter.width/2.;//53.34/2.;
			const double radius = parameter.holeRadius;//6.4135;
			const double half_notch_length = parameter.notchWidth/2.;//8.98/2.;
			const double notch_radius = parameter.ratio_x * radius; // 0.982

			const unsigned int n_additional_refinements = parameter.nbr_holeEdge_refinements;
			const unsigned int n_global_refinements = parameter.nbr_global_refinements;
			const int n_max_of_elements_in_the_coarse_area = 6;

			// The radius of the notch, e.g. the tool radius that was used to create the notch from the outside
			 const double R = ( half_notch_length*half_notch_length + (radius - notch_radius)*(radius - notch_radius) )
						      / ( 2.*(radius - notch_radius) );

			enum enum_coord_directions
			{
				x = 0, y = 1, z = 2
			};

			Assert(n_additional_refinements>0, ExcMessage("Rod<< Mesh not implemented for only 4 elements in total. Please increase the nbr_holeEdge_refinements to at least 1."));

			// Create in a first step the triangulation representing 1/8 of a cylinder
			 {
				// First we create a cylinder
				 Triangulation<dim> tria_full_cylinder;
				 GridGenerator::cylinder(tria_full_cylinder, radius, half_length);

				// Let's first refine the "cylinder" ones, because the initial mesh is a brick
				 tria_full_cylinder.refine_global( 1 );

				// We rotate the cylinder (oriented along x-axis by default) by 90° (=std::atan(1)*2 rad) around the z-axis
				 GridTools::rotate( std::atan(1)*2, z, tria_full_cylinder);

				// We only model 1/8 of the entire rod, hence we remove everything else
				std::set<typename Triangulation<dim>::active_cell_iterator > cells_to_remove;
				for (typename Triangulation<dim>::active_cell_iterator
					 cell = tria_full_cylinder.begin_active();
					 cell != tria_full_cylinder.end(); ++cell)
				{
					// Remove all cells that are not in the first quadrant.
					// The 1/8 shall reside in the positive x,y,z quadrant
					if (cell->center()[x] < 0.0 || cell->center()[y] < 0.0 || cell->center()[z] < 0.0 )
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
					  cell->face(face)->set_all_boundary_ids(0);
			}

			// Set boundary IDs and and manifolds
			for (typename Triangulation<dim>::active_cell_iterator
				 cell = triangulation.begin_active();
				 cell != triangulation.end(); ++cell)
			{
				for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
				 // Cells that describe the boundary can only describe the boundary when they possess a face that lies at the boundary:
				  if (cell->face(face)->at_boundary())
				  {
					// Cell at the x0-plane
					 if (std::abs(cell->face(face)->center()[x] - 0.0) < search_tolerance)
						cell->face(face)->set_boundary_id(parameters_internal.boundary_id_minus_X);
					// Cell at the y0-plane
					 else if (std::abs(cell->face(face)->center()[y] - 0.0) < search_tolerance)
						cell->face(face)->set_boundary_id(parameters_internal.boundary_id_minus_Y);
					// Cell at the z0-plane
					 else if (std::abs(cell->face(face)->center()[z] - 0.0) < search_tolerance)
						cell->face(face)->set_boundary_id(parameters_internal.boundary_id_minus_Z);
					// Cell at the other end of the rod
					 else if (std::abs(cell->face(face)->center()[y] - half_length) < search_tolerance)
						cell->face(face)->set_boundary_id(parameters_internal.boundary_id_plus_Y);
				  }
			}

			// Attach a manifold to the curved boundary
			// @note We can only guarantee that the vertices sit on the curve, so we must test with their position instead of the cell centre.
			for (typename Triangulation<dim>::active_cell_iterator
			   cell = triangulation.begin_active();
			   cell != triangulation.end(); ++cell)
			{
			  for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
				if (cell->face(face)->at_boundary())
				  for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_face; ++vertex)
				  {
					 // Compute the projected radius in the xz-plane, so the distance between the vertex and the y-axis
					  double distance_2d_xz = std::sqrt( cell->vertex(vertex)[x]*cell->vertex(vertex)[x] + cell->vertex(vertex)[z]*cell->vertex(vertex)[z] );
					 if ( std::abs(distance_2d_xz - radius) < search_tolerance )
					 {
						// This vertex lies on the outer surface, hence the face and the face belongs to the manifold
						// @note For some reason it is essential to use \æ set_all_manifold_ids() instead of just set_manifold_id
						cell->face(face)->set_all_manifold_ids(parameters_internal.manifold_id_surf);
						break;
					 }
				  }
			}

			// Create a cylindrical manifold to be put on the outer cylindrical surface
			 CylindricalManifold<dim> cylindrical_manifold_3d (y); // y-axis
			 triangulation.set_manifold( parameters_internal.manifold_id_surf, cylindrical_manifold_3d );

			// Global refinement of the mesh to get a better approximation of the contour:\n
			// Previous: 2 elements for quarter arc; After global refinement: 4 elements
			 triangulation.refine_global( 1 );

			// Add some local refinements:
			// Cells are cut in y-direction, so we simple get some more cells that will be rearranged subsequently
			// @note For some reason I cannot cut_y two cells that lie next to each other. Hence, I only refine the cell at the y0-plane
			 for (unsigned int refine_counter=0; refine_counter < n_additional_refinements; refine_counter++)
			 {
				for (typename Triangulation<dim>::active_cell_iterator
				   cell = triangulation.begin_active();
				   cell != triangulation.end(); ++cell)
				{
					for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
						if (cell->face(face)->at_boundary())
							if ( std::abs(cell->face(face)->center()[y])<search_tolerance )
							{
								cell->set_refine_flag(RefinementCase<dim>::cut_y); // refine only in the y-direction
								break;
							}
				}
				triangulation.execute_coarsening_and_refinement();
			 }

		  // Shift the refinement layers in y-direction:
		  // This is a bit tricky and can best be comprehended on paper for specific example values.
			double initial_pos, new_pos;

			const unsigned int nbr_of_y_cells = 4 + n_additional_refinements;
			unsigned int nbr_of_coarse_y_cells = std::min(int(std::ceil(nbr_of_y_cells/2.)),n_max_of_elements_in_the_coarse_area);
			const unsigned int nbr_of_fine_y_cells = nbr_of_y_cells - nbr_of_coarse_y_cells;

			// Shift the coarsest cells such that the coarser outer area is uniformly discretised
			 for ( unsigned int i=1; i<=3; i++ )
			 {
				initial_pos = half_length * (4-i)/4.;
				new_pos = (nbr_of_coarse_y_cells - i)/double(nbr_of_coarse_y_cells) * (half_length - half_notch_length) + half_notch_length;
				shift_vertex_layer( triangulation, initial_pos, new_pos, y );
			 }

			// We have to grab a few more cells from the local refinements in case we want more than 9 cells in y-direction
			 if ( nbr_of_coarse_y_cells>4 )
				for ( unsigned int i=3; i<=(nbr_of_coarse_y_cells-2); i++ )
				{
					initial_pos = half_length * 1./(std::pow(2,i));
					new_pos = (nbr_of_coarse_y_cells-1 - i)/double(nbr_of_coarse_y_cells) * (half_length - half_notch_length) + half_notch_length;
					shift_vertex_layer( triangulation, initial_pos, new_pos, y );
				}

			// A small trick to get this general framework to operate even for the two lowest refinements 1 and 2
			 if ( n_additional_refinements <= 2 )
				 nbr_of_coarse_y_cells = 4;

			// Now we are down to the notch length
			 for ( unsigned int i=(nbr_of_coarse_y_cells-1); i<=(n_additional_refinements+2); i++ )
			 {
				initial_pos = half_length * 1./(std::pow(2,i));
				new_pos = (nbr_of_y_cells-1 - i)/double(nbr_of_fine_y_cells)  * half_notch_length;
				shift_vertex_layer( triangulation, initial_pos, new_pos, y );
			 }

		    // Generate the notch
			 notch_body( triangulation, half_notch_length, radius, notch_radius, R );

			// Possibly some additional global isotropic refinements
			 triangulation.refine_global(n_global_refinements);	// ... Parameter.prm file

	//		// include the following two scopes to see directly how the variation of the input parameters changes the geometry of the grid
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
	
	
	// 2d grid
	template <int dim>
	void make_grid( Triangulation<2> &triangulation, const Parameter::GeneralParameters &parameter, std::vector<unsigned int> Vec_boundary_id_collection )
	{
		/*
		 * Input arguments:
		 * * boundary ids and manifold id
		 * * search tolerance
		 * * parameter.width,holeRadius,notchWidth,ratio_x,nbr_holeEdge_refinements,nbr_global_refinements
		 */

		parameterCollection parameters_internal ( Vec_boundary_id_collection );

		const double search_tolerance = parameters_internal.search_tolerance;

		const double half_length = parameter.width/2.;//53.34/2.;
		const double radius = parameter.holeRadius;//6.4135;
		const double half_notch_length = parameter.notchWidth/2.;//8.98/2.;
		const double notch_radius = parameter.ratio_x * radius; // 0.982

		const int n_max_of_elements_in_the_coarse_area = 6;

		// The radius of the notch, e.g. the tool radius that was used to create the notch from the outside
		 const double R = ( half_notch_length*half_notch_length + (radius - notch_radius)*(radius - notch_radius) )
						  / ( 2.*(radius - notch_radius) );

		enum enum_coord_directions
		{
			x = 0, y = 1
		};

		Assert(parameter.nbr_holeEdge_refinements>0, ExcMessage("Rod<< Mesh not implemented for only 4 elements in total. Please increase the nbr_holeEdge_refinements to at least 1."));

		// Create in a first step the triangulation representing 1/8 of a cylinder
		 {
			// First we create a cylinder
			 Point<dim> p1 (0,0);
			 Point<dim> p2 (radius,half_length);
			 GridGenerator::hyper_rectangle(triangulation, p1, p2);

			// Let's first refine the "cylinder" ones, because the initial mesh is a brick
			 triangulation.refine_global( 2 );
		 }

		// Clear boundary ID's
		for (typename Triangulation<dim>::active_cell_iterator
			 cell = triangulation.begin_active();
			 cell != triangulation.end(); ++cell)
		{
			for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
			  if (cell->face(face)->at_boundary())
				  cell->face(face)->set_all_boundary_ids(0);
		}

		// Set boundary IDs and and manifolds
		for (typename Triangulation<dim>::active_cell_iterator
			 cell = triangulation.begin_active();
			 cell != triangulation.end(); ++cell)
		{
			for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
			 // Cells that describe the boundary can only describe the boundary when they possess a face that lies at the boundary:
			  if (cell->face(face)->at_boundary())
			  {
				// Cell at the x0-plane
				 if (std::abs(cell->face(face)->center()[x] - 0.0) < search_tolerance)
					cell->face(face)->set_boundary_id(parameters_internal.boundary_id_minus_X);
				// Cell at the y0-plane
				 else if (std::abs(cell->face(face)->center()[y] - 0.0) < search_tolerance)
					cell->face(face)->set_boundary_id(parameters_internal.boundary_id_minus_Y);
				// Cell at the other end of the rod
				 else if (std::abs(cell->face(face)->center()[y] - half_length) < search_tolerance)
					cell->face(face)->set_boundary_id(parameters_internal.boundary_id_plus_Y);
			  }
		}

		// Add some local refinements:
		// Cells are cut in y-direction, so we simple get some more cells that will be rearranged subsequently
		// @note For some reason I cannot cut_y two cells that lie next to each other. Hence, I only refine the cell at the y0-plane
		 for (unsigned int refine_counter=0; refine_counter<parameter.nbr_holeEdge_refinements; refine_counter++)
		 {
			for (typename Triangulation<dim>::active_cell_iterator
			   cell = triangulation.begin_active();
			   cell != triangulation.end(); ++cell)
			{
				for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
					if (cell->face(face)->at_boundary())
						if ( std::abs(cell->face(face)->center()[y])<search_tolerance )
						{
							cell->set_refine_flag(RefinementCase<dim>::cut_y); // refine only in the y-direction
							break;
						}
			}
			triangulation.execute_coarsening_and_refinement();
		 }

	  // ToDo-optimize: Isn't this very similar to the 3D case? Maybe merge 2D and 3D, also the surrounding code seems familiar
	  // Shift the refinement layers in y-direction:
	  // This is a bit tricky and can best be comprehended on paper for specific example values.
		double initial_pos, new_pos;

		const unsigned int nbr_of_y_cells = 4 + parameter.nbr_holeEdge_refinements;
		unsigned int nbr_of_coarse_y_cells = std::min(int(std::ceil(nbr_of_y_cells/2.)),n_max_of_elements_in_the_coarse_area);
		const unsigned int nbr_of_fine_y_cells = nbr_of_y_cells - nbr_of_coarse_y_cells;

		// Shift the coarsest cells such that the coarser outer area is uniformly discretised
		 for ( unsigned int i=1; i<=3; i++ )
		 {
			initial_pos = half_length * (4-i)/4.;
			new_pos = (nbr_of_coarse_y_cells - i)/double(nbr_of_coarse_y_cells) * (half_length - half_notch_length) + half_notch_length;
			shift_vertex_layer( triangulation, initial_pos, new_pos, y );
		 }

		// We have to grab a few more cells from the local refinements in case we want more than 9 cells in y-direction
		 if ( nbr_of_coarse_y_cells>4 )
			for ( unsigned int i=3; i<=(nbr_of_coarse_y_cells-2); i++ )
			{
				initial_pos = half_length * 1./(std::pow(2,i));
				new_pos = (nbr_of_coarse_y_cells-1 - i)/double(nbr_of_coarse_y_cells) * (half_length - half_notch_length) + half_notch_length;
				shift_vertex_layer( triangulation, initial_pos, new_pos, y );
			}

		// A small trick to get this general framework to operate even for the two lowest refinements 1 and 2
		 if ( parameter.nbr_holeEdge_refinements<=2 )
			 nbr_of_coarse_y_cells = 4;

		// Now we are down to the notch length
		 for ( unsigned int i=(nbr_of_coarse_y_cells-1); i<=(parameter.nbr_holeEdge_refinements+2); i++ )
		 {
			initial_pos = half_length * 1./(std::pow(2,i));
			new_pos = (nbr_of_y_cells-1 - i)/double(nbr_of_fine_y_cells)  * half_notch_length;
			shift_vertex_layer( triangulation, initial_pos, new_pos, y );
		 }

		// Generate the notch
		 notch_body( triangulation, half_notch_length, radius, notch_radius, R );

		// Possibly some additional global isotropic refinements
		 triangulation.refine_global(parameter.nbr_global_refinements);	// ... Parameter.prm file

		// include the following two scopes to see directly how the variation of the input parameters changes the geometry of the grid
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

			if ( parameter.driver == enums::Dirichlet )
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
}
