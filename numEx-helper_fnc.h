#ifndef NUMEX_HELPERFNC
#define NUMEX_HELPERFNC

#include <deal.II/grid/grid_out.h>


#include <iostream>


namespace enums
{
	/**
	 * For the driver: prescribed displacement or force (for AL and normal disp-, force-control)
	 */
	 enum enum_driver
	 {
		Neumann =   1,
		Dirichlet = 2,
		Contact = 3
	 };

    enum enum_boundary_ids
	{
    	// @note Don't use id 0 for anything (it's the standard value or the reset value)
    	id_boundary_none = 999,
		id_boundary_fix = 99,

		id_boundary_xMinus = 2,
		id_boundary_xPlus =  3,

		id_boundary_yMinus = 4,
		id_boundary_yPlus =  5,

		id_boundary_zMinus = 6,
		id_boundary_zPlus =  7,

		id_boundary_load_default = 8,
		id_boundary_secondaryLoad_default = 9,

		id_body_dummy = 10,

		id_boundary_xMinus1 = 21,
		id_boundary_xMinus2 = 22,
		id_boundary_xPlus1 = 31,
		id_boundary_xPlus2 = 32,
		id_boundary_yPlus2 = 52
	};

    enum
	{
		id_primary_load = 0,
		id_secondary_load = 1
	};
    
    enum enum_loading_type
	{
    	standard = 0,
    	tension = 1,
		compression = 2,
		Brick_Seupel_etal_a = 3
	};

    enum enum_refine_special
	{
    	Mesh_refine_special_standard = 0,
		Rod_refine_special_uniform = 1,
		Mesh_refine_special_innermost = 2,
		Mesh_refine_special_Simo = 3,
		Mesh_refine_uniform = 4
	};
    
    enum enum_coord
	 {
		x = 0,
		y = 1,
		z = 2
		//r = 0,
		//theta = 2,
		//u = 0,
		//w = 1
	 };

   enum enum_special_QP
	{
		 tracked_QP = 4
	};
   
   enum enum_output_type
   {
	   output_eps = 0,
	   output_inp = 1
   };
   
   enum enum_notch_type
   {
	   notch_linear = 0,
	   notch_round = 1
   };
   
   /**
    * @todo Don't use x0 y0 etc.
    */
   enum enum_BC
   {
	   BC_none = 0,//!< BC_none
	   BC_sym = 1, //!< BC_sym
	   BC_fix = 2, //!< BC_fix
	   BC_x0 = 3,  //!< BC_x0
	   BC_x0_z0 = 4   //!< BC_x0_z0
   };
}

namespace numEx
{
	/**
	 * Apply the boundary condition on the given boundary id \a boundary_id for the component \a component.
	 * For symmetry constraints (zero displacement) we get a \load_increment of 0, so we use ZeroFunction, else we apply a ConstantFunciton.
	 */
	template<int dim>
	void BC_apply ( const enums::enum_boundary_ids boundary_id, const unsigned int component, const double load_increment, const bool &apply_dirichlet_bc,
					const DoFHandler<dim> &dof_handler,const FESystem<dim> &fe, AffineConstraints<double> &constraints )
	{	
		// @todo The component masks and displacement mask are HARDCODED and should depend on the order in the FESystem
		const FEValuesExtractors::Scalar displacement_component(component);
		const unsigned int n_components = fe.n_components();
		
		if (apply_dirichlet_bc == true )
		{
			// Apply the given load
			if ( load_increment!=0 )
			{
				VectorTools::interpolate_boundary_values(
															dof_handler,
															boundary_id,
															ConstantFunction<dim> (load_increment/*add only the increment*/, n_components),
															constraints,
															fe.component_mask(displacement_component)
														);
			}
			// Apply zero displacement BC
			else
			{
				VectorTools::interpolate_boundary_values(
															dof_handler,
															boundary_id,
															ZeroFunction<dim> ( n_components ),
															constraints,
															fe.component_mask(displacement_component)
														);
			}
		}
		else
		{
			VectorTools::interpolate_boundary_values(
														dof_handler,
														boundary_id,
														ZeroFunction<dim> ( n_components ),
														constraints,
														fe.component_mask(displacement_component)
													);
		}
	}
	
	template<int dim>
	void BC_apply_fix ( const enums::enum_boundary_ids boundary_id, const DoFHandler<dim> &dof_handler, 
						const FESystem<dim> &fe, AffineConstraints<double> &constraints )
	{
		FEValuesExtractors::Vector displacements(0);
		ComponentMask disp_mask = fe.component_mask (displacements);
		
		VectorTools::interpolate_boundary_values(
													dof_handler,
													boundary_id,
													ZeroFunction<dim> (fe.n_components()),
													constraints,
													disp_mask // all disp components
												);
	}
	

	/**
	 * Shift the given points \a points of the triangulation \a tria by the vectors in \a shift
	 */
	template<int dim>
	void shift_vertex_by_vector ( Triangulation<dim> &tria, const std::vector< Point<dim> > &points, const std::vector< Point<dim> > &shift, const std::string &numEx_name="" )
	{
		unsigned int shifted_node = 0;
		const unsigned int n_points = points.size();
		for (typename Triangulation<dim>::active_cell_iterator
		   cell = tria.begin_active();
		   cell != tria.end(); ++cell)
		{
		  for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_cell; ++vertex )
			  for ( unsigned int i=0; i < n_points; i++)
				  if ( cell->vertex(vertex).distance(points[i]) < 1e-12/*search_tolerance*/ )
				  {
					  cell->vertex(vertex) += shift[i];
					  shifted_node += 1; // -> We have shifted at least a single node
				  }
		}
		AssertThrow( shifted_node == n_points, ExcMessage(numEx_name+"<< Distortion, we only shifted "+std::to_string(shifted_node)+
														  " instead of "+std::to_string(n_points)+" vertices."));
	}

	
	template<int dim>
	Point<dim> extract_dim ( const Point<3> &point_3D )
	{
		Point<dim> point_dim;
		for ( unsigned int i=0; i<dim; i++)
			point_dim[i] = point_3D[i];
		
		return point_dim;
	}
	
	/**
	 * Clear boundary ID's
	 */
	template<int dim>
	void clear_boundary_IDs ( Triangulation<dim> &triangulation )
	{
		// Iterate over each cell of the triangulation
		for ( typename Triangulation<dim>::active_cell_iterator
			  cell = triangulation.begin_active();
			  cell != triangulation.end(); ++cell )
		{
			// Iterate over each face of this cell
			 for (unsigned int face=0; face < GeometryInfo<dim>::faces_per_cell; ++face)
				// Iterate only over faces that lay at the boundary
				 if ( cell->face(face)->at_boundary() )
					 // Set the the boundary ids to zero (standard|unassigned value)
					  cell->face(face)->set_all_boundary_ids(0);
		}
	}
	
	template<int dim>
	void output_triangulation ( const Triangulation<dim> &triangulation, const unsigned int output_type=enums::output_eps, const std::string numEx_name="numEx" )
	{
		std::cout << "numEx<< Writting triangulation to output ..." << std::endl;
		std::ostringstream filename;
		switch ( output_type )
		{
			case enums::output_eps:
			{
				filename << "grid-" << numEx_name << ".eps";
				std::ofstream out (filename.str().c_str()); // @todo-optimize That should be simpler?!
				GridOut grid_out;
				grid_out.write_eps (triangulation, out);
				break;
			}
			case enums::output_inp:
			{
				filename << "grid-" << numEx_name << ".inp";
				std::ofstream out_ucd(filename.str().c_str());
				GridOut grid_out;
				GridOutFlags::Ucd ucd_flags(true,true,true);
				grid_out.set_flags(ucd_flags);
				grid_out.write_ucd(triangulation, out_ucd);
				break;
			}
			default:
				AssertThrow(false, ExcMessage(numEx_name+" - output_triangulation<< You choose a not implemented output type, try eps or inp instead."));
		}
		std::cout << "numEx<< ... grid written to " << filename.str().c_str() << std::endl;
	}

	
	/**
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
		 AssertThrow( shifted_node==true, ExcMessage("shift_vertex_layer<< You haven't moved a single node. Please check the selection criterion initial_pos "
				 	 	 	 	 	 	 	 	 	 +std::to_string(initial_pos)+" vs your new_pos "+std::to_string(new_pos)+"."));
	}
	
	template <int dim>
	class BeamEnd : public Function<dim>
	{
	public:
		BeamEnd ( const double &length0, const double &height0, const double &lambda_n, const double &current_load_increment, const unsigned int &n_components )
		:
			l_0(length0),
			h_0(height0),
			theta_n(lambda_n),
			incr(current_load_increment),
			n_comp(n_components)
		{
		}
		virtual ~BeamEnd() {};

		double l_0, h_0;
		double theta_n;
		double incr;
		unsigned int n_comp;


		// return all components at one point
		virtual void vector_value(const Point<dim> &p, Vector<double>   &value) const
		{
			Tensor<1,dim> u_m_n;
			u_m_n[enums::x] = l_0 * ( std::sin(theta_n+1e-10)/(theta_n+1e-10) - 1.);
			u_m_n[enums::y] = l_0 * ( std::cos(theta_n+1e-10)/(theta_n+1e-10) - 1./(theta_n+1e-10));

			const double theta_n1 = theta_n + incr;
			Tensor<1,dim> u_m_n1;
			u_m_n1[enums::x] = l_0 * ( std::sin(theta_n1+1e-10)/(theta_n1+1e-10) - 1.);
			u_m_n1[enums::y] = l_0 * ( std::cos(theta_n1+1e-10)/(theta_n1+1e-10) - 1./(theta_n1+1e-10));

			Tensor<1,dim> u_n;
			u_n[enums::x] = 0 						+ u_m_n[enums::x] + (p[enums::y]-h_0/2.) * std::sin(theta_n);
			u_n[enums::y] = -(p[enums::y]-h_0/2.) 	+ u_m_n[enums::y] + (p[enums::y]-h_0/2.) * std::cos(theta_n);

			Tensor<1,dim> u_n1;
			u_n1[enums::x] = 0 						+ u_m_n1[enums::x] + (p[enums::y]-h_0/2.) * std::sin(theta_n1);
			u_n1[enums::y] = -(p[enums::y]-h_0/2.)	+ u_m_n1[enums::y] + (p[enums::y]-h_0/2.) * std::cos(theta_n1);

			for ( unsigned int i=0; i<dim; i++)
				value[i] = u_n1[i] - u_n[i];
		}

//		virtual double value (const Point<dim>   &p,
//							const unsigned int   = 0) const override final
//		{
//			double val = 1.;
//			SymmetricTensor<2,2> stress_tensor_dummy;
//			Point<dim> point_xy = p;
//			get_KIRSCH_stresses( point_xy, current_load, nu, stress_tensor_dummy, val);// 0929 stress_vonMises );
//			return val;
//		}
	};


	/**
	 * @todo-optimize Use deal.II object information to get normal and tangent for the face
	 */
	template <int dim>
	class NotchClass
	{
	public:
		// Linear notch
		NotchClass ( const unsigned int notch_type, const double &notch_length, const double &notch_depth, const Point<3> &notch_reference_point, 
					 const types::boundary_id &notch_face_BID, const Point<3> &face_normal_vector, const unsigned int notch_tangent_dir )
		:
		type(notch_type),
		length(notch_length),
		depth(notch_depth),
		ref_pos(notch_reference_point),
		face_BID(notch_face_BID),
		normal_vector(face_normal_vector),
		tangent_dir(notch_tangent_dir)
		{
		}
		
		// Round notch
		NotchClass ( const unsigned int notch_type, const double &notch_length, const double &notch_depth, const Point<3> &notch_reference_point,
					 const types::boundary_id &notch_face_BID, const Point<3> &face_normal_vector, const unsigned int notch_tangent_dir,
					 const types::manifold_id &notch_manifold_id )
		:
		type(notch_type),
		length(notch_length),
		depth(notch_depth),
		ref_pos(notch_reference_point),
		face_BID(notch_face_BID),
		normal_vector(face_normal_vector),
		tangent_dir(notch_tangent_dir),
		manifold_id(notch_manifold_id)
		{
			 radius = ( std::pow(notch_length/2.,2) + notch_depth*notch_depth ) / ( 2. * notch_depth );
			 cyl_center = notch_reference_point;
			 cyl_center += face_normal_vector * (radius - notch_depth);
		}
		
		// @todo Maybe also some constructors and functions in case we get the center point, etc.
		
		unsigned int type;
		double length;
		double depth;
		Point<3> ref_pos;
		types::boundary_id face_BID;
		Point<3> normal_vector;
		unsigned int tangent_dir;
		double radius = 0.;
		types::manifold_id manifold_id;
		Point<3> cyl_center;
	};
	
	double get_current_notch_radius( const unsigned int notch_type, double &y_coord, const double &half_notch_length, const double &radius, const double &notch_radius, const double &R  )
	{
		switch ( notch_type )
		{
			case enums::notch_round:
				return R + notch_radius - std::sqrt( R*R - y_coord*y_coord );
			case enums::notch_linear:
				return y_coord/half_notch_length * (radius-notch_radius) + notch_radius;
			default:
				AssertThrow(false, ExcMessage("numEx - << notch type not implemented"));
				return 0;
		}
	}
	
	template<int dim>
	double get_notching ( const NotchClass<dim> &notch, const double &delta_y )
	{
		switch ( notch.type )
		{
			case enums::notch_round:
				return std::sqrt( notch.radius*notch.radius - delta_y*delta_y ) - ( notch.radius - notch.depth );
			case enums::notch_linear:
				return notch.depth * ( 1. - delta_y / (notch.length/2.) );
			default:
				AssertThrow(false, ExcMessage("numEx - get_notching<< notch type not implemented"));
				return 0;
		}
	}

	
	template<int dim>
	void notch_body ( Triangulation<dim> &triangulation, const NotchClass<dim> &notch )
	{
		if ( /*soft notch*/true )
		{
			std::vector<unsigned int> index_list;

			// Generate the notch:
			// @note The following is far from trivial, because
			 for (typename Triangulation<dim>::active_cell_iterator
			   cell = triangulation.begin_active();
			   cell != triangulation.end(); ++cell)
			 {
				  for (unsigned int face=0; face < GeometryInfo<dim>::faces_per_cell; ++face)
					  if (cell->face(face)->at_boundary() && cell->face(face)->boundary_id()==notch.face_BID )
						  for (unsigned int vertex=0; vertex < GeometryInfo<dim>::vertices_per_face; ++vertex)
						  {
							  // Distance from the given vertex to the reference point \a POS of the notch. We can directly
							  // use the absolute value because the notch is symmetric
							   double distance_vertex2POS = std::abs( cell->face(face)->vertex(vertex)[notch.tangent_dir] - notch.ref_pos[notch.tangent_dir] );
							  // We can use "less than" without equal, because we woud do nothing for the corner points anyway.
							  // This might also avoid spreading of the manifold to an outer cell's face.
							   if ( distance_vertex2POS < notch.length/2. )
							   {
								  unsigned int index_vertex = cell->face(face)->vertex_index(vertex);
								  if ( std::find(index_list.begin(), index_list.end(), index_vertex) == index_list.end())
								  {
									  // Compute the absolute value we need to move the vertex inside
									   double notching = get_notching( notch, distance_vertex2POS );

									  // Shift the vertex inwards via the negative normal vector of the face
									   cell->face(face)->vertex(vertex) -= notching * extract_dim<dim>( notch.normal_vector );
									   
									  // Add the current vertex to the list of already shifted vertices
									   index_list.push_back(index_vertex);
								  } // end if(vertex not already shifted)
								  
								  // Assign cylindrical manifold for round notches. I guess we also have to do this
								  // even if we don't shift the vertex at hand, because it is possible that we already
								  // shifted the same vertex but from the neighbouring cell, but still need this cell's face to have
								  // the correct manifold id.
								   if ( notch.type == enums::notch_round )
										cell->face(face)->set_all_manifold_ids(notch.manifold_id);
							   } // end if(to be shifted)
						  } // end for(vertex)
			 } // end for(cell)
		}
	}

	
	/**
	 * @todo Think about using this fnc only for 2D and then using extrude_triangulation
	 */
	template<int dim>
	void prepare_tria_for_notching ( Triangulation<dim> &triangulation, const NotchClass<dim> &notch )
	{
		// First, we compute the corner points
		 Point<dim> point_corner_positive = extract_dim<dim>( notch.ref_pos ); 
		 point_corner_positive[notch.tangent_dir] += notch.length / 2.;
		
		 Point<dim> point_corner_negative = extract_dim<dim>( notch.ref_pos ); 
		 point_corner_negative[notch.tangent_dir] -= notch.length / 2.;
		 
		 // The point at the reference point shall be shifted into the body, so we look for
		 // vertex at this point.
		 Point<dim> point_corner_depth = extract_dim<dim>( notch.ref_pos );

		// Then, we find the nearest vertex and adapt them to the corner points either locally or globally
		// We limit the list of vertices only to the vertices at the boundary to reduce the number of vertices
		// to be tested and to avoid finding vertices inside the body.
		 const unsigned int vertexID_closest_positive = GridTools::find_closest_vertex( GridTools::get_all_vertices_at_boundary(triangulation), point_corner_positive );
		 const unsigned int vertexID_closest_negative = GridTools::find_closest_vertex( GridTools::get_all_vertices_at_boundary(triangulation), point_corner_negative );
		 const unsigned int vertexID_closest_depth = 	GridTools::find_closest_vertex( GridTools::get_all_vertices_at_boundary(triangulation), point_corner_depth );

		switch ( notch.type )
		{
			case enums::notch_round:
				if ( vertexID_closest_positive == vertexID_closest_negative )
					AssertThrow( false, ExcMessage( "prepare_tria_for_notching<< The two found points closest to the corner points of the notch are identical, "
													"but need to be different. Use a finer mesh.") );
				break;
			case enums::notch_linear:
				if ( vertexID_closest_positive == vertexID_closest_negative 
					||  vertexID_closest_depth == vertexID_closest_negative
					||  vertexID_closest_positive == vertexID_closest_depth )
					AssertThrow( false, ExcMessage( "prepare_tria_for_notching<< The three found points closest to the corner points of the notch coincident "
													"partly or fully, but need to be different. Use a finer mesh.") );
				break;
		}

		// @todo-extent The following limits the orientations of the notch to the xy-plane
		// Find the coordinates of the closest vertex.
		 Point<dim> vertex_closest_positive, vertex_closest_negative, vertex_closest_depth;
		 for (typename Triangulation<dim>::active_cell_iterator
		   cell = triangulation.begin_active();
		   cell != triangulation.end(); ++cell)
		 {
			  for (unsigned int vertex=0; vertex < GeometryInfo<dim>::vertices_per_cell; ++vertex)
			  {
				   unsigned int index_vertex = cell->vertex_index(vertex);
				   // @todo We should also check whether the point is reasonably far away.
				   if ( index_vertex == vertexID_closest_positive )
				   {
					   for ( unsigned int i=0; i<dim; i++ )
						   vertex_closest_positive[i] = cell->vertex(vertex)[i];
				   }
				   else if ( index_vertex == vertexID_closest_negative )
				   {
					   for ( unsigned int i=0; i<dim; i++ )
						   vertex_closest_negative[i] = cell->vertex(vertex)[i];
				   }
				   else if (notch.type==enums::notch_linear && index_vertex == vertexID_closest_depth )
				   {
					   for ( unsigned int i=0; i<dim; i++ )
						   vertex_closest_depth[i] = cell->vertex(vertex)[i];
				   }
			  }
		 }
		 
		// Loop over each vertex to find the above two and also the ones in the third dimension and shift them
		// @todo-optimize There should be a simpler way to find the vertex, when we know its index
		 for (typename Triangulation<dim>::active_cell_iterator
		   cell = triangulation.begin_active();
		   cell != triangulation.end(); ++cell)
		 {
			  for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_cell; ++vertex)
			  {
				  Point<dim> vector_to_positive = cell->vertex(vertex);
				  vector_to_positive -= vertex_closest_positive;
				  double distance2D_to_positive = std::sqrt( vector_to_positive[0]*vector_to_positive[0] + vector_to_positive[1]*vector_to_positive[1] );
				  Point<dim> vector_to_negative = cell->vertex(vertex);
				  vector_to_negative -= vertex_closest_negative;
				  double distance2D_to_negative = std::sqrt( vector_to_negative[0]*vector_to_negative[0] + vector_to_negative[1]*vector_to_negative[1] );
				  Point<dim> vector_to_depth = cell->vertex(vertex);
				  vector_to_depth -= vertex_closest_depth;
				  double distance2D_to_depth = std::sqrt( vector_to_depth[0]*vector_to_depth[0] + vector_to_depth[1]*vector_to_depth[1] );
				  
				   if ( distance2D_to_positive < 1e-8 )
				   {
					   cell->vertex(vertex)[0] = point_corner_positive[0];
					   cell->vertex(vertex)[1] = point_corner_positive[1];
				   }
				   else if ( distance2D_to_negative < 1e-8 )
				   {
					   cell->vertex(vertex)[0] = point_corner_negative[0];
					   cell->vertex(vertex)[1] = point_corner_negative[1];
				   }
				   else if ( notch.type==enums::notch_linear && distance2D_to_depth < 1e-8 )
				   {
					   cell->vertex(vertex)[0] = point_corner_depth[0];
					   cell->vertex(vertex)[1] = point_corner_depth[1];
				   }
			  }
		 }
	}
	
	
	/**
	 * @todo Use spherical manifold or sth on coarse mesh and only move a single node in, then do local refinements and all vertices will follow the curvature
	 * @param triangulation
	 * @param half_notch_length Half the length of the notch in y-direction. We only model 1/8 of the entire bar, hence only 1/2 of the notch length
	 * @param notch_radius The radius of the rod that is left at y=0 in the notch
	 * @param R The radius of the notch corresponds to the tool radius that could be used on a lathe to create the notch.
	 */
	template <int dim>
	void notch_body( Triangulation<dim> &triangulation, const double &half_notch_length, const double &radius, const double &notch_radius,
					 const double &R, const unsigned int notch_type, const bool geom_cylindrical=false )
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
						  if ( dim==2 || ( dim==3 && geom_cylindrical==false ) )
							  vertex_radius=x_coord;
						  else if ( dim==3 && geom_cylindrical )
						  {
							  z_coord = cell->vertex(vertex)[z];
							  vertex_radius = std::sqrt(x_coord*x_coord + z_coord*z_coord);
						  }
							  
						  // The radius of the leftover notched material describes an arc along the y-coord.
						  // Hence, the radius of the notch changes with the y-coord
						   double current_notch_radius = get_current_notch_radius( notch_type, y_coord, half_notch_length, radius, notch_radius, R );
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

	template <int dim>
	class EvalPointClass
	{
	public:
		EvalPointClass() = default;

		// Evaluation point and coordinate direction (x,y,z)
		 EvalPointClass ( Point<3> eval_point , unsigned int direction )
		 :
		 eval_point(eval_point),
		 direction(direction)
		 {
		 }
		// @todo-extent If needed we can also add a Point together with a direction vector for evaluation.
		// This flexibility is actually the reason this class was created and not replaced by a std::pair or similar.

		 Point<3> eval_point;
		 unsigned int direction = 99;

		 double extract_disp_component ( const Vector<double> &pt_solution )
		 {
			 if ( direction == 99 )
				 AssertThrow( false, ExcMessage("EvalPointClass<< You have not declared the evaluation point properly. "
						 	 	 	 	 	    "We require a direction (x,y,z) or a direction vector "
						 	 	 	 	 	    "(the latter is not yet implemented)."));
			 return pt_solution[direction];
		 }
	};
}


#endif // NUMEX_HELPERFNC
