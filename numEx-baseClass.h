#ifndef NUMEX_BASECLASS_H
#define NUMEX_BASECLASS_H

/**
 *
 */
template<int dim>
class numEx_class
{
  public:
    double search_tolerance = 1e-12;

	virtual std::string numEx_name() {
		return "numEx";
	};
	
	virtual unsigned int loading_direction() {
		return enums::coord_none;
	};
	
	virtual std::vector< enums::enum_boundary_ids > id_boundary_loads() {
		return std::vector< enums::enum_boundary_ids > ();
	}

	// @todo Check whether this virtual function is needed, because this is only called internally
	virtual std::vector< types::manifold_id > manifold_ids() {
		return std::vector< types::manifold_id > ();
	}

	// abstract function:
	// Is never executed, thus must be overwritten in each derived class
	 virtual void make_grid( /*input-> */ const Parameter::GeneralParameters &parameter,
							 /*output->*/ Triangulation<dim> &triangulation,
										  std::vector<double> &body_dimensions,
										  std::vector< numEx::EvalPointClass<3> > &eval_points_list,
										  const std::string relativePath ) =0;

	 virtual void make_constraints ( AffineConstraints<double> &constraints, const FESystem<dim> &fe, DoFHandler<dim> &dof_handler_ref,
	    							const bool &apply_dirichlet_bc, double &current_load_increment, const Parameter::GeneralParameters &parameter ) =0;
};

#endif //NUMEX_BASECLASS_H
