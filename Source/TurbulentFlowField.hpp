#pragma once

#include "FlowField.hpp"

class TurbulentFlowField: public FlowField {
private:

public:
  ScalarField viscosity_;
  ScalarField distance_;
  /** Constructs a field from parameters object
   *
   * Constructs a field from a parameters object, so that it dimensionality can be defined in
   * the configuration file.
   *
   * @param parameters Parameters object with geometric information
   */
  TurbulentFlowField(const Parameters& parameters);

  virtual ~TurbulentFlowField() = default;

#pragma omp declare target
  ScalarField& getViscosity();
  ScalarField& getDistance();

  void getDistanceAndViscosity(RealType& distance, RealType& viscosity, int i, int j);
  void getDistanceAndViscosity(RealType& distance, RealType& viscosity, int i, int j, int k);
#pragma omp end declare target
};
/* #pragma omp declare mapper(TurbulentFlowFieldMap: TurbulentFlowField f) \
  map(mapper(ScalarFieldMap), tofrom: f.pressure_, f.RHS_, f.distance_, f.viscosity_) \
  map(mapper(VectorFieldMap), tofrom: f.velocity_, f.FGH_) \
  map(mapper(IntScalarMap), tofrom: f.flags_)  */