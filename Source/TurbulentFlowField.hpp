#pragma once

#include "FlowField.hpp"

class TurbulentFlowField : FlowField {
private:
  ScalarField viscosity_;
  ScalarField distance_;

public:
  /** Constructs a field from parameters object
   *
   * Constructs a field from a parameters object, so that it dimensionality can be defined in
   * the configuration file.
   *
   * @param parameters Parameters object with geometric information
   */
  TurbulentFlowField(const Parameters& parameters);

  virtual ~TurbulentFlowField() override = default;

  ScalarField& getViscosity();
  ScalarField& getDistance();

  void getDistanceAndViscosity(RealType& distance, RealType& viscosity, int i, int j);
  void getDistanceAndViscosity(RealType& distance, RealType& viscosity, int i, int j, int k);
};