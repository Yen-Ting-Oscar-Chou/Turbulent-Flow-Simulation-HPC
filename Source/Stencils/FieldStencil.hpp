#pragma once

#include "Parameters.hpp"

namespace Stencils {

  /** Stencil class
   *
   * Abstract class for the definition of stencils and operations on the grids
   */
  template <class FlowFieldType>
  class FieldStencil {

  public:
    FieldStencil() = default;

    virtual ~FieldStencil() = default;

    /** Performs the operation in 2D in a given position
     * @param flowField Flow field data
     * @param parameters Parameters of the problem
     * @param i Position in the x direction
     * @param j Position in the y direction
     */
    virtual void apply(const Parameters& parameters, FlowFieldType& flowField, int i, int j) = 0;

    /** Performs the operation in 3D in a given position
     * @param flowField Flow field data
     * @param parameters Parameters of the problem
     * @param i Position in the x direction
     * @param j Position in the y direction
     * @param k Position in the z direction
     */
    virtual void apply(const Parameters& parameters, FlowFieldType& flowField, int i, int j, int k) = 0;
  };

} // namespace Stencils
