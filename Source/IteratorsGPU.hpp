#pragma once

#include "FlowField.hpp"
#include "Parameters.hpp"
#include "Stencils/StencilDelegate.hpp"
#include "TurbulentFlowField.hpp"

/** GPUIterator class
 *
 * Applies operations to a flow field
 */
template <class FlowFieldType>
class GPUIterator {
protected:
  const Parameters& parameters_;

public:
  FlowFieldType& flowField_;

  GPUIterator(FlowFieldType& flowfield, const Parameters& parameters):
    flowField_(flowfield),
    parameters_(parameters) {}

  virtual ~GPUIterator() = default;

  /** Perform the stencil operation on inner, non-ghost cells.
   */
  virtual void iterate() = 0;
};

template <class FlowFieldType>
class GPUFieldIterator: public GPUIterator<FlowFieldType> {
private:
  StencilDelegate& stencil_;

  //@brief Define the iteration domain to include more or less layers
  // Added since the ability to select the iteration domain provides more flexibility
  //@{
  const int lowOffset_;
  const int highOffset_;
  //@}

public:
  GPUFieldIterator(FlowFieldType& flowField, const Parameters& parameters, StencilDelegate& stencil, int lowOffset = 0, int highOffset = 0);

  virtual ~GPUFieldIterator() override = default;

  virtual void iterate() override;
};
