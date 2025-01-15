#include "StdAfx.hpp"

#include "GlobalBoundaryFactory.hpp"

GlobalBoundaryFactory::GlobalBoundaryFactory(Parameters& parameters):
  parameters_(parameters) {
  // The parameters will be modified, and therefore are not declared as constants.

  // All stencils are created, disregarding whether they will be used or not. This is less
  // complicated and doesn't seem that costly.

  periodic_[0] = new Stencils::PeriodicBoundaryVelocityStencil();
  periodic_[1] = new Stencils::PeriodicBoundaryFGHStencil();

  moving_[0] = new Stencils::MovingWallVelocityStencil();
  moving_[1] = new Stencils::MovingWallFGHStencil();

  outflow_[0] = new Stencils::NeumannVelocityBoundaryStencil();
  outflow_[1] = new Stencils::NeumannFGHBoundaryStencil();

  channelInput_[0] = new Stencils::BFInputVelocityStencil();
  channelInput_[1] = new Stencils::BFInputFGHStencil();

  if (parameters.simulation.scenario == CAVITY) {
    // Here, all is about setting the velocity at the boundaries.
    for (int i = 0; i < 6; i++) {
      velocityStencils_[i] = moving_[0];
      FGHStencils_[i]      = moving_[1];
    }
    parameters.walls.typeLeft   = DIRICHLET;
    parameters.walls.typeRight  = DIRICHLET;
    parameters.walls.typeBottom = DIRICHLET;
    parameters.walls.typeTop    = DIRICHLET;
    parameters.walls.typeFront  = DIRICHLET;
    parameters.walls.typeBack   = DIRICHLET;
  } else if (parameters.simulation.scenario == CHANNEL) {
    // To the left, we have the input
    velocityStencils_[0] = channelInput_[0];
    FGHStencils_[0]      = channelInput_[1];

    // To the right, there is an outflow boundary
    velocityStencils_[1] = outflow_[0];
    FGHStencils_[1]      = outflow_[1];

    // The other walls are moving walls
    for (int i = 2; i < 6; i++) {
      velocityStencils_[i] = moving_[0];
      FGHStencils_[i]      = moving_[1];
    }
    parameters.walls.typeLeft   = DIRICHLET;
    parameters.walls.typeRight  = NEUMANN;
    parameters.walls.typeBottom = DIRICHLET;
    parameters.walls.typeTop    = DIRICHLET;
    parameters.walls.typeFront  = DIRICHLET;
    parameters.walls.typeBack   = DIRICHLET;
  } else {
    throw std::runtime_error("Scenario not recognized");
  }
}

GlobalBoundaryFactory::~GlobalBoundaryFactory() {
  delete moving_[0];
  delete moving_[1];

  delete periodic_[0];
  delete periodic_[1];

  delete outflow_[0];
  delete outflow_[1];

  delete channelInput_[0];
  delete channelInput_[1];
}

GlobalBoundaryIterator<FlowField> GlobalBoundaryFactory::getGlobalBoundaryFGHIterator(FlowField& flowField) {
  if (parameters_.geometry.dim == 2) {
    return GlobalBoundaryIterator<FlowField>(flowField, parameters_, *(FGHStencils_[0]), *(FGHStencils_[1]), *(FGHStencils_[2]), *(FGHStencils_[3]), 1, 0);
  }
  return GlobalBoundaryIterator<FlowField>(
    flowField, parameters_, *(FGHStencils_[0]), *(FGHStencils_[1]), *(FGHStencils_[2]), *(FGHStencils_[3]), *(FGHStencils_[4]), *(FGHStencils_[5]), 1, 0
  );
}

GlobalBoundaryIterator<FlowField> GlobalBoundaryFactory::getGlobalBoundaryVelocityIterator(FlowField& flowField) {
  if (parameters_.geometry.dim == 2) {
    return GlobalBoundaryIterator<FlowField>(flowField, parameters_, *(velocityStencils_[0]), *(velocityStencils_[1]), *(velocityStencils_[2]), *(velocityStencils_[3]), 1, 0);
  }
  return GlobalBoundaryIterator<FlowField>(
    flowField, parameters_, *(velocityStencils_[0]), *(velocityStencils_[1]), *(velocityStencils_[2]), *(velocityStencils_[3]), *(velocityStencils_[4]), *(velocityStencils_[5]), 1, 0
  );
}
