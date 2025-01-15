#include "StdAfx.hpp"

#include "BFInputStencils.hpp"

RealType computeVelocity3D(int i, int j, int k, RealType stepSize, const Parameters& parameters) {
  const RealType posY = parameters.meshsize.getPosY(i, j, k);
  const RealType posZ = parameters.meshsize.getPosZ(i, j, k);
  const RealType dy   = parameters.meshsize.getDy(i, j, k);
  const RealType dz   = parameters.meshsize.getDz(i, j, k);

  if (posY + 0.5 * dy >= stepSize) {
    // Get the size of the inlet in Y. A 3 is subtracted because of the boundary cells.
    const RealType inletYSize = parameters.geometry.lengthY - stepSize;
    const RealType inletZSize = parameters.geometry.lengthZ;

    const RealType y = posY + 0.5 * dy - stepSize;
    const RealType z = posZ + 0.5 * dz;

    return 36.0 * parameters.walls.vectorLeft[0] / (inletZSize * inletZSize * inletYSize * inletYSize) * y * (y - inletYSize) * z * (z - inletZSize);
  } else {
    return 0.0;
  }
}

RealType computeVelocity2D(int i, int j, RealType stepSize, const Parameters& parameters) {
  const RealType posY = parameters.meshsize.getPosY(i, j);
  const RealType dy   = parameters.meshsize.getDy(i, j);

  if (posY + 0.5 * dy >= stepSize) {
    // Get the size of the inlet in Y. A 3 is subtracted because of the boundary cells.
    [[maybe_unused]] const RealType inletYSize = parameters.geometry.lengthY - stepSize;
    [[maybe_unused]] const RealType y          = posY + 0.5 * dy - stepSize;

    if (parameters.simulation.velocityProfile == PARABOLIC) {
      return 6.0 * parameters.walls.vectorLeft[0] / (inletYSize * inletYSize) * y * (inletYSize - y); // Parabolic inflow
    }
    return parameters.walls.vectorLeft[0]; // Uniform inflow
  } else {
    return 0.0;
  }
}

Stencils::BFInputVelocityStencil::BFInputVelocityStencil():
  stepSize_(-1) {}

void Stencils::BFInputVelocityStencil::initStepSize(const Parameters& parameters) {
  stepSize_ = parameters.bfStep.yRatio > 0.0 ? parameters.bfStep.yRatio * parameters.geometry.lengthY : 0.0;

  if (parameters.geometry.dim == 2) {
    RealType posY   = parameters.meshsize.getPosY(0, 0);
    RealType dy     = parameters.meshsize.getDy(0, 0);
    RealType nextDy = parameters.meshsize.getDy(0, 1);

    for (int j = 0; j < parameters.geometry.sizeY - 1; ++j) {
      posY   = parameters.meshsize.getPosY(0, j);
      dy     = parameters.meshsize.getDy(0, j);
      nextDy = parameters.meshsize.getDy(0, j + 1);

      // Check if stepSize is in this cell
      if (posY + 0.5 * dy < stepSize_ && stepSize_ <= posY + dy + 0.5 * nextDy) {
        stepSize_ = posY + dy;
        break;
      }
    }
  } else if (parameters.geometry.dim == 3) {
    RealType posY   = parameters.meshsize.getPosY(0, 0, 0);
    RealType dy     = parameters.meshsize.getDy(0, 0, 0);
    RealType nextDy = parameters.meshsize.getDy(0, 1, 0);

    for (int j = 0; j < parameters.geometry.sizeY - 1; ++j) {
      posY   = parameters.meshsize.getPosY(0, j, 0);
      dy     = parameters.meshsize.getDy(0, j, 0);
      nextDy = parameters.meshsize.getDy(0, j + 1, 0);

      if (posY + 0.5 * dy < stepSize_ && stepSize_ <= posY + dy + 0.5 * nextDy) {
        stepSize_ = posY + dy;
        break;
      }
    }
  }
}

void Stencils::BFInputVelocityStencil::applyLeftWall(const Parameters& parameters, FlowField& flowField, int i, int j) {
  if(stepSize_ == -1) {
    initStepSize(parameters);
  }
  flowField.getVelocity().getVectorElement(i, j, 0) = computeVelocity2D(i, j, stepSize_, parameters);
  flowField.getVelocity().getVectorElement(i, j, 1) = -flowField.getVelocity().getVectorElement(i + 1, j, 1);
}

// Most of the functions are empty, and they shouldn't be called, assuming that the input is always located at the left.
void Stencils::BFInputVelocityStencil::applyRightWall(
  [[maybe_unused]] const Parameters& parameters, [[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j
) {}
void Stencils::BFInputVelocityStencil::applyBottomWall(
  [[maybe_unused]] const Parameters& parameters, [[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j
) {}
void Stencils::BFInputVelocityStencil::applyTopWall(
  [[maybe_unused]] const Parameters& parameters, [[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j
) {}

void Stencils::BFInputVelocityStencil::applyLeftWall(const Parameters& parameters, FlowField& flowField, int i, int j, int k) {
  if(stepSize_ == -1) {
    initStepSize(parameters);
  }
  flowField.getVelocity().getVectorElement(i, j, k, 0) = computeVelocity3D(i, j, k, stepSize_, parameters);
  flowField.getVelocity().getVectorElement(i, j, k, 1) = -flowField.getVelocity().getVectorElement(i + 1, j, k, 1);
  flowField.getVelocity().getVectorElement(i, j, k, 2) = -flowField.getVelocity().getVectorElement(i + 1, j, k, 2);
}

void Stencils::BFInputVelocityStencil::applyRightWall(
  [[maybe_unused]] const Parameters& parameters, [[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k
) {}
void Stencils::BFInputVelocityStencil::applyBottomWall(
  [[maybe_unused]] const Parameters& parameters, [[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k
) {}
void Stencils::BFInputVelocityStencil::applyTopWall(
  [[maybe_unused]] const Parameters& parameters, [[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k
) {}
void Stencils::BFInputVelocityStencil::applyFrontWall(
  [[maybe_unused]] const Parameters& parameters, [[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k
) {}
void Stencils::BFInputVelocityStencil::applyBackWall(
  [[maybe_unused]] const Parameters& parameters, [[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k
) {}

Stencils::BFInputFGHStencil::BFInputFGHStencil():
  stepSize_(-1) {}

inline void Stencils::BFInputFGHStencil::initStepSize(const Parameters& parameters) {
  stepSize_ = parameters.bfStep.yRatio > 0.0 ? parameters.bfStep.yRatio * parameters.geometry.lengthY : 0.0;
}

void Stencils::BFInputFGHStencil::applyLeftWall(const Parameters& parameters, FlowField& flowField, int i, int j) {
  if(stepSize_ == -1) {
    initStepSize(parameters);
  }
  flowField.getFGH().getVectorElement(i, j, 0) = computeVelocity2D(i, j, stepSize_, parameters);
}

void Stencils::BFInputFGHStencil::applyRightWall(
  [[maybe_unused]] const Parameters& parameters, [[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j
) {}
void Stencils::BFInputFGHStencil::applyBottomWall(
  [[maybe_unused]] const Parameters& parameters, [[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j
) {}
void Stencils::BFInputFGHStencil::applyTopWall(
  [[maybe_unused]] const Parameters& parameters, [[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j
) {}

void Stencils::BFInputFGHStencil::applyLeftWall(const Parameters& parameters, FlowField& flowField, int i, int j, int k) {
  if(stepSize_ == -1) {
    initStepSize(parameters);
  }
  flowField.getFGH().getVectorElement(i, j, k, 0) = computeVelocity3D(i, j, k, stepSize_, parameters);
}

void Stencils::BFInputFGHStencil::applyRightWall(
  [[maybe_unused]] const Parameters& parameters, [[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k
) {}
void Stencils::BFInputFGHStencil::applyBottomWall(
  [[maybe_unused]] const Parameters& parameters, [[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k
) {}
void Stencils::BFInputFGHStencil::applyTopWall(
  [[maybe_unused]] const Parameters& parameters, [[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k
) {}
void Stencils::BFInputFGHStencil::applyFrontWall(
  [[maybe_unused]] const Parameters& parameters, [[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k
) {}
void Stencils::BFInputFGHStencil::applyBackWall(
  [[maybe_unused]] const Parameters& parameters, [[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k
) {}
