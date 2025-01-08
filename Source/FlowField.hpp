#pragma once

#include "DataStructures.hpp"
#include "Parameters.hpp"

/** Flow field
 *
 * Class intended to contain the state of the domain.
 */
class FlowField {
private:
  const int sizeX_; //! Size in the X direction
  const int sizeY_; //! Size in the Y direction
  const int sizeZ_; //! Size in the Z direction

  const int cellsX_;
  const int cellsY_;
  const int cellsZ_;

public:
  ScalarField pressure_; //! Scalar field representing the pressure
  VectorField velocity_; //! Multicomponent field representing velocity

  IntScalarField flags_; //! Integer field for the flags

  VectorField FGH_;
  ScalarField RHS_; //! Right hand side for the Poisson equation

  /** Constructor for the 2D flow field
   *
   * Constructor for the flow field. Allocates all the fields and sets
   * the sizes. Currently, this contructor is only used for testing purposes.
   *
   * @param Nx Size of the fuild domain (non-ghost cells), in the X direction
   * @param Ny Size of the fuild domain (non-ghost cells), in the Y direction
   */
  FlowField(int Nx, int Ny);

  /** Constructor for the 3D flow field
   *
   * Constructor for the flow field. Allocates all the fields and sets
   * the sizes. Currently, this contructor is only used for testing purposes.
   *
   * @param Nx Size of the fuild domain (non-ghost cells), in the X direction
   * @param Ny Size of the fuild domain (non-ghost cells), in the Y direction
   * @param Nz Size of the fuild domain (non-ghost cells), in the Z direction
   */
  FlowField(int Nx, int Ny, int Nz);

  /** Constructs a field from parameters object
   *
   * Constructs a field from a parameters object, so that it dimensionality can be defined in
   * the configuration file.
   *
   * @param parameters Parameters object with geometric information
   */
  FlowField(const Parameters& parameters);

  virtual ~FlowField() = default;

  FlowField(FlowField& flowField): 
    sizeX_(flowField.getNx()),
    sizeY_(flowField.getNy()),
    sizeZ_(flowField.getNz()),
    cellsX_(flowField.getCellsX()),
    cellsY_(flowField.getCellsY()),
    cellsZ_(flowField.getCellsZ()),
    pressure_(flowField.getPressure()),
    velocity_(flowField.getVelocity()),
    flags_(flowField.getFlags()),
    FGH_(flowField.getFGH()),
    RHS_(flowField.getRHS()) {}

#pragma omp declare target

  /** Obtain size in the X direction
   *
   * @return Number of cells in the X direction
   */
  int getNx() const;

  /** Obtain size in the Y direction
   *
   * @return Number of cells in the Y direction
   */
  int getNy() const;

  /** Obtain size in the Z direction
   *
   * @return Number of cells in the Z direction
   */
  int getNz() const;

  int getCellsX() const;
  int getCellsY() const;
  int getCellsZ() const;

  ScalarField& getPressure();
  VectorField& getVelocity();

  IntScalarField& getFlags();

  VectorField& getFGH();

  ScalarField& getRHS();

  void getPressureAndVelocity(RealType& pressure, RealType* const velocity, int i, int j);
  void getPressureAndVelocity(RealType& pressure, RealType* const velocity, int i, int j, int k);
#pragma omp end declare target
};
// #pragma omp declare mapper(FlowField f) \
  // map(f) \
  // map(f.pressure_, f.pressure_.data_[0:f.pressure_.size_]) \
  // map(f.RHS_, f.RHS_.data_[0:f.RHS_.size_]) \
  // map(f.velocity_, f.velocity_.data_[0:f.velocity_.size_]) \
  // map(f.FGH_, f.FGH_.data_[0:f.FGH_.size_]) \
  // map(f.flags_, f.flags_.data_[0:f.flags_.size_])
