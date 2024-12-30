#include <complex>
#include <memory>
#include <stdio.h>
#include <stdlib.h>
#include <vector>


#define N 100000
// Compile with llvm-14.0.0 and cuda-11.4.4
// clang++ -fopenmp-targets=nvptx64-nvidia-cuda --cuda-path=$CUDA_HOME -fopenmp mapper.cpp -o mapper.out

// GENERAL LEARNINGS:
// - no cpp functions in device code; hpp only

class Component {
private:
  int compNum_;

public:
  size_t  len_;
  double* data_;
  Component(int num, size_t len):
    compNum_(num),
    len_(len) {
    data_ = new double[len_];
  }

  int getComponentNumber() { return compNum_; }
};
// LEARNING this produces a compiler SegFault in conjunction with the other mapper
// #pragma omp declare mapper(Component c) map(c, c.data_[0 : c.len_])

class Aggregate {
private:
  int aggNum_;

public:
  size_t  len_;
  double* data_;
  Component comp_;
  Aggregate(int num, size_t len, Component& comp):
    aggNum_(num),
    len_(len),
    comp_(comp) {
    data_ = new double[len_];
  }

  int getAggregateNumber() { return aggNum_; }
};
#pragma omp declare mapper(Aggregate a) map(a, a.data_[0 : a.len_]) map(a.comp_, a.comp_.data_[0 : a.comp_.len_])

class Runner {
private:
  // LEARNING this procudes a memory error. References are not allowed
  // Aggregate& agg_;
  Component comp_;
  Aggregate agg_;

public:
  Runner(Aggregate& agg, Component& comp):
  agg_(agg),
  comp_(comp) {}

  // LEARNING this produces a memory error. Virtual functions not allowed
  // virtual void data_process(Aggregate& agg, int i) {
  void data_process(Aggregate& agg, int i) {
    agg.getAggregateNumber();
    agg.comp_.data_[i] = agg.comp_.getComponentNumber() * i + 1;
    return;
  }

  virtual void run() {
    #pragma omp target map(agg_)
    for (int i = 0; i < agg_.len_; i++) {
      data_process(agg_, i);
    }
    printf("agg_.comp_.data_[%d]=%lf\n", N - 1, agg_.comp_.data_[N - 1]);
  }
};

int main() {
  Component comp(5, N);
  Aggregate* agg = new Aggregate(0, N, comp);
  Runner* runner = new Runner(*agg, comp);
  runner->run();
}
