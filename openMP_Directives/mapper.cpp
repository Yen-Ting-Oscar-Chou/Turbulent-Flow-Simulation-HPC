#include <complex>
#include <memory>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <cmath>
#include <iostream>


#define N 100000
// Compile with llvm-14.0.0 and cuda-11.4.4
// clang++ -fopenmp-targets=nvptx64-nvidia-cuda --cuda-path=$CUDA_HOME -fopenmp mapper.cpp -o mapper.out

// GENERAL LEARNINGS:
// - no cpp functions in device code; hpp only (nvlink error) -> use #pragma omp declare target and #pragma omp end declare target

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

template <class T>
class Aggregate {
private:
  T aggNum_;

public:
  Component* ptr_ = NULL;
  size_t  len_;
  double* data_;
  Component comp_;

  Aggregate(T num, size_t len, Component& comp):
    aggNum_(num),
    len_(len),
    comp_(comp) {
    data_ = new double[len_];
    ptr_ = new Component(num, len);
  }

  int getAggregateNumber() { return aggNum_; }
};
#pragma omp declare mapper(Aggregate<int> a) map(a, a.data_[0 : a.len_], a.ptr_[0:1]) map(a.comp_, a.comp_.data_[0 : a.comp_.len_])

class Runner {
private:
  // LEARNING this procudes a memory error. References are not allowed
  // Aggregate& agg_;
  Component comp_;

public:
  Aggregate<int> agg_;
  
  Runner(Aggregate<int>& agg, Component& comp):
  agg_(agg),
  comp_(comp) {}

  // LEARNING this produces a memory error. Virtual functions not allowed
  // virtual void data_process(Aggregate& agg, int i) {
  void data_process(Aggregate<int>& agg, int num, int i) {
    agg.getAggregateNumber();
    agg.comp_.data_[i] = num;
    return;
  }

  virtual void run(int num) {
    #pragma omp target enter data map(to: agg_)
    //#pragma omp target map(agg_)
    {
      for (int i = 0; i < agg_.len_; i++) {
        data_process(agg_, num, i);
      }   
      printf("agg_.comp_.data_[%d]=%lf\n", N - 1, agg_.comp_.data_[N - 1]);
    }
    #pragma omp target update to(agg_)
    #pragma omp target exit data map(from: agg_)
  }
};

int main() {
  Component comp(5, N);
  Aggregate<int>* agg = new Aggregate<int>(0, N, comp);
  Runner* runner = new Runner(*agg, comp);
  runner->run(1);
  std::cout << "Run 1: " << runner->agg_.comp_.data_[N - 1] << std::endl;
  runner->run(2);
  std::cout << "Run 2: " << runner->agg_.comp_.data_[N - 1] << std::endl;
}
