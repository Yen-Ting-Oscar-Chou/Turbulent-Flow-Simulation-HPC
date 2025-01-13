#include <iostream>
#include <omp.h>

// Compile with llvm-14.0.0 and cuda-11.4.4
// clang++ -fopenmp-targets=nvptx64-nvidia-cuda --cuda-path=$CUDA_HOME -fopenmp mapper.cpp -o mapper.out
class SimpleField;

struct SimpleFieldGPUPtrs {
  SimpleField* gpu_field_ptr_;
  double* gpu_array_ptr_;
  SimpleFieldGPUPtrs(SimpleField* gpu_field_ptr, double* gpu_array_ptr) {
    gpu_field_ptr_ = gpu_field_ptr;
    gpu_array_ptr_ = gpu_array_ptr;
  }
};

class SimpleField {
private:
  double*      array_;
  int          size_;
  int          value_;

public:

  SimpleField(int size, int value) {
    size_  = size;
    value_ = value;
    array_ = new double[size];
  }

  ~SimpleField() { delete[] array_; }

  double* getArray() { return array_; }

#pragma omp begin declare target
  void fill() {
    // print the value of size_ and value_
    printf("Size: %d\n", size_);
    printf("Value: %d\n", value_);
    printf("A\n");
    for (int i = 0; i < size_; i++) {

      array_[i] = value_;
    }
    printf("B\n");
  }
#pragma omp end declare target

static SimpleFieldGPUPtrs map_to_gpu(int device, SimpleField& cpuField) {
    size_t obj_size        = sizeof(cpuField);
    size_t array_byte_size = cpuField.size_ * sizeof(double);

    // Allocate memory on the GPU
    auto gpu_field_ptr = static_cast<SimpleField*>(omp_target_alloc(obj_size, device));
    auto gpu_array_ptr = static_cast<double*>(omp_target_alloc(array_byte_size, device));

    // Associate the cpuField with the gpu_field_ptr
    bool associated_field = omp_target_associate_ptr(&cpuField, gpu_field_ptr, obj_size, 0, device) == 0;

    std::cout << "Successfully Associated Field: " << (associated_field ? "True" : "False") << std::endl;
    std::cout << gpu_field_ptr << std::endl;

    // Copy the SimpleField object
    bool copied_field = omp_target_memcpy(gpu_field_ptr, &cpuField, obj_size, 0, 0, device, omp_get_initial_device()) == 0;
    std::cout << "MEMCPY of Field to GPU Successful: " << (copied_field ? "True" : "False") << std::endl;

    // Copy the array to the GPU
    bool copied_array = omp_target_memcpy(gpu_array_ptr, &cpuField.array_, array_byte_size, 0, 0, device, omp_get_initial_device()) == 0;
    std::cout << "MEMCPY of Array to GPU Successful: " << (copied_array ? "True" : "False") << std::endl;

    // Update the device copy of SimpleField's array_ pointer
    omp_target_memcpy(&gpu_field_ptr->array_, &gpu_array_ptr, sizeof(double*), 0, 0, device, omp_get_initial_device());

    // Create and return SimpleFieldGPUPtrs
    SimpleFieldGPUPtrs gpuPtrs(gpu_field_ptr, gpu_array_ptr);
    return gpuPtrs;
  }

  static void map_to_cpu(int device, SimpleField& cpuField, SimpleFieldGPUPtrs ptrs) {
    size_t obj_size        = sizeof(cpuField);
    size_t array_byte_size = cpuField.size_ * sizeof(double);

    for (int i = 0; i < cpuField.size_; i++) {
      std::cout << cpuField.array_[i] << std::endl;
    }

    std::cout << "Copying Array from GPU to CPU" << std::endl;

    // Copy back from GPU to CPU
    omp_target_memcpy(cpuField.array_, ptrs.gpu_array_ptr_, array_byte_size, 0, 0, omp_get_initial_device(), device);
    std::cout << "Copied Array from GPU to CPU" << std::endl;
    for (int i = 0; i < cpuField.size_; i++) {
      std::cout << cpuField.array_[i] << std::endl;
    }

    // Deallocate memory and disassociate pointer
    omp_target_disassociate_ptr(&cpuField, device);
    omp_target_free(ptrs.gpu_array_ptr_, device);
    omp_target_free(ptrs.gpu_field_ptr_, device);
  }
};

int main() {
  // Initialize objects
  int device = omp_get_default_device();
  std::cout << "Default Device: " << device << std::endl;
  int         array_size = 10;
  SimpleField cpuField(array_size, 1);

  SimpleFieldGPUPtrs ptrs = SimpleField::map_to_gpu(device, cpuField);

  // Run on GPU
#pragma omp target device(device)
  { cpuField.fill(); }
  std::cout << "Returned from GPU" << std::endl;
  
  SimpleField::map_to_cpu(device, cpuField, ptrs);
}