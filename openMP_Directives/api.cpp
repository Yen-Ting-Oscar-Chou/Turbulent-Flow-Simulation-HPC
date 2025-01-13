#include <omp.h>
#include <iostream>

// Compile with llvm-14.0.0 and cuda-11.4.4
// clang++ -fopenmp-targets=nvptx64-nvidia-cuda --cuda-path=$CUDA_HOME -fopenmp mapper.cpp -o mapper.out

class SimpleField {
    public:
    double* array_;
    int size_;
    int value_;

    SimpleField(int size, int value) {
        size_ = size;
        value_ = value;
        array_ = new double[size];
    }

    ~SimpleField() {
        delete[] array_;
    }

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
};

int main() {
    // Initialize objects
    int device = omp_get_default_device();
    std::cout << "Default Device: " << device << std::endl;
    int array_size = 10;
    SimpleField cpuField(array_size, 1);
    
    // Define sizes
    size_t obj_size = sizeof(cpuField);
    size_t array_byte_size = array_size * sizeof(double);
    size_t combined_size = obj_size + array_byte_size;
    std::cout << "Combined Size: " << combined_size << std::endl;

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

    // Run on GPU
    #pragma omp target device(device)
    {
        cpuField.fill();
    }
    std::cout << "Returned from GPU" << std::endl;
    
    // Copy back from GPU to CPU
    omp_target_memcpy(cpuField.array_, gpu_array_ptr, array_byte_size, 0, 0, omp_get_initial_device(), device);
    std::cout << "Copied Array from GPU to CPU" << std::endl;
    for (int i = 0; i < array_size; i++) {
        std::cout << cpuField.array_[i] << std::endl;
    }

    // Deallocate memory and disassociate pointer
    omp_target_disassociate_ptr(&cpuField, device);
    omp_target_free(gpu_array_ptr, device);
    omp_target_free(gpu_field_ptr, device);
}