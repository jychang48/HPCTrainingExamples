Adding portable makefiles and cmake builds.

For ROCm environment:
   ```
   module load rocm
   module load cmake
   export CXX=${ROCM_PATH}/llvm/bin/clang++
   ```

For ROCm with `make`:
   ```
   make
   ```

For ROCm with `cmake`:
   ```
   mkdir build && cd build
   cmake ..
   make VERBOSE=1
   ./saxpy
   ctest
   ```

For CUDA environment:
   ```
   module load rocm
   module load CUDA/11.8
   module load cmake
   ```

For CUDA with `make`:
   ```
   HIP_PLATFORM=nvidia make
   ```

For CUDA with `cmake`:
   ```
   mkdir build && cd build
   cmake -DCMAKE_GPU_RUNTIME=CUDA ..
   make VERBOSE=1
   ./saxpy
   ctest
   ```
