cmake -B build \
    -DFETCHCONTENT_FULLY_DISCONNECTED=ON \
    -DCMAKE_BUILD_TYPE=Release \
    -DCUPDLP_GPU=ON \
    -DCUPDLP_FIND_CUDA=ON

if [ $? -ne 0 ]; then
    echo "Error: Failed to configure CMake"
    exit 1
fi

cmake --build build --parallel

if [ $? -ne 0 ]; then
    echo "Error: Failed to build the project"
    exit 1
fi


# -DCMAKE_CUDA_HOST_COMPILER=/usr/bin/g++-13 \
# -DCMAKE_C_COMPILER=/usr/bin/gcc-13 \
# -DCMAKE_CXX_COMPILER=/usr/bin/g++-13



# -DCUPDLP_GPU=ON \
    # -DCUPDLP_FIND_CUDA=ON
    # -DCMAKE_CUDA_HOST_COMPILER=/usr/bin/g++-13 \
    # -DCMAKE_C_COMPILER=/usr/bin/gcc-13 \
    # -DCMAKE_CXX_COMPILER=/usr/bin/g++-13