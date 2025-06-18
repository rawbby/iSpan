source /opt/intel/oneapi/setvars.sh
cmake -S . -B cmake-build-test -G Ninja -DBUILD_MPI=ON -DBUILD_TESTING=ON
cmake --build cmake-build-test
ctest --test-dir cmake-build-test --output-on-failure
