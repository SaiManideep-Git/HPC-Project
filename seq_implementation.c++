#include <iostream>
#include <vector>
#include <numeric>

// Define a 3D vector for easier use
using Cube = std::vector<std::vector<std::vector<float>>>;

// Function to initialize a 3D cube with some values
void initializeCube(Cube& data, int size) {
    data.assign(size, std::vector<std::vector<float>>(size, std::vector<float>(size)));
    for (int z = 0; z < size; ++z) {
        for (int y = 0; y < size; ++y) {
            for (int x = 0; x < size; ++x) {
                // Simple initialization for demonstration
                data[z][y][x] = static_cast<float>(x + y + z);
            }
        }
    }
}

// Function to print a slice of the cube for verification
void printSlice(const Cube& data, int z_slice) {
    if (z_slice >= data.size()) return;
    std::cout << "Slice Z = " << z_slice << ":" << std::endl;
    for (const auto& row : data[z_slice]) {
        for (float val : row) {
            std::cout << val << "\t";
        }
        std::cout << std::endl;
    }
}

int main() {
    const int CUBE_SIZE = 10; // A small size for easy verification
    const int STENCIL_RADIUS = 1; // A radius of 1 means a 3x3x3 stencil

    Cube sourceCube, resultCube;

    // Initialize the source data
    initializeCube(sourceCube, CUBE_SIZE);
    initializeCube(resultCube, CUBE_SIZE); // resultCube will store the output

    std::cout << "--- Source Data (Slice 0) ---" << std::endl;
    printSlice(sourceCube, 0);

    // --- Sequential 3D Mean Blur (Stencil Operation) ---
    // This is the core algorithm to be parallelized.
    // The logic mirrors the nested loops in Listing 2 of the paper.
    for (int z = 0; z < CUBE_SIZE; ++z) {
        for (int y = 0; y < CUBE_SIZE; ++y) {
            for (int x = 0; x < CUBE_SIZE; ++x) {
                
                float sum = 0.0f;
                int count = 0;

                // Iterate over the stencil area defined by the radius
                for (int mz = z - STENCIL_RADIUS; mz <= z + STENCIL_RADIUS; ++mz) {
                    for (int my = y - STENCIL_RADIUS; my <= y + STENCIL_RADIUS; ++my) {
                        for (int mx = x - STENCIL_RADIUS; mx <= x + STENCIL_RADIUS; ++mx) {
                            
                            // Boundary check: ensure we are within the cube's bounds
                            if (mz >= 0 && mz < CUBE_SIZE &&
                                my >= 0 && my < CUBE_SIZE &&
                                mx >= 0 && mx < CUBE_SIZE) {
                                
                                sum += sourceCube[mz][my][mx];
                                count++;
                            }
                        }
                    }
                }
                
                // Calculate the mean and store it in the result cube
                if (count > 0) {
                    resultCube[z][y][x] = sum / count;
                }
            }
        }
    }

    std::cout << "\n--- Result Data after Mean Blur (Slice 0) ---" << std::endl;
    printSlice(resultCube, 0);

    return 0;
}