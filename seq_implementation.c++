#include <iostream>
#include <vector>
#include <chrono>
#include <random>
#include <iomanip>
#include <cmath>
using namespace std;
using namespace std::chrono;
// 3D data structure using flat array for better memory layout
class Cube3D {
private:
    vector<float> data;
    int size_x, size_y, size_z;
    
public:
    Cube3D(int x, int y, int z) : size_x(x), size_y(y), size_z(z) {
        data.resize(x * y * z, 0.0f);
    }
    
    inline float& operator()(int x, int y, int z) {
        return data[z * size_x * size_y + y * size_x + x];
    }
    
    inline const float& operator()(int x, int y, int z) const {
        return data[z * size_x * size_y + y * size_x + x];
    }
    
    inline bool isValidIndex(int x, int y, int z) const {
        return x >= 0 && x < size_x && y >= 0 && y < size_y && z >= 0 && z < size_z;
    }
    
    int getSizeX() const { return size_x; }
    int getSizeY() const { return size_y; }
    int getSizeZ() const { return size_z; }
    size_t getTotalSize() const { return data.size(); }
};

// Initialize cube with test data
void initializeCube(Cube3D& cube, int pattern = 0) {
    random_device rd;
    mt19937 gen(42); // Fixed seed for reproducible results
    uniform_real_distribution<float> dis(0.0f, 100.0f);
    
    for (int z = 0; z < cube.getSizeZ(); ++z) {
        for (int y = 0; y < cube.getSizeY(); ++y) {
            for (int x = 0; x < cube.getSizeX(); ++x) {
                switch(pattern) {
                    case 0: // Random values
                        cube(x, y, z) = dis(gen);
                        break;
                    case 1: // Linear gradient
                        cube(x, y, z) = static_cast<float>(x + y + z);
                        break;
                    case 2: // Sine wave pattern
                        cube(x, y, z) = sin(x * 0.2f) * cos(y * 0.2f) * sin(z * 0.2f);
                        break;
                    default:
                        cube(x, y, z) = 1.0f;
                }
            }
        }
    }
}

// 3D Mean Blur Stencil Operation
void meanBlurStencil(const Cube3D& source, Cube3D& result, int radius) {
    const float neutral_value = 0.0f;
    
    for (int z = 0; z < source.getSizeZ(); ++z) {
        for (int y = 0; y < source.getSizeY(); ++y) {
            for (int x = 0; x < source.getSizeX(); ++x) {
                
                float sum = 0.0f;
                int count = 0;
                
                // Apply stencil centered at (x, y, z)
                for (int mz = z - radius; mz <= z + radius; ++mz) {
                    for (int my = y - radius; my <= y + radius; ++my) {
                        for (int mx = x - radius; mx <= x + radius; ++mx) {
                            
                            if (source.isValidIndex(mx, my, mz)) {
                                sum += source(mx, my, mz);
                                count++;
                            } else {
                                sum += neutral_value;
                                count++;
                            }
                        }
                    }
                }
                
                result(x, y, z) = (count > 0) ? sum / count : neutral_value;
            }
        }
    }
}

// Lattice Boltzmann Method (D3Q19 model)
class LBMSolver {
private:
    static const int Q = 19; // Number of velocity directions
    
    // D3Q19 velocity vectors
    const int velocities[Q][3] = {
        {0,0,0},   // Rest particle
        {1,0,0}, {-1,0,0}, {0,1,0}, {0,-1,0}, {0,0,1}, {0,0,-1},
        {1,1,0}, {-1,-1,0}, {1,-1,0}, {-1,1,0},
        {1,0,1}, {-1,0,-1}, {1,0,-1}, {-1,0,1},
        {0,1,1}, {0,-1,-1}, {0,1,-1}, {0,-1,1}
    };
    
    // Weights for D3Q19 model
    const float weights[Q] = {
        1.0f/3.0f,
        1.0f/18.0f, 1.0f/18.0f, 1.0f/18.0f, 1.0f/18.0f, 1.0f/18.0f, 1.0f/18.0f,
        1.0f/36.0f, 1.0f/36.0f, 1.0f/36.0f, 1.0f/36.0f,
        1.0f/36.0f, 1.0f/36.0f, 1.0f/36.0f, 1.0f/36.0f,
        1.0f/36.0f, 1.0f/36.0f, 1.0f/36.0f, 1.0f/36.0f
    };
    
    float tau;       // Relaxation parameter
    float delta_t;   // Time step
    
public:
    LBMSolver(float relaxation_time, float time_step) 
        : tau(relaxation_time), delta_t(time_step) {}
    
    // Single LBM time step (collision + streaming)
    void timeStep(const vector<Cube3D>& source, vector<Cube3D>& result) {
        int size_x = source[0].getSizeX();
        int size_y = source[0].getSizeY();
        int size_z = source[0].getSizeZ();
        
        for (int z = 0; z < size_z; ++z) {
            for (int y = 0; y < size_y; ++y) {
                for (int x = 0; x < size_x; ++x) {
                    
                    // Streaming step: gather distributions from neighbors
                    vector<float> f(Q);
                    for (int i = 0; i < Q; ++i) {
                        int source_x = x - velocities[i][0];
                        int source_y = y - velocities[i][1];
                        int source_z = z - velocities[i][2];
                        
                        if (source_x >= 0 && source_x < size_x &&
                            source_y >= 0 && source_y < size_y &&
                            source_z >= 0 && source_z < size_z) {
                            f[i] = source[i](source_x, source_y, source_z);
                        } else {
                            f[i] = 0.0f; // Boundary condition
                        }
                    }
                    
                    // Calculate macroscopic quantities
                    float density = 0.0f;
                    float velocity[3] = {0.0f, 0.0f, 0.0f};
                    
                    for (int i = 0; i < Q; ++i) {
                        density += f[i];
                        velocity[0] += velocities[i][0] * f[i];
                        velocity[1] += velocities[i][1] * f[i];
                        velocity[2] += velocities[i][2] * f[i];
                    }
                    
                    if (density > 1e-10f) {
                        velocity[0] /= density;
                        velocity[1] /= density;
                        velocity[2] /= density;
                    }
                    
                    // Collision step: BGK approximation
                    for (int i = 0; i < Q; ++i) {
                        float u_dot_ci = velocity[0] * velocities[i][0] + 
                                        velocity[1] * velocities[i][1] + 
                                        velocity[2] * velocities[i][2];
                        
                        float u_square = velocity[0] * velocity[0] + 
                                        velocity[1] * velocity[1] + 
                                        velocity[2] * velocity[2];
                        
                        // Equilibrium distribution function
                        float f_eq = weights[i] * density * 
                                    (1.0f + 3.0f * u_dot_ci + 
                                     4.5f * u_dot_ci * u_dot_ci - 
                                     1.5f * u_square);
                        
                        // BGK collision operator
                        result[i](x, y, z) = f[i] + delta_t / tau * (f_eq - f[i]);
                    }
                }
            }
        }
    }
    
    // Run multiple time steps
    void simulate(vector<Cube3D>& distributions, int num_steps) {
        vector<Cube3D> temp_distributions(Q, Cube3D(distributions[0].getSizeX(),
                                                    distributions[0].getSizeY(),
                                                    distributions[0].getSizeZ()));
        
        for (int step = 0; step < num_steps; ++step) {
            if (step % 2 == 0) {
                timeStep(distributions, temp_distributions);
            } else {
                timeStep(temp_distributions, distributions);
            }
        }
        
        // Ensure final result is in the original array
        if (num_steps % 2 == 1) {
            for (int i = 0; i < Q; ++i) {
                distributions[i] = temp_distributions[i];
            }
        }
    }
};

// Utility functions
void printSlice(const Cube3D& cube, int z_slice, int max_display = 8) {
    if (z_slice >= cube.getSizeZ()) return;
    
    cout << "Slice Z = " << z_slice << ":" << endl;
    int display_size = min(max_display, min(cube.getSizeX(), cube.getSizeY()));
    
    for (int y = 0; y < display_size; ++y) {
        for (int x = 0; x < display_size; ++x) {
            cout << fixed << setprecision(2) << setw(8) << cube(x, y, z_slice);
        }
        cout << endl;
    }
    cout << endl;
}

double calculateChecksum(const Cube3D& cube) {
    double sum = 0.0;
    for (int z = 0; z < cube.getSizeZ(); ++z) {
        for (int y = 0; y < cube.getSizeY(); ++y) {
            for (int x = 0; x < cube.getSizeX(); ++x) {
                sum += cube(x, y, z);
            }
        }
    }
    return sum;
}

int main() {
    cout << "3D Stencil Operations - Sequential Implementation" << endl;
    cout << "Based on: 'Optimizing Three-Dimensional Stencil-Operations on Heterogeneous Computing Environments'" << endl;
    cout << string(80, '=') << endl;
    
    // Test configurations
    const int CUBE_SIZE = 128;
    const int STENCIL_RADIUS = 2;
    const int LBM_SIZE = 64;
    const int LBM_STEPS = 100;
    
    // Test 1: Mean Blur Stencil Operation
    cout << "\nTest 1: Mean Blur Filter" << endl;
    cout << "Cube size: " << CUBE_SIZE << "³" << endl;
    cout << "Stencil radius: " << STENCIL_RADIUS << endl;
    
    Cube3D source_blur(CUBE_SIZE, CUBE_SIZE, CUBE_SIZE);
    Cube3D result_blur(CUBE_SIZE, CUBE_SIZE, CUBE_SIZE);
    
    initializeCube(source_blur, 1); // Linear gradient pattern
    
    cout << "\nSource data (first slice):" << endl;
    printSlice(source_blur, 0);
    
    auto start_time = high_resolution_clock::now();
    meanBlurStencil(source_blur, result_blur, STENCIL_RADIUS);
    auto end_time = high_resolution_clock::now();
    
    auto duration = duration_cast<milliseconds>(end_time - start_time);
    cout << "Mean blur completed in: " << duration.count() << " ms" << endl;
    
    cout << "\nResult after mean blur (first slice):" << endl;
    printSlice(result_blur, 0);
    
    cout << "Checksum - Source: " << fixed << setprecision(2) 
         << calculateChecksum(source_blur) << ", Result: " 
         << calculateChecksum(result_blur) << endl;
    
    // Test 2: Lattice Boltzmann Method
    cout << "\n" << string(80, '-') << endl;
    cout << "Test 2: Lattice Boltzmann Method (D3Q19)" << endl;
    cout << "Cube size: " << LBM_SIZE << "³" << endl;
    cout << "Time steps: " << LBM_STEPS << endl;
    
    vector<Cube3D> lbm_distributions;
    for (int i = 0; i < 19; ++i) {
        lbm_distributions.emplace_back(LBM_SIZE, LBM_SIZE, LBM_SIZE);
    }
    
    // Initialize with equilibrium distribution
    for (int i = 0; i < 19; ++i) {
        initializeCube(lbm_distributions[i], 0); // Random initialization
        // Scale down to reasonable values for LBM
        for (int z = 0; z < LBM_SIZE; ++z) {
            for (int y = 0; y < LBM_SIZE; ++y) {
                for (int x = 0; x < LBM_SIZE; ++x) {
                    lbm_distributions[i](x, y, z) *= 0.01f;
                }
            }
        }
    }
    
    LBMSolver solver(1.0f, 0.01f); // tau=1.0, dt=0.01
    
    start_time = high_resolution_clock::now();
    solver.simulate(lbm_distributions, LBM_STEPS);
    end_time = high_resolution_clock::now();
    
    duration = duration_cast<milliseconds>(end_time - start_time);
    cout << "LBM simulation completed in: " << duration.count() << " ms" << endl;
    
    // Display results
    cout << "\nFinal density distribution (component 0, first slice):" << endl;
    printSlice(lbm_distributions[0], 0);
    
    double total_mass = 0.0;
    for (const auto& dist : lbm_distributions) {
        total_mass += calculateChecksum(dist);
    }
    cout << "Total mass in system: " << fixed << setprecision(6) << total_mass << endl;
    
    // Performance summary
    cout << "\n" << string(80, '=') << endl;
    cout << "Performance Summary:" << endl;
    cout << "Mean Blur  - " << CUBE_SIZE << "³ elements, " << duration.count() << " ms" << endl;
    cout << "LBM        - " << LBM_SIZE << "³ x 19 distributions, " << LBM_STEPS << " steps, " 
         << duration.count() << " ms" << endl;
    
    return 0;
}