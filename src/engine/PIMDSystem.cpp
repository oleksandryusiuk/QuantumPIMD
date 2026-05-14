#include "engine/PIMDSystem.h"
#include <cmath>
#include <random>
#include <algorithm>

PIMDSystem::PIMDSystem(int N, int P, float boxSize, float temp) 
    : numAtoms(N), numBeads(P), L(boxSize), temperature(temp) {
    
    float hbar = 1.0f; 
    float kb = 1.0f;
    float beta = 1.0f / (kb * temperature);
    omegaP = std::sqrt((float)numBeads) / (beta * hbar);


    // 1. FCC Lattice Parameters
    // For Lennard-Jones, the ideal FCC lattice constant 'a' is approx 1.587
    // This perfectly spaces the nearest neighbors at the 1.122 minimum!
    float a = 1.587f; 
    
    // We place 4 atoms per unit cell
    int cellsPerSide = std::ceil(std::cbrt(numAtoms / 4.0f)); 
    
    float crystalSize = cellsPerSide * a;
    glm::vec3 offset = glm::vec3((L - crystalSize) / 2.0f);

    // The 4 local positions within an FCC unit cell
    glm::vec3 basis[4] = {
        glm::vec3(0.0f, 0.0f, 0.0f),
        glm::vec3(0.5f, 0.5f, 0.0f),
        glm::vec3(0.5f, 0.0f, 0.5f),
        glm::vec3(0.0f, 0.5f, 0.5f)
    };

    int placedAtoms = 0;

    for (int x = 0; x < cellsPerSide && placedAtoms < numAtoms; ++x) {
        for (int y = 0; y < cellsPerSide && placedAtoms < numAtoms; ++y) {
            for (int z = 0; z < cellsPerSide && placedAtoms < numAtoms; ++z) {
                
                // Loop through the 4 basis atoms in this specific cell
                for (int b = 0; b < 4 && placedAtoms < numAtoms; ++b) {
                    Atom atom;
                    
                    glm::vec3 cellCorner = offset + glm::vec3(x * a, y * a, z * a);
                    glm::vec3 center = cellCorner + (basis[b] * a);
                    
                    glm::vec3 atomVelocity(
                        ((float)rand() / RAND_MAX) - 0.5f,
                        ((float)rand() / RAND_MAX) - 0.5f,
                        ((float)rand() / RAND_MAX) - 0.5f
                    );
                    atomVelocity *= 0.01f; 

                    for (int k = 0; k < numBeads; ++k) {
                        Bead bead;
                        float u = (float)rand() / RAND_MAX;
                        float v = (float)rand() / RAND_MAX;
                        float theta = u * 2.0f * 3.14159f;
                        float phi = std::acos(2.0f * v - 1.0f);
                        float radius = 0.1f; 
                        
                        float bx = radius * std::sin(phi) * std::cos(theta);
                        float by = radius * std::sin(phi) * std::sin(theta);
                        float bz = radius * std::cos(phi);
                        
                        bead.pos = center + glm::vec3(bx, by, bz);
                        bead.vel = atomVelocity;
                        bead.force = glm::vec3(0.0f);
                        atom.beads.push_back(bead);
                    }
                    atoms.push_back(atom);
                    placedAtoms++;
                }
            }
        }
    }
}

void PIMDSystem::computeForces() {
    for (auto& atom : atoms) {
        for (auto& bead : atom.beads) bead.force = glm::vec3(0.0f);
    }

    float wallStiffness = 500.0f; 
    float wallDamping = 20.0f; 

    for (int i = 0; i < numAtoms; ++i) {
        float springK = atoms[i].mass * omegaP * omegaP;

        for (int k = 0; k < numBeads; ++k) {
            int prev = (k - 1 + numBeads) % numBeads;
            int next = (k + 1) % numBeads;
            
            glm::vec3 r_prev = atoms[i].beads[k].pos - atoms[i].beads[prev].pos;
            glm::vec3 r_next = atoms[i].beads[k].pos - atoms[i].beads[next].pos;
            atoms[i].beads[k].force -= springK * (r_prev + r_next);

            glm::vec3& pos = atoms[i].beads[k].pos;
            glm::vec3& vel = atoms[i].beads[k].vel; 

            if (pos.x < 0.0f) atoms[i].beads[k].force.x += wallStiffness * (0.0f - pos.x) - wallDamping * vel.x;
            else if (pos.x > L) atoms[i].beads[k].force.x -= wallStiffness * (pos.x - L) - wallDamping * vel.x;
            
            if (pos.y < 0.0f) atoms[i].beads[k].force.y += wallStiffness * (0.0f - pos.y) - wallDamping * vel.y;
            else if (pos.y > L) atoms[i].beads[k].force.y -= wallStiffness * (pos.y - L) - wallDamping * vel.y;
            
            if (pos.z < 0.0f) atoms[i].beads[k].force.z += wallStiffness * (0.0f - pos.z) - wallDamping * vel.z;
            else if (pos.z > L) atoms[i].beads[k].force.z -= wallStiffness * (pos.z - L) - wallDamping * vel.z;

            for (int j = i + 1; j < numAtoms; ++j) {
                glm::vec3 diff = atoms[i].beads[k].pos - atoms[j].beads[k].pos;
                float r = glm::length(diff);
                
                if (r < 0.8f) r = 0.8f; 
                
                float r6 = std::pow(r, 6);
                float forceMag = (48.0f / (r6 * r6 * r)) - (24.0f / (r6 * r)); 
                
                if (forceMag > 200.0f) forceMag = 200.0f;
                if (forceMag < -200.0f) forceMag = -200.0f;

                glm::vec3 f_dir = glm::normalize(diff);
                
                atoms[i].beads[k].force += f_dir * forceMag * (1.0f / numBeads);
                atoms[j].beads[k].force -= f_dir * forceMag * (1.0f / numBeads);
            }
        }
    }
}

void PIMDSystem::update(float dt) {
    for (auto& atom : atoms) {
        for (auto& b : atom.beads) {
            b.pos = b.pos + b.vel * dt + 0.5f * (b.force / atom.mass) * dt * dt;
            b.vel += 0.5f * (b.force / atom.mass) * dt;
        }
    }

    computeForces();

    float currentKineticEnergy = 0.0f;
    for (auto& atom : atoms) {
        for (auto& b : atom.beads) {
            b.vel += 0.5f * (b.force / atom.mass) * dt;
            currentKineticEnergy += 0.5f * atom.mass * glm::dot(b.vel, b.vel);
        }
    }

    float targetKineticEnergy = 1.5f * (numAtoms * numBeads) * temperature; 
    
    if (currentKineticEnergy > 0.0001f) {
        float scale = std::sqrt(targetKineticEnergy / currentKineticEnergy);
        // 3. More aggressive thermostat to ensure they keep moving in the huge box
        float coupling = 0.2f; 
        scale = 1.0f + coupling * (scale - 1.0f); 

        for (auto& atom : atoms) {
            for (auto& b : atom.beads) {
                b.vel *= scale;
            }
        }
    }
}