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

    // --- PHASE 3: 3D GRID INITIALIZATION ---
    int particlesPerSide = std::ceil(std::cbrt(numAtoms));
    float spacing = L / particlesPerSide;

    for (int i = 0; i < numAtoms; ++i) {
        Atom a;
        
        int xIndex = i % particlesPerSide;
        int yIndex = (i / particlesPerSide) % particlesPerSide;
        int zIndex = i / (particlesPerSide * particlesPerSide);

        glm::vec3 atomCenter(
            xIndex * spacing + spacing / 2.0f,
            yIndex * spacing + spacing / 2.0f,
            zIndex * spacing + spacing / 2.0f
        );

        // --- NEW: STRATIFIED ISOTOPE SETUP ---
        // If the atom is spawning in the bottom half of the box, make it a Boson.
        // If it's in the top half, make it a Fermion.
        if (atomCenter.y < L / 2.0f) {
            a.mass = 4.0f;     // Heavy Helium-4
            a.isBoson = true;  // Will form the superfluid web
        } else {
            a.mass = 3.0f;     // Lighter Helium-3
            a.isBoson = false; // Will remain isolated
        }

        glm::vec3 atomVelocity(0.0f); 

        for (int k = 0; k < numBeads; ++k) {
            Bead b;
            float u = (float)rand() / RAND_MAX;
            float v = (float)rand() / RAND_MAX;
            float theta = u * 2.0f * 3.14159f;
            float phi = std::acos(2.0f * v - 1.0f);
            
            float radius = 0.5f; 
            
            float bx = radius * std::sin(phi) * std::cos(theta);
            float by = radius * std::sin(phi) * std::sin(theta);
            float bz = radius * std::cos(phi);
            
            b.pos = atomCenter + glm::vec3(bx, by, bz);
            b.vel = atomVelocity;
            b.force = glm::vec3(0.0f);
            a.beads.push_back(b);
        }
        atoms.push_back(a);
    }
    // Make sure you delete the old "for (int i = 0; i < numAtoms; ++i)" 
    // isotope loop that used to be right below this!
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
    // 1. Recalculate omegaP so the quantum springs weaken dynamically as the system cools
    float hbar = 1.0f; 
    float kb = 1.0f;
    
    // Safety check to avoid division by zero near absolute zero
    if (temperature > 1e-6f) {
        float beta = 1.0f / (kb * temperature);
        omegaP = std::sqrt((float)numBeads) / (beta * hbar);
    } else {
        omegaP = 0.0f; // Springs die at complete delocalization
    }

    // 2. Velocity Verlet Step 1: Update positions and first half of velocities
    for (auto& atom : atoms) {
        for (auto& b : atom.beads) {
            b.pos = b.pos + b.vel * dt + 0.5f * (b.force / atom.mass) * dt * dt;
            b.vel += 0.5f * (b.force / atom.mass) * dt;
        }
    }
    
    // 3. Compute new forces based on updated positions
    computeForces();

    // 4. Velocity Verlet Step 2: Update second half of velocities (RESTORED)
    for (auto& atom : atoms) {
        for (auto& b : atom.beads) {
            b.vel += 0.5f * (b.force / atom.mass) * dt;
        }
    }

    // 5. LOCALIZED THERMOSTAT (Equipartition Enforcer)
    // Target kinetic energy for a single atom (all its beads)
    float targetAtomKE = 1.5f * numBeads * temperature; 

    for (auto& atom : atoms) {
        float atomKE = 0.0f;
        
        // Calculate the kinetic energy for THIS specific atom
        for (auto& b : atom.beads) {
            atomKE += 0.5f * atom.mass * glm::dot(b.vel, b.vel);
        }

        // Scale only this atom's velocities to match its mass-dependent thermal target
        if (atomKE > 0.0001f) {
            float scale = std::sqrt(targetAtomKE / atomKE);
            
            // Gentle coupling for stable quantum phase mixing
            float coupling = 0.1f; 
            scale = 1.0f + coupling * (scale - 1.0f); 

            for (auto& b : atom.beads) {
                b.vel *= scale;
            }
        }
    }
    
    // 6. Attempt multiple swaps per frame to allow the system to explore topologies
    int swapAttemptsPerFrame = numAtoms * 2; 
    for(int i = 0; i < swapAttemptsPerFrame; ++i) {
        attemptWormSwap();
    }
}

void PIMDSystem::attemptWormSwap() {
    if (numAtoms < 2) return;
        
    int a1 = rand() % numAtoms;
    int a2 = rand() % numAtoms;
    
    // Ensure we picked two different atoms
    if (a1 == a2) return; 

    // --- NEW: FERMIONIC EXCLUSION ---
    // Only Bosons can undergo quantum state exchange.
    if (!atoms[a1].isBoson || !atoms[a2].isBoson) {
        return; 
    }
    
    // Ensure identical isotopes (they must have the exact same mass)
    if (std::abs(atoms[a1].mass - atoms[a2].mass) > 0.01f) {
        return; 
    }

    swapAttempts++;

    // 2. Pick a random slice point (The Cut)
    int k = rand() % numBeads;
    int k_next = (k + 1) % numBeads;

    // 3. Calculate the change in Action (Spring Energy)
    float springK = atoms[a1].mass * omegaP * omegaP; 

    glm::vec3 rA_k = atoms[a1].beads[k].pos;
    glm::vec3 rA_next = atoms[a1].beads[k_next].pos;
    
    glm::vec3 rB_k = atoms[a2].beads[k].pos;
    glm::vec3 rB_next = atoms[a2].beads[k_next].pos;

    glm::vec3 diffA = rA_k - rA_next;
    glm::vec3 diffB = rB_k - rB_next;
    float distCurrent = glm::dot(diffA, diffA) + glm::dot(diffB, diffB);

    glm::vec3 diffCross1 = rA_k - rB_next;
    glm::vec3 diffCross2 = rB_k - rA_next;
    float distProposed = glm::dot(diffCross1, diffCross1) + glm::dot(diffCross2, diffCross2);

    float deltaAction = 0.5f * springK * (distProposed - distCurrent);

    // 4. The Metropolis-Hastings Acceptance Criterion
    bool accept = false;
    if (deltaAction < 0.0f) {
        accept = true; 
    } else {
        float prob = std::exp(-deltaAction);
        float randVal = (float)rand() / RAND_MAX;
        if (randVal < prob) accept = true; 
    }

    // 5. Execute the Topological Swap!
    if (accept) {
        swapAcceptances++;
        for (int j = k_next; j < numBeads; ++j) {
            Bead temp = atoms[a1].beads[j];
            atoms[a1].beads[j] = atoms[a2].beads[j];
            atoms[a2].beads[j] = temp;
        }
    }
}

float PIMDSystem::getAverageRadiusOfGyration(float targetMass) {
    float totalRg = 0.0f;
    int count = 0;

    for (const auto& atom : atoms) {
        if (std::abs(atom.mass - targetMass) > 0.1f) continue; 

        // 1. Calculate centroid using a reference bead to handle periodic boundaries
        glm::vec3 ref = atom.beads[0].pos;
        glm::vec3 offset(0.0f);
        
        for (const auto& b : atom.beads) {
            glm::vec3 diff = b.pos - ref;
            // Apply minimum image convention relative to reference bead
            if (diff.x > L * 0.5f) diff.x -= L;
            if (diff.x < -L * 0.5f) diff.x += L;
            if (diff.y > L * 0.5f) diff.y -= L;
            if (diff.y < -L * 0.5f) diff.y += L;
            if (diff.z > L * 0.5f) diff.z -= L;
            if (diff.z < -L * 0.5f) diff.z += L;
            
            offset += diff;
        }
        glm::vec3 centroid = ref + (offset / (float)numBeads);

        // 2. Calculate squared distance sum using the same periodic logic
        float rgSq = 0.0f;
        for (const auto& b : atom.beads) {
            glm::vec3 diff = b.pos - centroid;
            // Apply minimum image convention again for distance
            if (diff.x > L * 0.5f) diff.x -= L;
            if (diff.x < -L * 0.5f) diff.x += L;
            if (diff.y > L * 0.5f) diff.y -= L;
            if (diff.y < -L * 0.5f) diff.y += L;
            if (diff.z > L * 0.5f) diff.z -= L;
            if (diff.z < -L * 0.5f) diff.z += L;
            
            rgSq += glm::dot(diff, diff);
        }
        
        totalRg += std::sqrt(rgSq / (float)numBeads);
        count++;
    }

    return (count > 0) ? (totalRg / (float)count) : 0.0f;
}