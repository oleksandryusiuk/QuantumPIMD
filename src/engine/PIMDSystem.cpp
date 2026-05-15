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

    // --- PHASE 2: 2-ATOM SANDBOX ---
    glm::vec3 center(L / 2.0f);
    float spacing = 1.122f; // Lennard-Jones equilibrium resting distance

    for (int i = 0; i < numAtoms; ++i) {
        Atom a;
        // Place Atom 0 slightly left, Atom 1 slightly right
        glm::vec3 atomCenter = center;
        if (i == 0) atomCenter.x -= spacing / 2.0f;
        if (i == 1) atomCenter.x += spacing / 2.0f;

        glm::vec3 atomVelocity(0.0f); // Start perfectly calm

        for (int k = 0; k < numBeads; ++k) {
            Bead b;
            float u = (float)rand() / RAND_MAX;
            float v = (float)rand() / RAND_MAX;
            float theta = u * 2.0f * 3.14159f;
            float phi = std::acos(2.0f * v - 1.0f);
            
            // PRE-INFLATE the probability clouds so they instantly overlap
            float radius = 1.0f; 
            
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
        // Gentle thermostat for stable quantum phase mixing
        float coupling = 0.05f;
        scale = 1.0f + coupling * (scale - 1.0f); 

        for (auto& atom : atoms) {
            for (auto& b : atom.beads) {
                b.vel *= scale;
            }
        }
    }
}

void PIMDSystem::attemptWormSwap() {
    if (numAtoms < 2) return;
    swapAttempts++;

    // 1. Pick two random atoms
    int a1 = rand() % numAtoms;
    int a2 = rand() % numAtoms;
    while (a1 == a2) a2 = rand() % numAtoms; // Ensure they are different

    // 2. Pick a random slice point (The Cut)
    int k = rand() % numBeads;
    int k_next = (k + 1) % numBeads;

    // 3. Calculate the change in Action (Spring Energy)
    float springK = atoms[a1].mass * omegaP * omegaP; 

    glm::vec3 rA_k = atoms[a1].beads[k].pos;
    glm::vec3 rA_next = atoms[a1].beads[k_next].pos;
    
    glm::vec3 rB_k = atoms[a2].beads[k].pos;
    glm::vec3 rB_next = atoms[a2].beads[k_next].pos;

    // Current squared distances (Un-swapped)
    glm::vec3 diffA = rA_k - rA_next;
    glm::vec3 diffB = rB_k - rB_next;
    float distCurrent = glm::dot(diffA, diffA) + glm::dot(diffB, diffB);

    // Proposed squared distances (Cross-linked!)
    glm::vec3 diffCross1 = rA_k - rB_next;
    glm::vec3 diffCross2 = rB_k - rA_next;
    float distProposed = glm::dot(diffCross1, diffCross1) + glm::dot(diffCross2, diffCross2);

    float deltaAction = 0.5f * springK * (distProposed - distCurrent);

    // 4. The Metropolis-Hastings Acceptance Criterion
    bool accept = false;
    if (deltaAction < 0.0f) {
        accept = true; // Always accept if it relaxes the springs
    } else {
        float prob = std::exp(-deltaAction);
        float randVal = (float)rand() / RAND_MAX;
        if (randVal < prob) accept = true; // Accept based on quantum probability
    }

    // 5. Execute the Topological Swap!
    if (accept) {
        swapAcceptances++;
        // Physically trade all beads from the cut point to the end of the polymer
        for (int j = k_next; j < numBeads; ++j) {
            Bead temp = atoms[a1].beads[j];
            atoms[a1].beads[j] = atoms[a2].beads[j];
            atoms[a2].beads[j] = temp;
        }
    }
}