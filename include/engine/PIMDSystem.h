#ifndef PIMD_SYSTEM_H
#define PIMD_SYSTEM_H

#include <vector>
#include <glm/glm.hpp>

struct Bead {
    glm::vec3 pos;
    glm::vec3 vel;
    glm::vec3 force;
};

struct Atom {
    std::vector<Bead> beads;
    float mass = 1.0f;     
    float charge = 0.0f;
    bool isBoson; 
};

class PIMDSystem {
public:
    int numAtoms;
    int numBeads;
    float L;             // Box size
    float temperature;
    float omegaP;        

    std::vector<Atom> atoms;

    PIMDSystem(int N, int P, float boxSize, float temp);
    void update(float dt);
    // --- Phase 2: Worm Algorithm ---
    void attemptWormSwap();
    int swapAttempts = 0;
    int swapAcceptances = 0;

    int getSwapAttempts() const { return swapAttempts; }
    int getSwapAcceptances() const { return swapAcceptances; }
    
    float getSwapRate() const { 
        if (swapAttempts == 0) return 0.0f;
        return (float)swapAcceptances / swapAttempts; 
    }
    float getAverageRadiusOfGyration(float targetMass);
    
private:
    void computeForces();
    // Removed applyPBC and getMinimumImage!
};

#endif