#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <iostream>
#include "engine/PIMDSystem.h"
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <fstream>
#include <iostream>

int main() {
    if (!glfwInit()) return -1;
    
    GLFWwindow* window = glfwCreateWindow(800, 800, "Thermodynamics: Dynamic Phase Transition (Melting)", NULL, NULL);
    if (!window) { glfwTerminate(); return -1; }
    glfwMakeContextCurrent(window);
    
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) return -1;

    // --- Start in a Deep Freeze ---
    // Start at T = 0.05f to form the perfect crystal lattice
    // Mass = 20.0f (Good inertia), Temp = 0.05f (Cold, but NOT zero!)
    // 108 atoms for a perfect 3x3x3 FCC grid
    PIMDSystem system(4, 320, 8.0f, 0.05f);

    glEnable(GL_DEPTH_TEST);

    // --- CSV LOGGING SETUP ---
    std::ofstream csvFile("pimd_thermodynamics.csv");
    csvFile << "Time,Temperature,Rg_Heavy,Rg_Light,SwapRate\n";
    
    float timeElapsed = 0.0f;
    int frameCounter = 0;

    while (!glfwWindowShouldClose(window)) {
        float dt = 0.001f;

        system.update(dt);

        // --- DYNAMIC COOLING (Optional but highly recommended) ---
        // Slowly decrease the temperature by a tiny fraction every frame 
        // to sweep from classical gas down to quantum condensate.
        system.temperature *= 0.9995f; 

        // --- RECORD DATA (Once every 10 frames to avoid giant files) ---
        timeElapsed += dt;
        frameCounter++;

        if (frameCounter % 10 == 0) {
            float rgHeavy = system.getAverageRadiusOfGyration(1.0f);
            float rgLight = system.getAverageRadiusOfGyration(0.25f);
            float swapRate = system.getSwapRate();

            csvFile << timeElapsed << "," 
                    << system.temperature << "," 
                    << rgHeavy << "," 
                    << rgLight << "," 
                    << swapRate << "\n";
            
            // Print to console just so you know it's working
            std::cout << "T: " << system.temperature 
                      << " | Rg(H): " << rgHeavy 
                      << " | Rg(L): " << rgLight 
                      << " | Swaps: " << (swapRate * 100.0f) << "%" << std::endl;
        }

        // Calculate acceptance rate
        float swapRate = 0.0f;
        if (system.swapAttempts > 0) {
            swapRate = (float)system.swapAcceptances / system.swapAttempts * 100.0f;
        }

        // --- REAL-TIME TOPOLOGY UI ---
        std::string title = "PIMD Engine | Swaps: " + std::to_string(system.swapAcceptances) + 
                            " | Rate: " + std::to_string(swapRate).substr(0, 4) + "%"; 
        glfwSetWindowTitle(window, title.c_str());

        glClearColor(0.05f, 0.05f, 0.05f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        int width, height;
        glfwGetFramebufferSize(window, &width, &height);
        float ratio = width / (float)height;
        glViewport(0, 0, width, height);

        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        glm::mat4 proj = glm::perspective(glm::radians(45.0f), ratio, 0.1f, 100.0f);
        glLoadMatrixf(glm::value_ptr(proj));

        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();

        // Pushed camera in closer to see the two atoms interact
        float time = glfwGetTime();
        float camX = std::sin(time * 0.2f) * 10.0f; 
        float camZ = std::cos(time * 0.2f) * 10.0f;
        float center = system.L / 2.0f;
        glm::mat4 view = glm::lookAt(
            glm::vec3(camX + center, center + 2.0f, camZ + center), 
            glm::vec3(center, center, center),                      
            glm::vec3(0.0f, 1.0f, 0.0f)                             
        );
        glLoadMatrixf(glm::value_ptr(view));

        // Draw Boundary Box
        glColor3f(0.0f, 0.3f, 0.0f);
        glBegin(GL_LINES);
            float L = system.L;
            glVertex3f(0,0,0); glVertex3f(L,0,0); glVertex3f(L,0,0); glVertex3f(L,0,L);
            glVertex3f(L,0,L); glVertex3f(0,0,L); glVertex3f(0,0,L); glVertex3f(0,0,0);
            glVertex3f(0,L,0); glVertex3f(L,L,0); glVertex3f(L,L,0); glVertex3f(L,L,L);
            glVertex3f(L,L,L); glVertex3f(0,L,L); glVertex3f(0,L,L); glVertex3f(0,L,0);
            glVertex3f(0,0,0); glVertex3f(0,L,0); glVertex3f(L,0,0); glVertex3f(L,L,0);
            glVertex3f(L,0,L); glVertex3f(L,L,L); glVertex3f(0,0,L); glVertex3f(0,L,L);
        glEnd();

        // Draw Atoms (Modified to clearly show color swapping)
        for (int i = 0; i < system.numAtoms; ++i) {
            // Atom 0 is heavily Blue, Atom 1 is heavily Red
            float r = (i == 0) ? 0.2f : 1.0f;
            float g = 0.3f;
            float b = (i == 0) ? 1.0f : 0.2f;
            glColor3f(r, g, b);
            
            glBegin(GL_LINES); 
            for (int k = 0; k < system.numBeads; ++k) {
                glm::vec3 p1 = system.atoms[i].beads[k].pos;
                glm::vec3 p2 = system.atoms[i].beads[(k + 1) % system.numBeads].pos;
                glVertex3f(p1.x, p1.y, p1.z);
                glVertex3f(p2.x, p2.y, p2.z);
            }
            glEnd();

            glPointSize(4.0f); 
            glBegin(GL_POINTS);
            for (auto& bead : system.atoms[i].beads) glVertex3f(bead.pos.x, bead.pos.y, bead.pos.z);
            glEnd();
        }

        glfwSwapBuffers(window);
        glfwPollEvents();
    }
    csvFile.close();
    glfwTerminate();
    return 0;
}