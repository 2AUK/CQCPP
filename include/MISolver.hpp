#pragma once
#include "Molecule.hpp"

class MMDSolver{
    /* McMurchie-Davidson Method for evaluating molecular integrals
     * As is in Ch. 9 of Molecular Electronic-Structure Theory by Trygve Helgaker et al.
     */
    public:
        Molecule input_molecule; //Molecule object as input
        Eigen::ArrayXXf S(); //Overlap Matrix
        Eigen::ArrayXXf T(); //Kinetic Matrix
        Eigen::ArrayXXf V(); //Nuclear Attraction Matrix
        Eigen::ArrayXXf ERI(); //Electron Repulsion Interaction Matrix
    private:
        float E(); //Recursive definitions for expansion coeffiecients
        float p_overlap();
        float c_overlap();
        float p_kinetic();
        float c_kinetic();
        float R();
        float p_nuclear();
        float c_nuclear();
        float p_electron();
        float c_electron();
};

class OSSolver{
    /* Obara-Saika Method for evaluating molecular integrals
     */

};

class HGPSolver{
    /* Head-Gordon-Pople Method for evaluating molecular integrals
     */

}
