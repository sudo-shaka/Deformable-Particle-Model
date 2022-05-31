#include <vector>
#include <Cell.hpp>
#include <functional>
#include <string>
namespace DPM{
    class monolayer{
        public:
            //constructor
            monolayer(std::vector<DPM::Cell> Cells, double phi0);
            //vars
            double phi0;
            int NCELLS;
            int NSMALL;
            int VertDOF;
            double L;
            double Kc;
            double U;
            std::vector<DPM::Cell> Cells;
            std::vector<bool> overlaps;
            int cellidx = 0;
            std::string InteractingMethod;
            int InteractingMethodNum;
             //functions
            void disperse();
            void VertexFIRE(double alpha, double dt, int maxits, double Ftol);
            void UpdateEuler(int nsteps, double dt);
            void UpdateEuler(double dt);
            void UpdateVV(int nsteps, double dt);
            double GetPackingFraction();
            double GetVertexKineticEnergy();
            std::vector<bool> FindOverlaps(int ci, int cj);
            void RepulsiveForces(int ci);
            void AttractiveForces(int ci);
            void RetractingForceUpdate(int ci);
            void PinnedForces(int ci);
            void InteractingForceUpdate();
            void ShapeForceUpdate();
            void ResetForces();
            void DrivingForceUpdate();
            void UpdateDirectorDiffusion(double dt);
            void UpdateDirectorVicsek(double eta, double InteractingRadius);
            void CellDivision(int cellidx);
            void NearestNeighborUpdate();
            void FindNearestNeighbor(int ci);
            void SetInteractingMethod(std::string methodstring);
    };
}
