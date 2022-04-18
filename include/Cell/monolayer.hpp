#include <vector>
#include <Cell.hpp>
#include <functional>
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
             //functions
            void disperse();
            void VertexFIRE(std::function<void()> InteractingMethod, double alpha, double dt, int maxits, double Ftol);
            void UpdateEuler(std::function<void()> InteractingMethod,int nsteps, double dt);
            void UpdateEuler(double dt);
            void UpdateVV(std::function<void()> InteractingMethod, int nsteps, double dt);
            double GetPackingFraction();
            double GetVertexKineticEnergy();
            void FindOverlaps(int ci, int cj);
            void RepulsiveForces();
            void AttactiveForces();
            void RetractingForceUpdate();
            void PinnedForces();
            void MixedInteractingMethods(std::vector<bool> firstMethod,std::function<void()> Method1, std::function<void()> Method2);
            void ShapeForceUpdate();
            void ResetForces();
            void UpdateDirectorDiffusion(double dt);
            void UpdateDirectorVicsek(double eta, double InteractingRadius);
            void CellDivision(int cellidx);
            void NearestNeighborUpdate();
            void FindNearestNeighbor(int ci);
    };
}
