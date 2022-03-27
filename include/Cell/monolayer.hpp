#include <vector>
#include <Cell.hpp>
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
            double Ftol;
            double dt0;
            std::vector<DPM::Cell> Cells;
            std::vector<bool> overlaps;
            //functions
            void disperse();
            void VertexFIRE();
            void UpdateEuler(double dt);
            double GetPackingFraction();
            double GetVertexKineticEnergy();
            void FindOverlaps(int ci, int cj);
            void RepulsiveForceUpdate();
            void AttractiveForceUpdate();
            void ResetForces();
    };
}