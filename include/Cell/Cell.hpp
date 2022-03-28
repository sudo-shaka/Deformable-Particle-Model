#include <vector>
namespace DPM {
    class Cell
    {
        public:
            int idx;
            int NV;
            double calA0;
            double l0;
            double r0;
            std::vector<double> X;
            std::vector<double> Y;
            std::vector<double> vx;
            std::vector<double> vy;
            std::vector<double> Fx;
            std::vector<double> Fy;
            double v0;
            double Ka;
            double Kb;
            double Kl;
            std::vector<int> ip1;
            std::vector<int> im1;
            double vmin;
            double Dr;
            double Ds;
            double a0;
            double psi;
            std::vector<double> l1;
            std::vector<double> l2;
            std::vector<double> radii;


            Cell(int idx,
                double x1,double y1,
                double calA0,
                int NV, 
                double Kl, double Kb, double Ka, 
                double v0, double Dr, double Ds, 
                double a0, double psi);

            Cell(int NV);
            
            void ResetForces();
            void FindRadii();
            void PerimeterForceUpdate();
            void AreaForceUpdate();
            void BendingForceUpdate();
            void DrivingForceUpdate(double dt);
            void UpdateShapeForces();
            void UpdateEuler(double dt);
            void UpdateVV();
            void UpdateDirectorDiffusion(double dt);
            void UpdateDirectorPsi(double newpsi);
            void UpdateStickness(int vertex, double l1, double l2);
            double GetPerim();
            double GetArea();
            double GetCenterX();
            double GetCenterY();
            double GetShapeParam();
            bool isConvex();
    };
}
