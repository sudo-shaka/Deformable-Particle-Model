#include <cmath>
#include <Cell.hpp>
#include <vector>

namespace DPM{
    Cell::Cell(int idx1, double x1, double y1, double CalA01, int NV1, double Kl1, double Kb1, double Ka1, double v01, double Dr1, double Ds1, double a01, double psi1){
        idx = idx1;
        NV = NV1;
        calA0 = CalA01*(NV*tan(M_PI/NV)/M_PI);

        a0 = a01;
        r0 = sqrt((2.0*a0)/(NV*sin((2.0*M_PI)/NV)));
        l0 = 2.0*sqrt(M_PI*calA0*a0)/NV;
    
        Kl = Kl1;
        Kb = Kb1;
        Ka = Ka1;
        psi = psi1;
        Dr = Dr1;
        Ds = Ds1;
        v0 = v01;
        vmin = 1e-2*v0;

        X.resize(NV); Y.resize(NV); l1.resize(NV); l2.resize(NV);
        Fx.resize(NV); Fy.resize(NV); vx.resize(NV); vy.resize(NV);
        im1.resize(NV); ip1.resize(NV);

        for(int i=0; i<NV; i++){
            X[i] = r0*(cos(2.0*M_PI*(i+1)/NV)) + x1;
            Y[i] = r0*(sin(2.0*M_PI*(i+1)/NV)) + y1;
            l1[i] = 0.0;
            l2[i] = 0.0;
            vx[i] = 0.0;
            vy[i] = 0.0;
            im1[i] = i-1;
            ip1[i] = i+1;
        }
        im1[0] = NV-1; 
        ip1[NV-1] = 0;  
    }

    void Cell::PerimeterForceUpdate(){
        double lvx[NV], lvy[NV], ulvx[NV], ulvy[NV], length[NV], dli[NV], dlim1[NV];
        int i=0;
        for(i = 0;i<NV;i++){
            lvx[i] = X[ip1[i]] - X[i];
            lvy[i] = Y[ip1[i]] - Y[i];
        }
        for(i=0;i<NV;i++){
            length[i] = sqrt(lvx[i]*lvx[i] + lvy[i] * lvy[i]);   
        }
        for(i=0;i<NV;i++){
            ulvx[i] = lvx[i]/length[i];
            ulvy[i] = lvy[i]/length[i];
            dli[i] = (length[i]/l0) - 1;
            dlim1[i] = (length[im1[i]]/l0) - 1;
        }
        for(i=0;i<NV;i++){
            Fx[i] += (Kl*((sqrt(a0))/l0))*(dli[i]*ulvx[i] - dlim1[i]*ulvx[im1[i]]);
            Fy[i] += (Kl*((sqrt(a0))/l0))*(dli[i]*ulvy[i] - dlim1[i]*ulvy[im1[i]]);
        }
    }

    void Cell::AreaForceUpdate(){
        double area = GetArea();
        double areaStrain = (area/a0) - 1.0;
        for(int i=0; i<NV; i++){
            Fx[i] += (Ka/(sqrt(a0)))*0.5*areaStrain*(Y[im1[i]]-Y[ip1[i]]);
            Fy[i] += (Ka/(sqrt(a0)))*0.5*areaStrain*(X[im1[i]]-X[ip1[i]]);
        }
    }

    void Cell::BendingForceUpdate(){
        double rho0 = sqrt(a0);
        double fb = Kb*(rho0/(l0*l0));
        double lvx[NV], lvy[NV], six[NV], siy[NV];
        int i=0;
        for(i=0;i<NV;i++){
            lvx[i] = X[ip1[i]] - X[i];
            lvy[i] = Y[ip1[i]] - Y[i];
        }

        for(i=0;i<NV;i++){
            six[i] = lvx[i] - lvx[im1[i]];
            siy[i] = lvy[i] - lvy[im1[i]];
        }

        for(i=0;i<NV;i++){
            Fx[i] += 0.1*fb*(2.0*six[i] - six[im1[i]] - six[ip1[i]]);
            Fy[i] += 0.1*fb*(2.0*siy[i] - siy[im1[i]] - siy[ip1[i]]);
        }
    }

    void Cell::DrivingForceUpdate(double dt){
        if(v0 == 0.0){
            return;
        }
        double cx = GetCenterX();
        double cy = GetCenterY();
        double rx,ry, psiVi, v0tmp, rscale, dpsi;
        int i;
        for(i=0;i<NV;i++){
            rx = X[i] -cx;
            ry = Y[i] -cy;
            psiVi = atan2(rx,ry);
            dpsi = psiVi - psi;
            dpsi -= 2.0*M_PI*round(dpsi/(2.0*M_PI));
            v0tmp = v0*exp(-(dpsi*dpsi)/(2.0*Ds*Ds))+vmin;
            rscale = sqrt(rx*rx + ry*ry);
            Fx[i] += v0tmp*(rx/rscale);
            Fy[i] += v0tmp*(ry/rscale);
        }
    }

    void Cell::UpdateShapeForces(){
        for(int i=0;i<NV;i++){
            Fx[i] = 0.0;
            Fy[i] = 0.0;
        }
        AreaForceUpdate();
        PerimeterForceUpdate();
        BendingForceUpdate();
    }

    void Cell::UpdateEuler(double dt){
        UpdateShapeForces();
        DrivingForceUpdate(dt);
        for(int i=0; i<NV; i++){
            X[i] += Fx[i]*dt;
            Y[i] += Fy[i]*dt;
        }
    }

    void Cell::UpdateDirectorDiffusion(double dt){
        double r1, r2, grv;
        r1 = drand48();
        r2 = drand48();
        grv = sqrt(-2.0*log(r1))*cos(2.0*M_PI*r2);
        psi += sqrt(2.0*dt*Dr)*grv;
    }

    void Cell::UpdateStickness(int vertex, double newl1, double newl2){
        l1[vertex] = newl1;
        l2[vertex] = newl2;
    }

    double Cell::GetPerim(){
        double perimeter=0.0;
        double dx, dy; int i;
        for(i=0;i<NV;i++){
            dx = abs(X[im1[i]]-X[i]);
            dy = abs(Y[im1[i]]-Y[i]);
            perimeter += sqrt(dx*dx + dy*dy);
        }
        return perimeter;
    }

    double Cell::GetArea(){
	    double Area = 0.0;
        int j = NV-1;
	    for (int i = 0; i < NV; i++) {
	    	Area += 0.5 * ((X[j] + X[i]) * (Y[j] - Y[i]));
            j = i;
        }
        return abs(Area);
    }

    double Cell::GetShapeParam(){
        double perim = Cell::GetPerim();
        double area = Cell::GetArea();
        return (perim*perim)/(area*M_PI*4);
    }

    double Cell::GetCenterX(){
        double CenterX = 0.0;
        for(int i=0; i<NV; i++){
            CenterX += X[i];
        }
        return CenterX/NV;

    }

    double Cell::GetCenterY(){
        double CenterY = 0.0;
        for(int i=0; i<NV; i++){
            CenterY += Y[i];
        }
        return CenterY/NV;

    }
    bool Cell::isConvex(){
        bool neg = false;
        bool pos = false;
        for(int i=0;i<NV;i++){
            int a = i;
            int b = (i + 1) % NV;
            int c = (i + 2) % NV;
            double BAx = X[a] - X[b];
            double BAy = Y[a] - Y[b];
            double BCx = X[c] - X[b];
            double BCy = Y[c] - Y[b];
            double crossProduct = (BAx * BCy - BAy * BCx); 
            if(crossProduct < 0) neg = true;
            else if(crossProduct > 0) pos = true;
            if(neg && pos) return false;
        }
        return true;
    }
}
