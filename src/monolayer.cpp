#include <cmath>
#include <monolayer.hpp>
#include <vector>
#include <iostream>

using namespace std;

namespace DPM{
    monolayer::monolayer(std::vector<DPM::Cell> inputCells, double phi0input){
        Cells = inputCells;
        phi0 = phi0input;
        NCELLS = Cells.size();
        VertDOF = 0;
        double sumareas = 0.0;
        for(int i=0;i<NCELLS;i++){
            VertDOF += Cells[i].NV;
            sumareas += M_PI*(Cells[i].r0*Cells[i].r0);
            Cells[i].idx = i;
        }
        L = sqrt(sumareas)/phi0;
        Kc = 1.0;
        Ftol = 1e-2;
        dt0 = 0.005;
    }

    void monolayer::disperse(){
        std::vector<double> X,Y,Fx,Fy;
        X.resize(NCELLS);Y.resize(NCELLS);
        Fx.resize(NCELLS); Fy.resize(NCELLS);
        double ri,xi,yi,xj,yj,dx,dy,rj,dist;
        double ux,uy,ftmp,fx,fy;
        int i,j, count=0;
        for(i=0;i<NCELLS;i++){
            X[i] = drand48() * L;
            Y[i] = drand48() * L;
        }
        double oldU = 100, dU = 100, U;
        while(abs(dU) > 0.0001){
            U = 0;
            for(i=0;i<NCELLS;i++){
                Fx[i] = 0.0;
                Fy[i] = 0.0;
            }
            for(i=0;i<NCELLS;i++){
                xi = X[i];
                yi = Y[i];
                ri = Cells[i].r0;
                for(j=0;j<NCELLS;j++){
                    if(j != i){
                        xj = X[j];
                        yj = Y[j];
                        rj = Cells[j].r0;
                        dx = xj-xi;
                        dx -= L*round(dx/L);
                        dy = yj-yi;
                        dy -= L*round(dy/L);
                        dist = sqrt(dx*dx + dy*dy);
                        if(abs(dist) <= (ri+rj)){
                            ux = dx/dist;
                            uy = dy/dist;
                            ftmp = (1.0-dist/(ri+rj))/(ri+rj);
                            fx = ftmp*ux;
                            fy = ftmp*uy;
                            Fx[i] -= fx;
                            Fy[i] -= fy;
                            Fy[j] += fy;
                            Fx[j] += fx;
                            U += 0.5*(1-(dist/(ri+rj))*(1-dist/(ri+rj)));
                        }
                    }
                }
            }
            for(int i=0; i<NCELLS;i++){
                X[i] += 0.01*Fx[i];
                Y[i] += 0.01*Fy[i];
            }
            dU = U-oldU;
            oldU = U;
            count++;
            if(count > 1e5){
                break;
            }
        }
        for(int i=0;i<NCELLS;i++){
            for(j=0;j<Cells[i].NV;j++){
                Cells[i].X[j] = Cells[i].r0*(cos(2.0*M_PI*(j+1)/Cells[i].NV)) + X[i];
                Cells[i].Y[j] = Cells[i].r0*(sin(2.0*M_PI*(j+1)/Cells[i].NV)) + Y[i];
            }
        }
        Ftol = U/10;
    }
    void monolayer::VertexFIRE(){
        int ci,vi,i,NDELAY = 20, itmax = 5e4;

        double P, fnorm, fcheck, vnorm, alpha,alpha0, dtmax, dtmin;
        int npPos, npNeg, fireit;
        double dt = dt0;

        P=0;
        fnorm = 0;
        vnorm =0;
        alpha0 = 0.2;
        alpha = alpha0;

        dtmax = 10.0 * dt;
        dtmin = 1e-3 * dt;
        double finc = 1.01;
        double fdec = 0.99;
        double falpha = 0.99;
        int NNEGMAX = 1000;

        npPos = 0;
        npNeg = 0;

        fireit = 0;
        fcheck = 10*Ftol;

        ResetForces();

        while((fcheck > Ftol || fireit < NDELAY) && fireit < itmax){
            P = 0.0;

            for(i=0;i<NCELLS;i++){
                for(vi=0;vi<Cells[i].NV;i++){
                    P+= sqrt(Cells[i].vx[vi]*Cells[i].vy[vi]) + sqrt(Cells[i].Fx[vi]*Cells[i].Fy[vi]);
    
                }
            }

            if(P>0){
                npPos++;
                npNeg = 0;
                if (npPos > NDELAY){
                    if(dt * finc < dtmax)
                        dt *= finc;

                    alpha *= falpha;
                }
            }
            else{
                npPos = 0;
                npNeg++;

                if(npNeg > NNEGMAX){
                    cerr << "ERROR during FIRE, P<0 for too long..." << endl;
                    return;
                }

                for(ci=0;ci<NCELLS;ci++){
                    for(vi=0;vi<Cells[ci].NV;vi++){
                        Cells[ci].X[vi] -= 0.5 * dt * Cells[ci].vx[vi];
                        Cells[ci].Y[vi] -= 0.5 * dt * Cells[ci].vy[vi];
                        Cells[ci].vx[vi] = 0.0;
                        Cells[ci].vy[vi] = 0.0;
                    }
                }

                if(fireit > NDELAY){
                    if(dt*fdec > dtmin){
                        dt *= fdec;
                    }
                    alpha = alpha0;
                }
            }

            for(ci=0;ci<NCELLS;ci++){
                for(vi=0;vi<Cells[ci].NV;vi++){
                    Cells[ci].vx[vi] += 0.5 * dt * Cells[ci].Fx[vi];
                    Cells[ci].vy[vi] += 0.5 * dt * Cells[ci].Fy[vi];
                }
            }

            fnorm = 0.0;
            vnorm = 0.0;

            for(ci=0;ci<NCELLS;ci++){
                for(vi=0;vi<Cells[ci].NV;vi++){
                    fnorm += Cells[ci].Fx[vi] + Cells[ci].Fy[vi];
                    fnorm += Cells[ci].vx[vi] + Cells[ci].vy[vi];
                }
            }

            fnorm = sqrt(fnorm);
            vnorm = sqrt(vnorm);

            if(fnorm > 0){
                for(ci=0;ci<NCELLS;ci++){
                    for(vi=0;vi<Cells[ci].NV;vi++){
                        Cells[ci].vx[vi] = (1-alpha) * Cells[ci].vx[vi] + alpha * (Cells[ci].Fx[vi]/fnorm) *vnorm;
                        Cells[ci].vy[vi] = (1-alpha) * Cells[ci].vy[vi] + alpha * (Cells[ci].Fy[vi]/fnorm) *vnorm;
                    }
                }
            }

            for(ci=0;ci<NCELLS;ci++){
                for(vi=0;vi<Cells[ci].NV;vi++){
                    Cells[ci].X[vi] += dt*Cells[ci].vx[vi];
                    Cells[ci].Y[vi] += dt*Cells[ci].vy[vi];
                }
            }

            for(ci=0;ci<NCELLS;ci++){
                for(vi=0;vi<Cells[ci].NV;vi++){
                    Cells[ci].UpdateShapeForces();
                    RepulsiveForceUpdate();
                    AttractiveForceUpdate();
                }
            }

            for(ci=0;ci<NCELLS;ci++){
                for(vi=0;vi<Cells[ci].NV;vi++){  
                     Cells[ci].vx[vi] += 0.5 * dt * Cells[ci].Fx[vi];
                     Cells[ci].vy[vi] += 0.5 * dt * Cells[ci].Fy[vi];
                }
            }

            fcheck = 0.0;
            for(ci=0;ci<NCELLS;ci++){
                for(vi=0;vi<Cells[ci].NV;vi++){
                    fcheck += Cells[ci].Fx[vi] + Cells[ci].Fy[vi];
                }
            }
            fcheck = sqrt(fcheck)/VertDOF;

            fireit++;

        }
        if(fireit == itmax){
            cerr << "(!) Fire Minimization did not converge" << endl;
        }
        else{
            cout << "FIRE Minimization Converged!" << endl;
        }

    }
    void monolayer::UpdateEuler(double dt){
        ResetForces();
        for(int ci=0;ci<NCELLS;ci++){
            Cells[ci].AreaForceUpdate();
            Cells[ci].PerimeterForceUpdate();
            Cells[ci].BendingForceUpdate();
            Cells[ci].DrivingForceUpdate(dt);
        }
        RepulsiveForceUpdate();
        AttractiveForceUpdate();
        for(int ci=0;ci<NCELLS;ci++){
            for(int vi=0;vi<Cells[ci].NV;vi++){
                Cells[ci].X[vi] += dt*Cells[ci].Fx[vi];
                Cells[ci].Y[vi] += dt*Cells[ci].Fy[vi];
            }
        }
    }

    double monolayer::GetPackingFraction(){
        double sumCellArea=0, BoxArea = L*L;
        for(int i=0;i<NCELLS;i++){
            sumCellArea += Cells[i].GetArea();
        }

        return sumCellArea/BoxArea;
    }
    double monolayer::GetVertexKineticEnergy(){
        double K = 0.0;
        for(int ci=0;ci<NCELLS;ci++){
            for(int vi=0;vi<Cells[ci].NV;vi++){
                K += 0.5*abs(Cells[ci].vx[vi] + Cells[ci].vy[vi]);
            }
        }

        return K;
    }
    void monolayer::FindOverlaps(int ci, int cj){
        if(ci > NCELLS || cj > NCELLS){
            cout << "Error, cell index not found!" << endl;
            return;
        }
        int NVi = Cells[ci].NV, NVj = Cells[cj].NV, vi, vj;
        int nNeg=0, nPos=0;
        double x,y,test;
        overlaps.resize(NVi);
        double X[NVj],Y[NVj];
        for(vj=0;vj<NVj;vj++){
            X[vj] = Cells[cj].X[vj];
            Y[vj] = Cells[cj].Y[vj];
            X[vj] -= L * round((Cells[cj].X[0]-Cells[ci].X[0])/L);
            Y[vj] -= L * round((Cells[cj].Y[0]-Cells[ci].Y[0])/L);
        }
        if(Cells[cj].isConvex()){
            //use inclusion test
            for(vi=0;vi<NVi;vi++){
                nNeg = 0; nPos = 0;
                y = Cells[ci].Y[vi]; //point to test
                x = Cells[ci].X[vi]; 
                for(vj=0;vj<NVj;vj++){
                    //D = (x2 - x1) * (yp - y1) - (xp - x1) * (y2 - y1)
                    test = (Cells[cj].X[Cells[cj].ip1[vj]] - Cells[cj].X[vj]) * (y - Cells[cj].Y[vj]) - (x - Cells[cj].X[vj]) * (Cells[cj].Y[Cells[cj].ip1[vj]] - Cells[cj].Y[vj]);
                    if(test < 0)
                        nNeg++;
                    else
                        nPos++;
                    
                }
                if(nNeg == 0 || nPos == 0)
                    overlaps[vi] = true;
                else
                    overlaps[vi] = false;
            }
 
        }
        else{
            bool odd;
            //use crossings test
            for(vi=0;vi<NVi;vi++){
                x = Cells[ci].X[vi];
                y = Cells[ci].Y[vi];
                int i =0; int j = NVj-1;
                while(i<NVj){
                    i++;
                    if(((Y[i] > y) != (Y[j]>y)) && (x < ((X[j]-X[i])*(y-Y[i])/Y[j]-Y[i])+ X[i])){
                        bool odd = not odd;
                    }
                    j=i;
                }
                overlaps[vi] = odd;
            }
        }
        
    }
    void monolayer::RepulsiveForceUpdate(){
        int ci,cj,vi,vj;
        double cxi,cyi,cxj,cyj,dx,dy,sij,rij;
        double rho0 = sqrt(Cells[0].a0), ftmp;
        for(ci=0;ci<NCELLS;ci++){
            cxi = Cells[ci].GetCenterX();
            cyi = Cells[ci].GetCenterY();
            for(cj=0;cj<NCELLS;cj++){
                if(ci!=cj){
                    FindOverlaps(ci,cj);
                    cxj = Cells[cj].GetCenterX();
                    cyj = Cells[cj].GetCenterY();
                    dx = cxj-cxi;
                    dx -= L*round(dx/L);
                    dy  = cyj-cyi;
                    dy -= L*round(dx/L);
                    rij = sqrt(dx*dx + dy*dy);
                    for(vi=0;vi<Cells[ci].NV;vi++){
                        for(vj=0;vj<Cells[cj].NV;vj++){
                            sij = Cells[cj].r0 + Cells[ci].r0;
                            if(overlaps[vi]){
                                ftmp = Kc * (1-(rij/sij)) * (rho0/sij);
                                Cells[ci].Fx[vi] -= ftmp * (dx/rij);
                                Cells[ci].Fy[vi] -= ftmp * (dy/rij);
                                Cells[cj].Fx[vj] += ftmp * (dx/rij);
                                Cells[cj].Fy[vj] += ftmp * (dx/rij);
                            }
                        }
                    }
                } 
            }
        }
    }
    void monolayer::AttractiveForceUpdate(){

    }
    void monolayer::ResetForces(){
        for(int ci=0;ci<NCELLS;ci++){
            for(int vi=0;vi<Cells[ci].NV;vi++){
                Cells[ci].Fx[vi] = 0.0;
                Cells[ci].Fy[vi] = 0.0;
                Cells[ci].vx[vi] = 0.0;
                Cells[ci].vy[vi] = 0.0;

            }
        }
    }


}