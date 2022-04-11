#include <cmath>
#include <monolayer.hpp>
#include <vector>
#include <iostream>

namespace DPM{
    monolayer::monolayer(std::vector<DPM::Cell> inputCells, double phi0input){
        Cells = inputCells;
        phi0 = phi0input;
        NCELLS = Cells.size();
        VertDOF = 0;
        double sumareas = 0.0;
        for(int i=0;i<NCELLS;i++){
            VertDOF += Cells[i].NV;
            sumareas += Cells[i].GetArea();
        }
        L = sqrt(sumareas)/phi0;
        Kc = 1.0;
        Ftol = 1e-5;
        dt0 = 1e-3;
        U = 0.0;
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
        double oldU = 100, dU = 100;
        while(abs(dU) > 1e-6){
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
        dtmin = 1e-2 * dt;
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
                    P+= Cells[i].vx[vi]*Cells[i].Fx[vi] + Cells[i].vy[vi]*Cells[i].Fy[vi];
                }
            }

            if(P>0){
                npPos++;
                npNeg = 0;
                if (npPos > NDELAY && dt*finc < dtmax){
                    dt *= finc;
                }
                alpha *= falpha;
            }
            else{
                npPos = 0;
                npNeg++;

                if(npNeg > NNEGMAX){
                    std::cerr << "ERROR during FIRE, P<0 for too long..." << std::endl;
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

                if(fireit > NDELAY && dt*fdec > dtmin){
                    alpha = alpha0;
                    if(dt*fdec > dtmin){
                        dt *= fdec;
                    }
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
                    fnorm += pow(Cells[ci].Fx[vi],2.0) + pow(Cells[ci].Fy[vi],2.0);
                    vnorm += pow(Cells[ci].vx[vi],2.0) + pow(Cells[ci].vy[vi],2.0);
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
                }
            }

            InteractingForceUpdate();

            for(ci=0;ci<NCELLS;ci++){
                for(vi=0;vi<Cells[ci].NV;vi++){
                     Cells[ci].vx[vi] += 0.5 * dt * Cells[ci].Fx[vi];
                     Cells[ci].vy[vi] += 0.5 * dt * Cells[ci].Fy[vi];
                }
            }

            fcheck = 0.0;
            for(ci=0;ci<NCELLS;ci++){
                for(vi=0;vi<Cells[ci].NV;vi++){
                    fcheck += pow(Cells[ci].Fx[vi],2) + pow(Cells[ci].Fy[vi],2);
                }
            }
            fcheck = sqrt(fcheck)/VertDOF;
            fireit++;
            if(fireit % 1000 == 0){
                std::cout << "FIRE progress update:" << std::endl;
                std::cout << "	 iterations = " << fireit << std::endl;
                std::cout << "	 fcheck = " << fcheck << std::endl;
                std::cout << "	 dt = " << dt << std::endl;
                std::cout << "	 P = " << P << std::endl;
                std::cout << "	 alpha = " << alpha << std::endl;
            }
        }

        if(fireit == itmax){
            std::cerr << "(!) Fire Minimization did not converge" << std::endl;
        }
        else{
            std::cout << "FIRE Minimization Converged in "<< fireit << " timesteps" << std::endl;
            std::cout << "	 iterations = " << fireit << std::endl;
            std::cout << "	 fcheck = " << fcheck << std::endl;
            std::cout << "	 dt = " << dt << std::endl;
            std::cout << "	 P = " << P << std::endl;
            std::cout << "	 alpha = " << alpha << std::endl;
        }
    }
    void monolayer::UpdateEuler(int nsteps, double dt){
        int ci, vi;
        for(int step=0;step<nsteps;step++){
            ResetForces();
            //For all the cells update the forces based on intracellular shapes
            for(ci=0;ci<NCELLS;ci++){
                Cells[ci].UpdateShapeForces();
                Cells[ci].DrivingForceUpdate(dt);
            }
            //now update the forces they have on eachother
            InteractingForceUpdate();
            //Use that force to change the position of the vertices
            for(ci=0;ci<NCELLS;ci++){
                for(vi=0;vi<Cells[ci].NV;vi++){
                    Cells[ci].X[vi] += dt*Cells[ci].Fx[vi];
                    Cells[ci].Y[vi] += dt*Cells[ci].Fy[vi];
                    Cells[ci].vx[vi] = 0.5*dt*Cells[ci].Fx[vi];
                    Cells[ci].vy[vi] = 0.5*dt*Cells[ci].Fy[vi];
                }
            }
        }
    }

    void monolayer::UpdateVV(int nsteps, double dt){
      std::vector<double> oldFx;
      std::vector<double> oldFy;
      int ci, vi,step;
      for(step=0;step<nsteps;step++){
        for(ci=0;ci<NCELLS;ci++){
            oldFx.resize(Cells[ci].NV);
            oldFy.resize(Cells[ci].NV);
            for(vi=0;vi<Cells[ci].NV;vi++){
            oldFx[vi] = Cells[ci].Fx[vi];
            oldFy[vi] = Cells[ci].Fy[vi];
            //Position update
            Cells[ci].X[vi] += dt*Cells[ci].vx[vi] + 0.5*dt*dt*(Cells[ci].Fx[vi]);
            Cells[ci].Y[vi] += dt*Cells[ci].vy[vi] + 0.5*dt*dt*(Cells[ci].Fy[vi]);
            }
            Cells[ci].UpdateShapeForces();
            for(vi=0;vi<Cells[ci].NV;vi++){
            //Velocity update
            Cells[ci].vx[vi] += 0.5 * dt * (Cells[ci].Fx[vi]+oldFx[vi]);
            Cells[ci].vy[vi] += 0.5 * dt * (Cells[ci].Fy[vi]+oldFy[vi]);
            }
        }
        InteractingForceUpdate();
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
            std::cout << "Error, cell index not found!" << std::endl;
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
            //use inclusion test: is the point always to the right of the left of vector segements?
            for(vi=0;vi<NVi;vi++){
                nNeg = 0; nPos = 0;
                y = Cells[ci].Y[vi];
                x = Cells[ci].X[vi];
                for(vj=0;vj<NVj;vj++){
                    //D = (x2 - x1) * (yp - y1) - (xp - x1) * (y2 - y1)
                    test = (X[Cells[cj].ip1[vj]] - X[vj]) * (y - Y[vj]) - (x - X[vj]) * (Y[Cells[cj].ip1[vj]] - Y[vj]);
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
            //use crossing test: If you draw a line out from point, does it cross an odd number of times?
            int i, j;
            for(vi=0;vi<NVi;vi++){
                bool inside = false;
                y = Cells[ci].Y[vi];
                x = Cells[ci].X[vi];
                for (i = 0, j = NVj-1; i < NVj; j = i++) {
                    if ( ((Y[i]>y) != (Y[j]>y)) &&
                    (x < (X[j]-X[i]) * (y-Y[i]) / (Y[j]-Y[i]) + X[i]) )
                    {
                        inside = !inside;
                    }
                }
                overlaps[vi] = inside;
            }
        }
    }
    void monolayer::InteractingForceUpdate(){
        int ci,cj,vi,vj;
        double cxi,cyi,cxj,cyj,dx,dy,sij,rij;
        double shellij,cutij,xij,kint,ftmp;
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
                    dy -= L*round(dy/L);
                    rij = sqrt(dx*dx + dy*dy);
                    for(vi=0;vi<Cells[ci].NV;vi++){
                        for(vj=0;vj<Cells[cj].NV;vj++){
                            sij = sqrt(Cells[ci].a0) + sqrt(Cells[cj].a0);
                            xij = rij/sij;
                            ftmp = 0.0;
                            if(overlaps[vi]){
                                //if overlapping repell the particles away from eachothers centers
                                ftmp = Kc*(1-xij)/(sqrt(Cells[ci].a0)/sij);
                                U += 0.5*Kc * pow(1-xij,2);
                            }
                            else{
                                shellij = (1.0+Cells[ci].l2[vi])*sij;
                                cutij = (1.0+Cells[ci].l1[vi])*sij;
                                //if not overlapping, see if they are within attractive distance
                                if(rij >= cutij && rij <= shellij){
                                    kint = (Kc*Cells[ci].l1[vi])/(Cells[ci].l2[vi]-Cells[ci].l1[vi]);
                                    ftmp = kint * (xij - (1.0 - Cells[ci].l2[vi]))/sij;
                                    U += -0.5*kint*pow(1.0 + Cells[ci].l2[vi] - xij,2.0);
                                }
                            }
                            //force update
                            Cells[ci].Fx[vi] -= ftmp * (dx/rij);
                            Cells[ci].Fy[vi] -= ftmp * (dy/rij);
                            Cells[cj].Fx[vj] += ftmp * (dx/rij);
                            Cells[cj].Fy[vj] += ftmp * (dy/rij);
                        }
                    }
                }
            }
        }
    }
    void monolayer::ResetForces(){
        U = 0;
        for(int ci=0;ci<NCELLS;ci++){
            for(int vi=0;vi<Cells[ci].NV;vi++){
                Cells[ci].Fx[vi] = 0.0;
                Cells[ci].Fy[vi] = 0.0;
                Cells[ci].vx[vi] = 0.0;
                Cells[ci].vy[vi] = 0.0;
            }
        }
    }


    //Biological Processes
    void monolayer::CellDivision(int ci){
        //first do some checks
        if(ci > NCELLS-1){
            std::cout << "Error: Index out of range" << std::endl;
            return;
        }
        double packing = GetPackingFraction();
        double packingincrease = Cells[ci].GetArea() / (L*L);
        double newpacking = packingincrease + packing;
        if(newpacking > 0.95){
            std::cout << "Warning packing fraction is > 0.95\n";
        }
        else if(newpacking > 1.0){
            std::cout << "Cannot divide: Would make packing greater than 100%\n";
            return;
        }


        double movex, movey, radii;
        NCELLS ++;
        radii = Cells[ci].r0;
        movex = (radii * cos(Cells[ci].psi))/2.0;
        movey = (radii * sin(Cells[ci].psi))/2.0;
        Cells[ci].ResetForces();
        DPM::Cell DaugterCell = Cells[ci];
        for(int vi = 0; vi < Cells[ci].NV; vi++){
            Cells[ci].X[vi] += movex;
            Cells[ci].Y[vi] += movey;
            DaugterCell.X[vi] -= movex;
            DaugterCell.Y[vi] -= movey; 
        }

        Cells.push_back(DaugterCell);
    }

}
