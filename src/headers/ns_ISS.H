/*
 * This file is part of the Nano-Shell Simulation Project.
 * 
 * Copyright (C) 2025 Alessandro Veltri
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */
#include "extract.H"

std::complex<double> gammaD (std::complex<double>** A, double Ome){
    std::complex<double> gamdD = (A[1][1]-img*Ome)*(A[2][2]-img*Ome)-A[1][2]*A[2][1];
    return gamdD;
    }

std::complex<double> gamma1 (std::complex<double>** A, double Ome){
    std::complex<double> gam1 = (A[1][2]*A[2][0]-A[1][0]*(A[2][2]-img*Ome))/gammaD(A, Ome);
    return gam1;
    }

std::complex<double> gamma2 (std::complex<double>** A, double Ome){
    std::complex<double> gam2 = (A[1][0]*A[2][1]-A[2][0]*(A[1][1]-img*Ome))/gammaD(A, Ome);
    return gam2;
    }

    
std::complex<double> gimmeXi (std::complex<double>* p0, std::complex<double> gam1, std::complex<double> gam2, double tau1){
    std::complex<double> xi = tau1*(p0[0]+p0[1]*gam1+p0[2]*gam2);
    return xi;
    }
    
double* Dsqr(std::complex<double> p11, std::complex<double> p12, std::complex<double> p21, std::complex<double> p22, std::complex<double> ome, std::complex<double> gamd) {
    double* D;
    D = new double[5];

    D[0] = real(pow(ome, 8.) +
           (8. * pow(gamd, 2.) + 2. * p22 + 2. * p11) * pow(ome, 6.) +
           (16. * pow(gamd, 4.) + (8. * p22 + 8. * p11) * pow(gamd, 2.) + pow(p22, 2.) + 4. * p11 * p22 - 2. * p12 * p21 + pow(p11, 2.)) * pow(ome, 4.) +
           ((4. * pow(p22, 2.) + 8. * p12 * p21 + 4. * pow(p11, 2.)) * pow(gamd, 2.) + 2. * p11 * pow(p22, 2.) +
            (2. * pow(p11, 2.) - 2. * p12 * p21) * p22 - 2. * p11 * p12 * p21) * pow(ome, 2.) +
           pow(p11, 2.) * pow(p22, 2.) - 2. * p11 * p12 * p21 * p22 + pow(p12, 2.) * pow(p21, 2.));

    D[1] = real(-8. * pow(ome, 7.) +
           (-48. * pow(gamd, 2.) - 12. * p22 - 12. * p11) * pow(ome, 5.) +
           (-64. * pow(gamd, 4.) + (-32. * p22 - 32. * p11) * pow(gamd, 2.) - 4. * pow(p22, 2.) - 16. * p11 * p22 + 8. * p12 * p21 - 4. * pow(p11, 2.)) * pow(ome, 3.) +
           ((-8. * pow(p22, 2.) - 16. * p12 * p21 - 8. * pow(p11, 2.)) * pow(gamd, 2.) - 4. * p11 * pow(p22, 2.) +
            (4. * p12 * p21 - 4. * pow(p11, 2.)) * p22 + 4. * p11 * p12 * p21) * ome);

    D[2] = real(24. * pow(ome, 6.) +
           (104. * pow(gamd, 2.) + 24. * p22 + 24. * p11) * pow(ome, 4.) +
           (96. * pow(gamd, 4.) + (40. * p22 + 40. * p11) * pow(gamd, 2.) +
            4. * pow(p22, 2.) + 16. * p11 * p22 - 8. * p12 * p21 + 4. * pow(p11, 2.)) * pow(ome, 2.) +
           (4. * pow(p22, 2.) + 8. * p12 * p21 + 4. * pow(p11, 2.)) * pow(gamd, 2.));

    D[3] = real(-32. * pow(ome, 5.) +
           (-96. * pow(gamd, 2.) - 16. * p22 - 16. * p11) * pow(ome, 3.) +
           ((-16. * p22 - 16. * p11) * pow(gamd, 2.) - 64. * pow(gamd, 4.)) * ome);

    D[4] = real(16. * pow(ome, 4.) + 32. * pow(gamd, 2.) * pow(ome, 2.) + 16. * pow(gamd, 4.));

    return D;
}

double** Rnum(std::complex<double>* p1, std::complex<double>* p2, std::complex<double> ome, std::complex<double> gamd) {
    double** R;
    R = new double*[3];
    for (int i = 0; i < 3; i++)
        R[i] = new double[4];

    R[1][0] = real(-p1[0] * pow(ome, 6.) +
              (-4. * p1[0] * pow(gamd, 2.) - 2. * p1[0] * p2[2] + p1[2] * p2[0] - p1[0] * p1[1]) * pow(ome, 4.) +
              ((-4. * p1[2] * p2[0] - 4. * p1[0] * p1[1]) * pow(gamd, 2.) - p1[0] * pow(p2[2], 2.) +
               (p1[2] * p2[0] - 2. * p1[0] * p1[1]) * p2[2] + p1[0] * p1[2] * p2[1] + p1[1] * p1[2] * p2[0]) * pow(ome, 2.) -
              p1[0] * p1[1] * pow(p2[2], 2.) + (p1[0] * p1[2] * p2[1] + p1[1] * p1[2] * p2[0]) * p2[2] - pow(p1[2], 2.) * p2[0] * p2[1]);

    R[1][1] = real(6. * p1[0] * pow(ome, 5.) +
              (16. * p1[0] * pow(gamd, 2.) + 8. * p1[0] * p2[2] - 4. * p1[2] * p2[0] + 4. * p1[0] * p1[1]) * pow(ome, 3.) +
              ((8. * p1[2] * p2[0] + 8. * p1[0] * p1[1]) * pow(gamd, 2.) + 2. * p1[0] * pow(p2[2], 2.) +
               (4. * p1[0] * p1[1] - 2. * p1[2] * p2[0]) * p2[2] - 2. * p1[0] * p1[2] * p2[1] - 2. * p1[1] * p1[2] * p2[0]) * ome);

    R[1][2] = real(-12. * p1[0] * pow(ome, 4.) +
              (-20. * p1[0] * pow(gamd, 2.) - 8. * p1[0] * p2[2] + 4. * p1[2] * p2[0] - 4. * p1[0] * p1[1]) * pow(ome, 2.) +
              (-4. * p1[2] * p2[0] - 4. * p1[0] * p1[1]) * pow(gamd, 2.));

    R[1][3] = real(8. * p1[0] * pow(ome, 3.) + 8. * p1[0] * pow(gamd, 2.) * ome);

    R[2][0] = real(-p2[0] * pow(ome, 6.) +
              (-4. * p2[0] * pow(gamd, 2.) - p2[0] * p2[2] + p1[0] * p2[1] - 2. * p1[1] * p2[0]) * pow(ome, 4.) +
              ((-4. * p2[0] * p2[2] - 4. * p1[0] * p2[1]) * pow(gamd, 2.) + (p1[0] * p2[1] - 2. * p1[1] * p2[0]) * p2[2] +
               (p1[2] * p2[0] + p1[0] * p1[1]) * p2[1] - pow(p1[1], 2.) * p2[0]) * pow(ome, 2.) +
              (p1[0] * p1[1] * p2[1] - pow(p1[1], 2.) * p2[0]) * p2[2] - p1[0] * p1[2] * pow(p2[1], 2.) + p1[1] * p1[2] * p2[0] * p2[1]);

    R[2][1] = real(6. * p2[0] * pow(ome, 5.) +
              (16. * p2[0] * pow(gamd, 2.) + 4. * p2[0] * p2[2] - 4. * p1[0] * p2[1] + 8. * p1[1] * p2[0]) * pow(ome, 3.) +
              ((8. * p2[0] * p2[2] + 8. * p1[0] * p2[1]) * pow(gamd, 2.) + (4. * p1[1] * p2[0] - 2. * p1[0] * p2[1]) * p2[2] +
               (-2. * p1[2] * p2[0] - 2. * p1[0] * p1[1]) * p2[1] + 2. * pow(p1[1], 2.) * p2[0]) * ome);

    R[2][2] = real(-12. * p2[0] * pow(ome, 4.) +
              (-20. * p2[0] * pow(gamd, 2.) - 4. * p2[0] * p2[2] + 4. * p1[0] * p2[1] - 8. * p1[1] * p2[0]) * pow(ome, 2.) +
              (-4. * p2[0] * p2[2] - 4. * p1[0] * p2[1]) * pow(gamd, 2.));

    R[2][3] = real(8. * p2[0] * pow(ome, 3.) + 8. * p2[0] * pow(gamd, 2.) * ome);

    return R;
}

double** Inum(std::complex<double>* p1, std::complex<double>* p2, std::complex<double> ome, std::complex<double> gamd) {
    double** I;
    I = new double*[3];
    for (int i = 0; i < 3; i++)
        I[i] = new double[4];

    I[1][0] = real(2. * p1[0] * gamd * pow(ome, 5.) +
              (8. * p1[0] * pow(gamd, 3.) + (4. * p1[0] * p2[2] - 4. * p1[2] * p2[0]) * gamd) * pow(ome, 3.) +
              (2. * p1[0] * pow(p2[2], 2.) - 2. * p1[2] * p2[0] * p2[2] + 2. * p1[0] * p1[2] * p2[1] - 2. * p1[1] * p1[2] * p2[0]) * gamd * ome);

    I[1][1] = real(-10. * p1[0] * gamd * pow(ome, 4.) +
              ((12. * p1[2] * p2[0] - 12. * p1[0] * p2[2]) * gamd - 24. * p1[0] * pow(gamd, 3.)) * pow(ome, 2.) +
              (-2. * p1[0] * pow(p2[2], 2.) + 2. * p1[2] * p2[0] * p2[2] - 2. * p1[0] * p1[2] * p2[1] + 2. * p1[1] * p1[2] * p2[0]) * gamd);

    I[1][2] = real(16. * p1[0] * gamd * pow(ome, 3.) +
              (24. * p1[0] * pow(gamd, 3.) + (8. * p1[0] * p2[2] - 8. * p1[2] * p2[0]) * gamd) * ome);

    I[1][3] = real(-8. * p1[0] * gamd * pow(ome, 2.) - 8. * p1[0] * pow(gamd, 3.));

    I[2][0] = real(2. * p2[0] * gamd * pow(ome, 5.) +
              (8. * p2[0] * pow(gamd, 3.) + (4. * p1[1] * p2[0] - 4. * p1[0] * p2[1]) * gamd) * pow(ome, 3.) +
              (-2. * p1[0] * p2[1] * p2[2] + (2. * p1[2] * p2[0] - 2. * p1[0] * p1[1]) * p2[1] + 2. * pow(p1[1], 2.) * p2[0]) * gamd * ome);

    I[2][1] = real(-10. * p2[0] * gamd * pow(ome, 4.) +
              ((12. * p1[0] * p2[1] - 12. * p1[1] * p2[0]) * gamd - 24. * p2[0] * pow(gamd, 3.)) * pow(ome, 2.) +
              (2. * p1[0] * p2[1] * p2[2] + (2. * p1[0] * p1[1] - 2. * p1[2] * p2[0]) * p2[1] - 2. * pow(p1[1], 2.) * p2[0]) * gamd);

    I[2][2] = real(16. * p2[0] * gamd * pow(ome, 3.) +
              (24. * p2[0] * pow(gamd, 3.) + (8. * p1[1] * p2[0] - 8. * p1[0] * p2[1]) * gamd) * ome);

    I[2][3] = real(-8. * p2[0] * gamd * pow(ome, 2.) - 8. * p2[0] * pow(gamd, 3.));

    return I;
}

double *cffcOme (double *D, double **R, double **I, std::complex<double>* p0, double tau2, double ome, double ome_g){
    double *c;
    c = new double[5];

    c[0] = real(-(p0[2]*I[2][0]+p0[1]*I[1][0])*(ome-ome_g)*tau2+p0[2]*R[2][0]+p0[1]*R[1][0]+D[0]*p0[0]);
    c[1] = real(-(p0[2]*I[2][1]+p0[1]*I[1][1])*(ome-ome_g)*tau2+p0[2]*(I[2][0]*tau2+R[2][1])+
                                                          p0[1]*(I[1][0]*tau2+R[1][1])+p0[0]*D[1]);
    c[2] = real(-(p0[2]*I[2][2]+p0[1]*I[1][2])*(ome-ome_g)*tau2+p0[2]*(I[2][1]*tau2+R[2][2])+
                                                          p0[1]*(I[1][1]*tau2+R[1][2])+p0[0]*D[2]);
    c[3] = real(-(p0[2]*I[2][3]+p0[1]*I[1][3])*(ome-ome_g)*tau2+p0[2]*(I[2][2]*tau2+R[2][3])+
                                                          p0[1]*(I[1][2]*tau2+R[1][3])+p0[0]*D[3]);
    c[4] = real((p0[2]*I[2][3]+p0[1]*I[1][3])*tau2+p0[0]*D[4]);

    return c;
    }

std::vector<std::pair<double,std::complex<double>>> gimme_emi_kap(nanosphere ns, char* mdl, char* mtl, char* hst, double omemi, double omema, int omeN=1000, char* sol=0, double rho=0){
        std::vector<std::pair<double,std::complex<double>>> vkamp;
        int omi;
        std::complex<double> GG, OmeH, GamP, OmeP, **A, *p0, *p1, *p2, gam1, gam2, xi, *kap, *odv;
        
        double *nv, omeeV, dome  = (omema-omemi)/omeN, ome;
        double  eps_b, eps_s, tau2, gamd, ome_g;

        if ((sol != NULL)) {
            eps_s=ns.set_host(sol);
            } else eps_s=0;
        eps_b=ns.set_host(hst);

        A = new std::complex<double>*[3];
        for(int i = 0; i <= 2; i++)
            A[i] = new std::complex<double>[3];
        
        kap  = new std::complex<double>[3];
        nv   = new double[4];
        odv  = new std::complex<double> [3];

        // normalized variables
        nv = ns.normalized_variables ();
        tau2 = nv[1];
        gamd = nv[2];
        ome_g = nv[3];
        GG   = ns.set_GamG(ns.G, tau2);

        p0 = pcfc0(ns.eps_inf, eps_b, eps_s, rho);
        p1 = pcfc1(ns.eps_inf, eps_b, eps_s, rho);
        p2 = pcfc2(ns.eps_inf, eps_b, eps_s, rho);
        
        for (omi=0; omi<=omeN; omi++){
            omeeV = omemi + omi*dome;
        
            ome  = omeeV/ns.Ome_p;

            odv  = ns.set_ome_dep_vrbls(ome, ome_g, tau2, gamd);
            OmeH = odv[0];
            OmeP = odv[1];
            GamP = odv[2];
    
            A   = coefficients(OmeH, OmeP, GG, GamP, p0, p1, p2);
            kap  = eigenvalues(A);
            
            double maxVal = kap[0].real();
            int    imax = 0;

            for (int i = 1; i < 3; i++) {
                if (kap[i].real() > maxVal) {
                    maxVal = kap[i].real();
                    imax = i;
                    }
                }
            vkamp.push_back(std::make_pair(omeeV, kap[imax]));
            }
        return vkamp;
        }

double find_Omega(nanosphere ns, double omeeV, char* hst, char* sol=0, double rho=0){
    std::complex<double> GG, OmeH, GamP, OmeP, **A, *p0, *p1, *p2, gam1, gam2, xi, *vOme, *kap, *odv;
    double tau1, tau2, gamd, ome, ome_g, eps_b, eps_s, *D, **R, **I, *c;
    double Ome, equ1, equ2, kamp, dmem, *nv;
    
    eps_b=ns.set_host(hst);
    eps_s=ns.set_host(sol);
    
    A = new std::complex<double>*[3];
    for(int i = 0; i <= 2; i++)
        A[i] = new std::complex<double>[3];

    c = new double [5];

    D = new double [5];

    R = new double*[3];
    for(int i = 0; i < 3; i++)
        R[i] = new double[4];

    I = new double*[3];
    for(int i = 0; i < 3; i++)
        I[i] = new double[4];

    vOme = new std::complex<double>[4];
    kap  = new std::complex<double>[3];
    nv   = new double[4];
    odv  = new std::complex<double> [3];

    // normalized variables
    nv = ns.normalized_variables ();
    tau1 = nv[0];
    tau2 = nv[1];
    gamd = nv[2];
    ome_g = nv[3];
    GG   = ns.set_GamG(ns.G, tau2);

    ome  = omeeV/ns.Ome_p;

    odv  = ns.set_ome_dep_vrbls(ome, ome_g, tau2, gamd);
    OmeH = odv[0];
    OmeP = odv[1];
    GamP = odv[2];
    
    p0 = pcfc0(ns.eps_inf, eps_b, eps_s, rho);
    p1 = pcfc1(ns.eps_inf, eps_b, eps_s, rho);
    p2 = pcfc2(ns.eps_inf, eps_b, eps_s, rho);

    A   = coefficients(OmeH, OmeP, GG, GamP, p0, p1, p2);
        
    D    = Dsqr(p1[1], p1[2], p2[1], p2[2], ome, gamd);
    R    = Rnum(p1, p2, ome, gamd);
    I    = Inum(p1, p2, ome, gamd);
    c    = cffcOme (D, R, I, p0, tau2, ome, ome_g);
    

    
    vOme = solve4(c);
    kap  = eigenvalues(A);

    Ome=0.;

    kamp = 66.;
    dmem = 0.1;
    for (int i = 0; i < 3; i++) if (kap[i].real() > 0) kamp = kap[i].imag();
    if (kamp != 66){
        for (int i = 0; i < 5; i++){
            if (vOme[i].imag()==0){
                if (fabs(vOme[i].real()-kamp)<dmem) {
                    Ome = vOme[i].real();
                    dmem = fabs(vOme[i].real()-kamp);         
                    }
                }
            }
        }   
    gam1 = gamma1(A, Ome);
    gam2 = gamma2(A, Ome);
    xi   = gimmeXi(p0, gam1, gam2, tau1);
    equ1 = norm(gammaD(A, Ome))*(real(xi)+tau2*(Ome-ome+ome_g)*imag(xi))/tau1;
    equ2 = c[0]+c[1]*Ome+c[2]*pow(Ome,2)+c[3]*pow(Ome,3)+c[4]*pow(Ome,4);
    if (abs(equ1-equ2)>=1.e-18 && fabs(Ome)> 1.e-30){
        std::cout.precision(20);
        std::cout<<"> Fatal! The error on the coefficient equation is larger than "<<abs(equ1-equ2)<<std::endl;
        std::cout<<"> ome ="<<omeeV<<std::endl;
        exit(-66);
        }
    
    
    
    return Ome;
    }

double gOmega(double ome, double ome_g, double tau2, std::complex<double> xi, double Ome){
    double g = ome - ome_g - real(xi)/(tau2*imag(xi)) - Ome;
    return g;
    }
    
double find_Omega1(nanosphere ns, double omeeV, char* hst, char* sol, double rho){
    std::complex<double> GG, OmeH, GamP, OmeP, **A, *p0, *p1, *p2, cq0s, gam1, gam2, xi, *kap, *odv;
    double tau1, tau2, gamd, ome, ome_g, eps_b, eps_s;
    double Ome, Ome0, *nv, ddOme, gOme, gOmep, gOmem, DgOme, dOme;
    int ite=0;

    eps_b=ns.set_host(hst);
    eps_s=ns.set_host(sol);
    
    A = new std::complex<double>*[3];
    for(int i = 0; i <= 2; i++)
        A[i] = new std::complex<double>[3];

    kap  = new std::complex<double>[3];
    nv   = new double[4];
    odv  = new std::complex<double> [3];

    
    // normalized variables
    nv = ns.normalized_variables ();
    tau1 = nv[0];
    tau2 = nv[1];
    gamd = nv[2];
    ome_g = nv[3];

    GG   = ns.set_GamG(ns.G, tau2);
    
    ome  = omeeV/ns.Ome_p;

    odv  = ns.set_ome_dep_vrbls(ome, ome_g, tau2, gamd);
    OmeH = odv[0];
    OmeP = odv[1];
    GamP = odv[2];
    
    p0 = pcfc0(ns.ceps_inf, eps_b, eps_s, rho);
    p1 = pcfc1(ns.ceps_inf, eps_b, eps_s, rho);
    p2 = pcfc2(ns.ceps_inf, eps_b, eps_s, rho);
    
    A   = coefficients(OmeH, OmeP, GG, GamP, p0, p1, p2);

    kap  = eigenvalues(A);

    Ome0  = 0.;
    for (int i = 0; i < 3; i++) if (kap[i].real() > 0) Ome0 = kap[i].imag();
    
    cq0s.imag(666.);
    Ome = Ome0;
    ddOme = 1.e-4;

    while(fabs(cq0s.imag()) > 1.e-15 && fabs(Ome)>1.e-30){ // ideally should be 1.e-10 or lower
            ite++;
            gam1   = gamma1(A, Ome);
            gam2   = gamma2(A, Ome);
            xi     = gimmeXi(p0, gam1, gam2, tau1);
            gOme   = gOmega(ome, ome_g, tau2, xi, Ome);
            gam1   = gamma1(A, Ome + ddOme);
            gam2   = gamma2(A, Ome + ddOme);
            xi     = gimmeXi(p0, gam1, gam2, tau1);
            gOmep  = gOmega(ome, ome_g, tau2, xi, Ome + ddOme);
            
            gam1   = gamma1(A,Ome - ddOme);
            gam2   = gamma2(A,Ome - ddOme);
            xi     = gimmeXi(p0, gam1, gam2, tau1);
            gOmem  =  gOmega(ome, ome_g, tau2, xi, Ome - ddOme);
            
            DgOme  = .5*(gOmep-gOmem)/ddOme;
            dOme   = -gOme/DgOme;
            Ome = Ome+dOme;
            gam1   = gamma1(A,Ome);
            gam2   = gamma2(A,Ome);
            xi     = gimmeXi(p0, gam1, gam2, tau1);
            cq0s   = (GG*xi+tau1*(OmeH-img*Ome))/(GG*xi*xi.imag());
            }
    return Ome;
    }
    
double *ISS_results(nanosphere ns, double Ome, double ome1, double ome2, double omeeV, char* hst, char* sol=0, double rho=0){
    std::complex<double> GG, OmeH, GamP, OmeP, gam1, gam2, xi, **A, cq0s, *vec, *odv;
    double *res, q0s, p1s, p2s, p3s, N, tau1, tau2, gamd, ome, ome_g, *nv;
    std::complex<double>* p0, *p1, *p2, *p3;
    double eps_b, eps_s, *fro;
    eps_b=ns.set_host(hst);
    eps_s=ns.set_host(sol);

    res = new double[7];
    vec = new std::complex<double>[4];
    
    A = new std::complex<double>*[3];
    for(int i = 0; i <= 2; i++)
        A[i] = new std::complex<double>[3];
    
    nv  = new double[4];
    odv = new std::complex<double> [3];
    
    fro=ns.frohlich(2., 3.4, eps_b, eps_s, rho);
    // normalized variables
    nv = ns.normalized_variables ();
    tau1 = nv[0];
    tau2 = nv[1];
    gamd = nv[2];
    ome_g = nv[3];
    GG = ns.set_GamG(ns.G, tau2);
    
    ome  = omeeV/ns.Ome_p;

    odv  = ns.set_ome_dep_vrbls(ome, ome_g, tau2, gamd);
    OmeH = odv[0];
    OmeP = odv[1];
    GamP = odv[2];
    
    p0 = pcfc0(ns.eps_inf, eps_b, eps_s, rho);
    p1 = pcfc1(ns.eps_inf, eps_b, eps_s, rho);
    p2 = pcfc2(ns.eps_inf, eps_b, eps_s, rho);
    p3 = pcfc3(ns.eps_inf, eps_b, eps_s, rho);

    A   = coefficients(OmeH, OmeP, GG, GamP, p0, p1, p2);

    gam1 = gamma1(A, Ome);
    gam2 = gamma2(A, Ome);
    xi   = gimmeXi(p0, gam1, gam2, tau1);

    
    if (ome<ome1 || ome>ome2){
        cq0s.real(0);
        cq0s.imag(0);
        } else cq0s = 1./xi.imag()+tau1*(OmeH-img*Ome)/(GG*xi*xi.imag());
    
    if (fabs(cq0s.imag())>5.e-5 && fabs(ns.G) > fabs(fro[1])){
        std::cout.precision(20);
        std::cout<<"> Warning! The error on the imaginary part of |q1|² is larger than 5.e-5"
        <<std::endl<<"> err = "<<cq0s.imag()<<std::endl;
        std::cout<<"> ome ="<<omeeV<<", Ome = "<<Ome<<std::endl;
        } else if (fabs(cq0s.imag())>1. && fabs(ns.G) > fabs(fro[1])){
            std::cout.precision(20);
            std::cout<<"> Fatal! The error on the imaginary part of |q1|² is larger than 1"
            <<std::endl<<"> err = "<<cq0s.imag()<<std::endl;
            std::cout<<"> ome ="<<omeeV<<", Ome = "<<Ome<<std::endl;
            exit(-67); 
            } else q0s = cq0s.real();
    
    res[0] = q0s;
    res[1] = q0s*norm(gam1);
    res[2] = q0s*norm(gam2);
    
    vec[0] = {1,0};
    vec[1] = gam1;
    vec[2] = gam2;
    vec[3] = {0,0};
    
    p1s = norm(numerical_output(0, vec, p1))*q0s;
    res[3] = p1s;
    
    p2s = norm(numerical_output(0, vec, p2))*q0s;
    res[4] = p2s;
    
    p3s = norm(numerical_output(0, vec, p3))*q0s;
    res[5] = p3s;    

    N = 1 - xi.imag()*q0s;
    
    res[6] = N;
    
    return res;
    }


    double find_omeB(nanosphere ns, char* hst, char* sol, double rho, double omemi, double omema, int omeN=1000){
        double omeB, omeeV, dome  = (omema-omemi)/omeN, Ome, Omem=0.;
        double ome, Gmem=ns.G;
        ns.G=10.;
        int omi, ch=1, o1=0, o2=0, oB=0;
        for (omi=0; omi<=omeN; omi++){
            omeeV = omemi + omi*dome;
            Omem = Ome;
            Ome = find_Omega1(ns, omeeV, hst, sol, rho);
            ome  = omeeV/ns.Ome_p;
            if (o1==1 && Omem != 0) ch  = Omem*Ome/fabs(Omem*Ome);
            if (fabs(Ome)>1.e-30 && o1==0){
                o1=1;
                }
            if (fabs(Ome)<1.e-30 && o2==0 && o1==1 && oB==1){
                o2=1;
                }
            if (ch < 0 && o1==1 && oB==0) {
                omeB   = ome;
                oB=1;
                }
            }
        ns.G=Gmem;
        return omeB;
        }

    void find_ome1_ome2_unsafe(nanosphere &ns, double &ome1, double &ome2, double omemi, double omema, char* hst, char* sol, double rho){
        int omi, ch=1, o1=0, o2=0, oB=0;
        double *fro, omeeV, dome, Ome, Omem=0., ome;
        double eps3=ns.set_host(sol);
        double eps_b=ns.set_host(hst);

        fro=ns.frohlich(omemi, omema, eps_b, eps3, rho);
        dome  = (omema-omemi)/10000;
        for (omi=0; omi<=10000; omi++){
            omeeV = omemi + omi*dome;
            Omem = Ome;
            Ome = find_Omega(ns, omeeV, hst, sol, rho);
            ome  = omeeV/ns.Ome_p;
            if (o1==1 && Omem != 0) ch  = Omem*Ome/fabs(Omem*Ome);
            if (fabs(Ome)>1.e-30 && o1==0){
                ome1=ome;
                o1=1;
                }
            if (fabs(Ome)<1.e-30 && o2==0 && o1==1 && oB==1){
                ome2=ome-dome/ns.Ome_p;
                o2=1;
                }
            if (ch < 0 && o1==1 && oB==0) {
                ns.omeB     = ome;
                oB=1;
                }
            }
        if (fabs(ns.G)<fabs(fro[1])) {
            ome1=ns.omeB/ns.Ome_p+dome;
            ome2=ns.omeB/ns.Ome_p-dome;
            }
        if (ome1*ns.Ome_p==omemi || ome2 ==0){
            std::cout<<"Fatal! the frequency range is too small!"<<std::endl;
            exit(78);
            }
        }

    bool find_ome1_ome2(nanosphere &ns, double &ome1, double &ome2, double omemi, double omema, char* hst, char* sol, double rho){
        int omi, ch=1, o1=0, o2=0, oB=0;
        double *fro, omeeV, dome, Ome=0., Omem=0., ome;
        bool found1=false, found2=false, foundB=false;
        double eps3=ns.set_host(sol);
        double eps_b=ns.set_host(hst);

        fro=ns.frohlich(omemi, omema, eps_b, eps3, rho);
        dome  = (omema-omemi)/10000;
        for (omi=0; omi<=10000; omi++){
            omeeV = omemi + omi*dome;
            Omem = Ome;
            Ome = find_Omega(ns, omeeV, hst, sol, rho);
            ome  = omeeV/ns.Ome_p;
            if (o1==1 && Omem != 0) ch  = Omem*Ome/fabs(Omem*Ome);
            if (fabs(Ome)>1.e-30 && o1==0){
                ome1=ome;
                o1=1;
                found1=true;
                }
            if (fabs(Ome)<1.e-30 && o2==0 && o1==1 && oB==1){
                ome2=ome-dome/ns.Ome_p;
                o2=1;
                found2=true;
                }
            if (ch < 0 && o1==1 && oB==0) {
                ns.omeB     = ome;
                oB=1;
                foundB=true;
                }
            }
        if (fabs(ns.G)<fabs(fro[1])) {
            ome1=ns.omeB/ns.Ome_p+dome;
            ome2=ns.omeB/ns.Ome_p-dome;
            found1=true;
            found2=true;
            }
        if (!foundB) {
            ns.omeB = 0.5*(omemi+omema)/ns.Ome_p;
            }
        if (!(found1 && found2)) return false;
        if (ome1*ns.Ome_p==omemi || ome2 ==0) return false;

        return true;
        }
        
    std::vector<std::pair<double,double>> intensity_steady_state(nanosphere ns, char* mdl, char* mtl, char* hst, double omemi, double omema, char* sol, double rho, int omeN=1000){
        int omi;
        double *res, omeeV, dome  = (omema-omemi)/omeN, Ome, myOme, ome, shome;
        double ome1, ome2;
        std::vector<std::pair<double,double>> q0s_left;
        std::vector<std::pair<double,double>> q0s_rgth;
        std::vector<std::pair<double,double>> p1s_left;
        std::vector<std::pair<double,double>> p1s_rgth;
        std::vector<std::pair<double,double>> p3s_left;
        std::vector<std::pair<double,double>> p3s_rgth;
        
        res = new double[7];
                
        std::fstream stat, func, shft;
        stat.open("../data/output/ns_intensity_SS.dat", std::ios::out);
        shft.open("../data/output/shifted.dat", std::ios::out);
        
        stat<<"# * INTENSITY STEADY STATE * "<<std::endl;
        stat<<"#"<<std::endl;
        stat<<"# file generated by ** ns_ISS.H ** "<<std::endl;
        stat<<"#"<<std::endl;
        stat<<"# ------------- "<<std::endl;
        stat<<"# ome (eV)\tOme\t|p1|²\t|p2|²\t|p3|²\tN"<<std::endl;
        stat<<"#"<<std::endl;

        func.open("../data/output/ns_funct_ints_SS.dat", std::ios::out);
        
        func<<"# * INTENSITY STEADY STATE * "<<std::endl;
        func<<"#"<<std::endl;
        func<<"# file generated by ** ns_ISS.H ** "<<std::endl;
        func<<"#"<<std::endl;
        func<<"# ------------- "<<std::endl;
        func<<"# ome (eV)\tOme\t|q0|²\t|q1|²\t|q2|²\tN"<<std::endl;
        func<<"#"<<std::endl;

        /** STEADY STATE **/
        if (!find_ome1_ome2(ns, ome1, ome2, omemi, omema, hst, sol, rho)){
            std::cout<<"Warning! Strict ome1/ome2 scan failed, falling back to unsafe scan."<<std::endl;
            find_ome1_ome2_unsafe(ns, ome1, ome2, omemi, omema, hst, sol, rho);
            }
// ----------------------------------

        dome  = (omema-omemi)/omeN;
        for (omi=0; omi<=omeN; omi++){
            omeeV = omemi + omi*dome;
            Ome = find_Omega(ns, omeeV, hst, sol, rho);
            ome  = omeeV/ns.Ome_p;
            
            res = ISS_results(ns, Ome, ome1, ome2, omeeV, hst, sol, rho);

            if  (Ome<0) myOme=2.*(ome-ns.omeB)-Ome;
                else myOme=Ome;

            shome = omeeV-ns.Ome_p*Ome;
            if (fabs(Ome)>0 || (ome>ome1 && ome<ome2)){
                if (ome<ns.omeB){
                    q0s_rgth.push_back(std::make_pair(shome, res[0]));
                    p1s_rgth.push_back(std::make_pair(shome, res[3]));
                    p3s_rgth.push_back(std::make_pair(shome, res[5]));
                    }
                else{
                    q0s_left.push_back(std::make_pair(shome, res[0]));
                    p1s_left.push_back(std::make_pair(shome, res[3]));
                    p3s_left.push_back(std::make_pair(shome, res[5]));
                    }
                }
            
            stat<<"  "<<std::setw(8)<<std::setiosflags (std::ios::left)<<omeeV<<              //  1 ome
                  "\t"<<std::setw(11)<<std::setiosflags (std::ios::left)<<ns.Ome_p*Ome<<   //  2 Ome
                  "\t"<<std::setw(11)<<std::setiosflags (std::ios::left)<<res[3]<<           //  3 |p1|²
                  "\t"<<std::setw(11)<<std::setiosflags (std::ios::left)<<res[4]<<           //  4 |p2|²
                  "\t"<<std::setw(11)<<std::setiosflags (std::ios::left)<<res[5]<<           //  5 |p3|²
                  "\t"<<std::setw(11)<<std::setiosflags (std::ios::left)<<res[6]<<           //  6 N
                  "\t"<<std::setw(11)<<std::setiosflags (std::ios::left)<<ns.Ome_p*myOme<<   // 12 Ome
                    std::endl;
            func<<"  "<<std::setw(8)<<std::setiosflags (std::ios::left)<<omeeV<<              //  1 ome
                  "\t"<<std::setw(11)<<std::setiosflags (std::ios::left)<<ns.Ome_p*Ome<<   //  2 Ome
                  "\t"<<std::setw(11)<<std::setiosflags (std::ios::left)<<res[0]<<           //  3 |q0|²
                  "\t"<<std::setw(11)<<std::setiosflags (std::ios::left)<<res[1]<<           //  4 |q1|²
                  "\t"<<std::setw(11)<<std::setiosflags (std::ios::left)<<res[2]<<           //  5 |q2|²
                  "\t"<<std::setw(11)<<std::setiosflags (std::ios::left)<<res[6]<<           //  6 N
                    std::endl;
            }        
        if (q0s_rgth.size()==0){
            q0s_rgth.push_back(std::make_pair(ns.omeB, 0));
            p1s_rgth.push_back(std::make_pair(ns.omeB, 0));
            p3s_rgth.push_back(std::make_pair(ns.omeB, 0));
            }
        if (q0s_left.size()==0){
            q0s_left.push_back(std::make_pair(ns.omeB, 0));
            p1s_left.push_back(std::make_pair(ns.omeB, 0));
            p3s_left.push_back(std::make_pair(ns.omeB, 0));
            }
        int total=q0s_rgth.size()+ q0s_left.size();
        std::vector<std::pair<double,double>> q0s_shft(total);
        std::vector<std::pair<double,double>> p1s_shft(total);
        std::vector<std::pair<double,double>> p3s_shft(total);
        
        sort(q0s_left.begin(), q0s_left.end());
        sort(q0s_rgth.begin(), q0s_rgth.end());
        sort(p1s_left.begin(), p1s_left.end());
        sort(p1s_rgth.begin(), p1s_rgth.end());
        sort(p3s_left.begin(), p3s_left.end());
        sort(p3s_rgth.begin(), p3s_rgth.end());


        bool left_allZeros = std::all_of(q0s_left.begin(), q0s_left.end(), [](const auto& pair) {
            return pair.second == 0.0;
            });

        bool rgth_allZeros = std::all_of(q0s_rgth.begin(), q0s_rgth.end(), [](const auto& pair) {
            return pair.second == 0.0;
            });
        
        if (!left_allZeros && !rgth_allZeros){
            merge(q0s_left.begin(), q0s_left.end(), q0s_rgth.begin(), q0s_rgth.end(), q0s_shft.begin()); 
            merge(p1s_left.begin(), p1s_left.end(), p1s_rgth.begin(), p1s_rgth.end(), p1s_shft.begin());
            merge(p3s_left.begin(), p3s_left.end(), p3s_rgth.begin(), p3s_rgth.end(), p3s_shft.begin());
            }
        if (left_allZeros){
            q0s_shft=q0s_rgth;
            p1s_shft=p1s_rgth;
            p3s_shft=p3s_rgth;
            }
        if (rgth_allZeros){
            q0s_shft=q0s_left;
            p1s_shft=p1s_left;
            p3s_shft=p3s_left;
            }
            
        complete(q0s_rgth, omemi, omema, dome);
        complete(q0s_left, omemi, omema, dome);
        complete(p1s_rgth, omemi, omema, dome);
        complete(p1s_left, omemi, omema, dome);
        complete(p3s_rgth, omemi, omema, dome);
        complete(p3s_left, omemi, omema, dome);        
        complete(q0s_shft, omemi, omema, dome); 
        complete(p1s_shft, omemi, omema, dome); 
        complete(p3s_shft, omemi, omema, dome); 

        std::fstream wleft, wrght;
        wleft.open("../data/output/left.dat", std::ios::out);
        wrght.open("../data/output/right.dat", std::ios::out);
        
        for (int i=0; i <  int(q0s_left.size()); i++)
            wleft<<std::setw(20)<<std::setprecision(15)<<q0s_left[i].first<<
                                  "\t\t"<<std::setw(11)<<q0s_left[i].second<<
                                  "\t\t"<<std::setw(11)<<p1s_left[i].second<<
                                  "\t\t"<<std::setw(11)<<p3s_left[i].second<<std::endl;
        for (int i=0; i <  int(q0s_rgth.size()); i++)
            wrght<<std::setw(20)<<std::setprecision(15)<<q0s_rgth[i].first<<
                                  "\t\t"<<std::setw(11)<<q0s_rgth[i].second<<
                                  "\t\t"<<std::setw(11)<<p1s_rgth[i].second<<
                                  "\t\t"<<std::setw(11)<<p3s_rgth[i].second<<std::endl;

        for (int i=0; i <  int(q0s_shft.size()); i++)
            shft<<std::setw(20)<<std::setprecision(15)<<q0s_shft[i].first<<
                                 "\t\t"<<std::setw(11)<<q0s_shft[i].second<<
                                 "\t\t"<<std::setw(11)<<p1s_shft[i].second<<
                                 "\t\t"<<std::setw(11)<<p3s_shft[i].second<<std::endl;
        return p3s_shft;
        }
