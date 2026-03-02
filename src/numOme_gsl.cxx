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

#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <math.h>
#include <vector>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>
#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])

/** compila con 
 g++ -Iinclude ../src/numOme_gsl.cxx -o ../bin/nom -lgsl -lgslcblas -lm
**/ 

using namespace std;
double pi=acos(-1.);

double find_Ome_fourier(int nfft, vector<complex<double>> wave, double dt){
    double Ome, rOme, iOme, iMax=0, rMax=0;
    double data[2*nfft];
    gsl_fft_complex_wavetable * wavetable;
    gsl_fft_complex_workspace * workspace;

    wavetable = gsl_fft_complex_wavetable_alloc (nfft);
    workspace = gsl_fft_complex_workspace_alloc (nfft);
    ofstream  wve("../data/output/wave.dat");
    if (!wve) {
        cerr << "Error: Cannot open output file" << endl;
        return 1;
    }
    
    for (int i=0; i<nfft; i++) {
        REAL(data,i) = wave[i].real();
        IMAG(data,i) = wave[i].imag();
        wve<<REAL(data,i)<<" "<<IMAG(data,i)<<endl;
        }
    gsl_fft_complex_forward(data, 1, nfft, wavetable, workspace);

    for (int i=0; i<nfft; i++) {
        if (rMax<fabs(REAL(data,i))){
            rMax=fabs(REAL(data,i));
            rOme=i;
            }
        if (iMax<fabs(IMAG(data,i))){
            iMax=fabs(IMAG(data,i));
            iOme=i;
            }
        }
    cout<<rOme<<" "<<iOme<<endl;
    Ome=pi*(rOme+iOme)/(nfft*dt);
    gsl_fft_complex_wavetable_free (wavetable);
    gsl_fft_complex_workspace_free (workspace);
    return Ome;
    }

int main(){
    double Ome, nOme, t, dt=0.1;//, t0;
    Ome=3.9;
    int N=1100, pip=0;//, nmin, nmax;
    complex<double> img;
    img = complex<double>  (0., 1.);
    complex<double> value;
    vector<complex<double>> wave;
    for (int i=0; i<N; i++){
        t=i*dt;
        value=cos(Ome*t)+img*(sin(Ome*t));
        if (i>100 && i<=600){
            if (pip==0){
                // t0=t;
                pip=1;
                }
            wave.push_back(value);
            }
        }
    nOme=find_Ome_fourier(wave.size(), wave, dt);
    cout<<nOme<<" "<<Ome<<" err: "<<fabs(nOme-Ome)<<endl;
    return 0;
    }
