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

#include<iostream>
#include <stdio.h>
#include <stdlib.h>

#define c  299792458
/** c [m/s] Speed of light **/
#define h  6.626068e-34
/** h [m²kg/s] **/
#define j2eV 6.24150636309e18

/** Compila con
g++ -Iinclude ../src/eV2nm.cxx -o ../bin/eV2nm
**/
using namespace std;

int main(int argc, char** argv){
  double lam, ome;
   if (argc<2){
	cout<<endl<<"  Usage: "<<argv[0]<<" <wavelenght in eV>"<<endl<<endl;
	exit(0);
	}
  ome=atof(argv[1]);
  lam=h*j2eV*c/(ome*1.e-9);

  cout<<lam<<endl;
  return 0;
  } 
