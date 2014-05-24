/*
 *    This file is part of ACADO Toolkit.
 *
 *    ACADO Toolkit -- A Toolkit for Automatic Control and Dynamic Optimization.
 *    Copyright (C) 2008-2009 by Boris Houska and Hans Joachim Ferreau, K.U.Leuven.
 *    Developed within the Optimization in Engineering Center (OPTEC) under
 *    supervision of Moritz Diehl. All rights reserved.
 *
 *    ACADO Toolkit is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    ACADO Toolkit is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with ACADO Toolkit; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */


 /**
 *    \file interfaces/matlab/integrator/models/cstr.cpp
 *    \author Niels Haverbeke, Boris Houska, Hans Joachim Ferreau
 *    \date 2009
 *
 */

void tanks( DifferentialEquation *f ){

    // Define a Right-Hand-Side:
    // -------------------------

    
    const double a1=0.03;
    const double a2=0.02;

    const double A1=8.4;
    const double A2=7.1;
    const double A_u=1;

    const double g = 9.81;
    const double k = 0.4;

  DifferentialState cost;
  DifferentialState z1;
  DifferentialState z2;


  Control u;
    
	*f << dot(cost) ==  (pow((z1-2),2)+pow(u,2))/1000;
    *f << dot(z1) ==  -(a1/A1)*sqrt(2*g*z1)+(a2/A2)*sqrt(2*g*z2);
    *f << dot(z2) ==  -(a2/A2)*sqrt(2*g*z2)+(k/A_u)*u;
    
    
    
    
}


