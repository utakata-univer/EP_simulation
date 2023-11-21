/*			    -*- Mode: C -*-
// File: make_su2_mat.c
// Created: Thu Apr  6 1995
// Author: J. E. Hetrick <hetrick@corelli.physics.Arizona.EDU>
//
// Description: something like a constructor for su2_matrix
// Usage: su2_matrix make_su2_matrix(float a0, float a1, float a2, float a3);
//        returns the su2_matrix made from a0..a3
// $Id$
// $Log$
*/
#include "complex.h"
#include "globaldefs.h"
#include "su2.h"

su2_matrix make_su2_matrix(float a0, float a1, float a2, float a3) {
   su2_matrix a;

   a.e[0] = a0;
   a.e[1] = a1;
   a.e[2] = a2;
   a.e[3] = a3;

   return a;
}

