#ifndef UTILITY_FUNCTIONS_H_INCLUDED
#define UTILITY_FUNCTIONS_H_INCLUDED

#include "common.h"
#include "Rcpp.h"

Vertex nextNode(VisitMap&,graph&,Vertex);
void addOccs(MapOccurrences& ,Vertext,occ);

#endif // UTILITY_FUNCTIONS_H_INCLUDED
