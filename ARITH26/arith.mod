#---
#--- Modele optimisation des largeurs: arith-26
#--- Auteur: Hacene Ouzia, Sorbonne Universite (01-2019) 
#---

#---
#--- DIMENSIONS 
param p integer > 0;    #-- nbr de contraintes (p)
param np integer > 0;   #-- nbr de variables (n+p)

#---
#--- ENSEMBLES 
set I := 1 .. p;        #-- Indices contraintes
set J := 1 .. np;       #-- Indices variables

#---
#--- PARAMETRES 
param u{J} integer > 0;     #-- Valeur max var w

param E{I,J};                #-- Coefficients des contraintes ....
param eps{I};                #-- Second membre ...

param wtilde{J} integer > 0;
   #--
   #-- Il faut que 3 <= wtilde(j) < u(j), pr tt j ...
   check{j in J}:
      3 <= wtilde[j] and wtilde[j] < u[j];


param bsup{j in J} integer default u[j] - wtilde[j];    #-- val max de w_j - wtilde_j
param binf{j in J} integer default 2 - wtilde[j];       #-- val min de w_j - wtilde_j


#---
#--- VARIABLES 
var w{j in J} integer >= 2 <= u[j];     #-- word-length
var delta{J} binary;                    #-- delta


#--- 
#--- OBJECTIF 
minimize Longueur:     #-- Somme des largeurs ...
   sum{j in J} x[j]; 


#---
#--- CONTRAINTES 
subject to output_error{i in I}:        # output error
   sum{j in J} E[i,j]*2**(-w[j] + delta[j]) <= eps[i];

subject to deltaa{j in J}:
   x[j] - wtilde[j] <= (1 - delta[j])*bsup[j];

subject to Kondb{i in I, j in J}:
   x[j] - wtilde[j] >= (binf[j] - 1)*delta[j];
 
