/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2016 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "colvar/Colvar.h"
#include "colvar/ActionRegister.h"
#include "tools/Communicator.h"
#include "tools/Tools.h"
#include "tools/Angle.h"
#include "tools/IFile.h"

#include <string>
#include <math.h>

using namespace std;

namespace PLMD{
namespace colvar{

//+PLUMEDOC COLVAR MANY_ANGLE
/*

\plumedfile
LOAD FILE=ManyAngle2.cpp

# Define groups for the CV
INCLUDE FILE=centers.dat
C: GROUP ATOMS=1-864:8
O: GROUP ATOMS=2-864:8
N1: GROUP ATOMS=3-864:8
N2: GROUP ATOMS=6-864:8

MANY_ANGLE ...
 LABEL=ma
 CENTER=C
 START=N1
 END=N2
 RCUT=0.7
 UP_DOWN_SYMMETRY
... MANY_ANGLE

PRINT STRIDE=1  ARG=* FILE=COLVAR_MA
\endplumedfile

*/
//+ENDPLUMEDOC

class ManyAngle : public Colvar {
bool pbc, serial;
  vector<AtomNumber> center_lista, start_lista, end_lista;
  std::vector<PLMD::AtomNumber> atomsToRequest;
  double rcut, rcut2;
  // Up-down symmetry
  bool doUpDownSymmetry;
  // Low communication variant
  bool doLowComm;

public:
  explicit ManyAngle(const ActionOptions&);
  virtual void calculate();
  static void registerKeywords( Keywords& keys );
};

PLUMED_REGISTER_ACTION(ManyAngle,"MANY_ANGLE")

void ManyAngle::registerKeywords( Keywords& keys ){
  Colvar::registerKeywords(keys);
  keys.addFlag("SERIAL",false,"Perform the calculation in serial - for debug purpose");
  keys.add("atoms","CENTER","Center atoms");
  keys.add("atoms","START","Start point of vector defining orientation");
  keys.add("atoms","END","End point of vector defining orientation");
  keys.add("compulsory","RCUT","1","Maximum distance for being neighbor");
  keys.addFlag("UP_DOWN_SYMMETRY",false,"The symmetry is such that parallel and antiparallel vectors are not distinguished. The angle goes from 0 to pi/2 instead of from 0 to pi.");
  keys.addFlag("LOW_COMM",false,"Use an algorithm with less communication between processors");
}

ManyAngle::ManyAngle(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao),
pbc(true),
serial(false)
{

  parseFlag("SERIAL",serial);

  parseAtomList("CENTER",center_lista);
  parseAtomList("START",start_lista);
  parseAtomList("END",end_lista);
  if(center_lista.size()!=start_lista.size()) error("Number of atoms in START must be equal to the number of atoms in CENTER");
  if(center_lista.size()!=end_lista.size()) error("Number of atoms in END must be equal to the number of atoms in CENTER");
  
  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;

  if(pbc) log.printf("  using periodic boundary conditions\n");
  else    log.printf("  without periodic boundary conditions\n");

  addValueWithDerivatives(); setNotPeriodic();

  parse("RCUT",rcut);
  log.printf("  The neighbor is defined within 0 and %f \n", rcut);
  rcut2 = rcut*rcut;

  doUpDownSymmetry=false;
  parseFlag("UP_DOWN_SYMMETRY",doUpDownSymmetry);
  if (doUpDownSymmetry) log.printf("  The angle can take values between 0 and pi/2 due to the up down symmetry. \n");

  doLowComm=false;
  parseFlag("LOW_COMM",doLowComm);
  if (doLowComm) {
     log.printf("  Using the low communication variant of the algorithm");
  }

  checkRead();

  atomsToRequest.reserve ( center_lista.size() + start_lista.size() + end_lista.size() );
  atomsToRequest.insert (atomsToRequest.end(), center_lista.begin(), center_lista.end() );
  atomsToRequest.insert (atomsToRequest.end(), start_lista.begin(), start_lista.end() );
  atomsToRequest.insert (atomsToRequest.end(), end_lista.begin(), end_lista.end() );
  requestAtoms(atomsToRequest);

}

// calculator
void ManyAngle::calculate()
{
  double cv=0.0; // (T)
  unsigned nat=getNumberOfAtoms();
  vector<Vector> deriv(nat); // This will initialized vector with 0 modulo. (T)
  unsigned int num_cv=0; // (T)
  
  if(pbc) makeWhole();
  // Setup parallelization
  unsigned stride=comm.Get_size();
  unsigned rank=comm.Get_rank();
  if(serial){
    stride=1;
    rank=0;
  } else {
    stride=comm.Get_size();
    rank=comm.Get_rank(); // Somehow it is always rank=0. Does MPI really work? (T)
  }
  if (doLowComm) {
    for(unsigned int i=rank;i<center_lista.size();i+=stride) {
      unsigned atom1_mol1=i+center_lista.size(); // This start from index=i+center_lista, so it will actually read start_lista. (T)
      unsigned atom2_mol1=i+center_lista.size()+start_lista.size(); // This will actually read end_lista. (T)
      Vector dik=delta(getPosition(atom2_mol1),getPosition(atom1_mol1));
      for(unsigned int j=0;j<center_lista.size();j+=1) {
        double d2;
        Vector distance;
        if(getAbsoluteIndex(i)==getAbsoluteIndex(j)) continue;        
        distance=delta(getPosition(i),getPosition(j));
        if ( (d2=distance[0]*distance[0])<rcut2 && (d2+=distance[1]*distance[1])<rcut2 && (d2+=distance[2]*distance[2])<rcut2) {
          unsigned atom1_mol2=j+center_lista.size();
          unsigned atom2_mol2=j+center_lista.size()+start_lista.size();
          Vector dij=delta(getPosition(atom1_mol2),getPosition(atom2_mol2));

          Vector ddij,ddik; // (T)
          PLMD::Angle a; // (T)
          double angle=a.compute(dij,dik,ddij,ddik); // (T)
          
          double xi=tanh(20.0*(angle-1.85)); // (T)
          double pa=pi-2.0*angle; // (T)
          cv += 0.5*(pa*xi + pi); // (T)
          double dcv_dangle=10.0*pa*(1.0-xi*xi)-xi; // (T)
          
          deriv[atom1_mol1]+=dcv_dangle*ddik; // (T)
          deriv[atom2_mol1]+=-dcv_dangle*ddik; // (T)
          deriv[atom1_mol2]+=-dcv_dangle*ddij; // (T)
          deriv[atom2_mol2]+=dcv_dangle*ddij; // (T)

          num_cv++; // (T)
        }
      }
    }
  } else {
    for(unsigned int i=rank;i<(center_lista.size()-1);i+=stride) {
      unsigned atom1_mol1=i+center_lista.size();
      unsigned atom2_mol1=i+center_lista.size()+start_lista.size();
      Vector dik=delta(getPosition(atom2_mol1),getPosition(atom1_mol1));
      for(unsigned int j=i+1;j<center_lista.size();j+=1) {
        Vector distance;
        if(getAbsoluteIndex(i)==getAbsoluteIndex(j)) continue; // Only consider i != j. (T)        
        distance=delta(getPosition(i),getPosition(j));
        double d2;
        if ( (d2=distance[0]*distance[0])<rcut2 && (d2+=distance[1]*distance[1])<rcut2 && (d2+=distance[2]*distance[2])<rcut2) {
          unsigned atom1_mol2=j+center_lista.size();
          unsigned atom2_mol2=j+center_lista.size()+start_lista.size();
          Vector dij=delta(getPosition(atom1_mol2),getPosition(atom2_mol2));

          Vector ddij,ddik; // (T)
          PLMD::Angle a; // (T)
          double angle=a.compute(dij,dik,ddij,ddik); // (T)
          
          double xi=tanh(20.0*(angle-1.85)); // (T)
          double pa=pi-2.0*angle; // (T)
          cv += 0.5*(pa*xi + pi); // (T)
          double dcv_dangle=10.0*pa*(1.0-xi*xi)-xi; // (T)
          
          // Perform chain rule. d(cv)/d(angle)*d(angle)/d(dik)
          deriv[atom1_mol1]+=dcv_dangle*ddik; // (T)
          deriv[atom2_mol1]+=-dcv_dangle*ddik; // (T)
          deriv[atom1_mol2]+=-dcv_dangle*ddij; // (T)
          deriv[atom2_mol2]+=dcv_dangle*ddij; // (T)
          num_cv++; // (T)

          }
        }
      }
    }

  if(!serial){
    comm.Sum(cv);
    comm.Sum(num_cv);
    for(unsigned i=0;i<nat;++i) comm.Sum(deriv[i]);
  }
  
  // Assign output quantities
  Tensor virial; // (T)
  for(unsigned i=0;i<nat;++i){ 
    // If atom 3 is involved in calculation of 2 angles when there are 3 angles being actually calculated, 
    // This derivatives of atom 3 will be the sum of derivatives of atom 3 calculated in each angle divided
    // by 3 (angles).
    setAtomsDerivatives(i, deriv[i]/num_cv); // (T)
    // This line is calculated in the same way as what has been defined in Colvar::setBoxDerivativesNoPbc(Value* v).
    virial-=Tensor(getPosition(i), deriv[i]/num_cv); // (T)
  }
  setValue( cv/num_cv ); // (T)
  // If double dcv_dangle, the box derivatives will double. When there are [num_cv] angles being calculated,
  // the box derivatives will be the sum of box derivatives of each individual angles.
  setBoxDerivatives(virial); // (T)

}


}
}