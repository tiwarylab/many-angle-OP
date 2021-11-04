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
LOAD FILE=ManyAngle.cpp

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
  double cv=0.0;                         // (T)
  unsigned nat=getNumberOfAtoms();
  vector<Vector> deriv(nat);             // This will initialized vector with 0 modulo. (T)
  unsigned int num_cv=0;                 // (T)

  double angle(0.0);                     // (T)
  double xi(0.0);                        // (T)
  double pa(0.0);                        // (T)
  double dcv_dangle(0.0);                // (T)
  
  if(pbc) makeWhole();
  // Setup parallelization
  unsigned stride=comm.Get_size();
  unsigned rank=comm.Get_rank();

  if(serial){
    stride=1;
    rank=0;
  } else {
    stride=comm.Get_size();
    rank=comm.Get_rank();
  }

  for(unsigned int i=rank;i<(center_lista.size()-1);i+=stride) {
    unsigned atom1_mol1=i+center_lista.size();
    unsigned atom2_mol1=i+center_lista.size()+start_lista.size();
    Vector dik=pbcDistance(getPosition(atom2_mol1),getPosition(atom1_mol1));
    for(unsigned int j=i+1;j<center_lista.size();j+=1) {
      if(getAbsoluteIndex(i)==getAbsoluteIndex(j)) continue;
      Vector distance;
      if(pbc){
        distance=pbcDistance(getPosition(i),getPosition(j));
      } else {
        distance=delta(getPosition(i),getPosition(j));
      }

      double d2=0.0; // (T)
      if ( (d2=distance[0]*distance[0])<rcut2 && (d2+=distance[1]*distance[1])<rcut2 && (d2+=distance[2]*distance[2])<rcut2) {
        unsigned atom1_mol2=j+center_lista.size();
        unsigned atom2_mol2=j+center_lista.size()+start_lista.size();

        Vector dij=pbcDistance(getPosition(atom1_mol2),getPosition(atom2_mol2));

        Vector ddij,ddik;                            // (T)
        PLMD::Angle a;                               // (T)
        angle=a.compute(dij,dik,ddij,ddik);          // (T)
        
        xi=tanh(20.0*(angle-1.85));                  // (T)
        pa=pi-2.0*angle;                             // (T)
        cv += 0.5*(pa*xi + pi);                      // (T)
        dcv_dangle=10.0*pa*(1.0-xi*xi)-xi;           // (T)
        
        // Chain rule. d(cv)/d(angle)*d(angle)/d(dik)
        deriv[atom1_mol1]+=dcv_dangle*ddik;          // (T)
        deriv[atom2_mol1]+=-dcv_dangle*ddik;         // (T)
        deriv[atom1_mol2]+=-dcv_dangle*ddij;         // (T)
        deriv[atom2_mol2]+=dcv_dangle*ddij;          // (T)
        num_cv++;                                    // (T)

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
    setAtomsDerivatives(i, deriv[i]/num_cv);         // (T)
    // This line is calculated in the same way as what has been defined in Colvar::setBoxDerivativesNoPbc(Value* v).
    virial-=Tensor(getPosition(i), deriv[i]/num_cv); // (T)
  }
  setValue( cv/num_cv );                             // (T)
  // If double dcv_dangle, the box derivatives will double. When there are [num_cv] angles being calculated,
  // the box derivatives will be the sum of box derivatives of each individual angles.
  setBoxDerivatives(virial);                         // (T)

}


}
}