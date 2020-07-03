//
//  QFT-DDgenerator.cpp
//  QDAG
//
//  Created by Hassan Mallahzadeh on 2020-07-02.
//  Copyright © 2020 Hassan Mallahzadeh. All rights reserved.
//
#include "DDpackage.h"
#include "DDcomplex.h"
#include "QFT-DDgenerator.hpp"
using CN = dd::ComplexNumbers;
StateGenerator::StateGenerator() {
    dd = new dd::Package;
}


dd::Edge StateGenerator::dd_UniformState(int n) {
    dd::Edge e_state = dd->makeBasisState(n, 0);
    for (int i = 1; i < pow(2,n); ++i){
         e_state = dd->add(e_state, dd->makeBasisState(n, i));
    }
    dd::ComplexValue c{1/std::sqrt(pow(2,n)), 0.0 };
         dd::Complex cx = dd->cn.getTempCachedComplex(c.r,c.i);
             e_state.w =dd->cn.mulCached(e_state.w, cx);
   // e_state = dd->normalize(e_state, true);
   return e_state;
}

dd::Edge StateGenerator::dd_Sqrt3State(int n){//TODO: fix: function breaks after 9 steps (return negative for integer).
    fp sqrt3 = sqrt(3);
        dd::Edge e_state = dd->makeBasisState(n, 0);//3 gives 1.
    sqrt3 = sqrt3 * 10;
      for (int i = 1; i < pow(2,n); ++i){
          if((int)sqrt3 % 2)//bit is 1, so add
           e_state = dd->add(e_state, dd->makeBasisState(n, i));
           sqrt3 = sqrt3 * 10;
      }
      dd::ComplexValue c{1/std::sqrt(pow(2,n)), 0.0 };
           dd::Complex cx = dd->cn.getTempCachedComplex(c.r,c.i);
               e_state.w =dd->cn.mulCached(e_state.w, cx);
     return e_state;
}
dd::Edge StateGenerator::dd_RandomState(int n){
    dd::Edge e_state = dd->makeBasisState(n, 0);//assume first one is 1 (TODO: fix)
    srand (time(NULL));
    for (int i = 1; i < pow(2,n); ++i){
             if(rand() % 2)//bit is 1, so add
                 e_state = dd->add(e_state, dd->makeBasisState(n, i));
         }
    dd::ComplexValue c{1/std::sqrt(pow(2,n)), 0.0 };
              dd::Complex cx = dd->cn.getTempCachedComplex(c.r,c.i);
                  e_state.w =dd->cn.mulCached(e_state.w, cx);
  //  e_state = dd->normalize(e_state, true);
        return e_state;
}


