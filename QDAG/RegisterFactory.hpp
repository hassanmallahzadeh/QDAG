//
//  RegisterFactoryInterface.hpp
//  QDAG
//
//  Created by Hassan Mallahzadeh on 2020-08-21.
//  Copyright © 2020 Hassan Mallahzadeh. All rights reserved.
//

#ifndef RegisterFactoryInterface_hpp
#define RegisterFactoryInterface_hpp
#include "commonheaders.h"
#include "IIC-JKU/DDpackage.h"
#include "QFT-DDgenerator.hpp"
#include <stdio.h>
class RegisterFactory{
private:
    lli N = -1;
    lli a = -1;
    int n_in = -1;// input register size
    int n_out = -1;// output register size
    int n_car = -1;// auxiliary register (for quantum 'carry' values) size
    int n_addmod = -1;// auxiliary register (for classical 'fixed' value N) size in Adder modulo N.
    int n_mul = -1;// auxiliary register for modular multiplication
    int n_mod = -1;// auxiliary register (for quantum temp values) size in modular exponentiation.
    int n_tot = -1;// total number of qubits (input + output + auxiliary)
    dd::Edge state = dd::Package::DDnull;
    dd::Package *dd = nullptr;
    short* line = nullptr;//holds 'line'for setting basic gates.
    GateGenerator gg;
    vector<int> v_car;
    vector<int> v_mul_temp;
    vector<int> v_out;
    vector<int> v_in;
    vector<int> v_addmod_temp;
    vector<int> v_mod_temp;

    void DetermineRegSizes();
    void AllocateRegVectors();
    void MakeInitialState();
    void Carry(int c0,int a0, int b0, int c1);
    void CarryInv(int c0,int a0, int b0, int c1);
    void Sum(int c0,int a0, int b0);
    void RippleAdder(vector<int>, vector<int>);
    void RippleAdderInv(vector<int>, vector<int>);
    void ModuloNAdder(vector<int>, vector<int>);
    void ModuloNAdderInv(vector<int>, vector<int>);
    void TempNResetter();
    void extractedMapToTemp(vector<int>, int c1, int c2, lli cl);
    void CMulModN(int c, lli cfa);
    void CMulModNInv(int c, lli cfa);
public:
    RegisterFactory(lli N, lli a, dd::Package *dd);
    dd::Edge RippleAdderDebug(vector<int> num1, vector<int> num2);
    dd::Edge RippleAdder(vector<dd::ComplexValue> num1, vector<dd::ComplexValue> num2, int n);
    ~RegisterFactory();
    dd::Edge CExponentiation();
};
#endif /* RegisterFactoryInterface_hpp */
