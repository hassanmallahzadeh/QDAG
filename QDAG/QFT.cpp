//
//  QFT.cpp
//  cmakelearn
//
//  Created by Hassan Mallahzadeh on 2020-06-09.
//  Copyright © 2020 Hassan Mallahzadeh. All rights reserved.
//
#include "util.h"
#include "DDcomplex.h"
#include "QFT.hpp"
//test3
//#include "engine.h"
using namespace std::chrono;
using std::cout;
using std::endl;
/// Smatv1: Swap matrix construct using Pauli X, Y, Z, I matrices.
/// @param dd source package instance
/// @param n number of bits
/// @param b1 first bit in swap
/// @param b2 second bit in swap
dd::Edge Smatv1(dd::Package *dd, int n, int b1, int b2){
    short *line = new short[n]{};//HM: set 'line' for dd->makeGateDD
    for (int i = 0; i < n; ++i){
        line[i] = -1;
    }
    //a test
    line[b1] = 2;
    dd::Edge x1_gate = dd->makeGateDD(Xmat, n, line);
    dd::Edge y1_gate = dd->makeGateDD(Ymat, n, line);
    dd::Edge z1_gate = dd->makeGateDD(Zmat, n, line);
    line[b1] = -1;
    line[b2] = 2;
    dd::Edge x2_gate = dd->makeGateDD(Xmat, n, line);
    dd::Edge y2_gate = dd->makeGateDD(Ymat, n, line);
    dd::Edge z2_gate = dd->makeGateDD(Zmat, n, line);
    dd::Edge iMatnxn = dd->makeIdent(0, n-1);
    dd::Edge x_gate = dd->multiply(x1_gate, x2_gate);//make 2 qubit dd
    dd::Edge y_gate = dd->multiply(y1_gate, y2_gate);//make 2 qubit dd
    dd::Edge z_gate = dd->multiply(z1_gate, z2_gate);//make 2 qubit dd
    dd::Edge temp1 = dd->add(iMatnxn, x_gate);
    dd::Edge temp2 = dd->add(y_gate, z_gate);
    dd::Edge ret = dd->add(temp1, temp2);
    // ret = ret * 2; HM: code has to be added for this.
    delete[] line;
    return ret;
}
/// Smatv2: Swap matrix construct using C-Not gates.
/// @param dd source package instance
/// @param n number of bits
/// @param b1 first bit in swap
/// @param b2 second bit in swap
dd::Edge Smatv2(dd::Package *dd, int n, int b1, int b2){
    short *line = new short[n]{};//HM: set 'line' for dd->makeGateDD
    for (int i = 0; i < n; ++i){
        line[i] = -1;
    }
    line[b1] = 2;//HM: set target and control
    line[b2] = 1;
    dd::Edge c12_gate = dd->makeGateDD(Xmat, n, line);
    line[b1] = 1;
    line[b2] = 2;
    dd::Edge c21_gate = dd->makeGateDD(Xmat, n, line);
    dd::Edge temp_gate = dd->multiply(c21_gate, c12_gate);
    dd::Edge s_gate = dd->multiply(temp_gate, c21_gate);
    delete[] line;
    return s_gate;
}
///// Smatv3: Swap matrix construct using matrix construction.
///// @param dd source package instance
///// @param n number of bits
///// @param b1 first bit in swap
///// @param b2 second bit in swap
//dd::Edge Smatv3(dd::Package *dd, int n, int b1, int b2){
//    int ** mat = new int* [n]{};
//    for (int i = 0; i < n; ++i){
//        mat[i] = new int[n]{};
//    }
//    for (int i = 0; i < n; ++i)
//        for (int j = 0; j < n ; ++j){
//            if(i == j - n + 1)
//            mat[i][j] = 1;
//            else
//                mat[i][j] = 0;//HM: redunant I believe.
//        }
//    short *line = new short[n]{};//HM: set 'line' for dd->makeGateDD
//    for (int i = 0; i < n; ++i){
//        line[i] = -1;
//    }
//    line[b1] = 2;//HM: set target and control
//    line[b2] = 1;
//    dd::Edge c12_gate = dd->makeGateDD(Xmat, n, line);
//    line[b1] = 1;
//    line[b2] = 2;
//    dd::Edge c21_gate = dd->makeGateDD(Xmat, n, line);
//    dd::Edge temp_gate = dd->multiply(c12_gate, c21_gate);
//    dd::Edge s_gate = dd->multiply(temp_gate, c12_gate);
//    delete[] line;
//    return s_gate;
//}
////rotation gate
void RmatGenerator(dd::Matrix2x2 &m, int k){
    m[0][0] = { 1, 0 };
    m[0][1] = { 0, 0 };
    m[1][0] = { 0, 0 };
    fp angle = 2 * dd::PI/pow(2,k);
    m[1][1] = { cos(angle), sin(angle) };
}
void lineSet(short* line, int t,int c){
    if(t >= 0)
        line[t] = 2;
    if(c >= 0)
        line[c] = 1;
}
void lineReset(short* line, int t,int c){
    if(t >= 0)
        line[t] = -1;
    if(c >= 0)
        line[c] = -1;
}
dd::Edge permuteOperator(dd::Package *dd, int n){
    dd::Edge e_swap = dd->makeIdent(0, n-1);
    for(int i = 0; i < n/2; ++i){//apply n/2 swap gates.
        e_swap = dd->multiply(Smatv2(dd, n, i, n - i - 1), e_swap);
    }
    return e_swap;
}
dd::Edge dd_QFT(dd::Package *dd, int n, bool bApplySwap){
    // dd->reverseVariableOrder(n);
    dd::Edge e_ans = dd->makeIdent(0, n-1);//HM: this can be taken out, here for clarity
    
    short *line = new short[n]{};//HM: set 'line' for dd->makeGateDD
    for (int i = 0; i < n; ++i){
        line[i] = -1;
    }
    std::ofstream datafile;
    datafile.open ("timers.txt");
    time_point<system_clock> start, end;
    for (int i = 0; i < n; ++i){//HM: for each of n qubits do:
        start = system_clock::now();
        lineSet(line,i,-1);
        dd::Edge e_line = dd->makeGateDD(Hmat, n, line);//HM: start by hadamard as in circuit
        dd::Matrix2x2 Rmat;//HM: put rotation gates in place
        int j = 0;
        for (; j < n - i - 1; ++j){
            //  dd->export2Dot(e_line, "qft-l-"+std::to_string(i)+"-"+std::to_string(j)+".dot",false, true);
            lineSet(line, -1, i + j + 1);
            RmatGenerator(Rmat, j+2);
            //     cout<<"gate counter: "<<i<<", "<<j<<", "<<dd->size(e_line)<< endl;
            //    dd::Edge e_temp = dd->makeGateDD(Rmat, n, line);
            //     dd->export2Dot(e_temp, "qft-r-"+std::to_string(i)+"-"+std::to_string(j)+".dot",false, true);
            e_line = dd->multiply(dd->makeGateDD(Rmat, n, line), e_line);
            lineReset(line, -1, i + j + 1);
        }
        // cout<<"bit: "<<i<<", gate: "<< j << ", size: "<<dd->size(e_line) << endl;
        //   dd->export2Dot(e_line, "qft-l-"+std::to_string(i)+".dot",false, true);
        e_ans = dd->multiply(e_line, e_ans);//HM: apply current line to final operator
        //    dd->export2Dot(e_ans, "qft-a-"+std::to_string(i)+".dot",false, true);
        lineReset(line,i,n - i - 1);//HM: reset line array
        end = system_clock::now();
        duration<float> elapsed_seconds = end - start;
        datafile << i << " "<< elapsed_seconds.count() <<endl;
    }
    if(bApplySwap){
        e_ans = dd->multiply(e_ans, permuteOperator(dd,n));
        //   e_ans = dd->multiply(permuteOperator(dd,n),e_ans);//HM: temp, experimental, effect of reverse variable order.
    }
    //   dd->export2Dot(e_ans, "qft-f.dot",false, true);
    datafile.close();
    delete[] line;
    return e_ans;
}
dd::Edge dd_QFTV2(dd::Package *dd, int n, bool bApplySwap, dd::Edge state){//this one applies gates to state one by one, instead of making the overal gate then applying on state.
    // dd->reverseVariableOrder(n);
    dd::Edge e_ans = state;
    if(bApplySwap){
        e_ans = dd->multiply(permuteOperator(dd,n),e_ans);
    }
    short *line = new short[n]{};//HM: set 'line' for dd->makeGateDD
    for (int i = 0; i < n; ++i){
        line[i] = -1;
    }
    for (int i = 0; i < n; ++i){//HM: for each of n qubits do:
        lineSet(line,i,-1);
        dd::Edge e_line = dd->makeGateDD(Hmat, n, line);//HM: start by hadamard as in circuit
        e_ans = dd->multiply(e_line ,e_ans);
        dd::Matrix2x2 Rmat;//HM: put rotation gates in place
        for (int j = 0; j < n - i - 1; ++j){
            lineSet(line, -1, i + j + 1);
            RmatGenerator(Rmat, j+2);
            e_ans = dd->multiply(dd->makeGateDD(Rmat, n, line), e_ans);
            lineReset(line, -1, i + j + 1);
        }
        lineReset(line,i,n - i - 1);//HM: reset line array
    }
//    e_ans = dd->multiply(permuteOperator(dd,n),e_ans);//HM: temp, experimental, effect of variablre reverese.
    delete[] line;
    return e_ans;
}

void graphfunctuion(){
    
}
/// QFT in main. to be moved to function, class.
int main(){
    
    
    if(false){
        auto* dd = new dd::Package;
        int n = 5;
        
        //
        //        //    dd::Edge state3 = dd->makeBasisState(2, 3);
        //         dd::Edge state0 = dd->makeBasisState(n, 0);
        dd::Edge state1 = dd->makeBasisState(n, 1);
        //        dd::Edge state2 = dd->makeBasisState(n, 2);
        //        dd::Edge state3 = dd->makeBasisState(n, 3);
        cout << "Input vector:\n";
        dd::Edge state = state1;
        dd->export2Dot(state, "state5before.dot",true, true);
        
        //  dd->printVector(state);
        //  dd::Edge e_swap =permuteOperator(dd,n);//HM: this can be taken out, here for clarity
        state = dd_QFTV2(dd, n, true,state);
         dd->export2Dot(state, "state5after.dot",true, true);
        cout << "Output vector:\n";
        delete dd;
        //  dd->printVector(state);
        //  cout<<"dd->multiply(e_swap, e_qft): "<<endl;
        //  e_qft = dd->multiply(e_swap, e_qft);
        //        state = dd->multiply(e_qft, state);
        //        cout<<"size is: " << dd->size(e_qft)<<endl;
        //       // dd->printDD(e_qft, 100);
        //   dd->export2Dot(state, "qft3v2.dot",true, true);
        //  dd->export2Dot(e_swap, "swap3operator.dot",false, true);
        // dd::Edge state = dd->add(state3, state0);
        //  state = dd->normalize(state, false);
        //        short line[1]={2};
        ////
        //        dd::Edge test = dd->makeGateDD(Hone, 1, line);
        //        dd->export2Dot(test, "testnew.dot",false, true);
        //        dd::Edge testk2 = dd->kronecker(test, test);
        //       // dd->garbageCollect(true);
        //        dd->export2Dot(testk2, "testk2.dot",false, true);
    }
    
    if(true){//plot the graph (time and node count) as function of varible number.
        int N = 20;//points for graph
        int offset = 2;//starting num of qubits
        float a_et[2][N];//execution time
        int nodecounter[N];
        for (int i = 0; i < N; ++i){
            //Initialize package
            auto* dd = new dd::Package;
            time_point<system_clock> start, end;
            start = system_clock::now();
            // dd::Edge e_qft = dd_QFT(dd, offset + i, true);
            dd::Edge e_state = dd->makeBasisState(offset + i, offset + i - 1);
            cout<<dd->size(e_state)<<endl;
            e_state = dd_QFTV2(dd, offset + i, false, e_state);
            nodecounter[i] = dd->size(e_state);
            // state = dd->multiply(e_qft, state);
            end = system_clock::now();
            duration<float> elapsed_seconds = end - start;
            a_et[0][i] = offset+i;
            a_et[1][i] = elapsed_seconds.count();
            delete dd;
        }
        //  dd->printInformation();
        std::ofstream datafile;
        datafile.open ("timers.txt");
        if(datafile.is_open()){
            for (int i = 0; i < N; ++i){
                datafile << a_et[0][i] << " "<< a_et[1][i] <<endl;
            }
            datafile.close();
        }
        datafile.open ("counters.txt");
        if(datafile.is_open()){
            for (int i = 0; i < N; ++i){
                datafile << a_et[0][i] << " "<< nodecounter[i] <<endl;
            }
            datafile.close();
        }
    }
    //    Engine *pE;
    //    mxArray *T = NULL;
    //    pE = engOpen("/Applications/MATLAB_R2019a.app/bin/./matlab");
    //    if (pE == NULL)
    //    {
    //        cout << "Error: Matlab Engine could not be called." << endl;
    //        exit(1);
    //    }
    //    T = mxCreateDoubleMatrix(2, N, mxREAL);
    //    memcpy((void *)mxGetPr(T), (void *)a_et, sizeof(a_et));
    //    engPutVariable(pE, "T", T);
    //    engEvalString(pE, "plot(T(1,:),T(2,:));");
    //    printf("Hit return to continue\n\n");
    //    fgetc(stdin);
    //    mxDestroyArray(T);
    //    engEvalString(pE, "close;");
    return 0;
}
//int matlabCaller(){
//    Engine *ep;
//    mxArray *T = NULL, *result = NULL;
//    char buffer[BUFSIZE+1];
//    if (!(ep = engOpen("/Applications/MATLAB_R2019a.app/bin/./matlab"))) {
//        fprintf(stderr, "\nCan't start MATLAB engine\n");
//        return EXIT_FAILURE;
//}
//}
//swap gate
//n: num qubits. dd: Use existing package (making other package won't work I think)
//constexpr dd::Edge Smat(dd::Package *dd){
//     short line[2] = {2,-1};
//    dd::Edge x1_gate = dd->makeGateDD(Xmat, 2, line);
//    dd::Edge y1_gate = dd->makeGateDD(Ymat, 2, line);
//    dd::Edge z1_gate = dd->makeGateDD(Zmat, 2, line);
//    line[0] = -1;
//    line[1] = 2;
//    dd::Edge x2_gate = dd->makeGateDD(Xmat, 2, line);
//    dd::Edge y2_gate = dd->makeGateDD(Ymat, 2, line);
//    dd::Edge z2_gate = dd->makeGateDD(Zmat, 2, line);
//    dd::Edge iMat4x4 = dd->makeIdent(0, 1);
//    dd::Edge x_gate = dd->multiply(x1_gate, x2_gate);//make 2 qubit dd
//    dd::Edge y_gate = dd->multiply(y1_gate, y2_gate);//make 2 qubit dd
//    dd::Edge z_gate = dd->multiply(z1_gate, z2_gate);//make 2 qubit dd
//    dd::Edge temp1 = dd->add(iMat4x4, x_gate);
//    dd::Edge temp2 = dd->add(y_gate, z_gate);
//    return dd->add(temp1, temp2);
//    return iMat4x4;
//}
//swap gate
//main bin:
//   dd->printVector(state);
//    dd::Matrix2x2 Rmat;
//    RmatGenerator(Rmat,2);
//    short line[2] = {2,1};
//    dd::Edge r_gate = dd->makeGateDD(Rmat, 2, line);
//    line[0] = 2;
//    line[1] = -1;
//    dd::Edge h_gate = dd->makeGateDD(Hmat, 2, line);
//    dd::Edge bit1_edge = dd->multiply(h_gate, r_gate);
//    line[0] = -1;
//    line[1] = 2;
//    h_gate = dd->makeGateDD(Hmat, 2, line);
//    dd::Edge bit2_edge = h_gate;
//    dd::Edge multiplied = dd->multiply(bit1_edge, bit2_edge);
//  dd::Edge swapped = dd->multiply(Xmat, kronekered);
//    state = dd->multiply(multiplied, state);
//    state = dd->multiply(Smat(dd), state);
//    state = dd->normalize(state, false);
//    time_point<system_clock> start, end;
//    start = system_clock::now();
//
// state = dd->multiply(dd_QFT(dd, 13), state);
//    dd_QFT(dd, 5 );
//    end = system_clock::now();
//  state = dd->multiply(Smatv2(dd,2,0,1), state);
//    duration<double> elapsed_seconds = end - start;
//    cout << "elapsed time: "<< elapsed_seconds.count()<<endl;
//    cout << "Output vector:\n";
//   dd->printVector(state);
//    dd->export2Dot(state, "state.dot",true);
//   dd->printInformation();
//    dd::Edge zero_state = dd->makeZeroState(1);
//    dd::Edge one_state = dd->makeBasisState(1, 1);
//    dd::Edge input_state = dd->makeBasisState(2, 0);
//end main bin
void recordDDsize(const dd::Edge& e, std::ofstream& datafile, int v, int g){
    
}
