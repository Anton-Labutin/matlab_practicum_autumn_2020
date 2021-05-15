#include "mex.h"
#include <math.h>
#include "matrix.h"
#include <complex.h>
#include <complex>
#include <iostream>

using namespace std;

typedef struct Complex{
    double re;
    double im;
} Complex;


Complex
sum(const Complex& a, const Complex& b)
{
    Complex result;
    result.re = a.re + b.re;
    result.im = a.im + b.im;
    return result;
}


Complex 
diff(const Complex& a, const Complex& b)
{
    Complex result;
    result.re = a.re - b.re;
    result.im = a.im - b.im;
    return result;
}

Complex
div(const Complex& a, const Complex& b)
{
    Complex result;
    double tmp = b.re * b.re + b.im * b.im;
    result.re = (a.re * b.re + a.im * b.im) / tmp;
    result.im = (b.re * a.im - a.re * b.im) / tmp;
    
    return result;
}


Complex 
mul(const Complex& a, const Complex& b)
{
    Complex result;
    result.re = a.re * b.re - a.im * b.im;
    result.im = a.re * b.im + a.im * b.re;
    return result;
}


typedef struct Polar {
        double abs;
        double angle;
} Polar;
    
    
Polar 
get_polar(const Complex& c)
{
    Polar result;
    result.abs = pow(c.re * c.re + c.im * c.im, 0.5);
    
    double tmp = atan(c.im / c.re);;
    if (c.re > 0 && c.im >= 0)
        result.angle = tmp;
    else if (c.re > 0 && c.im < 0)
        result.angle = 2 * M_PI + tmp;
    else if (c.re < 0)
        result.angle = M_PI + tmp;
    else if (c.re == 0 && c.im > 0)
        result.angle = M_PI / 2.0;
    else
        result.angle = -(M_PI / 2.0);
    
    return result;
}


Complex*
root(const Complex& c, int power)
{
    Polar polar_c = get_polar(c);
           
    double rootabs = pow(polar_c.abs, 1 / (double) power);
    
    Complex *result = (Complex*) malloc(sizeof(Complex) * power);
    for(int i = 0; i < power; ++i) {
        (result[i]).re = rootabs * cos((polar_c.angle + i * 2 * M_PI) / power);
        (result[i]).im = rootabs * sin((polar_c.angle + i * 2 * M_PI) / power);
    }
    
    return result;
}


void 
cubesolve(const mxArray *inputA, const mxArray *inputB, const mxArray *inputC, 
        mxArray *inputX1, mxArray *inputX2, mxArray *inputX3, 
        mwSize clmCnt, mwSize rowCnt)
{    
    mxComplexDouble *A = mxGetComplexDoubles(inputA);
    mxComplexDouble *B = mxGetComplexDoubles(inputB);
    mxComplexDouble *C = mxGetComplexDoubles(inputC);

    mxComplexDouble *X1 = mxGetComplexDoubles(inputX1);
    mxComplexDouble *X2 = mxGetComplexDoubles(inputX2);
    mxComplexDouble *X3 = mxGetComplexDoubles(inputX3);
    
    Complex t1 = {0, 0};
    Complex t2 = {0, 0};
    Complex det;
    Complex p, q;
    Complex pc, qs;
    Complex detsq;
    Complex phi = {0, 0};
    
    mwSize AElemCnt = clmCnt * rowCnt;
    for(mwSize i = 0; i < AElemCnt; ++i) {
        Complex a, b;
        b.re = B[i].real;
        b.im = B[i].imag;
        a.re = A[i].real;
        a.im = A[i].imag;
        
        Complex result = div(b, a);
        
        p.re = result.re / 3.0;
        p.im = result.im / 3.0;
        
        a.re = C[i].real;
        a.im = C[i].imag;
        b.re = A[i].real;
        b.im = A[i].imag;
        
        result = div(a, b);
        
        q.re = result.re / 2.0;
        q.im = result.im / 2.0;

        result = mul(q, q);
        qs.re = result.re;
        qs.im = result.im;
        
        result = mul(p, p);
        pc.re = result.re;
        pc.im = result.im;
        
        result = mul(pc, p);
        pc.re = result.re;
        pc.im = result.im;

        result = sum(qs, pc);
        det.re = result.re;
        det.im = result.im;
        
        Complex *rootResult;
        rootResult = root(det, 2);
        detsq.re = (rootResult[0]).re;
        detsq.im = (rootResult[0]).im;
        free(rootResult);
      
        q.re *= -1;
        q.im *= -1;
        
        result = sum(q, detsq);
        t1.re = result.re;
        t1.im = result.im;

        result = diff(q, detsq);
        t2.re = result.re;
        t2.im = result.im;

        rootResult = root(t1, 3);
        X1[i].real = (rootResult[0]).re;
        X1[i].imag = (rootResult[0]).im;
        X2[i].real = (rootResult[1]).re;
        X2[i].imag = (rootResult[1]).im;
        X3[i].real = (rootResult[2]).re;
        X3[i].imag = (rootResult[2]).im;
        free(rootResult);
        
        rootResult = root(t2, 3);       
        if((p.re <= 0 && p.im <= 0 && q.re >= 0 && q.im > 0) || 
           (p.re >= 0 && p.im <=0 && q.re >= 0 && q.im >= 0) || 
           (p.re >= 0 && p.im <= 0 && q.re >= 0 && q.im <= 0) || 
           (p.re >= 0 && p.im < 0 && q.re <= 0 && q.im >= 0) || 
           (p.re >= 0 && p.im >= 0 && q.re >= 0 && q.im <= 0)
           ) {    
            X1[i].real += (rootResult[0]).re;
            X1[i].imag += (rootResult[0]).im;
            X2[i].real += (rootResult[2]).re;
            X2[i].imag += (rootResult[2]).im;
            X3[i].real += (rootResult[1]).re;
            X3[i].imag += (rootResult[1]).im;
        }
        else if((p.re <= 0 && p.im <= 0 && q.re >= 0 && q.im <= 0) || 
                (p.re <= 0 && p.im < 0 && q.re <= 0 && q.im >= 0) || 
                (p.re <= 0 && p.im <= 0 && q.re <= 0 && q.im <= 0) || 
                (p.re < 0 && p.im >= 0 && q.re >= 0 && q.im >= 0) || 
                (p.re == 0 && p.im >= 0 && q.re <= 0 && q.im >= 0)
                ) {
            X1[i].real += (rootResult[2]).re;
            X1[i].imag += (rootResult[2]).im;
            X2[i].real += (rootResult[1]).re;
            X2[i].imag += (rootResult[1]).im;
            X3[i].real += (rootResult[0]).re;
            X3[i].imag += (rootResult[0]).im;
        }
        else {
            X1[i].real += (rootResult[1]).re;
            X1[i].imag += (rootResult[1]).im;
            X2[i].real += (rootResult[0]).re;
            X2[i].imag += (rootResult[0]).im;
            X3[i].real += (rootResult[2]).re;
            X3[i].imag += (rootResult[2]).im;
        }
   
        free(rootResult);
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    size_t AClmCnt = mxGetM(prhs[0]);
    size_t ARowCnt = mxGetN(prhs[0]);
    size_t BClmCnt = mxGetM(prhs[1]);
    size_t BRowCnt = mxGetN(prhs[1]);
    size_t CClmCnt = mxGetM(prhs[2]);
    size_t CRowCnt = mxGetN(prhs[2]);
    
    if(AClmCnt != BClmCnt || BClmCnt != CClmCnt || ARowCnt != BRowCnt || BRowCnt != CRowCnt)
    {
        mexErrMsgIdAndTxt("rowCnt / clmCnt", "Dimensions must be the same.");
    }
    
    if(nrhs != 3)
    {
        mexErrMsgIdAndTxt("nrhs", "Three inputs required");
    }
    if(nlhs != 3)
    {
        mexErrMsgIdAndTxt("nlhs", "Three outputs reqiured");
    }
    
    plhs[0] = mxCreateDoubleMatrix((mwSize) ARowCnt, (mwSize) AClmCnt, mxCOMPLEX);
    plhs[1] = mxCreateDoubleMatrix((mwSize) ARowCnt, (mwSize) AClmCnt, mxCOMPLEX);
    plhs[2] = mxCreateDoubleMatrix((mwSize) ARowCnt, (mwSize) AClmCnt, mxCOMPLEX);

    cubesolve(prhs[0], prhs[1], prhs[2], plhs[0], plhs[1], plhs[2], (mwSize) AClmCnt, (mwSize) ARowCnt);
}