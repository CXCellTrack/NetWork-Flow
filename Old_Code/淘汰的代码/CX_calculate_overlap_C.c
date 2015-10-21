#include "mex.h"
#include<math.h>
#include<stdio.h>
#define pi 3.1415926535898
void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
    double *e1,*e2;
    double *overlap; 
    double a1,b1,x01,y01,alpha1,a2,b2,x02,y02,alpha2;
    double cc1,c11_x,c11_y,c12_x,c12_y,cc2,c21_x,c21_y,c22_x,c22_y;
    double x,y;
    double count = 0;
    double count1 = 0;
    double count2 = 0;
    
    if (nrhs!=2)
    {
        mexErrMsgTxt("two argument required: e1_4_c, e2_4_c");
    }
    e1 = mxGetPr(prhs[0]);
    e2 = mxGetPr(prhs[1]);
    plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
    overlap = mxGetPr(plhs[0]);
    
    a1 = e1[0];
    b1 = e1[1];
    x01 = e1[2];
    y01 = e1[3];
    alpha1 = e1[4]/180*pi;
    
    a2 = e2[0];
    b2 = e2[1];
    x02 = e2[2];
    y02 = e2[3];
    alpha2 = e2[4]/180*pi;
    
    
    cc1=sqrt(a1*a1-b1*b1);
    c11_x = x01 + cc1*cos(alpha1);
    c11_y = y01 + cc1*sin(alpha1);
    c12_x = x01 - cc1*cos(alpha1);
    c12_y = y01 - cc1*sin(alpha1);
//         for( x=x01-a1; x<=x01+a1; x++ )
//             for( y=y01-a1; y<=y01+a1; y++ )
//                 if (sqrt((x-c11_x)*(x-c11_x)+(y-c11_y)*(y-c11_y))+sqrt((x-c12_x)*(x-c12_x)+(y-c12_y)*(y-c12_y))<=2*a1)
//                     count1 = count1 + 1;
    count1 = pi*a1*b1;
    count2 = pi*a2*b2;//采用面积计算，加快速度很快

    cc2=sqrt(a2*a2-b2*b2);
    c21_x = x02 + cc2*cos(alpha2);
    c21_y = y02 + cc2*sin(alpha2);
    c22_x = x02 - cc2*cos(alpha2);
    c22_y = y02 + cc2*sin(alpha2);
    for( x=x02-a2; x<=x02+a2; x++ )
        for( y=y02-a2; y<=y02+a2; y++ )
            if (sqrt((x-c21_x)*(x-c21_x)+(y-c21_y)*(y-c21_y))+sqrt((x-c22_x)*(x-c22_x)+(y-c22_y)*(y-c22_y))<=2*a2)
            {
//                     count2 = count2 + 1;
                if (sqrt((x-c11_x)*(x-c11_x)+(y-c11_y)*(y-c11_y))+sqrt((x-c12_x)*(x-c12_x)+(y-c12_y)*(y-c12_y))<=2*a1)
                    count = count + 1; 
            }

    overlap[0] = count/(count1+count2-count);    
//         overlap[1] = count;
//         overlap[2] = count1;
//         overlap[3] = count2;
    
}    
    
    
    
    
    
    
    