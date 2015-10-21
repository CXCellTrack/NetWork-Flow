#include "mex.h"
#include "matrix.h"
#include <math.h>
#include <stdio.h>
#define pi 3.1415926
void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
//     char *para;
    double *e1; 
    double *im;
    double *output0;  
    
    
    // intensity feature
    int height, width;
    double ind_hist, max_intensity=0, min_intensity=1;
    double intensity_sum = 0;
    double intensity_mean = 0;
//     double intensity_devia = 0;
    
    double a1, b1, x01, y01, alpha1;
    double cc1, c11_x, c11_y, c12_x, c12_y;
    double x,y;
    double count = 0;
    double d2r = atan(1.0)/45.0;
    
//     
//     
//     if ( nrhs != 3 ) 
//     {
//         mexErrMsgTxt("three argument required: e_4_c, im, para");
//     }

    if (nrhs==2) //e_4_c, im, para
    {
        
//         e1 = mxCalloc(mxGetN(prhs[0]), sizeof(double));
        e1 = mxGetPr(prhs[0]);
        
//         im = mxCalloc(mxGetM(prhs[1])* mxGetN(prhs[1]), sizeof(double));
        im = mxGetPr(prhs[1]);
        
        height = mxGetM(prhs[1]);
        width = mxGetN(prhs[1]);
//         mexPrintf("%d %d",height,width);
//         para = mxGetPr(prhs[2]);
        
//         // ����Ĳ����ǡ�move��
//         if (*para == 'm' && mxGetN(prhs[2])==4)  
//         {
            plhs[0] = mxCreateDoubleMatrix(1, 11, mxREAL);
    //         plhs[1] = mxCreateDoubleMatrix(height, width, mxREAL);
            output0 = mxGetPr(plhs[0]);
    //         intensity_devia = mxGetPr(plhs[1]);

            a1 = e1[0];
            b1 = e1[1];
            x01 = e1[2];
            y01 = e1[3];
            alpha1 = e1[4]*d2r;

            cc1=sqrt(a1*a1-b1*b1);
            c11_x = x01 + cc1*cos(alpha1);
            c11_y = y01 + cc1*sin(alpha1);
            c12_x = x01 - cc1*cos(alpha1);
            c12_y = y01 - cc1*sin(alpha1);

            for( x=x01-a1; x<=x01+a1; x++ )
                for( y=y01-a1; y<=y01+a1; y++ )
                    // ����õ�����Բ�ڣ����������
                    if ( sqrt((x-c11_x)*(x-c11_x)+(y-c11_y)*(y-c11_y))+sqrt((x-c12_x)*(x-c12_x)+(y-c12_y)*(y-c12_y))<=2*a1 )
                    {
                        // ����ҶȲ�
                        if ( max_intensity < im[(int)x * height + (int)y] )
                            max_intensity = im[(int)x * height + (int)y];
                        if ( min_intensity > im[(int)x * height + (int)y] )
                            min_intensity = im[(int)x * height + (int)y];

                        // ����Ҷ��ܺ�
                        intensity_sum += im[(int)x * height + (int)y];
                        count = count + 1;
                        // ����Ҷ�ֱ��ͼ
                        for( ind_hist = 0; ind_hist<8; ind_hist++ )
                            if ( im[(int)x * height + (int)y] >= ind_hist/8 && im[(int)x * height + (int)y] < (ind_hist + 1.0)/8 )
                                output0[(int)ind_hist] += 1;  // ���ʹ�ûҶȼ�Ȩ��Ϊ im[(int)y * width + (int)x];
                        if ( im[(int)x * height + (int)y] == 1 )  //�Ҷ�Ϊ1���鵽��8��
                            output0[7] += 1;
                        //����Ҷ��ݶ�
    //                     intensity_devia[(int)x * height + (int)y] = sqrt( (im[(int)x * height + (int)y+1] - im[(int)x * height + (int)y-1]) * (im[(int)x * height + (int)y+1] - im[(int)x * height + (int)y-1]) + (im[(int)(x+1) * height + (int)y] - im[(int)(x-1) * height + (int)y]) * (im[(int)(x+1) * height + (int)y] - im[(int)(x-1) * height + (int)y]) );
                    }
            // ƽ���Ҷ�
            intensity_mean = intensity_sum/count;
        // output0�����ǰ8��Ϊ�Ҷ�ֱ��ͼ��һ��������ֱ�Ϊ�ҶȺ͡��ҶȾ�ֵ���ҶȲ�  
        for( ind_hist = 0; ind_hist<8; ind_hist++ )
            output0[(int)ind_hist] = output0[(int)ind_hist]/count;
        output0[8] = intensity_sum;
        output0[9] = intensity_mean;
        output0[10] = max_intensity - min_intensity;
        
        
//         mxFree(e1);
//         mxFree(im);
              
//     }
    } // end if nrhs==2
}

    
    
    
    
    
    
    