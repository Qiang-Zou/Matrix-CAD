#include<cstdio>
#include<algorithm>
#include<iostream>
#include<vector>
#include <thrust/device_vector.h>  
#include <thrust/host_vector.h>  
#include <thrust/transform.h>  
#include <thrust/functional.h>  
#include <sys/time.h>
#include "cuda_runtime.h"
#include "cublas_v2.h"
#include <Eigen/Dense>
#include <Eigen/Core>
// #include <unsupported/Eigen/Polynomials>

#include <Eigen/Eigen>

#include"lib/DataStructure.cuh"
#include"lib/MathModule.cuh"
#include"lib/InitModule.cuh"
#include"lib/GeometricModule.cuh"

// #define M_PI 3.14159265


double cpuSecond()
{
    struct timeval tp;
    gettimeofday(&tp,NULL);
    return ((double)tp.tv_sec+(double)tp.tv_usec*1.e-6);
}



int main()
{


    double *d_bigtransN_pointer, *d_inverse_B_pointer, *d_bigD_pointer,*d_big_Q_pointer;
    cudaMalloc((void**)&d_bigtransN_pointer, NUM * SIZE * SIZE * sizeof(double));
    cudaMalloc((void**)&d_inverse_B_pointer, SIZE * SIZE * sizeof(double));
    cudaMalloc((void**)&d_bigD_pointer, NUM * SIZE * SIZE * sizeof(double));
    cudaMalloc((void**)&d_big_Q_pointer, NUM * SIZE * 4 * sizeof(double));
    cudaMemcpyFromSymbol(d_inverse_B_pointer, d_inverse_B, SIZE * SIZE * sizeof(double), 0, cudaMemcpyDeviceToDevice);
    cudaMemcpy(inverse_B, d_inverse_B_pointer, SIZE * SIZE * sizeof(double), cudaMemcpyDeviceToHost);
    // cudaMemcpy(d_bigD_pointer, bigD, NUM * SIZE * SIZE * sizeof(double), cudaMemcpyHostToDevice);

    printf("1\n");

    system("nvidia-smi --query-compute-apps=pid,used_memory --format=csv >> mem_usage_before.txt");
    // init_degree4_2D_easy();
    // init_degree4_2D_normal();
    // init_degree4_2D_hard();
    // init_degree5_2D_easy();
    // init_degree5_2D_normal();
    // init_degree5_2D_hard();
    // init_degree6_2D_easy();
    // init_degree6_2D_normal();
    // init_degree6_2D_hard();
    // init_degree4_easy();
    // init_degree4_3D_normal();
    // init_degree4_complex();
    // init_degree5_easy();
    // init_degree5_normal();
    init_degree5_complex();
    // init_degree6_easy();
    // init_degree6();
    // init_degree6_complex();

    // init();
    // init_degree5();
    // init_degree6();
    // init_degree6_complex();
    for(int i=0;i<knots.size();i++)
    {
        printf("%f ",knots[i]);
    }
    printf("\n");
    for(int i=0;i<control_points.size();i++)
    {
        printf("%f %f %f %f\n",control_points[i].x,control_points[i].y,control_points[i].z,control_points[i].w);
    }
    printf("\n");

    initC<<<1,SIZE*2>>>(SIZE*2);
    cudaDeviceSynchronize();
    double iStart,iElaps;
    // insert_knot(2.5,3,3,0);

    // printf("refine\n");
    // X={1.0,1.0,2.0,2.0,3.0,3.0,4.0,4.0};
    // int r = X.size() - 1;
    // RefineKnotVectCurve(control_points.size()-1, p, knots, control_points, X, r, new_knots, new_control_points); // Provide the required arguments
    // for(int i=0;i<new_knots.size();i++)
    // {
    //     printf("%f ",new_knots[i]);
    // }
    // printf("\n");
    // for(int i=0;i<new_control_points.size();i++)
    // {
    //     printf("%f %f %f %f\n",new_control_points[i].x,new_control_points[i].y,new_control_points[i].z,new_control_points[i].w);
    // }
    // printf("\n");

    printf("decompose\n");
    iStart=cpuSecond();
    // DecomposeCurve(control_points.size()-1, p, knots, control_points, Qw);
    iElaps=cpuSecond()-iStart;
    double time1=iElaps;
 
    // for(int i=0;i<Qw.size();i++)
    // {
    //     printf("i = %d\n",i);
    //     for(int j=0;j<Qw[i].size();j++)
    //     {
    //         printf("%f %f %f %f\n",Qw[i][j].x,Qw[i][j].y,Qw[i][j].z,Qw[i][j].w);
    //     }
    //     printf("\n");
    // }

    int nbytes=NUM*NUM*SIZE*SIZE*sizeof(double);
    thrust::device_vector<Point> d_control_points(control_points.size());  
    thrust::copy(control_points.begin(), control_points.end(), d_control_points.begin());  
    Point* d_control_points_vector = thrust::raw_pointer_cast(d_control_points.data());  

    thrust::device_vector<double> d_knots(knots.size());  
    thrust::copy(knots.begin(), knots.end(), d_knots.begin());  
    double* d_knotsvector = thrust::raw_pointer_cast(d_knots.data());  
    dim3 numblocks=1;
    dim3 numthreads(SIZE,SIZE);


    iStart=cpuSecond();
    gpu_init2<<<(knots.size()-1+255)/256,256>>>(knots.size()-1,d_knotsvector);
    cudaDeviceSynchronize();

    pre_compute3<<<knots.size()-2,numthreads>>>(knots.size(),p,d_knotsvector);
    cudaDeviceSynchronize();
    iElaps=cpuSecond()-iStart;
    double time2=iElaps;
 
    cudaMemcpyFromSymbol(N,d_N,nbytes,0, cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();

    //对bezier段进行重参数化
    // int nbezierbytes=1005*4*4*4*sizeof(double);
    // bezier_compute<<<(knots.size()-1+255)/256,256>>>(knots.size()-1,p,d_knotsvector);
    // cudaDeviceSynchronize();
    // cudaMemcpyFromSymbol(B,d_B,nbezierbytes,0, cudaMemcpyDeviceToHost);
    // cudaDeviceSynchronize();

    //打印N矩阵
    printf("打印N矩阵\n");
    for(int j=0;j<=p;j++)
    {
        for(int i=0;i<knots.size()-j-1;i++)
        {
            printf("N %d %d\n",i,j);
            for(int k=0;k<knots.size()-1;k++)
            {
                printf("第 %d 段 ",k);
                for(int kk=0;kk<=p;kk++)
                {
                    printf("%f ",N[i][j][k][kk]);
                }
                printf("\n");
            }
        }
    }
    printf("\n");

    // for(int i=0;i<knots.size()-p-1;i++)
    // {
    //     printf("第%d个区间\n",i);
    //     for(int j=0;j<=p;j++)
    //     {
    //         printf("B %d %d: ",j,p);
    //         for(int k=0;k<=p;k++)
    //         {
    //             printf("%f ",B[i][j][p][k]);
    //         }
    //         printf("\n");
    //     }
    // }
    // printf("\n");

    thrust::device_vector<double> d_querys(querys.size());
    thrust::copy(querys.begin(),querys.end(),d_querys.begin());
    double* d_querysvector=thrust::raw_pointer_cast(d_querys.data());
    thrust::device_vector<Point> d_answers(querys.size());

    std::vector<Point> answers(querys.size());

    matrix_evaluate<<<(querys.size()-1+255)/256,256>>>(d_querysvector,thrust::raw_pointer_cast(&d_answers[0]),d_control_points_vector,p,d_knotsvector,knots.size());
    cudaDeviceSynchronize();
    thrust::copy(d_answers.begin(),d_answers.end(),answers.begin());

    for(int i=0;i<answers.size();i++)
        printf("%f %f %f %f\n",answers[i].x,answers[i].y,answers[i].z,answers[i].w);


    iStart=cpuSecond();
    //对N矩阵进行重参数化
   // transformation<<<(knots.size()+255)/256,256>>>(knots.size(),p,d_knotsvector);
    trans<<<(knots.size()+255)/256,256>>>(knots.size(),p,d_knotsvector);
    cudaDeviceSynchronize();
    iElaps=cpuSecond()-iStart;
    double time3=iElaps;
   
    int ntransbytes=NUM*SIZE*SIZE*SIZE*sizeof(double);
    cudaMemcpyFromSymbol(transform_N,d_transform_N,ntransbytes,0, cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();

    for(int i=0;i<knots.size()-p-1;i++)
    {
        printf("第%d个区间\n",i);
        for(int j=0;j<=p;j++)
        {
            printf("transform_N %d %d: ",j,p);
            for(int k=0;k<=p;k++)
            {
                printf("%f ",transform_N[i][j][p][k]);
            }
            printf("\n");
        }
    }
    printf("\n");

    int bigtransbytes=NUM*SIZE*SIZE*sizeof(double);
    cudaMemcpyFromSymbol(bigtransN,d_bigtransN,bigtransbytes,0, cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
    //打印bigtransN
    printf("打印bigtransN\n");
    for(int i=0;i<(knots.size()-p-1)*SIZE;i++)
    {
        printf("bigtransN %d: ",i);
        for(int j=0;j<SIZE;j++)
        {
            printf("%f ",bigtransN[(i*SIZE)+j]);
        }
        printf("\n");
    }

    printf("inverse B\n");
    for(int i=0;i<SIZE;i++)
    {
        for(int j=0;j<SIZE;j++)
        {
            printf("%f ",inverse_B[i*SIZE+j]);
        }
        printf("\n");
    }

    //利用cublas对d_bigtransN和d_inverse_B进行相乘
    cublasStatus_t status;
    cublasHandle_t handle;
    status = cublasCreate(&handle);
    if (status != CUBLAS_STATUS_SUCCESS)
    {
        if (status == CUBLAS_STATUS_NOT_INITIALIZED) {
            std::cout << "CUBLAS 对象实例化出错" << std::endl;
        }
        getchar ();
        return EXIT_FAILURE;
    }
    double alpha = 1.0;
    double beta = 0.0;
    cudaMemcpy(d_bigtransN_pointer, bigtransN, NUM * SIZE * SIZE * sizeof(double), cudaMemcpyHostToDevice);

    iStart=cpuSecond();

    getD(d_bigtransN_pointer,d_bigD_pointer,d_inverse_B_pointer,handle,alpha,beta);
    // cublasDgemm(
    //     handle,
    //     CUBLAS_OP_T,
    //     CUBLAS_OP_N,
    //     SIZE,
    //     NUM * SIZE,
    //     SIZE,
    //     &alpha, 
    //     d_inverse_B_pointer, 
    //     SIZE, 
    //     d_bigtransN_pointer,
    //     SIZE, 
    //     &beta, 
    //     d_bigD_pointer, 
    //     SIZE
    // );
    iElaps=cpuSecond()-iStart;
    
    cublasGetVector(
        NUM * SIZE *SIZE, 
        sizeof(double), 
        d_bigD_pointer, 
        1, 
        bigD, 
        1
    );

    double time6=iElaps;
    printf("time on tensor core matrix mul: %f\n",time6);
    // cublasDestroy(handle);
    printf("bigD\n");
    for(int i=0;i<(knots.size()-p-1) * SIZE;i++)
    {
        printf("bigD %d: ",i);
        for(int j=0;j<SIZE;j++)
        {
            printf("%f ",bigD[i*(SIZE)+j]);
        }
        printf("\n");
    }

    //计算B的逆与N相乘的中间矩阵D
    iStart=cpuSecond();
    // compute_D<<<(knots.size()+255)/256,256>>>(knots.size(),p,d_knotsvector);
    // cudaDeviceSynchronize();
    iElaps=cpuSecond()-iStart;
    double time4=iElaps;
  
    int nDbytes=NUM*SIZE*SIZE*SIZE*sizeof(double);
    cudaMemcpyFromSymbol(D,d_D,nDbytes,0, cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();

    // printf("中间矩阵D\n");
    // for(int i=0;i<knots.size()-p-1;i++)
    // {
    //     printf("第%d个区间\n",i);
    //     for(int j=0;j<=p;j++)
    //     {
    //         printf("D %d %d: ",j,p);
    //         for(int k=0;k<=p;k++)
    //         {
    //             printf("%f ",D[i][j][p][k]);
    //         }
    //         printf("\n");
    //     }
    // }
    // printf("\n");

    iStart=cpuSecond();
    //计算D与原控制点矩阵P相乘的控制点矩阵Q
    getQ<<<(knots.size()+255)/256,256>>>(knots.size(),p,d_bigD_pointer,d_control_points_vector,d_big_Q_pointer);
    // compute_Q<<<(knots.size()+255)/256,256>>>(knots.size(),p,d_control_points_vector);
    cudaDeviceSynchronize();
    iElaps=cpuSecond()-iStart;
    double time5=iElaps;

    int nQbytes=NUM*4*SIZE*sizeof(double);
    cudaMemcpyFromSymbol(Q,d_Q,nQbytes,0, cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();

    printf("Q\n");
    for(int i=0;i<knots.size()-p-1;i++)
    {
        printf("第%d个区间\n",i);
        for(int j=0;j<SIZE;j++)
        {
            printf("Q %d: ",j);
            for(int k=0;k<4;k++)
            {
                printf("%f ",Q[i][j][k]);
            }
            printf("\n");
        }
    }

    cudaMemcpy(big_Q, d_big_Q_pointer, NUM * SIZE * 4 * sizeof(double), cudaMemcpyDeviceToHost);
    printf("big Q\n");
    for(int i=0;i<(knots.size()-p-1);i++)
    {
        // printf("big Q %d: ",i);
        for(int j=0;j<4;j++)
        {
            for(int k=0;k<SIZE;k++)
            {
                printf("%f ",big_Q[i*SIZE*4+j*SIZE+k]);
            }
            printf("\n");
        }
        printf("\n");
    }

    initG_m_m<<<4,4>>>(3);//cubic bezier ,so p=3
    initG_m_n<<<SIZE,4>>>(3,p);
    cudaDeviceSynchronize();

    iStart=cpuSecond();
    get_delta<<<(knots.size()-p-1+255)/256,256>>>(3,p,d_big_Q_pointer,p,knots.size());
    cudaDeviceSynchronize();
    iElaps=cpuSecond()-iStart;
    double time7=iElaps;
    cudaMemcpyFromSymbol(cubic_Q,d_cubic_Q,NUM*4*4*sizeof(double),0, cudaMemcpyDeviceToHost);
    printf("cubic Q\n");
    for(int i=0;i<knots.size()-p-1;i++)
    {
        // printf("Q %d: ",i);
        for(int j=0;j<4;j++)
        {
            for(int k=0;k<4;k++)
            {
                printf("%f ",cubic_Q[i][j][k]);
            }
            printf("\n");
        }
        printf("\n");
    }

    CBC *d_bezier;
    cudaMalloc((void**)&d_bezier,SUB_NUM*sizeof(CBC));

    CBC *d_re_sub_bezier;
    cudaMalloc((void**)&d_re_sub_bezier,SUB_NUM*sizeof(CBC));

    CBC *d_after_sub_bezier;
    cudaMalloc((void**)&d_after_sub_bezier,SUB_NUM*sizeof(CBC));

    CBC *d_monotone_segment;
    cudaMalloc((void**)&d_monotone_segment,TOTAL_NUM*sizeof(CBC));


    int cnt_after_sub=0;
    int *d_cnt_after_sub;
    cudaMalloc((void **)&d_cnt_after_sub, sizeof(int));
    cudaMemcpy(d_cnt_after_sub,&cnt_after_sub, sizeof(int), cudaMemcpyHostToDevice);

    // int cnt_after_monotone[100];
    // for(int i=0;i<100;i++)
    // {
    //     cnt_after_monotone[i]=0;
    // }
    // int *d_cnt_after_monotone;
    // cudaMalloc((void **)&d_cnt_after_monotone, sizeof(int)*100);
    // cudaMemcpy(d_cnt_after_monotone,cnt_after_monotone, sizeof(int)*100, cudaMemcpyHostToDevice);

    int *d_cnt_after_monotone, *h_cnt_after_monotone;
    size_t size_test_point = cnt_point * sizeof(int);

    // 分配主机和设备内存
    h_cnt_after_monotone = (int *)malloc(size_test_point);
    cudaMalloc(&d_cnt_after_monotone, size_test_point);

    // 初始化数组
    for (int i = 0; i < cnt_point; i++) {
        h_cnt_after_monotone[i] = 0;
    }
    cudaMemcpy(d_cnt_after_monotone, h_cnt_after_monotone, size_test_point, cudaMemcpyHostToDevice);

    int *cnt_after_elimation=new int[cnt_point];
    for(int i=0;i<cnt_point;i++)
    {
        cnt_after_elimation[i]=0;
    }
    int *d_cnt_after_elimation;
    cudaMalloc((void **)&d_cnt_after_elimation, sizeof(int)*cnt_point);
    cudaMemcpy(d_cnt_after_elimation,cnt_after_elimation, sizeof(int)*cnt_point, cudaMemcpyHostToDevice);

    // calculate_error_kernel<<<(knots.size()-p-1+255)/256,256>>>(p,3,knots.size(),d_bezier);
    // cudaDeviceSynchronize();

    iStart=cpuSecond();

    first_calculate_kernel<<<(knots.size()-p-1+255)/256,256>>>(p,3,knots.size(),d_bezier,d_cnt_after_sub,d_after_sub_bezier);
    cudaDeviceSynchronize();

    iElaps=cpuSecond()-iStart;
    double time8=iElaps;
    printf("time on first_calculate_kernel = %f\n",time8);

    
    double *d_prefix_error_sum;
    int *d_seg;
    int *d_next_seg;
    cudaMalloc((void**)&d_seg,SUB_NUM*sizeof(int));
    cudaMalloc((void**)&d_prefix_error_sum,SUB_NUM*sizeof(double));
    cudaMalloc((void**)&d_next_seg,SUB_NUM*sizeof(int));
    cudaMemcpyFromSymbol(d_seg, seg, SUB_NUM * sizeof(int), 0, cudaMemcpyDeviceToDevice);
    cudaMemcpyFromSymbol(d_prefix_error_sum, prefix_error_sum, SUB_NUM * sizeof(double), 0, cudaMemcpyDeviceToDevice);

    // //利用thrust对seg数组求前缀和
    // thrust::inclusive_scan(d_seg , d_seg + NUM, d_seg);
    // thrust::inclusive_scan(d_prefix_error_sum , d_prefix_error_sum + NUM, d_prefix_error_sum);

    // 利用 Thrust 对 seg 和 prefix_error_sum 数组求前缀和
    thrust::device_ptr<int> d_seg_ptr(d_seg);
    thrust::device_ptr<double> d_prefix_error_sum_ptr(d_prefix_error_sum);
    thrust::inclusive_scan(d_seg_ptr, d_seg_ptr + SUB_NUM, d_seg_ptr);
    thrust::inclusive_scan(d_prefix_error_sum_ptr, d_prefix_error_sum_ptr + SUB_NUM, d_prefix_error_sum_ptr);

    thrust::device_ptr<int> d_next_seg_ptr(d_next_seg);
 
     // 将结果复制回主机并打印

    int* h_seg = new int[SUB_NUM];
    int* h_next_seg = new int[SUB_NUM];
    double *h_prefix_error_sum = new double[SUB_NUM];
    cudaMemcpy(h_seg, d_seg, SUB_NUM * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_prefix_error_sum, d_prefix_error_sum, SUB_NUM * sizeof(double), cudaMemcpyDeviceToHost);

    for(int i=0;i<knots.size()-p-1;i++)
    {
        printf("seg[%d]=%d\n",i,h_seg[i]);
    }
    for(int i=0;i<knots.size()-p-1;i++)
    {
        printf("prefix_error_sum[%d]=%f\n",i,h_prefix_error_sum[i]);
    }
    
    int num_seg=h_seg[knots.size()-p-2];

    iStart=cpuSecond();

    // first_re_subdivision<<<(num_seg+255)/256,256>>>(p,3,num_seg,d_bezier,d_re_sub_bezier,knots.size(),d_seg,d_next_seg);
    first_re_subdivision_at_max_value<<<(num_seg+255)/256,256>>>(p,3,num_seg,d_bezier,d_re_sub_bezier,knots.size(),d_seg,d_next_seg,d_cnt_after_sub,d_after_sub_bezier);
    cudaDeviceSynchronize();

    iElaps=cpuSecond()-iStart;
    double time9=iElaps;
    printf("time on first_re_subdivision = %f\n",time9);

    thrust::inclusive_scan(d_next_seg_ptr, d_next_seg_ptr + SUB_NUM, d_next_seg_ptr);
    cudaMemcpy(h_next_seg, d_next_seg, SUB_NUM * sizeof(int), cudaMemcpyDeviceToHost);
    for(int i=0;i<num_seg;i++)
    {
        printf("next_seg[%d]=%d\n",i,h_next_seg[i]);
    }
    int index_size=num_seg;
    num_seg=h_next_seg[num_seg-1];

    iStart=cpuSecond();
    // re_subdivision<<<(num_seg+255)/256,256>>>(p,3,num_seg,d_re_sub_bezier,d_bezier,index_size,d_next_seg,d_seg);
    re_subdivision_at_max_error<<<(num_seg+255)/256,256>>>(p,3,num_seg,d_re_sub_bezier,d_bezier,index_size,d_next_seg,d_seg,d_cnt_after_sub,d_after_sub_bezier);
    cudaDeviceSynchronize();

    int flag=0;
    for(int i=1;i<=10;i++)
    {
        index_size=num_seg;
        printf("the %d th subdivision\n",i+2);
        if(i%2==1)
        {
            thrust::inclusive_scan(d_seg_ptr, d_seg_ptr + SUB_NUM, d_seg_ptr);
            // num_seg=h_seg[num_seg-1];
            cudaMemcpy(&num_seg, &d_seg[index_size-1], sizeof(int), cudaMemcpyDeviceToHost);//只赋值一个位置的值
            if(num_seg==0)
            {
                printf("subdivision end\n");
                break;
            }
            printf("index_size=%d num_seg=%d\n",index_size,num_seg);
            // cudaMemcpy(h_seg, d_seg, NUM * sizeof(int), cudaMemcpyDeviceToHost);
            // for(int i=0;i<index_size;i++)
            // {
            //     printf("seg[%d]=%d\n",i,h_seg[i]);
            // }
            // re_subdivision<<<(num_seg+255)/256,256>>>(p,3,num_seg,d_bezier,d_re_sub_bezier,index_size,d_seg,d_next_seg);
            re_subdivision_at_max_error<<<(num_seg+255)/256,256>>>(p,3,num_seg,d_bezier,d_re_sub_bezier,index_size,d_seg,d_next_seg,d_cnt_after_sub,d_after_sub_bezier);
            cudaDeviceSynchronize();
        }
        else
        {
            thrust::inclusive_scan(d_next_seg_ptr, d_next_seg_ptr + SUB_NUM, d_next_seg_ptr);
            // num_seg=h_next_seg[num_seg-1];
            cudaMemcpy(&num_seg, &d_next_seg[index_size-1], sizeof(int), cudaMemcpyDeviceToHost);//只赋值一个位置的值
            if(num_seg==0)
            {
                flag=1;
                printf("subdivision end\n");
                break;
            }
            printf("index_size=%d num_seg=%d\n",index_size,num_seg);
            // cudaMemcpy(h_next_seg, d_next_seg, NUM * sizeof(int), cudaMemcpyDeviceToHost);
            // for(int i=0;i<index_size;i++)
            // {
            //     printf("next_seg[%d]=%d\n",i,h_next_seg[i]);
            // }
            // re_subdivision<<<(num_seg+255)/256,256>>>(p,3,num_seg,d_re_sub_bezier,d_bezier,index_size,d_next_seg,d_seg);
            re_subdivision_at_max_error<<<(num_seg+255)/256,256>>>(p,3,num_seg,d_re_sub_bezier,d_bezier,index_size,d_next_seg,d_seg,d_cnt_after_sub,d_after_sub_bezier);
            cudaDeviceSynchronize();
        }
    }
    iElaps=cpuSecond()-iStart;
    double time10=iElaps;
    printf("time on re-subdivision = %f\n",time10);
    printf("flag=%d\n",flag);

    cudaMemcpy(&cnt_after_sub, d_cnt_after_sub, sizeof(int), cudaMemcpyDeviceToHost);
    printf("cnt=%d\n",cnt_after_sub);

    CBC *h_after_sub_bezier=new CBC[cnt_after_sub];
    cudaMemcpy(h_after_sub_bezier,d_after_sub_bezier,cnt_after_sub*sizeof(CBC),cudaMemcpyDeviceToHost);

    for(int i=0;i<cnt_after_sub;i++)
    {
        printf("第%d个区间\n",i);
        printf("original bezier id=%d\n",h_after_sub_bezier[i].original_bezier_curve_id);
        printf("original left = %f right = %f\n",h_after_sub_bezier[i].original_left,h_after_sub_bezier[i].original_right);
        for(int j=0;j<4;j++)
        {
            printf("%f %f %f %f\n",h_after_sub_bezier[i].control_points[j].x,h_after_sub_bezier[i].control_points[j].y,h_after_sub_bezier[i].control_points[j].z,h_after_sub_bezier[i].control_points[j].w);
        }
        printf("\n");
    }

    Point test_point[cnt_point];
    // init_testpoint(test_point,cnt_point);
    // init_randompoint(test_point,cnt_point);
    init_randompoint_0to1(test_point,cnt_point);
    // init_randompoint_0to1_2D(test_point,cnt_point);
    Point *d_test_point; // 定义设备指针
    cudaMalloc((void**)&d_test_point, sizeof(Point) * cnt_point); // 分配设备内存
    cudaMemcpy(d_test_point, test_point, sizeof(Point) * cnt_point, cudaMemcpyHostToDevice); // 将主机内存中的数据拷贝到设备内存中



    double *h_minumum_distance=new double[cnt_point];
    double *h_ans=new double[cnt_point];
    double *h_parameters=new double[cnt_point];
    int *h_spline_id=new int[cnt_point];
    for(int i=0;i<cnt_point;i++)
    {
        h_minumum_distance[i]=1e19;
        h_ans[i]=0;
        h_parameters[i]=0.0;
        h_spline_id[i]=0;
    }
    double *d_minimum_distance;
    cudaMalloc((void**)&d_minimum_distance, cnt_point*sizeof(double));
    double *d_ans;
    cudaMalloc((void**)&d_ans, cnt_point*sizeof(double));
    double *d_parameters;
    cudaMalloc((void**)&d_parameters, cnt_point*sizeof(double));
    int *d_spline_id;
    cudaMalloc((void**)&d_spline_id, cnt_point*sizeof(int));
    cudaMemcpy(d_minimum_distance, h_minumum_distance, cnt_point*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_ans, h_ans, cnt_point*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_parameters, h_parameters, cnt_point*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_spline_id, h_spline_id, cnt_point*sizeof(int), cudaMemcpyHostToDevice);

    iStart=cpuSecond();

    // subdivide_into_monotonic_segments_using_newtwon<<<(cnt_after_sub+255)/256,256>>>(test_point,cnt_after_sub,d_after_sub_bezier,d_monotone_segment,d_cnt_after_monotone);
    // cudaDeviceSynchronize();
    subdivide_into_monotonic_segments_multiPoint<<<(cnt_after_sub*cnt_point+255)/256,256>>>(d_test_point,cnt_point,cnt_after_sub,d_after_sub_bezier,d_monotone_segment,d_cnt_after_monotone);
    cudaDeviceSynchronize();
    iElaps=cpuSecond()-iStart;
    double time11=iElaps;
    printf("time on subdivide_into_monotonic_segments = %f\n",time11);

    cudaMemcpy(h_cnt_after_monotone, d_cnt_after_monotone, size_test_point, cudaMemcpyDeviceToHost);
    printf("cnt_after_monotone\n");
    // for(int i=0;i<cnt_point;i++)
    // {
    //     printf("cnt_after_monotone[%d]=%d\n",i,h_cnt_after_monotone[i]);
    // }

    // printf("total_cnt_after_monotone=%d\n",h_cnt_after_monotone[cnt_point-1]);

    CBC *h_monotone_segment=new CBC[TOTAL_NUM];
    cudaMemcpy(h_monotone_segment,d_monotone_segment,TOTAL_NUM*sizeof(CBC),cudaMemcpyDeviceToHost);
    for(int i=0;i<3;i++)
    {
        printf("第%d个点\n",i);
        for(int j=0;j<h_cnt_after_monotone[i];j++)
        {
            printf("第%d个区间\n",j);
            printf("original bezier id=%d\n",h_monotone_segment[i*100+j].original_bezier_curve_id);
            printf("original left = %f right = %f\n",h_monotone_segment[i*100+j].original_left,h_monotone_segment[i*100+j].original_right);
            for(int k=0;k<4;k++)
            {
                printf("%f %f %f %f\n",h_monotone_segment[i*100+j].control_points[k].x,h_monotone_segment[i*100+j].control_points[k].y,h_monotone_segment[i*100+j].control_points[k].z,h_monotone_segment[i*100+j].control_points[k].w);
            }
            printf("\n");
        }

    }

    int *h_sum_of_cnt_after_monotone=new int[cnt_point];
    int *d_sum_of_cnt_after_monotone;
    cudaMalloc((void**)&d_sum_of_cnt_after_monotone,cnt_point*sizeof(int));
    thrust::device_ptr<int> d_sum_of_cnt_after_monotone_ptr(d_sum_of_cnt_after_monotone);
    thrust::device_ptr<int> d_cnt_after_monotone_ptr(d_cnt_after_monotone);

    iStart=cpuSecond();
    thrust::inclusive_scan(d_cnt_after_monotone_ptr, d_cnt_after_monotone_ptr + cnt_point, d_sum_of_cnt_after_monotone_ptr);
    iElaps=cpuSecond()-iStart;
    double time12=iElaps;
    printf("time on inclusive_scan = %f\n",time12);
    cudaMemcpy(h_sum_of_cnt_after_monotone,d_sum_of_cnt_after_monotone,cnt_point*sizeof(int),cudaMemcpyDeviceToHost);
    // for(int i=0;i<cnt_point;i++)
    // {
    //     printf("sum_of_cnt_after_monotone[%d]=%d\n",i,h_sum_of_cnt_after_monotone[i]);
    // }

    int total_cnt_after_monotone=h_sum_of_cnt_after_monotone[cnt_point-1];
    printf("total_cnt_after_monotone=%d\n",total_cnt_after_monotone);

    CNPB_5 *d_selected_segment;
    cudaMalloc((void**)&d_selected_segment,TOTAL_NUM*sizeof(CNPB_5));

    int *d_map_index;
    cudaMalloc((void**)&d_map_index,TOTAL_NUM*sizeof(int));

    int *h_sum_of_cnt_after_elimation=new int[cnt_point];
    int *d_sum_of_cnt_after_elimation;
    cudaMalloc((void**)&d_sum_of_cnt_after_elimation,cnt_point*sizeof(int));
    thrust::device_ptr<int> d_sum_of_cnt_after_elimation_ptr(d_sum_of_cnt_after_elimation);
    thrust::device_ptr<int> d_cnt_after_elimation_ptr(d_cnt_after_elimation);
    

    
    iStart=cpuSecond();
    elimation_criteria_multipoint<<<(total_cnt_after_monotone+255)/256,256>>>(total_cnt_after_monotone,cnt_point,d_monotone_segment,d_selected_segment,d_test_point,d_sum_of_cnt_after_monotone,d_cnt_after_elimation,d_map_index);
    cudaDeviceSynchronize();
    thrust::inclusive_scan(d_cnt_after_elimation_ptr, d_cnt_after_elimation_ptr + cnt_point, d_sum_of_cnt_after_elimation_ptr);
    iElaps=cpuSecond()-iStart;
    double time13=iElaps;
    printf("time on elimation_criteria = %f\n",time13);

    cudaMemcpy(cnt_after_elimation, d_cnt_after_elimation, cnt_point*sizeof(int), cudaMemcpyDeviceToHost);
    // for(int i=0;i<cnt_point;i++)
    // {
    //     printf("cnt_after_elimation[%d]=%d\n",i,cnt_after_elimation[i]);
    // }

    CNPB_5 *h_selected_segment=new CNPB_5[TOTAL_NUM];
    cudaMemcpy(h_selected_segment,d_selected_segment,TOTAL_NUM*sizeof(CNPB_5),cudaMemcpyDeviceToHost);
    // for(int i=0;i<10;i++)
    // {
    //     for(int j=0;j<cnt_after_elimation[i];j++)
    //     {
    //         printf("第%d个点的第%d个区间\n",i,j);
    //         for(int k=0;k<6;k++)
    //         {
    //             printf("%f ",h_selected_segment[i*cnt_point+j].coe[k]);
    //         }
    //         printf("\n");
    //     }
    // }

    cudaMemcpy(h_sum_of_cnt_after_elimation,d_sum_of_cnt_after_elimation,cnt_point*sizeof(int),cudaMemcpyDeviceToHost);
    // for(int i=0;i<cnt_point;i++)
    // {
    //     printf("sum_of_cnt_after_elimation[%d]=%d\n",i,h_sum_of_cnt_after_elimation[i]);
    // }

    int total_cnt_after_elimation=h_sum_of_cnt_after_elimation[cnt_point-1];
    printf("total_cnt_after_elimation=%d\n",total_cnt_after_elimation);

    CNPB_5 *d_non_parametric_Bezier;
    cudaMalloc((void**)&d_non_parametric_Bezier,total_cnt_after_elimation*sizeof(CNPB_5));
    CNPB_5 *h_non_parametric_Bezier=new CNPB_5[total_cnt_after_elimation];

    iStart=cpuSecond();
    Convert_polynomial_to_non_parametric_Bezier_multipoint<<<(total_cnt_after_elimation+255)/256,256>>>(total_cnt_after_elimation,cnt_point,d_test_point,d_selected_segment,d_non_parametric_Bezier,d_sum_of_cnt_after_elimation,d_map_index,d_monotone_segment);
    cudaDeviceSynchronize();
    iElaps=cpuSecond()-iStart;
    double time14=iElaps;
    printf("time on Convert_polynomial_to_non_parametric_Bezier = %f\n",time14);
    
    cudaMemcpy(h_non_parametric_Bezier,d_non_parametric_Bezier,total_cnt_after_elimation*sizeof(CNPB_5),cudaMemcpyDeviceToHost);
    // for(int i=0;i<10;i++)
    // {
    //     for(int j=0;j<cnt_after_elimation[i];j++)
    //     {
    //         printf("第%d个non_parametric_bezier的第%d个控制点\n",i,j);
    //         for(int k=0;k<6;k++)
    //         {
    //             if(i>0)
    //                 printf("%f ",h_non_parametric_Bezier[h_sum_of_cnt_after_elimation[i-1]+j].coe[k]);
    //             else
    //                 printf("%f ",h_non_parametric_Bezier[j].coe[k]);
    //         }
    //         printf("\n");
    //     }
    // }

    double *h_potential_t_value=new double[total_cnt_after_elimation];
    double *d_potential_t_value;
    cudaMalloc((void**)&d_potential_t_value,total_cnt_after_elimation*sizeof(double));

    
    iStart=cpuSecond();
    _3times_bezier_clipping_multipoint<<<(total_cnt_after_elimation+255)/256,256>>>(total_cnt_after_elimation,cnt_point,d_non_parametric_Bezier,d_potential_t_value,d_map_index,d_monotone_segment,d_test_point,d_sum_of_cnt_after_elimation);
    cudaDeviceSynchronize();
    iElaps=cpuSecond()-iStart;
    double time15=iElaps;
    printf("time on _3times_bezier_clipping = %f\n",time15);

    cudaMemcpy(h_potential_t_value,d_potential_t_value,total_cnt_after_elimation*sizeof(double),cudaMemcpyDeviceToHost);
    // for(int i=0;i<10;i++)
    // {
    //     for(int j=0;j<cnt_after_elimation[i];j++)
    //     {
    //         printf("第%d个点的第%d个potential_t_value=%f\n",i,j,h_potential_t_value[h_sum_of_cnt_after_elimation[i-1]+j]);
    //     }
    // }

    double iinitial=cpuSecond();
    initial_min_distance<<<(cnt_point+255)/256,256>>>(cnt_point,d_minimum_distance,d_test_point,control_points.size(),d_control_points_vector, d_parameters);
    cudaDeviceSynchronize();
    double afterinitial=cpuSecond()-iinitial;
    printf("time on initial_min_distance = %f\n",afterinitial);

    CBC *d_ans_bezier_curve;
    cudaMalloc((void**)&d_ans_bezier_curve,cnt_point*sizeof(CBC));
    CBC *h_ans_bezier_curve=new CBC[cnt_point];


    iStart=cpuSecond();
    get_minumum_distance<<<(total_cnt_after_elimation+255)/256,256>>>(total_cnt_after_elimation,cnt_point,d_sum_of_cnt_after_elimation,d_potential_t_value,d_minimum_distance,d_ans,d_map_index,d_monotone_segment,d_test_point,d_ans_bezier_curve,d_parameters,d_spline_id,d_knotsvector);
    // get_minumum_distance_without_atomic<<<(cnt_point+255)/256,256>>>(total_cnt_after_elimation,cnt_point,d_sum_of_cnt_after_elimation,d_potential_t_value,d_minimum_distance,d_ans,d_map_index,d_monotone_segment,d_test_point,d_ans_bezier_curve);
    cudaDeviceSynchronize();
    iElaps=cpuSecond()-iStart;
    double time16=iElaps;

    cudaMemcpy(h_minumum_distance, d_minimum_distance, cnt_point*sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_ans, d_ans, cnt_point*sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_ans_bezier_curve, d_ans_bezier_curve, cnt_point*sizeof(CBC), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_parameters, d_parameters, cnt_point*sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_spline_id, d_spline_id, cnt_point*sizeof(int), cudaMemcpyDeviceToHost);

    for(int i=0;i<10;i++)
    {
        printf("第%d个点的坐标 (%f, %f, %f, %f\t",i, test_point[i].x,test_point[i].y,test_point[i].z,test_point[i].w);
        printf("第%d个点的最小距离=%f\t",i,h_minumum_distance[i]);
        printf("第%d个点的最小距离对应的t值=%f\t",i,h_ans[i]);
        printf("第%d个点对应第%d个B-spline区间的piecewise bezier\t",i,h_spline_id[i]);
        printf("它的参数值=%f\t",h_parameters[i]);
        printf("第%d个点的最小距离对应的bezier曲线\n",i);
        for(int j=0;j<4;j++)
        {
            printf("%f %f %f %f\n",h_ans_bezier_curve[i].control_points[j].x,h_ans_bezier_curve[i].control_points[j].y,h_ans_bezier_curve[i].control_points[j].z,h_ans_bezier_curve[i].control_points[j].w);
        }

    }
    // cudaMemcpy(&cnt_after_monotone, d_cnt_after_monotone, sizeof(int), cudaMemcpyDeviceToHost);
    // printf("cnt_after_monotone=%d\n",cnt_after_monotone);

    // CBC *h_monotone_segment=new CBC[cnt_after_monotone];
    // cudaMemcpy(h_monotone_segment,d_monotone_segment,cnt_after_monotone*sizeof(CBC),cudaMemcpyDeviceToHost);

    // for(int i=0;i<cnt_after_monotone;i++)
    // {
    //     printf("第%d个区间\n",i);
    //     for(int j=0;j<4;j++)
    //     {
    //         printf("%f %f %f %f\n",h_monotone_segment[i].control_points[j].x,h_monotone_segment[i].control_points[j].y,h_monotone_segment[i].control_points[j].z,h_monotone_segment[i].control_points[j].w);
    //     }
    //     printf("\n");
    // }






    // CNPB_5 *d_non_parametric_Bezier;
    // cudaMalloc((void**)&d_non_parametric_Bezier,cnt_after_monotone*sizeof(CNPB_5));

    // int *d_map_index;
    // cudaMalloc((void**)&d_map_index,cnt_after_monotone*sizeof(int));

    // iStart=cpuSecond();

    // elimation_criteria<<<(cnt_after_monotone+255)/256,256>>>(cnt_after_monotone,d_monotone_segment,d_selected_segment,test_point,d_cnt_after_elimation,d_map_index);
    // cudaDeviceSynchronize();
    // iElaps=cpuSecond()-iStart;
    // double time12=iElaps;
    // printf("time on elimation_criteria = %f\n",time12);

    // cudaMemcpy(&cnt_after_elimation, d_cnt_after_elimation, sizeof(int), cudaMemcpyDeviceToHost);
    // printf("cnt_after_elimation=%d\n",cnt_after_elimation);

    // CNPB_5 *h_selected_segment=new CNPB_5[cnt_after_elimation];
    // cudaMemcpy(h_selected_segment,d_selected_segment,cnt_after_elimation*sizeof(CNPB_5),cudaMemcpyDeviceToHost);

    // for(int i=0;i<cnt_after_elimation;i++)
    // {
    //     printf("第%d个区间\n",i);
    //     for(int j=0;j<6;j++)
    //     {
    //         printf("%f ",h_selected_segment[i].coe[j]);
    //     }
    //     printf("\n");
    // }

    // iStart=cpuSecond();
    // Convert_polynomial_to_non_parametric_Bezier<<<(cnt_after_elimation+255)/256,256>>>(cnt_after_elimation,test_point,d_selected_segment,d_non_parametric_Bezier,d_map_index,d_monotone_segment);
    // cudaDeviceSynchronize();
    // iElaps=cpuSecond()-iStart;
    // double time13=iElaps;
    // printf("time on Convert_polynomial_to_non_parametric_Bezier = %f\n",time13);

    // CNPB_5 *h_non_parametric_Bezier=new CNPB_5[cnt_after_elimation];
    // cudaMemcpy(h_non_parametric_Bezier,d_non_parametric_Bezier,cnt_after_elimation*sizeof(CNPB_5),cudaMemcpyDeviceToHost);

    // for(int i=0;i<cnt_after_elimation;i++)
    // {
    //     printf("第%d个区间\n",i);
    //     for(int j=0;j<6;j++)
    //     {
    //         printf("%f ",h_non_parametric_Bezier[i].coe[j]);
    //     }
    //     printf("\n");
    // }

    // double *d_potential_projection_points_distance;
    // cudaMalloc((void**)&d_potential_projection_points_distance,cnt_after_elimation*sizeof(double));

    // iStart=cpuSecond();
    // _3times_bezier_clipping<<<(cnt_after_elimation+255)/256,256>>>(cnt_after_elimation,d_non_parametric_Bezier,d_potential_projection_points_distance,d_map_index,d_monotone_segment,test_point);
    // cudaDeviceSynchronize();
    // iElaps=cpuSecond()-iStart;
    // double time14=iElaps;
    // printf("time on _3times_bezier_clipping = %f\n",time14);

    // //把结果存进新的bezier段里
    // for(int i=0;i<knots.size()-p-1;i++)
    // {
    //     printf("第%d个区间\n",i);
    //     int judge=0;
    //     for(int j=0;j<=p;j++)
    //     {
    //         printf("Q %d: ",j);
    //         for(int k=0;k<4;k++)
    //         {
    //             printf("%f ",Q[i][j][k]);
    //             if(Q[i][j][k]!=0.0)
    //                 judge=1;
    //         }
    //         printf("\n");
    //     }
    //     if(judge==1)
    //     {
    //         BC temp;
    //         for(int j=0;j<=p;j++)
    //         {
                
    //             for(int k=0;k<4;k++)
    //             {
    //                 temp.control_points[j].x=Q[i][j][0];
    //                 temp.control_points[j].y=Q[i][j][1];
    //                 temp.control_points[j].z=Q[i][j][2];
    //                 temp.control_points[j].w=Q[i][j][3];
    //             }
    //         }
    //         bezier_curves.push_back(temp);
    //     }
    //     printf("\n");
    // }
    // printf("\n");

    // for(int i=0;i<bezier_curves.size();i++)
    // {
    //     printf("i = %d\n",i);
    //     for(int j=0;j<4;j++)
    //     {
    //         printf("%f %f %f %f\n",bezier_curves[i].control_points[j].x,bezier_curves[i].control_points[j].y,bezier_curves[i].control_points[j].z,bezier_curves[i].control_points[j].w);
    //     }
    //     printf("\n");
    // }

    // //初始化bezier曲线的subdivision相关的变量
    // thrust::device_vector<BC> d_bezier_curves(bezier_curves.size());
    // thrust::copy(bezier_curves.begin(),bezier_curves.end(),d_bezier_curves.begin());
    // BC* d_bezier_curves_vector=thrust::raw_pointer_cast(d_bezier_curves.data());
    
    // thrust::device_vector<BC> d_subdivision_result(bezier_curves.size()*2);
    // BC* d_subdivision_result_vector=thrust::raw_pointer_cast(d_subdivision_result.data());

    // //初始化subdivision的系数
    // initus();
    // thrust::device_vector<double> d_us(us.size());
    // thrust::copy(us.begin(),us.end(),d_us.begin());
    // double* d_us_vector=thrust::raw_pointer_cast(d_us.data());

    // //进行bezier曲线的subdivision
    // bezier_subdivision<<<(us.size()-1+255)/256,256>>>(d_bezier_curves_vector,d_subdivision_result_vector,us.size(),d_us_vector);
    // cudaDeviceSynchronize();

    // std::vector<BC> subdivision_result;
    // subdivision_result.resize(bezier_curves.size()*2);
    // thrust::copy(d_subdivision_result.begin(),d_subdivision_result.end(),subdivision_result.begin());

    // printf("subdivision\n");
    // for(int i=0;i<subdivision_result.size();i++)
    // {
    //     printf("i = %d\n",i);
    //     for(int j=0;j<4;j++)
    //     {
    //         printf("%f %f %f %f\n",subdivision_result[i].control_points[j].x,subdivision_result[i].control_points[j].y,subdivision_result[i].control_points[j].z,subdivision_result[i].control_points[j].w);
    //     }
    //     printf("\n");
    // }

    system("nvidia-smi --query-compute-apps=pid,used_memory --format=csv >> mem_usage_after.txt");

    printf("time on cpu=%f\n",time1);
    printf("time on pre compute N=%f\n",time2);
    printf("time on reparametration by recurrence relation=%f\n",time3);
    printf("time on compute D using tensor core=%f\n",time6);
    printf("time on compute Q by Eigen=%f\n",time5);
    printf("time on cubic bezier approximation =%f\n",time7);
    printf("time on total re-subdivision = %f\n",time8+time9+time10);
    printf("time on subdivide_into_monotonic_segments = %f\n",time11+time12);
    printf("time on elimation_criteria = %f\n",time13);
    printf("time on Convert_polynomial_to_non_parametric_Bezier = %f\n",time14);
    printf("time on _3times_bezier_clipping = %f\n",time15);
    printf("time on get_minumum_distance = %f\n",time16+afterinitial);
    double total_time=time2+time3+time6+time5+time7+time8+time9+time10+time11+time12+time13+time14+time15+time16+afterinitial;
    printf("time on GPU = %f\n",time2+time3+time6+time5+time7+time8+time9+time10+time11+time12+time13+time14+time15+time16+afterinitial);
    printf("time on total decomposition = %f\n",time2+time3+time5+time6+time7+time8+time9+time10+time11+time12);
    printf("time on projection = %f\n",time13+time14+time15+time16+afterinitial);
    printf("average time=%f\n",total_time/cnt_point*1e6);

    // printtestG<<<1,1>>>(3,4);
    // cudaDeviceSynchronize();

    printf("Size Of CBC is %ld\n",sizeof(CBC));

    printf("time on decomposition = %f\n",time2+time3+time5+time6);
    printf("time on cubic decomposition = %f\n",time7+time8+time9+time10);
    printf("time on error controlled decomposition = %f\n",time2+time3+time5+time6+time7+time8+time9+time10);
    printf("time on monotonic decomposition = %f\n",time11+time12);
    printf("time on GPU-optimized projection = %f\n",time13+time14+time15+time16+afterinitial);

    //试下求解四次方程
    double equ_a = 1;
    double equ_b = -10;
    double equ_c = 35;
    double equ_d = -50;
    double equ_e = 24; // 示例系数
    solveQuartic(equ_a,equ_b,equ_c,equ_d,equ_e);

    solvePolynomial(equ_a,equ_b,equ_c,equ_d,equ_e);

    double equ_x[4], equ_p[5];
	int rootCount;
	UINT32 errNo = ERR_NO_ERROR;

	//一元四次方程测试
	// （1）(x - 1) * (x - 2) * (x^2 + 1) = 0 (x^4 - 3*x^3 + 3*x^2 - 3*x + 2 = 0)

	//(2) (x - 1)^2 * (x^2 + 1) = 0 (x^4 - 2*x^3 + 2*x^2 - 2*x + 1 = 0)
	//a = 1;
	//b = -2;
	//c = 2;
	//d = -2;
	//e = 1;

	//(3) (x - 1) * (x - 2) * (x - 3) * (x - 4) = 0 (x^4 - 10*x^3 + 35*x^2 - 50*x + 24 = 0)
	//a = 1;
	//b = -10;
	//c = 35;
	//d = -50;
	//e = 24;

	//(4) (x - 1)^2 * (x - 2)^2 = 0 (x^4 - 6*x^3 + 13*x^2 - 12*x + 4 = 0)
	//a = 1;
	//b = -6;
	//c = 13;
	//d = -12;
	//e = 4;

	//(5) 0*x^4 + x^3 - 3*x^2 + 3*x - 1 = 0
	//a = 0;
	//b = 1;
	//c = -3;
	//d = 3;
	//e = -1;

	//(6) 0*x^4 + 0*x^3 + x^2 - 2*x + 1 = 0
	//a = 0;
	//b = 0;
	//c = 1;
	//d = -2;
	//e = 1;

	//(7) 0*x^4 + 0*x^3 + 0*x^2 - 2*x + 1 = 0
	//a = 0;
	//b = 0;
	//c = 0;
	//d = -2;
	//e = 1;

	equ_p[0] = equ_e;
	equ_p[1] = equ_d;
	equ_p[2] = equ_c;
	equ_p[3] = equ_b;
	equ_p[4] = equ_a;
	errNo = solve_quartic_equation(equ_p, equ_x, &rootCount);

    // printf("费尔曼方法\n");
    // printf("rootCount = %u\n",rootCount);
    // for(int i=0;i<rootCount;i++)
    // printf("%f ",equ_x[i]);
    // printf("\n");

    // printf("新方法\n");
    // std::vector<double> coefficients = {equ_a, equ_b, equ_c, equ_d, equ_e};  // 例：x^4 - 4x^3 + 6x^2 - 4x + 1
    // auto result = cal_quartic_ik(coefficients);
    
    // std::cout << "Number of real roots: " << result.first << std::endl;
    // for (const auto& root : result.second) {
    //     std::cout << "Root: " << root << std::endl;
    // }



    // compute_pre_G<<<1,1>>>();

    //print C
    // for(int i=0;i<SIZE;i++)
    // {
    //     for(int j=0;j<SIZE;j++)
    //     {
    //         printf("%d ",C[i][j]);
    //     }
    //     printf("\n");
    // }
    // printC<<<1,32>>>();
    printG<<<SIZE,4>>>(3,p);
    cudaDeviceSynchronize();

    // init_degree4_complex();
    for(int i=0;i<control_points.size();i++)
    {
        printf("%f %f %f %f\n",control_points[i].x,control_points[i].y,control_points[i].z,control_points[i].w);
    }


    //printinverseB
    // cudaMemcpyFromSymbol(inverse_B,d_inverse_B,SIZE*SIZE*sizeof(double),0, cudaMemcpyDeviceToHost);
    // cudaDeviceSynchronize();
    // printf("inverse B\n");
    // for(int i=0;i<SIZE;i++)
    // {
    //     for(int j=0;j<SIZE;j++)
    //     {
    //         printf("%f ",inverse_B[i][j]);
    //     }
    //     printf("\n");
    // }
}