#pragma once
#include"DataStructure.cuh"
#include"MathModule.cuh"

__device__ double calculate_point_distance(Point a,Point b)
{
    return sqrt((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y)+(a.z-b.z)*(a.z-b.z));
}

__device__ Point evaluate_cubic_bezier(const Point& p0, const Point& p1, const Point& p2, const Point& p3, double t) 
{
    double u = 1 - t;
    double tt = t * t;
    double uu = u * u;
    double uuu = uu * u;
    double ttt = tt * t;

    Point p =  p0 * uuu;
    p =  (p1 * 3 * uu * t) + p;
    p =  (p2 * 3 * u * tt) + p;
    p =  (p3 * ttt) + p;

    return p;
}
// __device__ struct d_cubic_bezier_curve
// {
//     Point control_points[4];
//     int dep;
//     int seg;
//     int id;
//     d_cubic_bezier_curve *child=NULL;
// };
// typedef d_cubic_bezier_curve d_CBC;

// std::vector<double> new_weights;

// void insert_knot_once(double u,int p)
// {
//     int n = control_points.size();
//     int m = knots.size();
//     int k = std::distance(knots.begin(), std::upper_bound(knots.begin(), knots.end(), u));
//     new_control_points.resize(n+1);
//     new_weights.resize(n+1);
//     new_knots.resize(m+1);
//     for(int i=0; i<=k-p; i++)
//     {
//         new_control_points[i] = control_points[i];
//         new_weights[i] = weights[i];
//     }
//     for(int i=k-p+1; i<=k; i++)
//     {
//         double alpha = (u-knots[i])/(knots[i+p]-knots[i]);
//         new_control_points[i] = control_points[i-1]*alpha + control_points[i]*(1-alpha);
//         new_weights[i] = weights[i-1]*alpha + weights[i]*(1-alpha);
//     }
//     for(int i=k+1; i<=n; i++)
//     {
//         new_control_points[i] = control_points[i];
//         new_weights[i] = weights[i];
//     }
//     for(int i=0; i<=k; i++)
//     {
//         new_knots[i] = knots[i];
//     }
//     new_knots[k+1] = u;
//     for(int i=k+1; i<=m; i++)
//     {
//         new_knots[i+1] = knots[i];
//     }
// }

// __global__ void print(int x,int n)
// {
//     int id=blockDim.x*blockIdx.x+threadIdx.x;
//     if(id<n)
//     {
//         printf("x=%d\n",x);
//         printf("x*10=%d\n",x*10);
//         printf("x*100=%d\n",x*100);
//     }
// }

// void initC() 
// {
//     for (int i = 0; i < SIZE; ++i) 
//         for (int j = 0; j <= i; ++j)
//             if (!j) 
//                 C[i][j] = 1;
//             else 
//                 C[i][j] = (C[i-1][j] + C[i-1][j-1]);
// }

__device__ Point evaluate_bezier_curve(double u, Point *control_points, int p)
{
    Point result;
    double left[SIZE];
    double right[SIZE];
    left[0] = 1;
    right[0] = 1;
    for (int i = 1; i <= p; i++)
    {
        left[i] = left[i-1] * (1 - u);
        right[i] = right[i-1] * u;
    }
    for (int i = 0; i <= p; i++)
    {
        result = result + control_points[i] * left[p-i] * right[i] * C[p][i];
    }
    return result;
}

// 二分查找确定 `u` 所属的区间
int findSpan(double u, int p, int n, const std::vector<double>& knots) {
    if (u >= knots[n]) return n - 1; // 处理 u == knots[n] 的情况

    int low = p;
    int high = n;
    int mid;

    while (high - low > 1) {
        mid = (high + low) >> 1;
        if (u < knots[mid] - eps) high = mid;
        else low = mid;
    }

    return low;
}

// 计算 B 样条基函数
std::vector<double> computeBasis(int span, double u, int p, const std::vector<double>& knots) {
    std::vector<double> basis(p + 1, 0.0);
    std::vector<double> left(p + 1, 0.0);
    std::vector<double> right(p + 1, 0.0);

    basis[0] = 1.0;

    for (int j = 1; j <= p; ++j) {
        left[j] = u - knots[span + 1 - j];
        right[j] = knots[span + j] - u;

        double saved = 0.0;
        for (int r = 0; r < j; ++r) {
            double basisMultiplier = basis[r] / (right[r + 1] + left[j - r]);
            basis[r] = saved + right[r + 1] * basisMultiplier;
            saved = left[j - r] * basisMultiplier;
        }
        basis[j] = saved;
    }

    return basis;
}

// 计算 B 样条曲线上 `u` 对应的点
Point evaluate_Bspline(const std::vector<double>& knots, const std::vector<Point>& control_points, int p, int n, double u) {
    if (u >= knots[n]) return control_points[n - 1]; // 确保 u 在合法范围内

    int span = findSpan(u, p, n, knots);
    std::vector<double> basis = computeBasis(span, u, p, knots);
    
    Point result;
    for (int i = 0; i <= p; ++i) {
        double basisMultiplier = basis[i];
        result = result + (control_points[span - p + i] * basisMultiplier);
    }

    return result;
}

void insert_knot(double u, int p, int r, int s) //r = p - s
{  
    int np = control_points.size()-1;
    int mp = np + p + 1;  
    int nq = np + r;  
    int k = std::distance(knots.begin(), std::upper_bound(knots.begin(), knots.end(), u)) - 1;
    printf("k = %d\n",k);
    int L;
    new_knots.resize(mp+r+1);
    new_control_points.resize(nq+1);
  
    // Load new knot vector  
    for (int i = 0; i <= k; i++) 
    {  
        new_knots[i] = knots[i];
    }  
    for (int i = 1; i <= r; i++) 
    {  
        new_knots[k + i] = u;  
    }  
    for (int i = k + 1; i <= mp; i++) 
    {  
        new_knots[i + r] = knots[i];  
    }  
  
    // Save unaltered control points  
    for (int i = 0; i <= k - p; i++) 
    {  
        new_control_points[i] = control_points[i];  
    }  
    for (int i = k - s; i <= np; i++) 
    {  
        new_control_points[i + r] = control_points[i];  
    }  

    Point Rw[p+1]; 
  
    for (int i = 0; i <= p - s; i++) 
    {  
        Rw[i] = control_points[k - p + i];  
    }  
  
    for (int j = 1; j <= r; j++) 
    { // Insert the knot r times  
        L = k - p + j;  
        for (int i = 0; i <= p - j - s; i++) 
        {  
            // Perform the alpha calculation  
            printf("u = %f knots[%d] = %f knots[%d] = %f\t",u,L+i,knots[L+i],i+k+1,knots[i+k+1]);
            printf("Rw[%d] = %f %f %f %f\t",i,Rw[i].x,Rw[i].y,Rw[i].z,Rw[i].w);
            printf("Rw[%d] = %f %f %f %f\n",i+1,Rw[i+1].x,Rw[i+1].y,Rw[i+1].z,Rw[i+1].w);
            double alpha = (u - knots[L + i]) / (knots[i + k + 1] - knots[L + i]);  
            Rw[i] = Rw[i + 1] * alpha + Rw[i] * (1.0 - alpha);  
            printf("alpha = %f\t",alpha);
            printf("Rw[%d] = %f %f %f %f\n",i,Rw[i].x,Rw[i].y,Rw[i].z,Rw[i].w);
        }  
        new_control_points[L] = Rw[0];  
        new_control_points[k + r - j - s] = Rw[p - j - s];  
    }  
  
    for (int i = L + 1; i < k - s; i++) 
    {  
        new_control_points[i] = Rw[i - L];  
    }  
}  

int FindSpan(int n, int p, double u, const std::vector<double>& knots)
{
    if (u >= knots[n+1])
        return n;
    
    int low = p;
    int high = n + 1;
    int mid = (low + high) / 2;
    
    while (u < knots[mid] || u >= knots[mid+1])
    {
        if (u < knots[mid])
            high = mid;
        else
            low = mid;
        
        mid = (low + high) / 2;
    }
    
    return mid;
}

void RefineKnotVectCurve(int n, int p, std::vector<double>& knots, std::vector<Point>& control_points, std::vector<double>& X, int r, std::vector<double>& new_knots, std::vector<Point>& new_control_points) 
{  
    int m = n + p + 1;  //n = control_points.size()-1
    int a = FindSpan(n, p, X[0], knots);  
    int b = FindSpan(n, p, X[r], knots) + 1;  //r = X.size()-1;

    new_knots.resize(m+r+2);
    new_control_points.resize(n+2+r);
      
    for (int j = 0; j <= a - p; j++) 
    {  
        new_control_points[j] = control_points[j]; 
    }  
      
    for (int j = b - 1; j <= n; j++) 
    {  
        new_control_points[j+r+1] = control_points[j];  
    }  
      
    for (int j = 0; j <= a; j++) 
    {  
        new_knots[j] = knots[j];  
    }  
      
    for (int j = b + p; j <= m; j++) 
    {  
        new_knots[j+r+1] = knots[j];  
    }  
      
    int i = b + p - 1;  
    int k = b + p + r;  
      
    for (int j = r; j >= 0; j--) 
    {  
        printf("j = %d\n",j);
        while (X[j] <= knots[i] && i > a) 
        {  
            new_control_points[k - p - 1] = control_points[i - p - 1];  
            new_knots[k] = knots[i];  
            k = k - 1;  
            i = i - 1;  
        }  
        printf("i = %d\t",i);
        printf("k = %d\n",k);
          
        new_control_points[k - p - 1] = new_control_points[k - p];  
          
        for (int l = 1; l <= p; l++) 
        {  
            int ind = k - p + l;  
            double alfa = new_knots[k + l] - X[j];  
            if (std::abs(alfa) == 0.0) 
            {  
                new_control_points[ind - 1] = new_control_points[ind];  
            } 
            else 
            {  
                alfa = alfa / (new_knots[k + l] - knots[i - p + l]);  
                new_control_points[ind - 1] =  new_control_points[ind - 1] * alfa + new_control_points[ind] * (1.0 - alfa);  
            }  
        }
        new_knots[k] = X[j];  
        k = k - 1;    
    }  
}

void DecomposeCurve(int n, int p, std::vector<double>& U, std::vector<Point>& Pw, std::vector<std::vector<Point>>& Qw) 
{  
    int m = n + p + 1;  
    int a = p;  
    int b = p + 1;  
    int nb = 0;  
    std::vector<double> alphas(p+1, 0.0);  
    Qw.resize(m-2*p);
    for(int i=0;i<Qw.size();i++)
    {
        Qw[i].resize(p+1);
    }
      
    for (int i = 0; i <= p; i++) 
    {  
        Qw[nb][i] = Pw[i];  
    }  
      
    while (b < m) 
    {  
        int i = b;  
        while (b < m && U[b + 1] == U[b]) 
        {  
            b++;  
        }  
        int mult = b - i + 1;  
          
        if (mult < p) 
        {  
            double numer = U[b] - U[a];  
            for (int j = p; j > mult; j--) 
            {  
                alphas[j - mult - 1] = numer / (U[a + j] - U[a]);  
            }  
            int r = p - mult;  
            for (int j = 1; j <= r; j++) 
            {  
                int save = r - j;  
                int s = mult + j;  
                  
                for (int k = p; k >= s; k--) 
                {  
                    double alpha = alphas[k - s];  
                    Qw[nb][k] = Qw[nb][k] * alpha + Qw[nb] [k-1] * (1.0-alpha); 
                }  
                if (b < m) 
                {  
                    Qw[nb+1][save] = Qw[nb][p];  
                }  
            }  
        }
        nb++;  
        if (b < m) 
        {  
            for (int i = p - mult; i <= p; i++) 
            {  
                Qw[nb][i] = Pw[b - p + i];  
            }  
            a = b;  
            b++;  
        }
    }  

}  

__global__ void pre_compute3(int n,int p,double *knot)
{
   // double iStart=cpuSecond();
  //  gpu_init(n,knot);
    int ix=threadIdx.x;
    int iy=threadIdx.y;
    int j=blockIdx.x;
    for(int i=1;i<=p;i++)
    {
        if(j+i>=n-1)
            return;
        double leftunder=knot[j+i]-knot[j];
        double rightunder=knot[j+i+1]-knot[j+1];
        double co_left_var,co_left_con,co_right_var,co_right_con;
        if(leftunder<eps)
        {
            co_left_con=0.0;
            co_left_var=0.0;
        }
        else
        {
          //  co_left_con=-knot[j]/leftunder;
            co_left_var=1.0/leftunder;
            co_left_con=-knot[j]*co_left_var;
        }
        if(rightunder<eps)
        {
            co_right_con=0.0;
            co_right_var=0.0;
        }
        else
        {
            // co_right_con=knot[j+i+1]/rightunder;
            co_right_var=-1.0/rightunder;
            co_right_con=-knot[j+i+1]*co_right_var;
        }
        d_N[j][i][j+ix][iy]=d_N[j][i-1][j+ix][iy]*co_left_con+d_N[j+1][i-1][j+ix][iy]*co_right_con;
        if(iy>0)
            d_N[j][i][j+ix][iy]+=d_N[j][i-1][j+ix][iy-1]*co_left_var+d_N[j+1][i-1][j+ix][iy-1]*co_right_var;
        __syncthreads();
        //printf("N %d %d %d %d : %f\n",j,i,j+ix,iy,d_N[j][i][j+ix][iy]);
    }
}

__global__ void bezier_compute(int n,int p,double *knot)
{
    int id=blockDim.x*blockIdx.x+threadIdx.x;
    if(id<n)
    {
        double u_i=knot[id];
        double u_i1=knot[id+1];
        if(u_i1-u_i<eps)
        {
            return;
        }
        double b=u_i1;
        double a=u_i;
        d_B[id][0][2][0]=1+a*a/(b-a)/(b-a)+2*a/(b-a);
        d_B[id][0][2][1]=(2+2*a)/(a-b);
        d_B[id][0][2][2]=1/(b-a)/(b-a);
        d_B[id][1][2][0]=2*a/(a-b)-2*a*a/(b-a)/(b-a);
        d_B[id][1][2][1]=2/(b-a)+4*a/(b-a)/(b-a);
        d_B[id][1][2][2]=-2/(b-a)/(b-a);
        d_B[id][2][2][0]=a*a/(b-a)/(b-a);
        d_B[id][2][2][1]=-2*a/(b-a)/(b-a);
        d_B[id][2][2][2]=1/(b-a)/(b-a);
        // printf("b=%f a=%f ",b,a);
        // printf("d_B= %f\n",d_B[id][2][2][2]);
    }
}

__global__ void trans(int m,int p,double *knots)
{
    int id=blockDim.x*blockIdx.x+threadIdx.x;
    if(id<p||id>m-p-1)
        return;
    double a=knots[id];
    double b=knots[id+1];
    if(b-a>eps)
    {
        double left[8];//p+1
        double right[8];//p+1
        left[0]=1;
        right[0]=1;
        for(int i=1;i<=p;i++)
        {
            left[i]=left[i-1]*(b-a);
            right[i]=right[i-1]*a;
        }
        for(int i=0;i<=p;i++)
        {
            for(int j=0;j<=p;j++)
            {
                double sum=0;
                for(int k=j;k<=p;k++)
                {
                    sum+=C[k][j]*right[k-j]*d_N[id-p+i][p][id][k];
                }
                d_transform_N[id][i][p][j]=left[j]*sum;
                // d_bigtransN[id*SIZE+i][j]=d_transform_N[id][i][p][j];
                d_bigtransN[id*SIZE*SIZE+i*SIZE+j]=d_transform_N[id][i][p][j];
            }
        }
    }
}


// __global__ void transformation(int m,int p,double *knots)
// {
//     int id=blockDim.x*blockIdx.x+threadIdx.x;
//     if(id<p||id>m-p-1)
//         return;
//     if(p==2)
//     {
//         double a=knots[id];
//         double b=knots[id+1];
//         if(b-a>eps)
//         {
//             d_transform_N[id][0][p][0] = d_N[id-2][p][id][0]    +   d_N[id-2][p][id][1]*a        +   d_N[id-2][p][id][2]*a*a;
//             d_transform_N[id][0][p][1] =                            d_N[id-2][p][id][1]*(b-a)    +   d_N[id-2][p][id][2]*a*(b-a)*2;
//             d_transform_N[id][0][p][2] =                                                             d_N[id-2][p][id][2]*(b-a)*(b-a);

//             d_transform_N[id][1][p][0] = d_N[id-1][p][id][0]    +   d_N[id-1][p][id][1]*a        +   d_N[id-1][p][id][2]*a*a;
//             d_transform_N[id][1][p][1] =                            d_N[id-1][p][id][1]*(b-a)    +   d_N[id-1][p][id][2]*a*(b-a)*2;
//             d_transform_N[id][1][p][2] =                                                             d_N[id-1][p][id][2]*(b-a)*(b-a);

//             d_transform_N[id][2][p][0] = d_N[id][p][id][0]    +     d_N[id][p][id][1]*a        +     d_N[id][p][id][2]*a*a;
//             d_transform_N[id][2][p][1] =                            d_N[id][p][id][1]*(b-a)    +     d_N[id][p][id][2]*a*(b-a)*2;
//             d_transform_N[id][2][p][2] =                                                             d_N[id][p][id][2]*(b-a)*(b-a);

//         }
//     }
//     if(p==3)
//     {
//         double a=knots[id];
//         double b=knots[id+1];
//         if(b-a>eps)
//         {
//             d_transform_N[id][0][p][0] = d_N[id-3][p][id][0]    +   d_N[id-3][p][id][1]*a        +   d_N[id-3][p][id][2]*a*a            +   d_N[id-3][p][id][3]*a*a*a;
//             d_transform_N[id][0][p][1] =                            d_N[id-3][p][id][1]*(b-a)    +   d_N[id-3][p][id][2]*a*(b-a)*2      +   d_N[id-3][p][id][3]*a*a*(b-a)*3;
//             d_transform_N[id][0][p][2] =                                                             d_N[id-3][p][id][2]*(b-a)*(b-a)    +   d_N[id-3][p][id][3]*a*(b-a)*(b-a)*3;
//             d_transform_N[id][0][p][3] =                                                                                                    d_N[id-3][p][id][3]*(b-a)*(b-a)*(b-a);

//             d_transform_N[id][1][p][0] = d_N[id-2][p][id][0]    +   d_N[id-2][p][id][1]*a        +   d_N[id-2][p][id][2]*a*a            +   d_N[id-2][p][id][3]*a*a*a;
//             d_transform_N[id][1][p][1] =                            d_N[id-2][p][id][1]*(b-a)    +   d_N[id-2][p][id][2]*a*(b-a)*2      +   d_N[id-2][p][id][3]*a*a*(b-a)*3;
//             d_transform_N[id][1][p][2] =                                                             d_N[id-2][p][id][2]*(b-a)*(b-a)    +   d_N[id-2][p][id][3]*a*(b-a)*(b-a)*3;
//             d_transform_N[id][1][p][3] =                                                                                                    d_N[id-2][p][id][3]*(b-a)*(b-a)*(b-a);

//             d_transform_N[id][2][p][0] = d_N[id-1][p][id][0]    +   d_N[id-1][p][id][1]*a        +   d_N[id-1][p][id][2]*a*a            +   d_N[id-1][p][id][3]*a*a*a;
//             d_transform_N[id][2][p][1] =                            d_N[id-1][p][id][1]*(b-a)    +   d_N[id-1][p][id][2]*a*(b-a)*2      +   d_N[id-1][p][id][3]*a*a*(b-a)*3;
//             d_transform_N[id][2][p][2] =                                                             d_N[id-1][p][id][2]*(b-a)*(b-a)    +   d_N[id-1][p][id][3]*a*(b-a)*(b-a)*3;
//             d_transform_N[id][2][p][3] =                                                                                                    d_N[id-1][p][id][3]*(b-a)*(b-a)*(b-a);

//             d_transform_N[id][3][p][0] = d_N[id][p][id][0]      +   d_N[id][p][id][1]*a          +   d_N[id][p][id][2]*a*a              +   d_N[id][p][id][3]*a*a*a;
//             d_transform_N[id][3][p][1] =                            d_N[id][p][id][1]*(b-a)      +   d_N[id][p][id][2]*a*(b-a)*2        +   d_N[id][p][id][3]*a*a*(b-a)*3;
//             d_transform_N[id][3][p][2] =                                                             d_N[id][p][id][2]*(b-a)*(b-a)      +   d_N[id][p][id][3]*a*(b-a)*(b-a)*3;
//             d_transform_N[id][3][p][3] =                                                                                                    d_N[id][p][id][3]*(b-a)*(b-a)*(b-a);
//         }
//     }

// }

// __global__ void compute_D(int m,int p,double *knots)
// {
//     int id=blockDim.x*blockIdx.x+threadIdx.x;
//     if(id<p||id>m-p-1)
//         return;
//     if(p==2)
//     {
//         d_D[id][0][p][0] = d_transform_N[id][0][p][0]; 
//         d_D[id][0][p][1] = d_transform_N[id][0][p][0] + 0.5 * d_transform_N[id][0][p][1];
//         d_D[id][0][p][2] = d_transform_N[id][0][p][0] +       d_transform_N[id][0][p][1] + d_transform_N[id][0][p][2];

//         d_D[id][1][p][0] = d_transform_N[id][1][p][0];
//         d_D[id][1][p][1] = d_transform_N[id][1][p][0] + 0.5 * d_transform_N[id][1][p][1];
//         d_D[id][1][p][2] = d_transform_N[id][1][p][0] +       d_transform_N[id][1][p][1] + d_transform_N[id][1][p][2];

//         d_D[id][2][p][0] = d_transform_N[id][2][p][0];
//         d_D[id][2][p][1] = d_transform_N[id][2][p][0] + 0.5 * d_transform_N[id][2][p][1];
//         d_D[id][2][p][2] = d_transform_N[id][2][p][0] +       d_transform_N[id][2][p][1] + d_transform_N[id][2][p][2];
//     }
//     if(p==3)
//     {
//         d_D[id][0][p][0] = d_transform_N[id][0][p][0];
//         d_D[id][0][p][1] = d_transform_N[id][0][p][0] +       d_transform_N[id][0][p][1] / 3.0;
//         d_D[id][0][p][2] = d_transform_N[id][0][p][0] + 2.0 * d_transform_N[id][0][p][1] / 3.0 + d_transform_N[id][0][p][2] / 3.0;
//         d_D[id][0][p][3] = d_transform_N[id][0][p][0] +       d_transform_N[id][0][p][1]       + d_transform_N[id][0][p][2]       + d_transform_N[id][0][p][3];

//         d_D[id][1][p][0] = d_transform_N[id][1][p][0];
//         d_D[id][1][p][1] = d_transform_N[id][1][p][0] +       d_transform_N[id][1][p][1] / 3.0;
//         d_D[id][1][p][2] = d_transform_N[id][1][p][0] + 2.0 * d_transform_N[id][1][p][1] / 3.0 + d_transform_N[id][1][p][2] / 3.0;
//         d_D[id][1][p][3] = d_transform_N[id][1][p][0] +       d_transform_N[id][1][p][1]       + d_transform_N[id][1][p][2]       + d_transform_N[id][1][p][3];

//         d_D[id][2][p][0] = d_transform_N[id][2][p][0];
//         d_D[id][2][p][1] = d_transform_N[id][2][p][0] +       d_transform_N[id][2][p][1] / 3.0;
//         d_D[id][2][p][2] = d_transform_N[id][2][p][0] + 2.0 * d_transform_N[id][2][p][1] / 3.0 + d_transform_N[id][2][p][2] / 3.0;
//         d_D[id][2][p][3] = d_transform_N[id][2][p][0] +       d_transform_N[id][2][p][1]       + d_transform_N[id][2][p][2]       + d_transform_N[id][2][p][3];

//         d_D[id][3][p][0] = d_transform_N[id][3][p][0];
//         d_D[id][3][p][1] = d_transform_N[id][3][p][0] +       d_transform_N[id][3][p][1] / 3.0;
//         d_D[id][3][p][2] = d_transform_N[id][3][p][0] + 2.0 * d_transform_N[id][3][p][1] / 3.0 + d_transform_N[id][3][p][2] / 3.0;
//         d_D[id][3][p][3] = d_transform_N[id][3][p][0] +       d_transform_N[id][3][p][1]       + d_transform_N[id][3][p][2]       + d_transform_N[id][3][p][3];
//     }
// }



__global__ void gpu_init2(int n,double *knots)//曲线
{
    int id=blockDim.x*blockIdx.x+threadIdx.x;
    if(id<n)
    {
        double u_i=knots[id];
        double u_i1=knots[id+1];
        if(u_i1-u_i>eps)
        {
            d_N[id][0][id][0]=1.0;
        }
    }
}

void getD(double *d_bigtransN_pointer,double *d_bigD_pointer,double *d_inverse_B_pointer,cublasHandle_t handle,double alpha,double beta)
{
    cublasDgemm(
        handle,
        CUBLAS_OP_T,
        CUBLAS_OP_N,
        SIZE,
        NUM * SIZE,
        SIZE,
        &alpha, 
        d_inverse_B_pointer, 
        SIZE, 
        d_bigtransN_pointer,
        SIZE, 
        &beta, 
        d_bigD_pointer, 
        SIZE
    );
}

__global__ void getQ(int m,int p,double *D,Point* control_points,double *d_big_Q)
{
    int id=blockDim.x*blockIdx.x+threadIdx.x;
    if(id<p||id>m-p-1)
        return;
    // double a=knots[id];
    // double b=knots[id+1];
    Eigen::Matrix <double,SIZE,SIZE> A;
    Eigen::Matrix <double,SIZE,4> B;
    Eigen::Matrix <double,SIZE,4> C;
    for(int i=0;i<SIZE;i++)
    {
        for(int j=0;j<SIZE;j++)
        {
            A(j,i)=D[id*SIZE*SIZE+i*SIZE+j];
        }
    }
    for(int i=0;i<=p;i++)
    {
        B(i,0)=control_points[id+i-p].x;
        B(i,1)=control_points[id+i-p].y;
        B(i,2)=control_points[id+i-p].z;
        B(i,3)=control_points[id+i-p].w;
    }

    C=A*B;//等价于A(i,j).transpose()
    for(int i=0;i<SIZE;i++)
    {
        for(int j=0;j<4;j++)
        {
            // if(id==5)
            // printf("C %d %d : %f\n",i,j,C(i,j));
            d_Q[id][i][j]=C(i,j);
            d_big_Q[id*SIZE*4+j*SIZE+i]=C(i,j);
        }
    }

    // if(id==5)
    // {
    //     printf("A\n");
    //     for(int i=0;i<SIZE;i++)
    //     {
    //         for(int j=0;j<SIZE;j++)
    //         {
    //             printf("%f ",A(i,j));
    //         }
    //         printf("\n");
    //     }
    //     printf("B\n");
    //     for(int i=0;i<SIZE;i++)
    //     {
    //         for(int j=0;j<4;j++)
    //         {
    //             printf("%f ",B(i,j));
    //         }
    //         printf("\n");
    //     }
    
    // }
}
// __global__ void get_delta_v1(int m,int n,double *d_big_Q,int p,int knot_size)
// {
//     int id=blockDim.x*blockIdx.x+threadIdx.x;
//     if(id<p||id>knot_size-p-1)
//         return;
//     Eigen::Matrix <double,1,SIZE> A;
//     Eigen::Matrix <double,1,SIZE> B;
//     Eigen::Matrix <double,1,SIZE> C;
//     Eigen::Matrix <double,4,SIZE> P;
//     Eigen::Matrix <double,4,1> P0;
//     Eigen::Matrix <double,4,1> P3;
//     Eigen::Matrix <double,4,1> deltaPn_1;
//     Eigen::Matrix <double,1,4> Co1;
//     Eigen::Matrix <double,1,4> Co2;
//     Eigen::Matrix <double,1,4> Co3;
//     Eigen::Matrix <double,1,4> Co4;
//     for(int i=0;i<SIZE;i++)
//     {
//         A(0,i)=G_m_n[2][i];
//         B(0,i)=G_m_n[1][i];
//     }
//     C=A-(B*G_m_m[2][1]/G_m_m[1][1]);
//     for(int i=0;i<4;i++)
//     {
//         for(int j=0;j<SIZE;j++)
//         {
//             P(i,j)=d_big_Q[id*SIZE*4+i*SIZE+j];
//         }
//     }
//     for(int i=0;i<4;i++)
//     {
//         P0(i,0)=P(i,0);
//         P3(i,0)=P(i,p);
//         deltaPn_1(i,0)=P(i,p)-P(i,p-1);
//     }
//     Co1=C*P.transpose();
//     Co2=(G_m_m[2][1]/G_m_m[1][1]*(G_m_m[1][0]+G_m_m[1][1])-(G_m_m[2][0]+G_m_m[2][1]))*P0;
//     Co3=(G_m_m[2][1]/G_m_m[1][1]*(G_m_m[1][2]+G_m_m[1][3])-(G_m_m[2][2]+G_m_m[2][3]))*P3;
//     Co4=Co1+Co2+Co3;
//     double coe=1.0*m/n/G_m_m[2][2]/(1-G_m_m[2][1]*G_m_m[1][2]/G_m_m[1][1]/G_m_m[2][2]);
//     d_delta[id][1][0]=coe*deltaPn_1(0,0);
//     for(int i=0;i<4;i++)
//     {
//         // printf("fabs(Co4(0,i))=%f\n",fabs(Co4(0,i)));
//         if(fabs(Co4(0,i))>eps)
//         {
//             // printf("yes\n");
//             d_delta[id][1][i]=coe*Co4(0,i)/deltaPn_1(i,0);
//         }
//         else
//             d_delta[id][1][i]=0;
//     }
//     if(id==5)
//     {
//         printf("C\n");
//         for(int i=0;i<SIZE;i++)
//         {
//             printf("%f ",C(0,i));
//         }
//         printf("\n");
//         printf("P\n");
//         for(int i=0;i<4;i++)
//         {
//             for(int j=0;j<SIZE;j++)
//             {
//                 printf("%f ",P(i,j));
//             }
//             printf("\n");
//         }
//         printf("P0\n");
//         for(int i=0;i<4;i++)
//         {
//             printf("%f ",P0(i,0));
//         }
//         printf("\n");
//         printf("P3\n");
//         for(int i=0;i<4;i++)
//         {
//             printf("%f ",P3(i,0));
//         }
//         printf("\n");
//         printf("Co1\n");
//         for(int i=0;i<4;i++)
//         {
//             printf("%f ",Co1(0,i));
//         }
//         printf("\n");
//         printf("Co2\n");
//         for(int i=0;i<4;i++)
//         {
//             printf("%f ",Co2(0,i));
//         }
//         printf("\n");
//         printf("Co3\n");
//         for(int i=0;i<4;i++)
//         {
//             printf("%f ",Co3(0,i));
//         }
//         printf("\n");
//         printf("Co4\n");
//         for(int i=0;i<4;i++)
//         {
//             printf("%f ",Co4(0,i));
//         }
//         printf("\n");
//         printf("deltaPn_1\n");
//         for(int i=0;i<4;i++)
//         {
//             printf("%f ",deltaPn_1(i,0));
//         }
//         printf("\n");
//         printf("coe=%f\n",coe);
//         printf("delta\n");
//         for(int i=0;i<4;i++)
//         {
//             printf("%f ",d_delta[id][1][i]);
//         }
//         printf("\n");
        
//     }

// }

__global__ void get_delta(int m,int n,double *d_big_Q,int p,int knot_size)
{
    int id=blockDim.x*blockIdx.x+threadIdx.x;
    if(id<p||id>knot_size-p-1)
        return;
    Eigen::Matrix <double,1,SIZE> A;
    Eigen::Matrix <double,1,SIZE> B;
    Eigen::Matrix <double,1,SIZE> C;
    Eigen::Matrix <double,4,SIZE> P;
    Eigen::Matrix <double,4,1> P0;
    Eigen::Matrix <double,4,1> P3;
    Eigen::Matrix <double,4,1> deltaPn_1;
    Eigen::Matrix <double,4,1> deltaPn_0;
    Eigen::Matrix <double,1,4> Co1;
    Eigen::Matrix <double,1,4> Co2;
    Eigen::Matrix <double,1,4> Co3;
    Eigen::Matrix <double,1,4> Co4;
    Eigen::Matrix <double,1,4> Do1;
    Eigen::Matrix <double,1,4> Do2;
    Eigen::Matrix <double,1,4> Do3;
    Eigen::Matrix <double,1,4> Do4;
    Eigen::Matrix <double,1,4> Do5;
    for(int i=0;i<SIZE;i++)
    {
        A(0,i)=G_m_n[2][i];
        B(0,i)=G_m_n[1][i];
    }
    C=A-(B*G_m_m[2][1]/G_m_m[1][1]);
    for(int i=0;i<4;i++)
    {
        for(int j=0;j<SIZE;j++)
        {
            P(i,j)=d_big_Q[id*SIZE*4+i*SIZE+j];
        }
    }
    for(int i=0;i<4;i++)
    {
        P0(i,0)=P(i,0);
        P3(i,0)=P(i,p);
        deltaPn_1(i,0)=P(i,p)-P(i,p-1);
        deltaPn_0(i,0)=P(i,1)-P(i,0);
    }
    Co1=C*P.transpose();
    Co2=(G_m_m[2][1]/G_m_m[1][1]*(G_m_m[1][0]+G_m_m[1][1])-(G_m_m[2][0]+G_m_m[2][1]))*P0;
    Co3=(G_m_m[2][1]/G_m_m[1][1]*(G_m_m[1][2]+G_m_m[1][3])-(G_m_m[2][2]+G_m_m[2][3]))*P3;
    Co4=Co1+Co2+Co3;
    double coe=1.0*m/n/G_m_m[2][2]/(G_m_m[2][1]*G_m_m[1][2]/G_m_m[1][1]/G_m_m[2][2]-1);
    d_delta[id][1][0]=coe*deltaPn_1(0,0);

    Do1=B*P.transpose();
    Do2=(G_m_m[1][0]+G_m_m[1][1])*P0;
    Do3=(G_m_m[1][2]+G_m_m[1][3])*P3;
    Do4=G_m_m[1][2]*n/m*deltaPn_1.transpose();
    for(int i=0;i<4;i++)
    {
        // printf("fabs(Co4(0,i))=%f\n",fabs(Co4(0,i)));
        if(fabs(deltaPn_1(i,0))>eps)
        {
            // printf("yes\n");
            d_delta[id][1][i]=coe*Co4(0,i)/deltaPn_1(i,0);
        }
        else
            d_delta[id][1][i]=1;
        Do5(0,i)=Do1(0,i)-Do2(0,i)-Do3(0,i)+Do4(0,i)*d_delta[id][1][i];
    }
    for(int i=0;i<4;i++)
    {
        if(fabs(deltaPn_0(i,0))>eps)
            d_delta[id][0][i]=Do5(0,i)*m/n/G_m_m[1][1]/deltaPn_0(i,0);
        else
            d_delta[id][0][i]=1;
    }
    for(int i=0;i<4;i++)
    {
        d_cubic_Q[id][0][i]=P(i,0);
        d_cubic_Q[id][3][i]=P(i,p);
        d_cubic_Q[id][1][i]=P(i,0)+1.0*n/m*deltaPn_0(i,0)*d_delta[id][0][i];
        d_cubic_Q[id][2][i]=P(i,p)-1.0*n/m*deltaPn_1(i,0)*d_delta[id][1][i];
        // printf("d_cubic_Q[%d][%d][%d]=%f\n",id,0,i,d_cubic_Q[id][0][i]);
        // printf("d_cubic_Q[%d][%d][%d]=%f\n",id,1,i,d_cubic_Q[id][1][i]);
        // printf("d_cubic_Q[%d][%d][%d]=%f\n",id,2,i,d_cubic_Q[id][2][i]);
        // printf("d_cubic_Q[%d][%d][%d]=%f\n",id,3,i,d_cubic_Q[id][3][i]);
    }

    // if(id==5)
    // {
    //     // printf("C\n");
    //     // for(int i=0;i<SIZE;i++)
    //     // {
    //     //     printf("%f ",C(0,i));
    //     // }
    //     // printf("\n");
    //     // printf("P\n");
    //     // for(int i=0;i<4;i++)
    //     // {
    //     //     for(int j=0;j<SIZE;j++)
    //     //     {
    //     //         printf("%f ",P(i,j));
    //     //     }
    //     //     printf("\n");
    //     // }
    //     // printf("P0\n");
    //     // for(int i=0;i<4;i++)
    //     // {
    //     //     printf("%f ",P0(i,0));
    //     // }
    //     // printf("\n");
    //     // printf("P3\n");
    //     // for(int i=0;i<4;i++)
    //     // {
    //     //     printf("%f ",P3(i,0));
    //     // }
    //     // printf("\n");
    //     // printf("Co1\n");
    //     // for(int i=0;i<4;i++)
    //     // {
    //     //     printf("%f ",Co1(0,i));
    //     // }
    //     // printf("\n");
    //     // printf("Co2\n");
    //     // for(int i=0;i<4;i++)
    //     // {
    //     //     printf("%f ",Co2(0,i));
    //     // }
    //     // printf("\n");
    //     // printf("Co3\n");
    //     // for(int i=0;i<4;i++)
    //     // {
    //     //     printf("%f ",Co3(0,i));
    //     // }
    //     // printf("\n");
    //     // printf("Co4\n");
    //     // for(int i=0;i<4;i++)
    //     // {
    //     //     printf("%f ",Co4(0,i));
    //     // }
    //     // printf("\n");
    //     // printf("deltaPn_1\n");
    //     // for(int i=0;i<4;i++)
    //     // {
    //     //     printf("%f ",deltaPn_1(i,0));
    //     // }
    //     // printf("\n");
    //     // printf("deltaPn_0\n");
    //     // for(int i=0;i<4;i++)
    //     // {
    //     //     printf("%f ",deltaPn_0(i,0));
    //     // }
    //     // printf("\n");
    //     // printf("coe=%f\n",coe);
    //     printf("delta 1\n");
    //     for(int i=0;i<4;i++)
    //     {
    //         printf("%f ",d_delta[id][1][i]);
    //     }
    //     printf("\n");
    //     printf("delta 0\n");
    //     for(int i=0;i<4;i++)
    //     {
    //         printf("%f ",d_delta[id][0][i]);
    //     }
    //     printf("\n");
    // }

}
// __global__ void compute_Q(int m,int p,Point *control_points)
// {
//     int id=blockDim.x*blockIdx.x+threadIdx.x;
//     if(id<p||id>m-p-1)
//         return;
//     if(p==2)
//     {
//         d_Q[id][0][0]=d_D[id][0][2][0]*control_points[id-2].x+d_D[id][1][2][0]*control_points[id-1].x+d_D[id][2][2][0]*control_points[id].x;
//         d_Q[id][0][1]=d_D[id][0][2][0]*control_points[id-2].y+d_D[id][1][2][0]*control_points[id-1].y+d_D[id][2][2][0]*control_points[id].y;
//         d_Q[id][0][2]=d_D[id][0][2][0]*control_points[id-2].z+d_D[id][1][2][0]*control_points[id-1].z+d_D[id][2][2][0]*control_points[id].z;
//         d_Q[id][0][3]=d_D[id][0][2][0]*control_points[id-2].w+d_D[id][1][2][0]*control_points[id-1].w+d_D[id][2][2][0]*control_points[id].w;

//         d_Q[id][1][0]=d_D[id][0][2][1]*control_points[id-2].x+d_D[id][1][2][1]*control_points[id-1].x+d_D[id][2][2][1]*control_points[id].x;
//         d_Q[id][1][1]=d_D[id][0][2][1]*control_points[id-2].y+d_D[id][1][2][1]*control_points[id-1].y+d_D[id][2][2][1]*control_points[id].y;
//         d_Q[id][1][2]=d_D[id][0][2][1]*control_points[id-2].z+d_D[id][1][2][1]*control_points[id-1].z+d_D[id][2][2][1]*control_points[id].z;
//         d_Q[id][1][3]=d_D[id][0][2][1]*control_points[id-2].w+d_D[id][1][2][1]*control_points[id-1].w+d_D[id][2][2][1]*control_points[id].w;

//         d_Q[id][2][0]=d_D[id][0][2][2]*control_points[id-2].x+d_D[id][1][2][2]*control_points[id-1].x+d_D[id][2][2][2]*control_points[id].x;
//         d_Q[id][2][1]=d_D[id][0][2][2]*control_points[id-2].y+d_D[id][1][2][2]*control_points[id-1].y+d_D[id][2][2][2]*control_points[id].y;
//         d_Q[id][2][2]=d_D[id][0][2][2]*control_points[id-2].z+d_D[id][1][2][2]*control_points[id-1].z+d_D[id][2][2][2]*control_points[id].z;
//         d_Q[id][2][3]=d_D[id][0][2][2]*control_points[id-2].w+d_D[id][1][2][2]*control_points[id-1].w+d_D[id][2][2][2]*control_points[id].w;
//     }
//     if(p==3)
//     {
//         d_Q[id][0][0] = d_D[id][0][3][0] * control_points[id-3].x + d_D[id][1][3][0] * control_points[id-2].x + d_D[id][2][3][0] * control_points[id-1].x + d_D[id][3][3][0] * control_points[id].x;
//         d_Q[id][0][1] = d_D[id][0][3][0] * control_points[id-3].y + d_D[id][1][3][0] * control_points[id-2].y + d_D[id][2][3][0] * control_points[id-1].y + d_D[id][3][3][0] * control_points[id].y;
//         d_Q[id][0][2] = d_D[id][0][3][0] * control_points[id-3].z + d_D[id][1][3][0] * control_points[id-2].z + d_D[id][2][3][0] * control_points[id-1].z + d_D[id][3][3][0] * control_points[id].z;
//         d_Q[id][0][3] = d_D[id][0][3][0] * control_points[id-3].w + d_D[id][1][3][0] * control_points[id-2].w + d_D[id][2][3][0] * control_points[id-1].w + d_D[id][3][3][0] * control_points[id].w;

//         d_Q[id][1][0] = d_D[id][0][3][1] * control_points[id-3].x + d_D[id][1][3][1] * control_points[id-2].x + d_D[id][2][3][1] * control_points[id-1].x + d_D[id][3][3][1] * control_points[id].x;
//         d_Q[id][1][1] = d_D[id][0][3][1] * control_points[id-3].y + d_D[id][1][3][1] * control_points[id-2].y + d_D[id][2][3][1] * control_points[id-1].y + d_D[id][3][3][1] * control_points[id].y;
//         d_Q[id][1][2] = d_D[id][0][3][1] * control_points[id-3].z + d_D[id][1][3][1] * control_points[id-2].z + d_D[id][2][3][1] * control_points[id-1].z + d_D[id][3][3][1] * control_points[id].z;
//         d_Q[id][1][3] = d_D[id][0][3][1] * control_points[id-3].w + d_D[id][1][3][1] * control_points[id-2].w + d_D[id][2][3][1] * control_points[id-1].w + d_D[id][3][3][1] * control_points[id].w;
        
//         d_Q[id][2][0] = d_D[id][0][3][2] * control_points[id-3].x + d_D[id][1][3][2] * control_points[id-2].x + d_D[id][2][3][2] * control_points[id-1].x + d_D[id][3][3][2] * control_points[id].x;
//         d_Q[id][2][1] = d_D[id][0][3][2] * control_points[id-3].y + d_D[id][1][3][2] * control_points[id-2].y + d_D[id][2][3][2] * control_points[id-1].y + d_D[id][3][3][2] * control_points[id].y;
//         d_Q[id][2][2] = d_D[id][0][3][2] * control_points[id-3].z + d_D[id][1][3][2] * control_points[id-2].z + d_D[id][2][3][2] * control_points[id-1].z + d_D[id][3][3][2] * control_points[id].z;
//         d_Q[id][2][3] = d_D[id][0][3][2] * control_points[id-3].w + d_D[id][1][3][2] * control_points[id-2].w + d_D[id][2][3][2] * control_points[id-1].w + d_D[id][3][3][2] * control_points[id].w;

//         d_Q[id][3][0] = d_D[id][0][3][3] * control_points[id-3].x + d_D[id][1][3][3] * control_points[id-2].x + d_D[id][2][3][3] * control_points[id-1].x + d_D[id][3][3][3] * control_points[id].x;
//         d_Q[id][3][1] = d_D[id][0][3][3] * control_points[id-3].y + d_D[id][1][3][3] * control_points[id-2].y + d_D[id][2][3][3] * control_points[id-1].y + d_D[id][3][3][3] * control_points[id].y;
//         d_Q[id][3][2] = d_D[id][0][3][3] * control_points[id-3].z + d_D[id][1][3][3] * control_points[id-2].z + d_D[id][2][3][3] * control_points[id-1].z + d_D[id][3][3][3] * control_points[id].z;
//         d_Q[id][3][3] = d_D[id][0][3][3] * control_points[id-3].w + d_D[id][1][3][3] * control_points[id-2].w + d_D[id][2][3][3] * control_points[id-1].w + d_D[id][3][3][3] * control_points[id].w;
//     }
// }

__device__ double calculate_point_error(Point a,Point b)
{
    return ((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y)+(a.z-b.z)*(a.z-b.z)+(a.w-b.w)*(a.w-b.w));
}

__device__ double calculate_error(int n, int m, Point *P_n, Point *P_m,int id)
{
    double error=0;
    double N_n[SIZE][4];
    double N_m[4][4];
    double N_minus[SIZE][4];
    double N_new[2*SIZE-1][4];
    double coe=1;
    double temp[4];
    for(int j=0;j<=n;j++)
    {
        temp[0]=0;
        temp[1]=0;
        temp[2]=0;
        temp[3]=0;
        for(int i=0;i<=j;i++)
        {
            if((j-i)%2==0)
                coe=1;
            else
                coe=-1;
            temp[0]+=coe*C[n][i]*C[n-i][j-i]*P_n[i].x;
            temp[1]+=coe*C[n][i]*C[n-i][j-i]*P_n[i].y;
            temp[2]+=coe*C[n][i]*C[n-i][j-i]*P_n[i].z;
            temp[3]+=coe*C[n][i]*C[n-i][j-i]*P_n[i].w;
        }
        N_n[j][0]=temp[0];
        N_n[j][1]=temp[1];
        N_n[j][2]=temp[2];
        N_n[j][3]=temp[3];
    }
    for(int j=0;j<=m;j++)
    {
        temp[0]=0;
        temp[1]=0;
        temp[2]=0;
        temp[3]=0;
        for(int i=0;i<=j;i++)
        {
            if((j-i)%2==0)
                coe=1;
            else
                coe=-1;
            temp[0]+=coe*C[m][i]*C[m-i][j-i]*P_m[i].x;
            temp[1]+=coe*C[m][i]*C[m-i][j-i]*P_m[i].y;
            temp[2]+=coe*C[m][i]*C[m-i][j-i]*P_m[i].z;
            temp[3]+=coe*C[m][i]*C[m-i][j-i]*P_m[i].w;
        }
        N_m[j][0]=temp[0];
        N_m[j][1]=temp[1];
        N_m[j][2]=temp[2];
        N_m[j][3]=temp[3];
    }
    for(int i=0;i<=m;i++)
    {
        N_minus[i][0]=N_n[i][0]-N_m[i][0];
        N_minus[i][1]=N_n[i][1]-N_m[i][1];
        N_minus[i][2]=N_n[i][2]-N_m[i][2];
        N_minus[i][3]=N_n[i][3]-N_m[i][3];
        // printf("N_minus[%d][0]=%f\n",i,N_minus[i][0]);
    }
    for(int i=m+1;i<=n;i++)
    {
        N_minus[i][0]=N_n[i][0];
        N_minus[i][1]=N_n[i][1];
        N_minus[i][2]=N_n[i][2];
        N_minus[i][3]=N_n[i][3];
    }
    for(int i=0;i<=2*n;i++)
    {
        double temx=0;
        double temy=0;
        double temz=0;
        double temw=0;
        for(int j=max(i-n,0);j<=min(n,i);j++)
        {
            temx+=N_minus[j][0]*N_minus[i-j][0];
            temy+=N_minus[j][1]*N_minus[i-j][1];
            temz+=N_minus[j][2]*N_minus[i-j][2];
            temw+=N_minus[j][3]*N_minus[i-j][3];
        }
        N_new[i][0]=temx;
        N_new[i][1]=temy;
        N_new[i][2]=temz;
        N_new[i][3]=temw;
    }
    for(int i=0;i<=2*n;i++)
    {
        error+=(N_new[i][0]+N_new[i][1]+N_new[i][2]+N_new[i][3])/(i+1);
    }
    // if(id==7)
    // {
    //     for(int i=0;i<=n;i++)
    //     {
    //         printf("N_minus[%d][0]=%f\n",i,N_minus[i][0]);
    //         printf("N_minus[%d][1]=%f\n",i,N_minus[i][1]);
    //         printf("N_minus[%d][2]=%f\n",i,N_minus[i][2]);
    //         printf("N_minus[%d][3]=%f\n",i,N_minus[i][3]);
    //     }
    //     for(int i=0;i<=2*n;i++)
    //     {
    //         printf("N_new[%d][0]=%f\n",i,N_new[i][0]);
    //         printf("N_new[%d][1]=%f\n",i,N_new[i][1]);
    //         printf("N_new[%d][2]=%f\n",i,N_new[i][2]);
    //         printf("N_new[%d][3]=%f\n",i,N_new[i][3]);
    //     }
    //     printf("error=%f\n",error);
    //     printf("%f\n",N_minus[4][0]*N_minus[4][0]);
    // }
    return error;
}


__device__ double calculate_L1_error(int n, int m, Point *P_n, Point *P_m,int id,double t)
{
    double error=0;
    double N_n[SIZE][4];
    double N_m[4][4];
    double N_minus[SIZE][4];
    double N_new[2*SIZE-1][4];
    double coe=1;
    double temp[4];
    for(int j=0;j<=n;j++)
    {
        temp[0]=0;
        temp[1]=0;
        temp[2]=0;
        temp[3]=0;
        for(int i=0;i<=j;i++)
        {
            if((j-i)%2==0)
                coe=1;
            else
                coe=-1;
            temp[0]+=coe*C[n][i]*C[n-i][j-i]*P_n[i].x;
            temp[1]+=coe*C[n][i]*C[n-i][j-i]*P_n[i].y;
            temp[2]+=coe*C[n][i]*C[n-i][j-i]*P_n[i].z;
            temp[3]+=coe*C[n][i]*C[n-i][j-i]*P_n[i].w;
        }
        N_n[j][0]=temp[0];
        N_n[j][1]=temp[1];
        N_n[j][2]=temp[2];
        N_n[j][3]=temp[3];
    }
    for(int j=0;j<=m;j++)
    {
        temp[0]=0;
        temp[1]=0;
        temp[2]=0;
        temp[3]=0;
        for(int i=0;i<=j;i++)
        {
            if((j-i)%2==0)
                coe=1;
            else
                coe=-1;
            temp[0]+=coe*C[m][i]*C[m-i][j-i]*P_m[i].x;
            temp[1]+=coe*C[m][i]*C[m-i][j-i]*P_m[i].y;
            temp[2]+=coe*C[m][i]*C[m-i][j-i]*P_m[i].z;
            temp[3]+=coe*C[m][i]*C[m-i][j-i]*P_m[i].w;
        }
        N_m[j][0]=temp[0];
        N_m[j][1]=temp[1];
        N_m[j][2]=temp[2];
        N_m[j][3]=temp[3];
    }
    for(int i=0;i<=m;i++)
    {
        N_minus[i][0]=N_n[i][0]-N_m[i][0];
        N_minus[i][1]=N_n[i][1]-N_m[i][1];
        N_minus[i][2]=N_n[i][2]-N_m[i][2];
        N_minus[i][3]=N_n[i][3]-N_m[i][3];
        // printf("N_minus[%d][0]=%f\n",i,N_minus[i][0]);
    }
    for(int i=m+1;i<=n;i++)
    {
        N_minus[i][0]=N_n[i][0];
        N_minus[i][1]=N_n[i][1];
        N_minus[i][2]=N_n[i][2];
        N_minus[i][3]=N_n[i][3];
    }
    for(int i=0;i<=2*n;i++)
    {
        double temx=0;
        double temy=0;
        double temz=0;
        double temw=0;
        for(int j=max(i-n,0);j<=min(n,i);j++)
        {
            temx+=N_minus[j][0]*N_minus[i-j][0];
            temy+=N_minus[j][1]*N_minus[i-j][1];
            temz+=N_minus[j][2]*N_minus[i-j][2];
            temw+=N_minus[j][3]*N_minus[i-j][3];
        }
        N_new[i][0]=temx;
        N_new[i][1]=temy;
        N_new[i][2]=temz;
        N_new[i][3]=temw;
    }
    double lr=1;
    for(int i=0;i<=2*n;i++)
    {
        error+=(N_new[i][0]+N_new[i][1]+N_new[i][2]+N_new[i][3])*lr;
        lr*=t;
    }
    double coef[SIZE*2-1];
    for(int i=0;i<SIZE*2-1;i++)
    {
        coef[i]=N_new[i][0]+N_new[i][1]+N_new[i][2]+N_new[i][3];
    }
    Eigen::VectorXd coeffs = definePolynomial(2*SIZE-1,coef);
    Eigen::VectorXd derivativeCoeffs = derivative(coeffs);
    double initial_guess = 0.5;
    // 使用牛顿迭代法求解多项式根
    double root = newtonIteration(derivativeCoeffs, initial_guess);
    double max_l1_error=evaluatePolynomial(coeffs,root);
    printf("id=%d root=%f max_l1_error=%f\n",id,root,max_l1_error);

    // // 求导数的根（极值点）
    // Eigen::PolynomialSolver<double, Eigen::Dynamic> solver;
    // solver.compute(derivativeCoeffs);
    // Eigen::VectorXd roots = solver.roots().real();

    // 定义区间 [a, b]
    double a = 0.0;
    double b = 1.0;

    // // 计算边界值
    // double f_a = evaluatePolynomial(coeffs, a);
    // double f_b = evaluatePolynomial(coeffs, b);

    // // 找到落在区间内的极值点并计算对应的多项式值
    // std::vector<double> candidates = {f_a, f_b};
    // for (int i = 0; i < roots.size(); ++i) {
    //     if (roots[i] >= a && roots[i] <= b) {
    //         candidates.push_back(evaluatePolynomial(coeffs, roots[i]));
    //     }
    // }
    // if(id==7)
    // {
    //     for(int i=0;i<=n;i++)
    //     {
    //         printf("N_minus[%d][0]=%f\n",i,N_minus[i][0]);
    //         printf("N_minus[%d][1]=%f\n",i,N_minus[i][1]);
    //         printf("N_minus[%d][2]=%f\n",i,N_minus[i][2]);
    //         printf("N_minus[%d][3]=%f\n",i,N_minus[i][3]);
    //     }
    //     for(int i=0;i<=2*n;i++)
    //     {
    //         printf("N_new[%d][0]=%f\n",i,N_new[i][0]);
    //         printf("N_new[%d][1]=%f\n",i,N_new[i][1]);
    //         printf("N_new[%d][2]=%f\n",i,N_new[i][2]);
    //         printf("N_new[%d][3]=%f\n",i,N_new[i][3]);
    //     }
    //     printf("error=%f\n",error);
    //     printf("%f\n",N_minus[4][0]*N_minus[4][0]);
    // }
    return error;
}

__device__ void calculate_max_L1_error(int n, int m, Point *P_n, Point *P_m,int id,double *ans)
{
    double error=0;
    double N_n[SIZE][4];
    double N_m[4][4];
    double N_minus[SIZE][4];
    double N_new[2*SIZE-1][4];
    double coe=1;
    double temp[4];
    for(int j=0;j<=n;j++)
    {
        temp[0]=0;
        temp[1]=0;
        temp[2]=0;
        temp[3]=0;
        for(int i=0;i<=j;i++)
        {
            if((j-i)%2==0)
                coe=1;
            else
                coe=-1;
            temp[0]+=coe*C[n][i]*C[n-i][j-i]*P_n[i].x;
            temp[1]+=coe*C[n][i]*C[n-i][j-i]*P_n[i].y;
            temp[2]+=coe*C[n][i]*C[n-i][j-i]*P_n[i].z;
            temp[3]+=coe*C[n][i]*C[n-i][j-i]*P_n[i].w;
        }
        N_n[j][0]=temp[0];
        N_n[j][1]=temp[1];
        N_n[j][2]=temp[2];
        N_n[j][3]=temp[3];
    }
    for(int j=0;j<=m;j++)
    {
        temp[0]=0;
        temp[1]=0;
        temp[2]=0;
        temp[3]=0;
        for(int i=0;i<=j;i++)
        {
            if((j-i)%2==0)
                coe=1;
            else
                coe=-1;
            temp[0]+=coe*C[m][i]*C[m-i][j-i]*P_m[i].x;
            temp[1]+=coe*C[m][i]*C[m-i][j-i]*P_m[i].y;
            temp[2]+=coe*C[m][i]*C[m-i][j-i]*P_m[i].z;
            temp[3]+=coe*C[m][i]*C[m-i][j-i]*P_m[i].w;
        }
        N_m[j][0]=temp[0];
        N_m[j][1]=temp[1];
        N_m[j][2]=temp[2];
        N_m[j][3]=temp[3];
    }
    for(int i=0;i<=m;i++)
    {
        N_minus[i][0]=N_n[i][0]-N_m[i][0];
        N_minus[i][1]=N_n[i][1]-N_m[i][1];
        N_minus[i][2]=N_n[i][2]-N_m[i][2];
        N_minus[i][3]=N_n[i][3]-N_m[i][3];
        // printf("N_minus[%d][0]=%f\n",i,N_minus[i][0]);
    }
    for(int i=m+1;i<=n;i++)
    {
        N_minus[i][0]=N_n[i][0];
        N_minus[i][1]=N_n[i][1];
        N_minus[i][2]=N_n[i][2];
        N_minus[i][3]=N_n[i][3];
    }
    for(int i=0;i<=2*n;i++)
    {
        double temx=0;
        double temy=0;
        double temz=0;
        double temw=0;
        for(int j=max(i-n,0);j<=min(n,i);j++)
        {
            temx+=N_minus[j][0]*N_minus[i-j][0];
            temy+=N_minus[j][1]*N_minus[i-j][1];
            temz+=N_minus[j][2]*N_minus[i-j][2];
            temw+=N_minus[j][3]*N_minus[i-j][3];
        }
        N_new[i][0]=temx;
        N_new[i][1]=temy;
        N_new[i][2]=temz;
        N_new[i][3]=temw;
    }
    double coef[SIZE*2-1];
    for(int i=0;i<SIZE*2-1;i++)
    {
        coef[i]=N_new[i][0]+N_new[i][1]+N_new[i][2]+N_new[i][3];
    }
    Eigen::VectorXd coeffs = definePolynomial(2*SIZE-1,coef);
    Eigen::VectorXd derivativeCoeffs = derivative(coeffs);
    
    //利用牛顿迭代求方程的根
    double initial_guess = 0.1;
    double root=0;
    double max_l1_error=0;
    for(int i=1;i<=9;i++)
    {
        double temp_root=newtonIteration(derivativeCoeffs, initial_guess*i);
        double temp_l1_error=evaluatePolynomial(coeffs,temp_root);
        if(temp_l1_error>max_l1_error)
        {
            max_l1_error=temp_l1_error;
            root=temp_root;
        }
    }
    // double root = newtonIteration(derivativeCoeffs, initial_guess);
    // double max_l1_error=evaluatePolynomial(coeffs,root);
    ans[0]=root;
    ans[1]=max_l1_error;
}

__device__ double calculate_trans_L1_error(int n, int m, Point *P_n, Point *P_m,double a,double b, int id,double u)
{
    double error=0;
    double N_n[SIZE][4];
    double N_tra[SIZE][4];
    double N_m[4][4];
    double N_minus[SIZE][4];
    double N_new[2*SIZE-1][4];
    double coe=1;
    double temp[4];
    double M[SIZE][SIZE];
    double left[SIZE],right[SIZE];
    left[0]=1;
    right[0]=1;
    for(int i=1;i<=n;i++)
    {
        left[i]=left[i-1]*a;   
        right[i]=right[i-1]*(b-a);
    }
    for(int i=0;i<=n;i++)
    {
        for(int j=0;j<=n;j++)
        {
            if(j<i)
            {
                M[i][j]=0;
            }
            else
            {
                M[i][j]=C[j][i]*left[j-i]*right[i];
            }
        }
    }
    // if(id==8)
    // {
    //     for(int i=0;i<=n;i++)
    //     {
    //         for(int j=0;j<=n;j++)
    //         {
    //             printf("M[%d][%d]=%f ",i,j,M[i][j]);
    //         }
    //         printf("\n");
    //     }
    // }
    for(int j=0;j<=n;j++)
    {
        temp[0]=0;
        temp[1]=0;
        temp[2]=0;
        temp[3]=0;
        for(int i=0;i<=j;i++)
        {
            if((j-i)%2==0)
                coe=1;
            else
                coe=-1;
            temp[0]+=coe*C[n][i]*C[n-i][j-i]*P_n[i].x;
            temp[1]+=coe*C[n][i]*C[n-i][j-i]*P_n[i].y;
            temp[2]+=coe*C[n][i]*C[n-i][j-i]*P_n[i].z;
            temp[3]+=coe*C[n][i]*C[n-i][j-i]*P_n[i].w;
        }
        N_n[j][0]=temp[0];
        N_n[j][1]=temp[1];
        N_n[j][2]=temp[2];
        N_n[j][3]=temp[3];
    }
    for(int j=0;j<=m;j++)
    {
        temp[0]=0;
        temp[1]=0;
        temp[2]=0;
        temp[3]=0;
        for(int i=0;i<=j;i++)
        {
            if((j-i)%2==0)
                coe=1;
            else
                coe=-1;
            temp[0]+=coe*C[m][i]*C[m-i][j-i]*P_m[i].x;
            temp[1]+=coe*C[m][i]*C[m-i][j-i]*P_m[i].y;
            temp[2]+=coe*C[m][i]*C[m-i][j-i]*P_m[i].z;
            temp[3]+=coe*C[m][i]*C[m-i][j-i]*P_m[i].w;
        }
        N_m[j][0]=temp[0];
        N_m[j][1]=temp[1];
        N_m[j][2]=temp[2];
        N_m[j][3]=temp[3];
    }

    for(int i=0;i<=n;i++)
    {
        for(int j=0;j<=3;j++)
        {
            double temp=0;
            for(int k=0;k<=n;k++)
            {
                temp+=M[i][k]*N_n[k][j];
            }
            N_tra[i][j]=temp;
        }
    }
    // if(id==8)
    // {
    //     for(int j=0;j<=n;j++)
    //     {
    //         for(int i=0;i<=j;i++)
    //         {
    //             if((j-i)%2==0)
    //                 coe=1;
    //             else
    //                 coe=-1;
    //             printf("coe*C[n][i]*C[n-i][j-i] = %f ",coe*C[n][i]*C[n-i][j-i]);
    //         }
    //         printf("\n");
    //     }
    //     for(int i=0;i<=n;i++)
    //     {
    //         printf("N_tra[%d][0]=%f ",i,N_tra[i][0]);
    //         printf("N_tra[%d][1]=%f ",i,N_tra[i][1]);
    //         printf("N_tra[%d][2]=%f ",i,N_tra[i][2]);
    //         printf("N_tra[%d][3]=%f\n",i,N_tra[i][3]);
    //     }
    //     for(int i=0;i<=n;i++)
    //     {
    //         printf("N_n[%d][0]=%f ",i,N_n[i][0]);
    //         printf("N_n[%d][1]=%f ",i,N_n[i][1]);
    //         printf("N_n[%d][2]=%f ",i,N_n[i][2]);
    //         printf("N_n[%d][3]=%f\n",i,N_n[i][3]);
    //     }
    //     for(int i=0;i<=m;i++)
    //     {
    //         printf("N_m[%d][0]=%f ",i,N_m[i][0]);
    //         printf("N_m[%d][1]=%f ",i,N_m[i][1]);
    //         printf("N_m[%d][2]=%f ",i,N_m[i][2]);
    //         printf("N_m[%d][3]=%f\n",i,N_m[i][3]);
    //     }
    // }

    for(int i=0;i<=m;i++)
    {
        N_minus[i][0]=N_tra[i][0]-N_m[i][0];
        N_minus[i][1]=N_tra[i][1]-N_m[i][1];
        N_minus[i][2]=N_tra[i][2]-N_m[i][2];
        N_minus[i][3]=N_tra[i][3]-N_m[i][3];
        // printf("N_minus[%d][0]=%f\n",i,N_minus[i][0]);
    }
    for(int i=m+1;i<=n;i++)
    {
        N_minus[i][0]=N_tra[i][0];
        N_minus[i][1]=N_tra[i][1];
        N_minus[i][2]=N_tra[i][2];
        N_minus[i][3]=N_tra[i][3];
    }
    for(int i=0;i<=2*n;i++)
    {
        double temx=0;
        double temy=0;
        double temz=0;
        double temw=0;
        for(int j=max(i-n,0);j<=min(n,i);j++)
        {
            temx+=N_minus[j][0]*N_minus[i-j][0];
            temy+=N_minus[j][1]*N_minus[i-j][1];
            temz+=N_minus[j][2]*N_minus[i-j][2];
            temw+=N_minus[j][3]*N_minus[i-j][3];
        }
        N_new[i][0]=temx;
        N_new[i][1]=temy;
        N_new[i][2]=temz;
        N_new[i][3]=temw;
    }
    // if(id==8)
    // {
    //     for(int i=0;i<=n;i++)
    //     {
    //         printf("N_minus[%d][0]=%f ",i,N_minus[i][0]);
    //         printf("N_minus[%d][1]=%f ",i,N_minus[i][1]);
    //         printf("N_minus[%d][2]=%f ",i,N_minus[i][2]);
    //         printf("N_minus[%d][3]=%f\n",i,N_minus[i][3]);
    //     }
    //     for(int i=0;i<=2*n;i++)
    //     {
    //         printf("N_new[%d][0]=%f ",i,N_new[i][0]);
    //         printf("N_new[%d][1]=%f ",i,N_new[i][1]);
    //         printf("N_new[%d][2]=%f ",i,N_new[i][2]);
    //         printf("N_new[%d][3]=%f\n",i,N_new[i][3]);
    //     }
    // }
    // for(int i=0;i<=2*n;i++)
    // {
    //     error+=( (N_new[i][0]+N_new[i][1]+N_new[i][2]+N_new[i][3]) / (i+1) );
    // }
    double lr=1;
    for(int i=0;i<=2*n;i++)
    {
        error+=(N_new[i][0]+N_new[i][1]+N_new[i][2]+N_new[i][3])*lr;
        lr*=u;
    }
    return error;
}


__device__ double calculate_trans_L1_max_error(int n, int m, Point *P_n, Point *P_m,double a,double b, int id)
{
    double error=0;
    double N_n[SIZE][4];
    double N_tra[SIZE][4];
    double N_m[4][4];
    double N_minus[SIZE][4];
    double N_new[2*SIZE-1][4];
    double coe=1;
    double temp[4];
    double M[SIZE][SIZE];
    double left[SIZE],right[SIZE];
    left[0]=1;
    right[0]=1;
    for(int i=1;i<=n;i++)
    {
        left[i]=left[i-1]*a;   
        right[i]=right[i-1]*(b-a);
    }
    for(int i=0;i<=n;i++)
    {
        for(int j=0;j<=n;j++)
        {
            if(j<i)
            {
                M[i][j]=0;
            }
            else
            {
                M[i][j]=C[j][i]*left[j-i]*right[i];
            }
        }
    }

    for(int j=0;j<=n;j++)
    {
        temp[0]=0;
        temp[1]=0;
        temp[2]=0;
        temp[3]=0;
        for(int i=0;i<=j;i++)
        {
            if((j-i)%2==0)
                coe=1;
            else
                coe=-1;
            temp[0]+=coe*C[n][i]*C[n-i][j-i]*P_n[i].x;
            temp[1]+=coe*C[n][i]*C[n-i][j-i]*P_n[i].y;
            temp[2]+=coe*C[n][i]*C[n-i][j-i]*P_n[i].z;
            temp[3]+=coe*C[n][i]*C[n-i][j-i]*P_n[i].w;
        }
        N_n[j][0]=temp[0];
        N_n[j][1]=temp[1];
        N_n[j][2]=temp[2];
        N_n[j][3]=temp[3];
    }
    for(int j=0;j<=m;j++)
    {
        temp[0]=0;
        temp[1]=0;
        temp[2]=0;
        temp[3]=0;
        for(int i=0;i<=j;i++)
        {
            if((j-i)%2==0)
                coe=1;
            else
                coe=-1;
            temp[0]+=coe*C[m][i]*C[m-i][j-i]*P_m[i].x;
            temp[1]+=coe*C[m][i]*C[m-i][j-i]*P_m[i].y;
            temp[2]+=coe*C[m][i]*C[m-i][j-i]*P_m[i].z;
            temp[3]+=coe*C[m][i]*C[m-i][j-i]*P_m[i].w;
        }
        N_m[j][0]=temp[0];
        N_m[j][1]=temp[1];
        N_m[j][2]=temp[2];
        N_m[j][3]=temp[3];
    }

    for(int i=0;i<=n;i++)
    {
        for(int j=0;j<=3;j++)
        {
            double temp=0;
            for(int k=0;k<=n;k++)
            {
                temp+=M[i][k]*N_n[k][j];
            }
            N_tra[i][j]=temp;
        }
    }


    for(int i=0;i<=m;i++)
    {
        N_minus[i][0]=N_tra[i][0]-N_m[i][0];
        N_minus[i][1]=N_tra[i][1]-N_m[i][1];
        N_minus[i][2]=N_tra[i][2]-N_m[i][2];
        N_minus[i][3]=N_tra[i][3]-N_m[i][3];
        // printf("N_minus[%d][0]=%f\n",i,N_minus[i][0]);
    }
    for(int i=m+1;i<=n;i++)
    {
        N_minus[i][0]=N_tra[i][0];
        N_minus[i][1]=N_tra[i][1];
        N_minus[i][2]=N_tra[i][2];
        N_minus[i][3]=N_tra[i][3];
    }
    for(int i=0;i<=2*n;i++)
    {
        double temx=0;
        double temy=0;
        double temz=0;
        double temw=0;
        for(int j=max(i-n,0);j<=min(n,i);j++)
        {
            temx+=N_minus[j][0]*N_minus[i-j][0];
            temy+=N_minus[j][1]*N_minus[i-j][1];
            temz+=N_minus[j][2]*N_minus[i-j][2];
            temw+=N_minus[j][3]*N_minus[i-j][3];
        }
        N_new[i][0]=temx;
        N_new[i][1]=temy;
        N_new[i][2]=temz;
        N_new[i][3]=temw;
    }

    // if(id==8)
    // {
    //     for(int j=0;j<=n;j++)
    //     {
    //         for(int i=0;i<=j;i++)
    //         {
    //             if((j-i)%2==0)
    //                 coe=1;
    //             else
    //                 coe=-1;
    //             printf("coe*C[n][i]*C[n-i][j-i] = %f ",coe*C[n][i]*C[n-i][j-i]);
    //         }
    //         printf("\n");
    //     }
    //     for(int i=0;i<=n;i++)
    //     {
    //         printf("N_tra[%d][0]=%f ",i,N_tra[i][0]);
    //         printf("N_tra[%d][1]=%f ",i,N_tra[i][1]);
    //         printf("N_tra[%d][2]=%f ",i,N_tra[i][2]);
    //         printf("N_tra[%d][3]=%f\n",i,N_tra[i][3]);
    //     }
    //     for(int i=0;i<=n;i++)
    //     {
    //         printf("N_n[%d][0]=%f ",i,N_n[i][0]);
    //         printf("N_n[%d][1]=%f ",i,N_n[i][1]);
    //         printf("N_n[%d][2]=%f ",i,N_n[i][2]);
    //         printf("N_n[%d][3]=%f\n",i,N_n[i][3]);
    //     }
    //     for(int i=0;i<=m;i++)
    //     {
    //         printf("N_m[%d][0]=%f ",i,N_m[i][0]);
    //         printf("N_m[%d][1]=%f ",i,N_m[i][1]);
    //         printf("N_m[%d][2]=%f ",i,N_m[i][2]);
    //         printf("N_m[%d][3]=%f\n",i,N_m[i][3]);
    //     }
    //     for(int i=0;i<=n;i++)
    //     {
    //         printf("N_minus[%d][0]=%f ",i,N_minus[i][0]);
    //         printf("N_minus[%d][1]=%f ",i,N_minus[i][1]);
    //         printf("N_minus[%d][2]=%f ",i,N_minus[i][2]);
    //         printf("N_minus[%d][3]=%f\n",i,N_minus[i][3]);
    //     }
    // }

    double coef[SIZE*2-1];
    for(int i=0;i<SIZE*2-1;i++)
    {
        coef[i]=N_new[i][0]+N_new[i][1]+N_new[i][2]+N_new[i][3];
    }
    Eigen::VectorXd coeffs = definePolynomial(2*SIZE-1,coef);
    Eigen::VectorXd derivativeCoeffs = derivative(coeffs);
    
    //利用牛顿迭代求方程的根
    double initial_guess = 0.1;
    double max_l1_error=0;
    double root;
    for(int i=1;i<=9;i++)
    {
        double temp_root=newtonIteration(derivativeCoeffs, initial_guess*i);
        double temp_error=evaluatePolynomial(coeffs,temp_root);
        if(temp_error>max_l1_error)
        {
            max_l1_error=temp_error;
            root=temp_root;
        }
        // max_l1_error=max(max_l1_error,evaluatePolynomial(coeffs,root[i]));
    }
    // double root;
    // root = newtonIteration(derivativeCoeffs, initial_guess);
    // max_l1_error=evaluatePolynomial(coeffs,root);

    // printf("a= %f b=%f\n",a,b);
    // for(int i=0;i<=2*n;i++)
    // {
    //     printf("coef[%d]=%f\n",i,coef[i]);
    // }
    // printf("root = %f\n",root);

    // double max_l1_error=evaluatePolynomial(coeffs,root);
    return max_l1_error;
}

__device__ void calculate_trans_L1_max_error_and_root(int n, int m, Point *P_n, Point *P_m,double a,double b, int id,double *ans)
{
    double error=0;
    double N_n[SIZE][4];
    double N_tra[SIZE][4];
    double N_m[4][4];
    double N_minus[SIZE][4];
    double N_new[2*SIZE-1][4];
    double coe=1;
    double temp[4];
    double M[SIZE][SIZE];
    double left[SIZE],right[SIZE];
    left[0]=1;
    right[0]=1;
    for(int i=1;i<=n;i++)
    {
        left[i]=left[i-1]*a;   
        right[i]=right[i-1]*(b-a);
    }
    for(int i=0;i<=n;i++)
    {
        for(int j=0;j<=n;j++)
        {
            if(j<i)
            {
                M[i][j]=0;
            }
            else
            {
                M[i][j]=C[j][i]*left[j-i]*right[i];
            }
        }
    }

    for(int j=0;j<=n;j++)
    {
        temp[0]=0;
        temp[1]=0;
        temp[2]=0;
        temp[3]=0;
        for(int i=0;i<=j;i++)
        {
            if((j-i)%2==0)
                coe=1;
            else
                coe=-1;
            temp[0]+=coe*C[n][i]*C[n-i][j-i]*P_n[i].x;
            temp[1]+=coe*C[n][i]*C[n-i][j-i]*P_n[i].y;
            temp[2]+=coe*C[n][i]*C[n-i][j-i]*P_n[i].z;
            temp[3]+=coe*C[n][i]*C[n-i][j-i]*P_n[i].w;
        }
        N_n[j][0]=temp[0];
        N_n[j][1]=temp[1];
        N_n[j][2]=temp[2];
        N_n[j][3]=temp[3];
    }
    for(int j=0;j<=m;j++)
    {
        temp[0]=0;
        temp[1]=0;
        temp[2]=0;
        temp[3]=0;
        for(int i=0;i<=j;i++)
        {
            if((j-i)%2==0)
                coe=1;
            else
                coe=-1;
            temp[0]+=coe*C[m][i]*C[m-i][j-i]*P_m[i].x;
            temp[1]+=coe*C[m][i]*C[m-i][j-i]*P_m[i].y;
            temp[2]+=coe*C[m][i]*C[m-i][j-i]*P_m[i].z;
            temp[3]+=coe*C[m][i]*C[m-i][j-i]*P_m[i].w;
        }
        N_m[j][0]=temp[0];
        N_m[j][1]=temp[1];
        N_m[j][2]=temp[2];
        N_m[j][3]=temp[3];
    }

    for(int i=0;i<=n;i++)
    {
        for(int j=0;j<=3;j++)
        {
            double temp=0;
            for(int k=0;k<=n;k++)
            {
                temp+=M[i][k]*N_n[k][j];
            }
            N_tra[i][j]=temp;
        }
    }


    for(int i=0;i<=m;i++)
    {
        N_minus[i][0]=N_tra[i][0]-N_m[i][0];
        N_minus[i][1]=N_tra[i][1]-N_m[i][1];
        N_minus[i][2]=N_tra[i][2]-N_m[i][2];
        N_minus[i][3]=N_tra[i][3]-N_m[i][3];
        // printf("N_minus[%d][0]=%f\n",i,N_minus[i][0]);
    }
    for(int i=m+1;i<=n;i++)
    {
        N_minus[i][0]=N_tra[i][0];
        N_minus[i][1]=N_tra[i][1];
        N_minus[i][2]=N_tra[i][2];
        N_minus[i][3]=N_tra[i][3];
    }
    for(int i=0;i<=2*n;i++)
    {
        double temx=0;
        double temy=0;
        double temz=0;
        double temw=0;
        for(int j=max(i-n,0);j<=min(n,i);j++)
        {
            temx+=N_minus[j][0]*N_minus[i-j][0];
            temy+=N_minus[j][1]*N_minus[i-j][1];
            temz+=N_minus[j][2]*N_minus[i-j][2];
            temw+=N_minus[j][3]*N_minus[i-j][3];
        }
        N_new[i][0]=temx;
        N_new[i][1]=temy;
        N_new[i][2]=temz;
        N_new[i][3]=temw;
    }


    double coef[SIZE*2-1];
    for(int i=0;i<SIZE*2-1;i++)
    {
        coef[i]=N_new[i][0]+N_new[i][1]+N_new[i][2]+N_new[i][3];
    }
    Eigen::VectorXd coeffs = definePolynomial(2*SIZE-1,coef);
    Eigen::VectorXd derivativeCoeffs = derivative(coeffs);
    
    //利用牛顿迭代求方程的根
    double initial_guess = 0.1;
    double max_l1_error=0;
    double root;
    for(int i=1;i<=9;i++)
    {
        double temp_root=newtonIteration(derivativeCoeffs, initial_guess*i);
        double temp_error=evaluatePolynomial(coeffs,temp_root);
        if(temp_error>max_l1_error)
        {
            max_l1_error=temp_error;
            root=temp_root;
        }
        // max_l1_error=max(max_l1_error,evaluatePolynomial(coeffs,root[i]));
    }
    ans[0]=root;
    ans[1]=max_l1_error;
}


__device__ void calculate_tt(int n, int m, Point *P_n, Point *P_m,int id,int num_seg,double *tt)
{
    double error=0;
    double N_n[SIZE][4];
    double N_m[4][4];
    double N_minus[SIZE][4];
    double N_new[2*SIZE-1][4];
    double coe=1;
    double temp[4];
    for(int j=0;j<=n;j++)
    {
        temp[0]=0;
        temp[1]=0;
        temp[2]=0;
        temp[3]=0;
        for(int i=0;i<=j;i++)
        {
            if((j-i)%2==0)
                coe=1;
            else
                coe=-1;
            temp[0]+=coe*C[n][i]*C[n-i][j-i]*P_n[i].x;
            temp[1]+=coe*C[n][i]*C[n-i][j-i]*P_n[i].y;
            temp[2]+=coe*C[n][i]*C[n-i][j-i]*P_n[i].z;
            temp[3]+=coe*C[n][i]*C[n-i][j-i]*P_n[i].w;
        }
        N_n[j][0]=temp[0];
        N_n[j][1]=temp[1];
        N_n[j][2]=temp[2];
        N_n[j][3]=temp[3];
    }
    for(int j=0;j<=m;j++)
    {
        temp[0]=0;
        temp[1]=0;
        temp[2]=0;
        temp[3]=0;
        for(int i=0;i<=j;i++)
        {
            if((j-i)%2==0)
                coe=1;
            else
                coe=-1;
            temp[0]+=coe*C[m][i]*C[m-i][j-i]*P_m[i].x;
            temp[1]+=coe*C[m][i]*C[m-i][j-i]*P_m[i].y;
            temp[2]+=coe*C[m][i]*C[m-i][j-i]*P_m[i].z;
            temp[3]+=coe*C[m][i]*C[m-i][j-i]*P_m[i].w;
        }
        N_m[j][0]=temp[0];
        N_m[j][1]=temp[1];
        N_m[j][2]=temp[2];
        N_m[j][3]=temp[3];
    }
    for(int i=0;i<=m;i++)
    {
        N_minus[i][0]=N_n[i][0]-N_m[i][0];
        N_minus[i][1]=N_n[i][1]-N_m[i][1];
        N_minus[i][2]=N_n[i][2]-N_m[i][2];
        N_minus[i][3]=N_n[i][3]-N_m[i][3];
        // printf("N_minus[%d][0]=%f\n",i,N_minus[i][0]);
    }
    for(int i=m+1;i<=n;i++)
    {
        N_minus[i][0]=N_n[i][0];
        N_minus[i][1]=N_n[i][1];
        N_minus[i][2]=N_n[i][2];
        N_minus[i][3]=N_n[i][3];
    }
    for(int i=0;i<=2*n;i++)
    {
        double temx=0;
        double temy=0;
        double temz=0;
        double temw=0;
        for(int j=max(i-n,0);j<=min(n,i);j++)
        {
            temx+=N_minus[j][0]*N_minus[i-j][0];
            temy+=N_minus[j][1]*N_minus[i-j][1];
            temz+=N_minus[j][2]*N_minus[i-j][2];
            temw+=N_minus[j][3]*N_minus[i-j][3];
        }
        N_new[i][0]=temx;
        N_new[i][1]=temy;
        N_new[i][2]=temz;
        N_new[i][3]=temw;
    }
    for(int i=0;i<=2*n;i++)
    {
        error+=(N_new[i][0]+N_new[i][1]+N_new[i][2]+N_new[i][3])/(i+1);
    }
    for(int i=1;i<=num_seg-1;i++)
    {
        double err=error*i/num_seg;
        printf("err=%f\n",err);
        double left=0;
        double right=1;
        double ans=0;
        while(right-left>1e-4)
        {
            double mid=(left+right)/2;
            double temp=0;
            double co=mid;
            for(int j=0;j<=2*n;j++)
            {
                temp+=(N_new[j][0]+N_new[j][1]+N_new[j][2]+N_new[j][3])/(j+1)*co;
                co*=mid;
            }
            if(temp>err)
            {
                right=mid;
            }
            else
            {
                left=mid;
            }
            ans=temp;
        }
        tt[i-1]=left;
        printf("tt[%d]=%f temp =%f\n",i-1,tt[i-1],ans);
    }
    // if(id==7)
    // {
    //     for(int i=0;i<=n;i++)
    //     {
    //         printf("N_minus[%d][0]=%f\n",i,N_minus[i][0]);
    //         printf("N_minus[%d][1]=%f\n",i,N_minus[i][1]);
    //         printf("N_minus[%d][2]=%f\n",i,N_minus[i][2]);
    //         printf("N_minus[%d][3]=%f\n",i,N_minus[i][3]);
    //     }
    //     for(int i=0;i<=2*n;i++)
    //     {
    //         printf("N_new[%d][0]=%f\n",i,N_new[i][0]);
    //         printf("N_new[%d][1]=%f\n",i,N_new[i][1]);
    //         printf("N_new[%d][2]=%f\n",i,N_new[i][2]);
    //         printf("N_new[%d][3]=%f\n",i,N_new[i][3]);
    //     }
    //     printf("error=%f\n",error);
    //     printf("%f\n",N_minus[4][0]*N_minus[4][0]);
    // }
}

__device__ double calculate_trans_error(int n, int m, Point *P_n, Point *P_m,double a,double b, int id)
{
    double error=0;
    double N_n[SIZE][4];
    double N_tra[SIZE][4];
    double N_m[4][4];
    double N_minus[SIZE][4];
    double N_new[2*SIZE-1][4];
    double coe=1;
    double temp[4];
    double M[SIZE][SIZE];
    double left[SIZE],right[SIZE];
    left[0]=1;
    right[0]=1;
    for(int i=1;i<=n;i++)
    {
        left[i]=left[i-1]*a;   
        right[i]=right[i-1]*(b-a);
    }
    for(int i=0;i<=n;i++)
    {
        for(int j=0;j<=n;j++)
        {
            if(j<i)
            {
                M[i][j]=0;
            }
            else
            {
                M[i][j]=C[j][i]*left[j-i]*right[i];
            }
        }
    }
    if(id==8)
    {
        for(int i=0;i<=n;i++)
        {
            for(int j=0;j<=n;j++)
            {
                printf("M[%d][%d]=%f ",i,j,M[i][j]);
            }
            printf("\n");
        }
    }
    for(int j=0;j<=n;j++)
    {
        temp[0]=0;
        temp[1]=0;
        temp[2]=0;
        temp[3]=0;
        for(int i=0;i<=j;i++)
        {
            if((j-i)%2==0)
                coe=1;
            else
                coe=-1;
            temp[0]+=coe*C[n][i]*C[n-i][j-i]*P_n[i].x;
            temp[1]+=coe*C[n][i]*C[n-i][j-i]*P_n[i].y;
            temp[2]+=coe*C[n][i]*C[n-i][j-i]*P_n[i].z;
            temp[3]+=coe*C[n][i]*C[n-i][j-i]*P_n[i].w;
        }
        N_n[j][0]=temp[0];
        N_n[j][1]=temp[1];
        N_n[j][2]=temp[2];
        N_n[j][3]=temp[3];
    }
    for(int j=0;j<=m;j++)
    {
        temp[0]=0;
        temp[1]=0;
        temp[2]=0;
        temp[3]=0;
        for(int i=0;i<=j;i++)
        {
            if((j-i)%2==0)
                coe=1;
            else
                coe=-1;
            temp[0]+=coe*C[m][i]*C[m-i][j-i]*P_m[i].x;
            temp[1]+=coe*C[m][i]*C[m-i][j-i]*P_m[i].y;
            temp[2]+=coe*C[m][i]*C[m-i][j-i]*P_m[i].z;
            temp[3]+=coe*C[m][i]*C[m-i][j-i]*P_m[i].w;
        }
        N_m[j][0]=temp[0];
        N_m[j][1]=temp[1];
        N_m[j][2]=temp[2];
        N_m[j][3]=temp[3];
    }

    for(int i=0;i<=n;i++)
    {
        for(int j=0;j<=3;j++)
        {
            double temp=0;
            for(int k=0;k<=n;k++)
            {
                temp+=M[i][k]*N_n[k][j];
            }
            N_tra[i][j]=temp;
        }
    }
    if(id==8)
    {
        for(int j=0;j<=n;j++)
        {
            for(int i=0;i<=j;i++)
            {
                if((j-i)%2==0)
                    coe=1;
                else
                    coe=-1;
                printf("coe*C[n][i]*C[n-i][j-i] = %f ",coe*C[n][i]*C[n-i][j-i]);
            }
            printf("\n");
        }
        for(int i=0;i<=n;i++)
        {
            printf("N_tra[%d][0]=%f ",i,N_tra[i][0]);
            printf("N_tra[%d][1]=%f ",i,N_tra[i][1]);
            printf("N_tra[%d][2]=%f ",i,N_tra[i][2]);
            printf("N_tra[%d][3]=%f\n",i,N_tra[i][3]);
        }
        for(int i=0;i<=n;i++)
        {
            printf("N_n[%d][0]=%f ",i,N_n[i][0]);
            printf("N_n[%d][1]=%f ",i,N_n[i][1]);
            printf("N_n[%d][2]=%f ",i,N_n[i][2]);
            printf("N_n[%d][3]=%f\n",i,N_n[i][3]);
        }
        for(int i=0;i<=m;i++)
        {
            printf("N_m[%d][0]=%f ",i,N_m[i][0]);
            printf("N_m[%d][1]=%f ",i,N_m[i][1]);
            printf("N_m[%d][2]=%f ",i,N_m[i][2]);
            printf("N_m[%d][3]=%f\n",i,N_m[i][3]);
        }
    }

    for(int i=0;i<=m;i++)
    {
        N_minus[i][0]=N_tra[i][0]-N_m[i][0];
        N_minus[i][1]=N_tra[i][1]-N_m[i][1];
        N_minus[i][2]=N_tra[i][2]-N_m[i][2];
        N_minus[i][3]=N_tra[i][3]-N_m[i][3];
        // printf("N_minus[%d][0]=%f\n",i,N_minus[i][0]);
    }
    for(int i=m+1;i<=n;i++)
    {
        N_minus[i][0]=N_tra[i][0];
        N_minus[i][1]=N_tra[i][1];
        N_minus[i][2]=N_tra[i][2];
        N_minus[i][3]=N_tra[i][3];
    }
    for(int i=0;i<=2*n;i++)
    {
        double temx=0;
        double temy=0;
        double temz=0;
        double temw=0;
        for(int j=max(i-n,0);j<=min(n,i);j++)
        {
            temx+=N_minus[j][0]*N_minus[i-j][0];
            temy+=N_minus[j][1]*N_minus[i-j][1];
            temz+=N_minus[j][2]*N_minus[i-j][2];
            temw+=N_minus[j][3]*N_minus[i-j][3];
        }
        N_new[i][0]=temx;
        N_new[i][1]=temy;
        N_new[i][2]=temz;
        N_new[i][3]=temw;
    }
    if(id==8)
    {
        for(int i=0;i<=n;i++)
        {
            printf("N_minus[%d][0]=%f ",i,N_minus[i][0]);
            printf("N_minus[%d][1]=%f ",i,N_minus[i][1]);
            printf("N_minus[%d][2]=%f ",i,N_minus[i][2]);
            printf("N_minus[%d][3]=%f\n",i,N_minus[i][3]);
        }
        for(int i=0;i<=2*n;i++)
        {
            printf("N_new[%d][0]=%f ",i,N_new[i][0]);
            printf("N_new[%d][1]=%f ",i,N_new[i][1]);
            printf("N_new[%d][2]=%f ",i,N_new[i][2]);
            printf("N_new[%d][3]=%f\n",i,N_new[i][3]);
        }
    }
    for(int i=0;i<=2*n;i++)
    {
        error+=( (N_new[i][0]+N_new[i][1]+N_new[i][2]+N_new[i][3]) / (i+1) );
    }
    return error;
}

__device__ void add_child(CBC *parent)
{
    // cudaMalloc((void**)&parent->child, sizeof(CBC*) * parent->seg);
    // for (int i = 0; i < parent->seg; i++) 
    // {
    //     cudaMalloc((void**)&parent->child[i], sizeof(CBC));
    // }
    parent->child=(CBC**)malloc(sizeof(CBC*)*parent->seg);
    for(int i=0;i<parent->seg;i++)
    {
        parent->child[i]=(CBC*)malloc(sizeof(CBC));
    }

}

__device__ void set_child(CBC *parent, int id, double error, int dep,double le,double ri,Point *P)
{
    if (id < 0 || id >= parent->seg) return;  // 确保id在合法范围内
    parent->child[id]->error = error;
    parent->child[id]->id = id;
    parent->child[id]->seg = 1;
    parent->child[id]->dep = dep;
    parent->child[id]->left=le;
    parent->child[id]->right=ri;

    for(int i=0;i<4;i++)
    {
        parent->child[id]->control_points[i]=P[i];
    }
}

__device__ void print_child(CBC *parent)
{
    for(int i=0;i<parent->seg;i++)
    {
        printf("child[%d].error=%f\n",i,parent->child[i]->error);
        printf("child[%d].id=%d\n",i,parent->child[i]->id);
        printf("child[%d].seg=%d\n",i,parent->child[i]->seg);
        printf("child[%d].dep=%d\n",i,parent->child[i]->dep);
        printf("child[%d].left=%f\n",i,parent->child[i]->left);
        printf("child[%d].right=%f\n",i,parent->child[i]->right);
        printf("child[%d].control_points:\n",i);
        for(int j=0;j<4;j++)
        {
            printf("child[%d].P[%d].x=%f\t",i,j,parent->child[i]->control_points[j].x);
            printf("child[%d].P[%d].y=%f\t",i,j,parent->child[i]->control_points[j].y);
            printf("child[%d].P[%d].z=%f\t",i,j,parent->child[i]->control_points[j].z);
            printf("child[%d].P[%d].w=%f\n",i,j,parent->child[i]->control_points[j].w);
        }
    }
}

__device__ void get_control_points(Point *P, Point *P_left,Point *P_right,double cu)
{
    P_left[0]=P[0];
    P_left[1]=P[1]*cu+P[0]*(1-cu);
    P_left[2]=P[2]*cu*cu +P[1]*2*cu*(1-cu)+P[0]*(1-cu)*(1-cu);
    P_left[3]=P[3]*cu*cu*cu+P[2]*3*cu*cu*(1-cu)+P[1]*3*cu*(1-cu)*(1-cu)+P[0]*(1-cu)*(1-cu)*(1-cu);

    P_right[3]=P[3];
    P_right[2]=P[3]*cu+P[2]*(1-cu);
    P_right[1]=P[3]*cu*cu +P[2]*2*cu*(1-cu)+P[1]*(1-cu)*(1-cu);
    P_right[0]=P[3]*cu*cu*cu+P[2]*3*cu*cu*(1-cu)+P[1]*3*cu*(1-cu)*(1-cu)+P[0]*(1-cu)*(1-cu)*(1-cu);

    // for(int i=0;i<4;i++)
    // {
    //     printf("P_left[%d].x=%f\t",i,P_left[i].x);
    //     printf("P_left[%d].y=%f\t",i,P_left[i].y);
    //     printf("P_left[%d].z=%f\t",i,P_left[i].z);
    //     printf("P_left[%d].w=%f\n",i,P_left[i].w);
    // }
    // for(int i=0;i<4;i++)
    // {
    //     printf("P_right[%d].x=%f\t",i,P_right[i].x);
    //     printf("P_right[%d].y=%f\t",i,P_right[i].y);
    //     printf("P_right[%d].z=%f\t",i,P_right[i].z);
    //     printf("P_right[%d].w=%f\n",i,P_right[i].w);
    // }
}

__device__ void point_re_subdivision(Point *former,double le,double ri, Point *now, Point *P_n, int p_n)
{
    Point a[4];
    Point b[4];
    get_control_points(former,a,b,le);
    get_control_points(b,now,a,(ri-le)/(1.0-le));
    Point x,y;
    x=evaluate_bezier_curve(le,P_n,p_n);
    y=evaluate_bezier_curve(ri,P_n,p_n);
    now[0]=x;
    now[3]=y;
}

__device__ void point_re_subdivision_v2(Point *former,double le,double ri, Point *now, Point *P_n, int p_n,double original_left,double original_right)
{
    Point a[4];
    Point b[4];
    get_control_points(former,a,b,le);
    get_control_points(b,now,a,(ri-le)/(1.0-le));
    Point x,y;
    x=evaluate_bezier_curve(original_left,P_n,p_n);
    y=evaluate_bezier_curve(original_right,P_n,p_n);
    now[0]=x;
    now[3]=y;
}



__global__ void calculate_error_kernel(int p_n, int p_m,int knot_size, CBC * d_b)
{
    int id=blockDim.x*blockIdx.x+threadIdx.x;
    double result;
    double result_L1;
    if(id>=knot_size-p_n-1||id<p_n)
        return;
    Point P_n[SIZE];
    Point P_m[4];
    for(int i=0;i<SIZE;i++)
    {
        // P_n[i].x=d_big_Q[id*SIZE*4+i];
        // P_n[i].y=d_big_Q[id*SIZE*4+SIZE+i];
        // P_n[i].z=d_big_Q[id*SIZE*4+2*SIZE+i];
        // P_n[i].w=d_big_Q[id*SIZE*4+3*SIZE+i];
        P_n[i].x=d_Q[id][i][0];
        P_n[i].y=d_Q[id][i][1];
        P_n[i].z=d_Q[id][i][2];
        P_n[i].w=d_Q[id][i][3];
    }
    for(int i=0;i<4;i++)
    {
        P_m[i].x=d_cubic_Q[id][i][0];
        P_m[i].y=d_cubic_Q[id][i][1];
        P_m[i].z=d_cubic_Q[id][i][2];
        P_m[i].w=d_cubic_Q[id][i][3];
    }
    result=calculate_error(p_n,p_m,P_n,P_m,id);
    double t_l1[2];
    calculate_max_L1_error(p_n,p_m,P_n,P_m,id,t_l1);
    // result_L1=calculate_L1_error(p_n,p_m,P_n,P_m,id,0.5);

    // printf("result %d = %f\n",id,result);
    // printf("result_L1 %d = %.9f\n",id,result_L1);
    printf("id = %d t = %f max l1 error = %f\n",id,t_l1[0],t_l1[1]);
    seg[id]=ceil(t_l1[1]/epsl);
    printf("seg %d = %d\n",id,seg[id]);
    d_b[id].error=t_l1[1];
    d_b[id].id=id;
    d_b[id].seg=min(seg[id],4);
    d_b[id].dep=1;

    if(seg[id]>1)
    {
        add_child(&d_b[id]);
        // set_child(&d_b[id],0,result/seg[id],d_b[id].dep);
    }
    // if(id==7)
    // {
    //     double u=0.2;
    //     Point left[4];
    //     Point right[4];
    //     get_control_points(P_m,left,right,u);
    //     double error=calculate_trans_error(p_n,p_m,P_n,left,0,u,id);
    //     printf("seg[%d].error when get point =%f   u=%f\n",id,error,u);
    // }
    if(id==8)
    {
        double tt[3];
        calculate_tt(p_n,p_m,P_n,P_m,id,4,tt);
        Point left[4];
        Point right[4];
        double sum=0;
        Point temp[4];
        Point max_pp=evaluate_bezier_curve(0.154903,P_n,p_n);
        Point max_ppp=evaluate_bezier_curve(0.154903,P_m,p_m);
        double max_p_2_p_error=calculate_point_error(max_pp,max_ppp);
        printf("max_p_2_p_error=%f\n",max_p_2_p_error);
        double cal_error=calculate_L1_error(p_n,p_m,P_n,P_m,id,0.154903);
        printf("cal_error=%f\n",cal_error);
        for(int i=0;i<d_b[id].seg;i++)
        {
            double u;
            if(i<d_b[id].seg-1)
                u=(tt[i]-sum)/(1-sum);
            if(i==0)
            {
                // get_control_points(P_m,left,right,u);
                sum=tt[i];

                point_re_subdivision(P_m,0,tt[i],left,P_n,p_n);


                // double error=calculate_trans_L1_error(p_n,p_m,P_n,left,0,u,id,1.0);
                double error=calculate_trans_L1_max_error(p_n,p_m,P_n,left,0,u,id);
                Point pp=evaluate_bezier_curve(u,P_n,p_n);
                Point ppp=evaluate_bezier_curve(1,left,3);
                double p_2_p_error=calculate_point_error(pp,ppp);
                printf("child[%d].error when get point =%f   u=%f\n",i,error,u);
                printf("pp[%d] = %f %f %f %f\n",i,pp.x,pp.y,pp.z,pp.w);
                printf("ppp[%d] = %f %f %f %f\n",i,ppp.x,ppp.y,ppp.z,ppp.w);
                printf("p2p_error=%f\n",p_2_p_error);

                set_child(&d_b[id],i,error,d_b[id].dep + 1,0,tt[i],left);
            }
            else if(i<d_b[id].seg-1)
            {
                // get_control_points(temp,left,right,u);
                sum=tt[i];

                point_re_subdivision(P_m,tt[i-1],tt[i],left,P_n,p_n);


                double error=calculate_trans_L1_max_error(p_n,p_m,P_n,left,tt[i-1],tt[i],id);
                // double error=calculate_trans_L1_error(p_n,p_m,P_n,left,tt[i-1],tt[i],id,1.0);
                Point pp=evaluate_bezier_curve(tt[i],P_n,p_n);
                Point ppp=evaluate_bezier_curve(1.0,left,3);
                double p_2_p_error=calculate_point_error(pp,ppp);
                printf("child[%d].error when get point =%f   u=%f\n",i,error,u);
                printf("pp[%d] = %f %f %f %f\n",i,pp.x,pp.y,pp.z,pp.w);
                printf("ppp[%d] = %f %f %f %f\n",i,ppp.x,ppp.y,ppp.z,ppp.w);
                printf("p2p_error=%f\n",p_2_p_error);

                set_child(&d_b[id],i,error,d_b[id].dep + 1,tt[i-1],tt[i],left);
            }
            else
            {

                point_re_subdivision(P_m,tt[i-1],1.0,right,P_n,p_n);

                double error=calculate_trans_L1_max_error(p_n,p_m,P_n,right,tt[i-1],1.0,id);
                // double error=calculate_trans_L1_error(p_n,p_m,P_n,right,tt[i-1],1.0,id,0.0);
                Point pp=evaluate_bezier_curve(tt[i-1],P_n,p_n);
                Point ppp=evaluate_bezier_curve(0.0,right,3);
                double p_2_p_error=calculate_point_error(pp,ppp);
                printf("child[%d].error when get point =%f   u=%f\n",i,error,u);
                printf("pp[%d] = %f %f %f %f\n",i,pp.x,pp.y,pp.z,pp.w);
                printf("ppp[%d] = %f %f %f %f\n",i,ppp.x,ppp.y,ppp.z,ppp.w);
                printf("p2p_error=%f\n",p_2_p_error);

                set_child(&d_b[id],i,error,d_b[id].dep+1,tt[i-1],1,right);
            }
            // if(i==0)
            // {
            //     set_child(&d_b[id],i,0,1,0,tt[i],left);
            // }
            // else if(i==d_b[id].seg-1)
            // {
            //     set_child(&d_b[id],i,result*i/d_b[id].seg,d_b[id].dep,tt[i-1],1,right);
            // }
            // else
            // {
            //     set_child(&d_b[id],i,result*i/d_b[id].seg,d_b[id].dep,tt[i-1],tt[i],left);
            // }
            temp[0]=right[0];
            temp[1]=right[1];
            temp[2]=right[2];
            temp[3]=right[3];
            // set_child(&d_b[id],i,result*i/d_b[id].seg,d_b[id].dep);
        }
        for(int i=0;i<3;i++)
        {
            printf("tt[%d]=%f\n",i,tt[i]);
        }
        printf("db[id].error=%f\n",d_b[id].error);
        printf("db[id].id=%d\n",d_b[id].id);
        printf("db[id].seg=%d\n",d_b[id].seg);
        printf("db[id].dep=%d\n",d_b[id].dep);

        print_child(&d_b[id]);

    }
}

__global__ void first_calculate_kernel(int p_n, int p_m,int knot_size, CBC * d_b,int *cnt,CBC *d_final_ans)
{
    int id=blockDim.x*blockIdx.x+threadIdx.x;
    double result_L1;
    if(id>=knot_size-p_n-1||id<p_n)
        return;
    Point P_n[SIZE];
    Point P_m[4];
    for(int i=0;i<SIZE;i++)
    {
        P_n[i].x=d_Q[id][i][0];
        P_n[i].y=d_Q[id][i][1];
        P_n[i].z=d_Q[id][i][2];
        P_n[i].w=d_Q[id][i][3];
    }
    for(int i=0;i<4;i++)
    {
        P_m[i].x=d_cubic_Q[id][i][0];
        P_m[i].y=d_cubic_Q[id][i][1];
        P_m[i].z=d_cubic_Q[id][i][2];
        P_m[i].w=d_cubic_Q[id][i][3];
    }
    double t_l1[2];
    calculate_max_L1_error(p_n,p_m,P_n,P_m,id,t_l1);
    // printf("id = %d t = %f max l1 error = %f\n",id,t_l1[0],t_l1[1]);
    seg[id]=ceil(t_l1[1]/epsl);
    prefix_error_sum[id]=t_l1[1];
    if(seg[id]==1)
    {
        seg[id]=0;
    }
    else
    {
        // seg[id]*=1;
        seg[id]=2;
    }
    // printf("seg %d = %d\n",id,seg[id]);
    d_b[id].error=t_l1[1];
    d_b[id].id=id;
    d_b[id].seg=seg[id];
    d_b[id].dep=1;
    d_b[id].control_points[0]=P_m[0];
    d_b[id].control_points[1]=P_m[1];
    d_b[id].control_points[2]=P_m[2];
    d_b[id].control_points[3]=P_m[3];
    d_b[id].left=0;
    d_b[id].right=1;
    d_b[id].original_left=0;
    d_b[id].original_right=1;
    d_b[id].original_bezier_curve_id=id;
    d_b[id].sub_value=t_l1[0];

    if(d_b[id].seg==0)
    {
        auto old=atomicAdd(cnt,1);
        // printf("cnt=%d\n",old);
        d_final_ans[old]=d_b[id];
    }

}

__device__ int find_index(int *arr, int x, int size)
{
    int l = 0;
    int r = size - 1;
    int result = -1;
    while (l <= r) 
    {
        int mid = l + (r - l) / 2;
        if (arr[mid] > x) 
        {
            result = mid;
            r = mid - 1;
        } 
        else 
        {
            l = mid + 1;
        }
    }
    return result;
}

__global__ void first_re_subdivision(int p_n,int p_m,int num_seg, CBC *d_parent,CBC *d_child,int knot_size,int *seg_sum,int *next_seg)
{
    int id=blockDim.x*blockIdx.x+threadIdx.x;
    if(id>=num_seg)
        return;
    int index_size=knot_size-p_n;
    d_child[id].id=id;
    d_child[id].parent_id=find_index(seg_sum,id,knot_size);
    d_child[id].original_bezier_curve_id=d_child[id].parent_id;

    // printf("d_child[%d].parent_id=%d\n",id,d_child[id].parent_id);
    int parent_id=d_child[id].parent_id;
    int start=seg_sum[parent_id-1];
    int end=seg_sum[parent_id]-1;
    d_child[id].left=(double)(id-start)/(end-start+1);
    d_child[id].right=(double)(id-start+1)/(end-start+1);
    // printf("d_child[%d].left= %f right= %f\n",id,d_child[id].left,d_child[id].right);

    Point P_father[4];
    P_father[0]=d_parent[parent_id].control_points[0];
    P_father[1]=d_parent[parent_id].control_points[1];
    P_father[2]=d_parent[parent_id].control_points[2];
    P_father[3]=d_parent[parent_id].control_points[3];

    Point P_child[4];
    Point P_n[SIZE];
    for(int i=0;i<SIZE;i++)
    {
        P_n[i].x=d_Q[d_child[id].original_bezier_curve_id][i][0];
        P_n[i].y=d_Q[d_child[id].original_bezier_curve_id][i][1];
        P_n[i].z=d_Q[d_child[id].original_bezier_curve_id][i][2];
        P_n[i].w=d_Q[d_child[id].original_bezier_curve_id][i][3];
    }

    double original_left=d_child[id].left * (d_parent[parent_id].right-d_parent[parent_id].left)+d_parent[parent_id].left;
    double original_right=d_child[id].right * (d_parent[parent_id].right-d_parent[parent_id].left)+d_parent[parent_id].left;

    // point_re_subdivision(P_father,d_child[id].left,d_child[id].right,P_child,P_n,p_n);
    point_re_subdivision_v2(P_father,d_child[id].left,d_child[id].right,P_child,P_n,p_n,original_left,original_right);
    double error=calculate_trans_L1_max_error(p_n,p_m,P_n,P_child,d_child[id].left,d_child[id].right,parent_id);
    d_child[id].error=error;
    d_child[id].dep=d_parent[parent_id].dep+1;
    d_child[id].seg=ceil(error/epsl);
    if(d_child[id].seg==1)
    {
        d_child[id].seg=0;
    }
    else
    {
        d_child[id].seg=2;
    }
    next_seg[id]=d_child[id].seg;
    // printf("next_seg[%d]=%d\n",id,next_seg[id]);
    d_child[id].control_points[0]=P_child[0];
    d_child[id].control_points[1]=P_child[1];
    d_child[id].control_points[2]=P_child[2];
    d_child[id].control_points[3]=P_child[3];

    d_child[id].original_left=original_left;
    d_child[id].original_right=original_right;
    printf("d_child[%d].error=%f d_father.error=%f\n",id,d_child[id].error,d_parent[parent_id].error);

    // if(id==3)
    // {
    //     printf("d_child[%d].original_bezier_id = %d\n",id,d_child[id].original_bezier_curve_id);
    //     printf("P_n:\n");
    //     for(int i=0;i<SIZE;i++)
    //     {
    //         printf("P_n[%d].x=%f\t",i,P_n[i].x);
    //         printf("P_n[%d].y=%f\t",i,P_n[i].y);
    //         printf("P_n[%d].z=%f\t",i,P_n[i].z);
    //         printf("P_n[%d].w=%f\n",i,P_n[i].w);
    //     }
    //     printf("P_father:\n");
    //     for(int i=0;i<4;i++)
    //     {
    //         printf("P_father[%d].x=%f\t",i,P_father[i].x);
    //         printf("P_father[%d].y=%f\t",i,P_father[i].y);
    //         printf("P_father[%d].z=%f\t",i,P_father[i].z);
    //         printf("P_father[%d].w=%f\n",i,P_father[i].w);
    //     }
    //     printf("P_child:\n");
    //     for(int i=0;i<4;i++)
    //     {
    //         printf("P_child[%d].x=%f\t",i,P_child[i].x);
    //         printf("P_child[%d].y=%f\t",i,P_child[i].y);
    //         printf("P_child[%d].z=%f\t",i,P_child[i].z);
    //         printf("P_child[%d].w=%f\n",i,P_child[i].w);
    //     }
    //     Point p1=evaluate_bezier_curve(original_left,P_n,p_n);
    //     Point p2=evaluate_bezier_curve(original_right,P_n,p_n);
    //     printf("p1.x=%f p1.y=%f p1.z=%f p1.w=%f\n",p1.x,p1.y,p1.z,p1.w);
    //     printf("p2.x=%f p2.y=%f p2.z=%f p2.w=%f\n",p2.x,p2.y,p2.z,p2.w);
    // }
}

__global__ void first_re_subdivision_at_max_value(int p_n,int p_m,int num_seg, CBC *d_parent,CBC *d_child,int knot_size,int *seg_sum,int *next_seg,int *cnt,CBC *d_final_ans)
{
    int id=blockDim.x*blockIdx.x+threadIdx.x;
    if(id>=num_seg)
        return;
    int index_size=knot_size-p_n;
    d_child[id].id=id;
    d_child[id].parent_id=find_index(seg_sum,id,knot_size);
    d_child[id].original_bezier_curve_id=d_child[id].parent_id;

    // printf("d_child[%d].parent_id=%d\n",id,d_child[id].parent_id);
    int parent_id=d_child[id].parent_id;
    int start=seg_sum[parent_id-1];
    int end=seg_sum[parent_id]-1;
    d_child[id].left=(double)(id-start)/(end-start+1);
    d_child[id].right=(double)(id-start+1)/(end-start+1);
    // printf("d_child[%d].left= %f right= %f\n",id,d_child[id].left,d_child[id].right);

    if(d_child[id].left==0)
    {
        d_child[id].right=d_parent[parent_id].sub_value;
    }
    else
    {
        d_child[id].left=d_parent[parent_id].sub_value;
    }

    Point P_father[4];
    P_father[0]=d_parent[parent_id].control_points[0];
    P_father[1]=d_parent[parent_id].control_points[1];
    P_father[2]=d_parent[parent_id].control_points[2];
    P_father[3]=d_parent[parent_id].control_points[3];

    Point P_child[4];
    Point P_n[SIZE];
    for(int i=0;i<SIZE;i++)
    {
        P_n[i].x=d_Q[d_child[id].original_bezier_curve_id][i][0];
        P_n[i].y=d_Q[d_child[id].original_bezier_curve_id][i][1];
        P_n[i].z=d_Q[d_child[id].original_bezier_curve_id][i][2];
        P_n[i].w=d_Q[d_child[id].original_bezier_curve_id][i][3];
    }

    double original_left=d_child[id].left * (d_parent[parent_id].right-d_parent[parent_id].left)+d_parent[parent_id].left;
    double original_right=d_child[id].right * (d_parent[parent_id].right-d_parent[parent_id].left)+d_parent[parent_id].left;

    // point_re_subdivision(P_father,d_child[id].left,d_child[id].right,P_child,P_n,p_n);
    point_re_subdivision_v2(P_father,d_child[id].left,d_child[id].right,P_child,P_n,p_n,original_left,original_right);
    double ans[2];
    // double error=calculate_trans_L1_max_error(p_n,p_m,P_n,P_child,d_child[id].left,d_child[id].right,parent_id);
    calculate_trans_L1_max_error_and_root(p_n,p_m,P_n,P_child,d_child[id].left,d_child[id].right,parent_id,ans);
    d_child[id].sub_value=ans[0];
    d_child[id].error=ans[1];
    d_child[id].dep=d_parent[parent_id].dep+1;
    d_child[id].seg=ceil(ans[1]/epsl);
    if(d_child[id].seg==1)
    {
        d_child[id].seg=0;
    }
    else
    {
        d_child[id].seg=2;
    }
    next_seg[id]=d_child[id].seg;
    // printf("next_seg[%d]=%d\n",id,next_seg[id]);
    d_child[id].control_points[0]=P_child[0];
    d_child[id].control_points[1]=P_child[1];
    d_child[id].control_points[2]=P_child[2];
    d_child[id].control_points[3]=P_child[3];

    d_child[id].original_left=original_left;
    d_child[id].original_right=original_right;
    // printf("d_child[%d].error=%f d_father.error=%f\n",id,d_child[id].error,d_parent[parent_id].error);

    if(d_child[id].seg==0)
    {
        auto old=atomicAdd(cnt,1);
        // printf("cnt=%d\n",old);
        d_final_ans[old]=d_child[id];
    }

}

__global__ void re_subdivision(int p_n,int p_m,int num_seg, CBC *d_parent,CBC *d_child,int seg_size,int *seg_sum,int *next_seg)
{
    int id=blockDim.x*blockIdx.x+threadIdx.x;
    if(id>=num_seg)
        return;
    int index_size=seg_size;
    d_child[id].id=id;
    d_child[id].parent_id=find_index(seg_sum,id,seg_size);
    d_child[id].original_bezier_curve_id=d_parent[d_child[id].parent_id].original_bezier_curve_id;

    // printf("d_child[%d].parent_id=%d\n",id,d_child[id].parent_id);
    int parent_id=d_child[id].parent_id;
    int start=seg_sum[parent_id-1];
    int end=seg_sum[parent_id]-1;
    d_child[id].left=(double)(id-start)/(end-start+1);
    d_child[id].right=(double)(id-start+1)/(end-start+1);
    // printf("d_child[%d].left= %f right= %f\n",id,d_child[id].left,d_child[id].right);

    Point P_father[4];
    P_father[0]=d_parent[parent_id].control_points[0];
    P_father[1]=d_parent[parent_id].control_points[1];
    P_father[2]=d_parent[parent_id].control_points[2];
    P_father[3]=d_parent[parent_id].control_points[3];

    Point P_child[4];
    Point P_n[SIZE];
    for(int i=0;i<SIZE;i++)
    {
        P_n[i].x=d_Q[d_child[id].original_bezier_curve_id][i][0];
        P_n[i].y=d_Q[d_child[id].original_bezier_curve_id][i][1];
        P_n[i].z=d_Q[d_child[id].original_bezier_curve_id][i][2];
        P_n[i].w=d_Q[d_child[id].original_bezier_curve_id][i][3];
    }

    // double original_left=d_child[id].left * (d_parent[parent_id].right-d_parent[parent_id].left)+d_parent[parent_id].left;
    // double original_right=d_child[id].right * (d_parent[parent_id].right-d_parent[parent_id].left)+d_parent[parent_id].left;

    double original_left=d_child[id].left * (d_parent[parent_id].original_right-d_parent[parent_id].original_left)+d_parent[parent_id].original_left;
    double original_right=d_child[id].right * (d_parent[parent_id].original_right-d_parent[parent_id].original_left)+d_parent[parent_id].original_left;

    

    // point_re_subdivision(P_father,d_child[id].left,d_child[id].right,P_child,P_n,p_n);
    point_re_subdivision_v2(P_father,d_child[id].left,d_child[id].right,P_child,P_n,p_n,original_left,original_right);
    double error=calculate_trans_L1_max_error(p_n,p_m,P_n,P_child,original_left,original_right,parent_id);
    d_child[id].error=error;
    d_child[id].dep=d_parent[parent_id].dep+1;
    d_child[id].seg=ceil(error/epsl);
    if(d_child[id].seg==1)
    {
        d_child[id].seg=0;
    }
    else
    {
        d_child[id].seg=2;
    }
    next_seg[id]=d_child[id].seg;
    // printf("next_seg[%d]=%d\n",id,next_seg[id]);
    d_child[id].control_points[0]=P_child[0];
    d_child[id].control_points[1]=P_child[1];
    d_child[id].control_points[2]=P_child[2];
    d_child[id].control_points[3]=P_child[3];

    d_child[id].original_left=original_left;
    d_child[id].original_right=original_right;
    printf("d_child[%d].error=%f d_father.error=%f\n",id,d_child[id].error,d_parent[parent_id].error);

    // if(id==20)
    // {
    //     printf("d_parent[parent_id].right=%f d_parent[parent_id].left=%f\n",d_parent[parent_id].right,d_parent[parent_id].left);
    //     printf("d_parent[parent_id].original_right=%f d_parent[parent_id].original_left=%f\n",d_parent[parent_id].original_right,d_parent[parent_id].original_left);
    //     printf("d_child[%d].original_bezier_id = %d\n",id,d_child[id].original_bezier_curve_id);
    //     printf("P_n:\n");
    //     printf("original left=%f original right=%f\n",original_left,original_right);
    //     printf("left=%f right=%f\n",d_child[id].left,d_child[id].right);
    //     for(int i=0;i<SIZE;i++)
    //     {
    //         printf("P_n[%d].x=%f\t",i,P_n[i].x);
    //         printf("P_n[%d].y=%f\t",i,P_n[i].y);
    //         printf("P_n[%d].z=%f\t",i,P_n[i].z);
    //         printf("P_n[%d].w=%f\n",i,P_n[i].w);
    //     }
    //     printf("P_father:\n");
    //     for(int i=0;i<4;i++)
    //     {
    //         printf("P_father[%d].x=%f\t",i,P_father[i].x);
    //         printf("P_father[%d].y=%f\t",i,P_father[i].y);
    //         printf("P_father[%d].z=%f\t",i,P_father[i].z);
    //         printf("P_father[%d].w=%f\n",i,P_father[i].w);
    //     }
    //     printf("P_child:\n");
    //     for(int i=0;i<4;i++)
    //     {
    //         printf("P_child[%d].x=%f\t",i,P_child[i].x);
    //         printf("P_child[%d].y=%f\t",i,P_child[i].y);
    //         printf("P_child[%d].z=%f\t",i,P_child[i].z);
    //         printf("P_child[%d].w=%f\n",i,P_child[i].w);
    //     }
    //     Point p1=evaluate_bezier_curve(original_left,P_n,p_n);
    //     Point p2=evaluate_bezier_curve(original_right,P_n,p_n);
    //     printf("p1.x=%f p1.y=%f p1.z=%f p1.w=%f\n",p1.x,p1.y,p1.z,p1.w);
    //     printf("p2.x=%f p2.y=%f p2.z=%f p2.w=%f\n",p2.x,p2.y,p2.z,p2.w);
    // }
}

__global__ void re_subdivision_at_max_error(int p_n,int p_m,int num_seg, CBC *d_parent,CBC *d_child,int seg_size,int *seg_sum,int *next_seg,int *cnt,CBC *d_final_ans)
{
    int id=blockDim.x*blockIdx.x+threadIdx.x;
    if(id>=num_seg)
        return;
    int index_size=seg_size;
    d_child[id].id=id;
    d_child[id].parent_id=find_index(seg_sum,id,seg_size);
    d_child[id].original_bezier_curve_id=d_parent[d_child[id].parent_id].original_bezier_curve_id;

    // printf("d_child[%d].parent_id=%d\n",id,d_child[id].parent_id);
    int parent_id=d_child[id].parent_id;
    int start=seg_sum[parent_id-1];
    int end=seg_sum[parent_id]-1;
    d_child[id].left=(double)(id-start)/(end-start+1);
    d_child[id].right=(double)(id-start+1)/(end-start+1);
    if(d_child[id].left==0)
    {
        d_child[id].right=d_parent[parent_id].sub_value;
    }
    else
    {
        d_child[id].left=d_parent[parent_id].sub_value;
    }
    // printf("d_child[%d].left= %f right= %f\n",id,d_child[id].left,d_child[id].right);

    Point P_father[4];
    P_father[0]=d_parent[parent_id].control_points[0];
    P_father[1]=d_parent[parent_id].control_points[1];
    P_father[2]=d_parent[parent_id].control_points[2];
    P_father[3]=d_parent[parent_id].control_points[3];

    Point P_child[4];
    Point P_n[SIZE];
    for(int i=0;i<SIZE;i++)
    {
        P_n[i].x=d_Q[d_child[id].original_bezier_curve_id][i][0];
        P_n[i].y=d_Q[d_child[id].original_bezier_curve_id][i][1];
        P_n[i].z=d_Q[d_child[id].original_bezier_curve_id][i][2];
        P_n[i].w=d_Q[d_child[id].original_bezier_curve_id][i][3];
    }

    // double original_left=d_child[id].left * (d_parent[parent_id].right-d_parent[parent_id].left)+d_parent[parent_id].left;
    // double original_right=d_child[id].right * (d_parent[parent_id].right-d_parent[parent_id].left)+d_parent[parent_id].left;

    double original_left=d_child[id].left * (d_parent[parent_id].original_right-d_parent[parent_id].original_left)+d_parent[parent_id].original_left;
    double original_right=d_child[id].right * (d_parent[parent_id].original_right-d_parent[parent_id].original_left)+d_parent[parent_id].original_left;

    

    // point_re_subdivision(P_father,d_child[id].left,d_child[id].right,P_child,P_n,p_n);
    point_re_subdivision_v2(P_father,d_child[id].left,d_child[id].right,P_child,P_n,p_n,original_left,original_right);
    // double error=calculate_trans_L1_max_error(p_n,p_m,P_n,P_child,original_left,original_right,parent_id);
    double ans[2];
    calculate_trans_L1_max_error_and_root(p_n,p_m,P_n,P_child,original_left,original_right,parent_id,ans);
    d_child[id].sub_value=ans[0];
    d_child[id].error=ans[1];
    d_child[id].dep=d_parent[parent_id].dep+1;
    d_child[id].seg=ceil(ans[1]/epsl);
    if(d_child[id].seg==1)
    {
        d_child[id].seg=0;
    }
    else
    {
        d_child[id].seg=2;
    }
    next_seg[id]=d_child[id].seg;
    // printf("next_seg[%d]=%d\n",id,next_seg[id]);
    d_child[id].control_points[0]=P_child[0];
    d_child[id].control_points[1]=P_child[1];
    d_child[id].control_points[2]=P_child[2];
    d_child[id].control_points[3]=P_child[3];

    d_child[id].original_left=original_left;
    d_child[id].original_right=original_right;
    // printf("d_child[%d].error=%f d_father.error=%f\n",id,d_child[id].error,d_parent[parent_id].error);
    if(d_child[id].seg==0)
    {
        auto old=atomicAdd(cnt,1);
        // printf("cnt=%d\n",old);
        d_final_ans[old]=d_child[id];
    }

}
// __device__ void re_subdivision(CBC* parent)
// {

// }

// __device__ void transB(int m,int p,double a,double b)
// {
//     int id=blockDim.x*blockIdx.x+threadIdx.x;
//     if(id<p||id>m-p-1)
//         return;
//     Eigen::Matrix <double,SIZE,SIZE> B_new;
//     if(b-a>eps)
//     {
//         double left[8];//p+1
//         double right[8];//p+1
//         left[0]=1;
//         right[0]=1;
//         for(int i=1;i<=p;i++)
//         {
//             left[i]=left[i-1]*(b-a);
//             right[i]=right[i-1]*a;
//         }
//         for(int i=0;i<=p;i++)
//         {
//             for(int j=0;j<=p;j++)
//             {
//                 double sum=0;
//                 for(int k=j;k<=p;k++)
//                 {
//                     sum+=C[k][j]*right[k-j]*B_N[i][k];
//                 }
//                 B_new(i,j)=left[j]*sum;
//             }
//         }
//     }
// }

__global__ void print_after_sub_bezier(int num,CBC *d_after_sub_bezier)
{
    int id=blockDim.x*blockIdx.x+threadIdx.x;
    if(id>=num)
        return;
    printf("id=%d\n",id);
    for(int i=0;i<3;i++)
    {
        printf("d_b[%d].control_points[%d].x=%f  ",id,i,d_after_sub_bezier[id].control_points[i].x);
        printf("d_b[%d].control_points[%d].y=%f  ",id,i,d_after_sub_bezier[id].control_points[i].y);
        printf("d_b[%d].control_points[%d].z=%f  ",id,i,d_after_sub_bezier[id].control_points[i].z);
        printf("d_b[%d].control_points[%d].w=%f\n",id,i,d_after_sub_bezier[id].control_points[i].w);
    }
}

__global__ void bezier_subdivision(BC *input,BC *result,int n,double *u)
{
    int id=blockDim.x*blockIdx.x+threadIdx.x;
    if(id>=n)
        return;
    double cu=u[id];
    BC temp=input[id];
    double a[4][4];
    result[id*2].control_points[0]=temp.control_points[0];
    result[id*2].control_points[1]=temp.control_points[0] * (1.0-cu) + temp.control_points[1] * cu;
    result[id*2].control_points[2]=temp.control_points[0] * (1.0-cu) * (1.0-cu) + temp.control_points[1] * 2.0 * cu * (1.0-cu) + temp.control_points[2] * cu * cu;
    result[id*2].control_points[3]=temp.control_points[0] * (1.0-cu) * (1.0-cu) * (1.0-cu) + temp.control_points[1] * 3.0 * cu * (1.0-cu) * (1.0-cu) + temp.control_points[2] * 3.0 * cu * cu * (1.0-cu) + temp.control_points[3] * cu * cu * cu;

    result[id*2+1].control_points[0]=temp.control_points[0] * (1.0-cu) * (1.0-cu) * (1.0-cu) + temp.control_points[1] * 3.0 * cu * (1.0-cu) * (1.0-cu) + temp.control_points[2] * 3.0 * cu * cu * (1.0-cu) + temp.control_points[3] * cu * cu * cu;
    result[id*2+1].control_points[1]=temp.control_points[1] * (1.0-cu) * (1.0-cu) + temp.control_points[2] * 2.0 * cu * (1.0-cu) + temp.control_points[3] * cu * cu;
    result[id*2+1].control_points[2]=temp.control_points[2] * (1.0-cu) + temp.control_points[3] * cu;
    result[id*2+1].control_points[3]=temp.control_points[3];
}

__device__ int binary_search(double var, double *knots,int m,int p)
{
    int l=p;
    int r=m-p-1;
    int mid=0;
    // __shared__ double innerknots[1000];
    // innerknots[id]=knots[id];
    // __syncthreads();
    while(r>l+1)
    {
        mid=(l+r)>>1;
        if(knots[mid]>var+eps)
            r=mid;
        else
            l=mid;
    }
    return l;
}

__device__ int binary_search_int(int var, int *arr, int n)
{
    int l=0;
    int r=n-1;
    int mid=0;
    while(l<=r)
    {
        mid=(l+r)>>1;
        if(arr[mid]>=var)
            r=mid-1;
        else
            l=mid+1;
    }
    return l;
}


__global__ void matrix_evaluate(double *numbers, Point* results,Point* control_points,int p,double *knots,int m)
{
    int id=blockDim.x*blockIdx.x+threadIdx.x;
    double t=numbers[id];
    int st=binary_search(t,knots,m,p);
    Point ans{0,0,0,0};
    for(int i=p;i>=0;i--)
    {
        double weight=0;
        double accumulator=1;
        for(int j=0;j<=p;j++)
        {
            weight+=d_N[st-i][p][st][j]*accumulator;
            accumulator*=t;
        }
        Point temp=control_points[st-i]*weight;
        ans=ans+temp;
    }
    results[id]=ans;
}

__global__ void compute_pre_G()
{
    Eigen::Matrix <double,2,4> A;
    Eigen::Matrix <double,4,4> A2;
    Eigen::Matrix <double,4,4> A3;
    Eigen::Matrix <double,4,2> B;
    Eigen::Matrix <double,2,SIZE> C;
    Eigen::Matrix <double,4,SIZE> res;
    for(int i=0;i<2;i++)
    {
        for(int j=0;j<4;j++)
        {
            A(i,j)=G_m_m[i+1][j];
        }
    }
    printf("A\n");
    for(int i=0;i<2;i++)
    {
        for(int j=0;j<4;j++)
            printf("%f ",A(i,j));
        printf("\n");
    }
    A2=A.transpose()*A;
    A3=A2.inverse();
    B=A3*A.transpose();
    for(int i=0;i<4;i++)
    {
        for(int j=0;j<2;j++)
            d_pre_G[i][j]=B(i,j);
    }
    printf("A2\n");
    for(int i=0;i<4;i++)
    {
        for(int j=0;j<4;j++)
            printf("%f ",A2(i,j));
        printf("\n");
    }
    printf("A3\n");
    for(int i=0;i<4;i++)
    {
        for(int j=0;j<4;j++)
            printf("%f ",A3(i,j));
        printf("\n");
    }
    printf("pre G\n");
    for(int i=0;i<4;i++)
    {
        for(int j=0;j<2;j++)
            printf("%f ",d_pre_G[i][j]);
        printf("\n");
    }
    for(int i=0;i<2;i++)
    {
        for(int j=0;j<SIZE;j++)
            C(i,j)=G_m_n[i+1][j];
    }
    res=B*C;
    printf("res\n");
    for(int i=0;i<4;i++)
    {
        for(int j=0;j<SIZE;j++)
            printf("%f ",res(i,j));
        printf("\n");
    }
}

__global__ void print_prefix()
{
    int id=blockDim.x*blockIdx.x+threadIdx.x;
    printf("prefix_seg[%d] = %d\n",id,seg[id]);
}


struct Result {
    int numRoots;
    double roots[4]; // 存储最多4个实根
};

__device__ void subdivide_to_Monotone_segment(CBC *former,CBC *now,Result result,int *cnt,int id)
{
    Point a[4];
    Point b[4];
    Point former_points[4];
    for(int i=0;i<4;i++)
    {
        former_points[i]=former->control_points[i];
    }
    double sum=1.0;
    double last=0.0;
    // if(id==31)
    // {
    //     printf("former_points:\n");
    //     for(int i=0;i<4;i++)
    //     {
    //         printf("former_points[%d].x=%f\t",i,former_points[i].x);
    //         printf("former_points[%d].y=%f\t",i,former_points[i].y);
    //         printf("former_points[%d].z=%f\t",i,former_points[i].z);
    //         printf("former_points[%d].w=%f\n",i,former_points[i].w);
    //     }
    // }
    for(int i=0;i<result.numRoots;i++)
    {
        auto old=atomicAdd(cnt,1);
        double t=result.roots[i];
        get_control_points(former_points,a,b,(t-last)/sum);
        former_points[0]=b[0];
        former_points[1]=b[1];
        former_points[2]=b[2];
        former_points[3]=b[3];
        now[old].control_points[0]=a[0];
        now[old].control_points[1]=a[1];
        now[old].control_points[2]=a[2];
        now[old].control_points[3]=a[3];


        // if(id==31)
        // {
        //     printf("a:  sub_t=%f\n",(t-last)/sum);
        //     for(int i=0;i<4;i++)
        //     {
        //         printf("a[%d].x=%f\t",i,a[i].x);
        //         printf("a[%d].y=%f\t",i,a[i].y);
        //         printf("a[%d].z=%f\t",i,a[i].z);
        //         printf("a[%d].w=%f\n",i,a[i].w);
        //     }
        // }


        sum=1-t;
        last=t;
    }
    auto old=atomicAdd(cnt,1);
    now[old].control_points[0]=b[0];
    now[old].control_points[1]=b[1];
    now[old].control_points[2]=b[2];
    now[old].control_points[3]=b[3];


    // if(id==31)
    // {
    //     printf("b:\n");
    //     for(int i=0;i<4;i++)
    //     {
    //         printf("b[%d].x=%f\t",i,b[i].x);
    //         printf("b[%d].y=%f\t",i,b[i].y);
    //         printf("b[%d].z=%f\t",i,b[i].z);
    //         printf("b[%d].w=%f\n",i,b[i].w);
    //     }
    // }



    // get_control_points(former_points,a,b,le);
    // get_control_points(b,now,a,(ri-le)/(1.0-le));
}

__device__ void subdivide_to_Monotone_segment_multiPoint(CBC *former,CBC *d_monotone_segments,Result result,int *cnt,int id)
{
    Point a[4];
    Point b[4];
    Point former_points[4];
    for(int i=0;i<4;i++)
    {
        former_points[i]=former->control_points[i];
    }
    double sum=1.0;
    double last=0.0;

    for(int i=0;i<result.numRoots;i++)
    {
        auto old=atomicAdd(&cnt[id],1);
        double t=result.roots[i];
        get_control_points(former_points,a,b,(t-last)/sum);
        former_points[0]=b[0];
        former_points[1]=b[1];
        former_points[2]=b[2];
        former_points[3]=b[3];
        // now[old].control_points[0]=a[0];
        // now[old].control_points[1]=a[1];
        // now[old].control_points[2]=a[2];
        // now[old].control_points[3]=a[3];

        d_monotone_segments[id*100+old].control_points[0]=a[0];
        d_monotone_segments[id*100+old].control_points[1]=a[1];
        d_monotone_segments[id*100+old].control_points[2]=a[2];
        d_monotone_segments[id*100+old].control_points[3]=a[3];
        d_monotone_segments[id*100+old].original_bezier_curve_id=former->original_bezier_curve_id;
        d_monotone_segments[id*100+old].original_left = (former->original_right- former->original_left) * last + former->original_left;
        d_monotone_segments[id*100+old].original_right = (former->original_right- former->original_left) * t + former->original_left;

        sum=1-t;
        last=t;
    }
    auto old=atomicAdd(&cnt[id],1);
    // now[old].control_points[0]=b[0];
    // now[old].control_points[1]=b[1];
    // now[old].control_points[2]=b[2];
    // now[old].control_points[3]=b[3];
    d_monotone_segments[id*100+old].control_points[0]=b[0];
    d_monotone_segments[id*100+old].control_points[1]=b[1];
    d_monotone_segments[id*100+old].control_points[2]=b[2];
    d_monotone_segments[id*100+old].control_points[3]=b[3];
    d_monotone_segments[id*100+old].original_bezier_curve_id=former->original_bezier_curve_id;
    d_monotone_segments[id*100+old].original_left = (former->original_right- former->original_left) * last + former->original_left;
    d_monotone_segments[id*100+old].original_right = former->original_right;


}
// 设备端函数：计算四次方程根
__device__ Result cal_quartic_ik_device(const double* args) 
{
    double a = args[0], b = args[1], c = args[2], d = args[3], e = args[4];

    double D = 3 * pow(b, 2) - 8 * a * c;
    double E = -pow(b, 3) + 4 * a * b * c - 8 * pow(a, 2) * d;
    double F = 3 * pow(b, 4) + 16 * pow(a, 2) * pow(c, 2) - 16 * a * pow(b, 2) * c + 16 * pow(a, 2) * b * d - 64 * pow(a, 3) * e;

    double A = D * D - 3 * F;
    double B = D * F - 9 * pow(E, 2);
    double C = F * F - 3 * D * pow(E, 2);

    double delta = B * B - 4 * A * C;  // 总判别式

    Result result = {0, {0, 0, 0, 0}};

    if (D == 0 && E == 0 && F == 0) {
        // 四重实根
        result.numRoots = 1;
        result.roots[0] = -b / (4 * a);
        return result;
    }
    if (A == 0 && B == 0 && C == 0 && D * E * F != 0) {
        // 两个实根，其中一个三重实根
        result.numRoots = 2;
        result.roots[0] = (-b * D + 9 * E) / (4 * a * D);
        result.roots[1] = (-b * D - 3 * E) / (4 * a * D);
        return result;
    }
    if (E == 0 && F == 0 && D != 0) {
        // 一对二重根
        if (D > 0) {  // 根为实数
            result.numRoots = 2;
            result.roots[0] = (-b + sqrt(D)) / (4 * a);
            result.roots[1] = (-b - sqrt(D)) / (4 * a);
        }
        return result;
    }
    if (A * B * C != 0 && delta == 0) {
        // 一对二重实根
        double x3 = (-b - copysign(1.0, A * B * E) * sqrt(D - B / A)) / (4 * a);
        double x4 = x3; // 重根
        if (A * B > 0) {  // 其余两根为不等实根
            result.numRoots = 4;
            result.roots[0] = (-b + copysign(1.0, A * B * E) * sqrt(D - B / A) + sqrt(2 * B / A)) / (4 * a);
            result.roots[1] = (-b + copysign(1.0, A * B * E) * sqrt(D - B / A) - sqrt(2 * B / A)) / (4 * a);
            result.roots[2] = x3;
            result.roots[3] = x4;
        } else {  // 其余两根为共轭虚根
            result.numRoots = 2;
            result.roots[0] = x3;
            result.roots[1] = x4;
        }
        return result;
    }
    if (delta > 0) {
        // 两个不等实根和一对共轭虚根
        double z1 = A * D + 3 * ((-B + sqrt(delta)) / 2.0);
        double z2 = A * D + 3 * ((-B - sqrt(delta)) / 2.0);

        double z = D * D - D * (copysign(1.0, z1) * pow(fabs(z1), 1.0 / 3.0) + copysign(1.0, z2) * pow(fabs(z2), 1.0 / 3.0)) +
            (copysign(1.0, z1) * pow(fabs(z1), 1.0 / 3.0) + copysign(1.0, z2) * pow(fabs(z2), 1.0 / 3.0)) * 
            (copysign(1.0, z1) * pow(fabs(z1), 1.0 / 3.0) + copysign(1.0, z2) * pow(fabs(z2), 1.0 / 3.0)) - 3 * A;

        result.numRoots = 2;
        result.roots[0] = (-b + copysign(1.0, E) * sqrt((D + copysign(1.0, z1) * pow(fabs(z1), 1.0 / 3.0) + copysign(1.0, z2) * pow(fabs(z2), 1.0 / 3.0)) / 3.0) +
            sqrt((2 * D - copysign(1.0, z1) * pow(fabs(z1), 1.0 / 3.0) - copysign(1.0, z2) * pow(fabs(z2), 1.0 / 3.0) + 2 * sqrt(z)) / 3.0)) / (4 * a);
        result.roots[1] = (-b + copysign(1.0, E) * sqrt((D + copysign(1.0, z1) * pow(fabs(z1), 1.0 / 3.0) + copysign(1.0, z2) * pow(fabs(z2), 1.0 / 3.0)) / 3.0) -
            sqrt((2 * D - copysign(1.0, z1) * pow(fabs(z1), 1.0 / 3.0) - copysign(1.0, z2) * pow(fabs(z2), 1.0 / 3.0) + 2 * sqrt(z)) / 3.0)) / (4 * a);
        return result;
    }
    if (delta < 0 && E == 0 && D > 0 && F > 0) {
        // 四个不等实根
        result.numRoots = 4;
        result.roots[0] = (-b + sqrt(D + 2 * sqrt(F))) / (4 * a);
        result.roots[1] = (-b - sqrt(D + 2 * sqrt(F))) / (4 * a);
        result.roots[2] = (-b + sqrt(D - 2 * sqrt(F))) / (4 * a);
        result.roots[3] = (-b - sqrt(D - 2 * sqrt(F))) / (4 * a);
        return result;
    }
    return result; // 默认返回空结果
}

__device__ Result filter_roots_in_range(const Result& inputResult, double lower, double upper) {
    Result outputResult = {0, {0, 0, 0, 0}}; // 初始化输出结果

    // 遍历输入结果中的解
    for (int i = 0; i < inputResult.numRoots; ++i) {
        double root = inputResult.roots[i];

        // 判断是否在 [lower, upper] 区间内
        if (root > lower && root < upper) {
            // 添加到输出结果中
            outputResult.roots[outputResult.numRoots] = root;
            outputResult.numRoots++;
        }
    }

    return outputResult;
}

__global__ void subdivide_into_monotonic_segments_using_newtwon(Point Q,int num,CBC *d_final_ans, CBC *d_monotone_segments, int *cnt)
{
    int id=blockDim.x*blockIdx.x+threadIdx.x;
    if(id>=num)
        return;
    CBC temp=d_final_ans[id];
    Point P[4];
    P[0]=temp.control_points[0];
    P[1]=temp.control_points[1];
    P[2]=temp.control_points[2];
    P[3]=temp.control_points[3];

    Point coe[4];
    coe[0]=P[0]+Q*(-1.0);
    coe[1]=P[0]*(-3)+P[1]*3;
    coe[2]=P[0]*3+P[1]*(-6)+P[2]*3;
    coe[3]=P[0]*(-1)+P[1]*3+P[2]*(-3)+P[3];

    double final_coe[7];
    for(int i=0;i<7;i++)
    {
        final_coe[i]=0;
    }
    for(int i=0;i<4;i++)
    {
        for(int j=0;j<4;j++)
        {
            // if(id==20)
            // {
            //     printf("coe[%d].dot(coe[%d])=%f\n",i,j,coe[i].dot(coe[j]));
            // }
            final_coe[i+j]+=(coe[i].dot(coe[j]));
        }
    }

    double one_differentiate[6];
    double two_differentiate[5];

    for(int i=0;i<6;i++)
    {
        one_differentiate[i] = final_coe[i + 1] * (i + 1);
    }
    for (int i = 0; i < 5; i++)
    {
        two_differentiate[i] = one_differentiate[i + 1] * (i + 1);
    }

    double reverse_two_differentiate[5];
    for(int i=0;i<5;i++)
    {
        reverse_two_differentiate[i]=two_differentiate[4-i];
    }

    Eigen::VectorXd coeffs = definePolynomial(5,reverse_two_differentiate);
    double initial_guess=0.5;
    double root = newtonIteration(coeffs, initial_guess);



    // Result res=cal_quartic_ik_device(reverse_two_differentiate);
    // Result final_res=filter_roots_in_range(res,0.0,1.0);


    // if(final_res.numRoots==0)
    // {
    //     auto old=atomicAdd(cnt,1);
    //     d_monotone_segments[old]=temp;
    // }
    // else if(final_res.numRoots>0)
    // {
    //     // printf("id=%d   final_res.numRoots=%d\n",id,final_res.numRoots);  
    //     // for(int i=0;i<final_res.numRoots;i++)
    //     // {
    //     //     printf("final_res.roots[%d]=%f\n",i,final_res.roots[i]);
    //     // }

    //     subdivide_to_Monotone_segment(&temp,d_monotone_segments,final_res,cnt,id);

    //     // printf("\n");
    // }



}

__global__ void subdivide_into_monotonic_segments(Point Q,int num,CBC *d_final_ans, CBC *d_monotone_segments, int *cnt)
{
    int id=blockDim.x*blockIdx.x+threadIdx.x;
    if(id>=num)
        return;
    CBC temp=d_final_ans[id];
    Point P[4];
    P[0]=temp.control_points[0];
    P[1]=temp.control_points[1];
    P[2]=temp.control_points[2];
    P[3]=temp.control_points[3];

    Point coe[4];
    coe[0]=P[0]+Q*(-1.0);
    coe[1]=P[0]*(-3)+P[1]*3;
    coe[2]=P[0]*3+P[1]*(-6)+P[2]*3;
    coe[3]=P[0]*(-1)+P[1]*3+P[2]*(-3)+P[3];

    double final_coe[7];
    for(int i=0;i<7;i++)
    {
        final_coe[i]=0;
    }
    for(int i=0;i<4;i++)
    {
        for(int j=0;j<4;j++)
        {
            // if(id==20)
            // {
            //     printf("coe[%d].dot(coe[%d])=%f\n",i,j,coe[i].dot(coe[j]));
            // }
            final_coe[i+j]+=(coe[i].dot(coe[j]));
        }
    }

    double one_differentiate[6];
    double two_differentiate[5];

    for(int i=0;i<6;i++)
    {
        one_differentiate[i] = final_coe[i + 1] * (i + 1);
    }
    for (int i = 0; i < 5; i++)
    {
        two_differentiate[i] = one_differentiate[i + 1] * (i + 1);
    }

    double reverse_two_differentiate[5];
    for(int i=0;i<5;i++)
    {
        reverse_two_differentiate[i]=two_differentiate[4-i];
    }

    Result res=cal_quartic_ik_device(reverse_two_differentiate);
    Result final_res=filter_roots_in_range(res,0.0,1.0);

    // if(id==31)
    // {
    //     printf("coe\n");
    //     for(int i=0;i<4;i++)
    //     {
    //         printf("coe[%d].x=%f coe[%d].y=%f coe[%d].z=%f coe[%d].w=%f\n",i,coe[i].x,i,coe[i].y,i,coe[i].z,i,coe[i].w);
    //     }
    //     printf("final_coe\n");
    //     for(int i=0;i<7;i++)
    //     {
    //         printf("final_coe[%d]=%f\n",i,final_coe[i]);
    //     }

    //     printf("one_differentiate\n");
    //     for(int i=0;i<6;i++)
    //     {
    //         printf("one_differentiate[%d]=%f\n",i,one_differentiate[i]);
    //     }

    //     printf("two_differentiate\n");
    //     for(int i=0;i<5;i++)
    //     {
    //         printf("two_differentiate[%d]=%f\n",i,two_differentiate[i]);
    //     }

    //     printf("res.numRoots=%d\n",res.numRoots);
    //     for(int i=0;i<res.numRoots;i++)
    //     {
    //         printf("res.roots[%d]=%f\n",i,res.roots[i]);
    //     }

    //     printf("final_res.numRoots=%d\n",final_res.numRoots);
    //     for(int i=0;i<final_res.numRoots;i++)
    //     {
    //         printf("final_res.roots[%d]=%f\n",i,final_res.roots[i]);
    //     }
    // }

    if(final_res.numRoots==0)
    {
        auto old=atomicAdd(cnt,1);
        d_monotone_segments[old]=temp;
    }
    else if(final_res.numRoots>0)
    {
        // printf("id=%d   final_res.numRoots=%d\n",id,final_res.numRoots);  
        // for(int i=0;i<final_res.numRoots;i++)
        // {
        //     printf("final_res.roots[%d]=%f\n",i,final_res.roots[i]);
        // }

        subdivide_to_Monotone_segment(&temp,d_monotone_segments,final_res,cnt,id);

        // printf("\n");
    }



}

__global__ void subdivide_into_monotonic_segments_multiPoint(Point *Q,int num_point, int num_seg, CBC *d_final_ans,CBC *d_monotone_segments, int *cnt_after_monotone)
{
    int num=num_point*num_seg;
    int id=blockDim.x*blockIdx.x+threadIdx.x;
    if(id>=num)
        return;
    int idx=id/num_point;
    int idy=id%num_point;
    if(idx>=num_seg||idy>=num_point)
        return;

    // auto old=atomicAdd(&cnt_after_monotone[idy],1);

    // printf("id=%d\n",id);

    Point temp_Q=Q[idy];
    CBC temp=d_final_ans[idx];

    Point P[4];
    P[0]=temp.control_points[0];
    P[1]=temp.control_points[1];
    P[2]=temp.control_points[2];
    P[3]=temp.control_points[3];

    Point coe[4];
    coe[0]=P[0]+temp_Q*(-1.0);
    coe[1]=P[0]*(-3)+P[1]*3;
    coe[2]=P[0]*3+P[1]*(-6)+P[2]*3;
    coe[3]=P[0]*(-1)+P[1]*3+P[2]*(-3)+P[3];

    double final_coe[7];
    for(int i=0;i<7;i++)
    {
        final_coe[i]=0;
    }
    for(int i=0;i<4;i++)
    {
        for(int j=0;j<4;j++)
        {
            // if(id==20)
            // {
            //     printf("coe[%d].dot(coe[%d])=%f\n",i,j,coe[i].dot(coe[j]));
            // }
            final_coe[i+j]+=(coe[i].dot(coe[j]));
        }
    }

    double one_differentiate[6];
    double two_differentiate[5];

    for(int i=0;i<6;i++)
    {
        one_differentiate[i] = final_coe[i + 1] * (i + 1);
    }
    for (int i = 0; i < 5; i++)
    {
        two_differentiate[i] = one_differentiate[i + 1] * (i + 1);
    }

    double reverse_two_differentiate[5];
    for(int i=0;i<5;i++)
    {
        reverse_two_differentiate[i]=two_differentiate[4-i];
    }

    Result res=cal_quartic_ik_device(reverse_two_differentiate);
    Result final_res=filter_roots_in_range(res,0.0,1.0);

    if(final_res.numRoots==0)
    {
        auto old=atomicAdd(&cnt_after_monotone[idy],1);
        d_monotone_segments[idy*100+old]=temp;
    }
    else if(final_res.numRoots>0)
    {

        subdivide_to_Monotone_segment_multiPoint(&temp,d_monotone_segments,final_res,cnt_after_monotone,idy);
        // printf("\n");
    }



}


__device__ double evaluate_5degree_polynomial(double *coe,double x)
{
    double ans=0;
    for(int i=0;i<6;i++)
    {
        ans+=coe[i]*pow(x,i);
    }
    return ans;
}

__device__ double evaluate_5degree_non_parameterc_bezier(double *coe,double x)
{
    double ans=0;
    ans=coe[0]*pow((1-x),5)+5*coe[1]*pow((1-x),4)*x+10*coe[2]*pow((1-x),3)*pow(x,2)+10*coe[3]*pow((1-x),2)*pow(x,3)+5*coe[4]*(1-x)*pow(x,4)+coe[5]*pow(x,5);
    return ans;
}


__global__ void elimation_criteria(int num,CBC *monotone_segment,CNPB_5 *seleted_segment,Point Q,int *cnt,int *map_seleted_monotone)
{
    int id=blockDim.x*blockIdx.x+threadIdx.x;
    if(id>=num)
        return;
    CBC temp=monotone_segment[id];
    Point P[4];
    P[0]=temp.control_points[0];
    P[1]=temp.control_points[1];
    P[2]=temp.control_points[2];
    P[3]=temp.control_points[3];

    Point coe[4];
    coe[0]=P[0]+Q*(-1.0);
    coe[1]=P[0]*(-3)+P[1]*3;
    coe[2]=P[0]*3+P[1]*(-6)+P[2]*3;
    coe[3]=P[0]*(-1)+P[1]*3+P[2]*(-3)+P[3];

    double final_coe[7];
    for(int i=0;i<7;i++)
    {
        final_coe[i]=0;
    }
    for(int i=0;i<4;i++)
    {
        for(int j=0;j<4;j++)
        {
            final_coe[i+j]+=(coe[i].dot(coe[j]));
        }
    }

    double one_differentiate[6];

    for(int i=0;i<6;i++)
    {
        one_differentiate[i] = final_coe[i + 1] * (i + 1);
    }
    double t_at_0=evaluate_5degree_polynomial(one_differentiate,0.0);
    double t_at_1=evaluate_5degree_polynomial(one_differentiate,1.0);
    if(t_at_0*t_at_1<0)
    {
        auto old=atomicAdd(cnt,1);
        seleted_segment[old].coe[0]=one_differentiate[0];
        seleted_segment[old].coe[1]=one_differentiate[1];
        seleted_segment[old].coe[2]=one_differentiate[2];
        seleted_segment[old].coe[3]=one_differentiate[3];
        seleted_segment[old].coe[4]=one_differentiate[4];
        seleted_segment[old].coe[5]=one_differentiate[5];
        // printf("%d needs clipping\n",id);
        map_seleted_monotone[old]=id;
    }
}

__global__ void elimation_criteria_multipoint(int total_cnt_after_monotone,int num_point,CBC *monotone_segment,CNPB_5 *seleted_segment,Point *Q,int *sum_of_cnt_after_monotone,int *cnt_after_elimation,int *map_seleted_monotone)
{
    int id=blockDim.x*blockIdx.x+threadIdx.x;
    if(id>=total_cnt_after_monotone)
        return;
    int idy=binary_search_int(id,sum_of_cnt_after_monotone,num_point);
    int idx;
    if(id==sum_of_cnt_after_monotone[idy])
        idy++;
    if(idy==0)
        idx=id;
    else
        idx=id-sum_of_cnt_after_monotone[idy-1];
    if(id==10||id==55||id==56||id==57||id==110||id==111||id==112||id==5500)
    {
        printf("id=%d idy=%d idx=%d\n",id,idy,idx);
    }

    Point temp_Q=Q[idy];
    CBC temp=monotone_segment[idy*100+idx];
    Point P[4];
    P[0]=temp.control_points[0];
    P[1]=temp.control_points[1];
    P[2]=temp.control_points[2];
    P[3]=temp.control_points[3];

    Point coe[4];
    coe[0]=P[0]+temp_Q*(-1.0);
    coe[1]=P[0]*(-3)+P[1]*3;
    coe[2]=P[0]*3+P[1]*(-6)+P[2]*3;
    coe[3]=P[0]*(-1)+P[1]*3+P[2]*(-3)+P[3];

    double final_coe[7];
    for(int i=0;i<7;i++)
    {
        final_coe[i]=0;
    }
    for(int i=0;i<4;i++)
    {
        for(int j=0;j<4;j++)
        {
            final_coe[i+j]+=(coe[i].dot(coe[j]));
        }
    }

    double one_differentiate[6];

    for(int i=0;i<6;i++)
    {
        one_differentiate[i] = final_coe[i + 1] * (i + 1);
    }
    double t_at_0=evaluate_5degree_polynomial(one_differentiate,0.0);
    double t_at_1=evaluate_5degree_polynomial(one_differentiate,1.0);
    if(t_at_0 * t_at_1 <= numerial_eps && t_at_0 < numerial_eps)
    {
        auto old=atomicAdd(&cnt_after_elimation[idy],1);
        seleted_segment[idy*100+old].coe[0]=one_differentiate[0];
        seleted_segment[idy*100+old].coe[1]=one_differentiate[1];
        seleted_segment[idy*100+old].coe[2]=one_differentiate[2];
        seleted_segment[idy*100+old].coe[3]=one_differentiate[3];
        seleted_segment[idy*100+old].coe[4]=one_differentiate[4];
        seleted_segment[idy*100+old].coe[5]=one_differentiate[5];


        // printf("%d needs clipping\n",id);
        map_seleted_monotone[idy*100+old]=idy*100+idx;
    }
}


__global__ void Convert_polynomial_to_non_parametric_Bezier(int num,Point Q,CNPB_5 *seleted_segment,CNPB_5 *non_parametric_Bezier,int *map_seleted_monotone, CBC *monotone_segment)
{
    int id=blockDim.x*blockIdx.x+threadIdx.x;
    if(id>=num)
        return;
    CBC cubic_bezier=monotone_segment[map_seleted_monotone[id]];

    CNPB_5 temp=seleted_segment[id];
    non_parametric_Bezier[id].coe[0]=temp.coe[0];
    non_parametric_Bezier[id].coe[1]=temp.coe[0]+temp.coe[1]*0.2;
    non_parametric_Bezier[id].coe[2]=temp.coe[0]+temp.coe[1]*0.4+temp.coe[2]*0.1;
    non_parametric_Bezier[id].coe[3]=temp.coe[0]+temp.coe[1]*0.6+temp.coe[2]*0.3+temp.coe[3]*0.1;
    non_parametric_Bezier[id].coe[4]=temp.coe[0]+temp.coe[1]*0.8+temp.coe[2]*0.6+temp.coe[3]*0.4+temp.coe[4]*0.2;
    non_parametric_Bezier[id].coe[5]=temp.coe[0]+temp.coe[1]+    temp.coe[2]+    temp.coe[3]+    temp.coe[4]+temp.coe[5];

    double ans1=evaluate_5degree_non_parameterc_bezier(non_parametric_Bezier[id].coe,0.5);
    double ans2=evaluate_5degree_polynomial(temp.coe,0.5);
    // double ans3=calculate_point_distance(evaluate_bezier_curve(0.5,cubic_bezier.control_points,3),Q);

    double t0=evaluate_5degree_non_parameterc_bezier(non_parametric_Bezier[id].coe,0.0);
    double t1=evaluate_5degree_non_parameterc_bezier(non_parametric_Bezier[id].coe,1.0);
    printf("id=%d ans1=%f ans2=%f t0=%f t1=%f\n",id,ans1,ans2,t0,t1);
    // printf("id=%d ans1=%f ans2=%f\n",id,ans1,ans2);

}

__global__ void Convert_polynomial_to_non_parametric_Bezier_multipoint(int total_elimated_num,int num_point,Point *Q,CNPB_5 *seleted_segment,CNPB_5 *non_parametric_Bezier,int *sum_of_cnt_after_elimation,int *map_seleted_monotone, CBC *monotone_segment)
{
    int id=blockDim.x*blockIdx.x+threadIdx.x;
    if(id>=total_elimated_num)
        return;
    int idy=binary_search_int(id,sum_of_cnt_after_elimation,num_point);
    int idx;
    if(id==sum_of_cnt_after_elimation[idy])
        idy++;
    if(idy==0)
        idx=id;
    else
        idx=id-sum_of_cnt_after_elimation[idy-1];
    if(id==10||id==0||id==1||id==2||id==3||id==99)
    {
        printf("id=%d idy=%d idx=%d\n",id,idy,idx);
    }

    CBC cubic_bezier=monotone_segment[map_seleted_monotone[idy*100+idx]];
    // if(id==10||id==0||id==1||id==2||id==3||id==99)
    // {
    //     printf("cubic bezier:\n");
    //     for(int i=0;i<4;i++)
    //     {
    //         printf("cubic_bezier.control_points[%d].x=%f cubic_bezier.control_points[%d].y=%f cubic_bezier.control_points[%d].z=%f cubic_bezier.control_points[%d].w=%f\n",i,cubic_bezier.control_points[i].x,i,cubic_bezier.control_points[i].y,i,cubic_bezier.control_points[i].z,i,cubic_bezier.control_points[i].w);
    //     }
    // }

    CNPB_5 temp=seleted_segment[idy*100+idx];
    non_parametric_Bezier[id].coe[0]=temp.coe[0];
    non_parametric_Bezier[id].coe[1]=temp.coe[0]+temp.coe[1]*0.2;
    non_parametric_Bezier[id].coe[2]=temp.coe[0]+temp.coe[1]*0.4+temp.coe[2]*0.1;
    non_parametric_Bezier[id].coe[3]=temp.coe[0]+temp.coe[1]*0.6+temp.coe[2]*0.3+temp.coe[3]*0.1;
    non_parametric_Bezier[id].coe[4]=temp.coe[0]+temp.coe[1]*0.8+temp.coe[2]*0.6+temp.coe[3]*0.4+temp.coe[4]*0.2;
    non_parametric_Bezier[id].coe[5]=temp.coe[0]+temp.coe[1]+    temp.coe[2]+    temp.coe[3]+    temp.coe[4]+temp.coe[5];

    // double ans1=evaluate_5degree_non_parameterc_bezier(non_parametric_Bezier[id].coe,0.5);
    // double ans2=evaluate_5degree_polynomial(temp.coe,0.5);
    // // double ans3=calculate_point_distance(evaluate_bezier_curve(0.5,cubic_bezier.control_points,3),Q);

    // double t0=evaluate_5degree_non_parameterc_bezier(non_parametric_Bezier[id].coe,0.0);
    // double t1=evaluate_5degree_non_parameterc_bezier(non_parametric_Bezier[id].coe,1.0);
    // printf("id=%d ans1=%f ans2=%f t0=%f t1=%f\n",id,ans1,ans2,t0,t1);
    // // printf("id=%d ans1=%f ans2=%f\n",id,ans1,ans2);

}

// 设备端点结构体
struct DevicePoint {
    double x, y;
};

// 计算叉积的设备端函数
__device__ double device_cross(const DevicePoint& O, const DevicePoint& A, const DevicePoint& B) {
    return (A.x - O.x)*(B.y - O.y) - (A.y - O.y)*(B.x - O.x);
}

// 设备端凸包构建函数（Andrew's Monotone Chain Algorithm）
__device__ void device_convex_hull(DevicePoint* points, int n, DevicePoint* hull, int* hull_size) {
    int k = 0;
    // Lower hull
    for(int i = 0; i < n; i++) {
        while(k >= 2 && device_cross(hull[k-2], hull[k-1], points[i]) < 0) {
            k--;
        }
        hull[k++] = points[i];
    }
    
    // Upper hull
    int t = k+1;
    for(int i = n-2; i >=0; i--){
        while(k >= t && device_cross(hull[k-2], hull[k-1], points[i]) <= 0){
            k--;
        }
        hull[k++] = points[i];
    }
    
    // Remove the last point because it is the same as the first point
    if(k > 1) k--;
    *hull_size = k;
}

// 设备端交点计算函数
__device__ void device_find_intersections(DevicePoint* hull, int hull_size, double* intersections, int* count) {
    int cnt = 0;
    for(int i = 0; i < hull_size; i++) {
        int j = (i + 1) % hull_size;
        double y1 = hull[i].y;
        double y2 = hull[j].y;
        if( (y1 > 0 && y2 < 0) || (y1 < 0 && y2 > 0) ){
            // 线性插值求交点
            double t = -y1 / (y2 - y1);
            double x = hull[i].x + t*(hull[j].x - hull[i].x);
            intersections[cnt++] = x;
        }
        // 处理点在x轴上的情况
        else if(y1 == 0 && y2 !=0){
            intersections[cnt++] = hull[i].x;
        }
    }
    *count = cnt;
}

__device__ bool findIntersection_segment_line(const DevicePoint& A, const DevicePoint& B, DevicePoint& intersection){
    // 检查是否有端点在x轴上
    if(A.y == 0 && B.y == 0){
        // 线段在x轴上，交点为整个线段
        // 这里选择输出两个端点之一作为交点
        intersection = A; // 或者选择B，根据需求调整
        return true;
    }
    else if(A.y == 0){
        // 点A在x轴上
        intersection = A;
        return true;
    }
    else if(B.y == 0){
        // 点B在x轴上
        intersection = B;
        return true;
    }
    else{
        // 检查y1和y2是否异号
        if(A.y * B.y < 0){
            // 线段跨越x轴，计算交点
            // 使用线性插值法计算t，使得y = y1 + t*(y2 - y1) = 0
            double t = -A.y / (B.y - A.y);
            // 确保t在0到1之间
            if(t < 0.0 || t > 1.0){
                // 理论上不会发生，因为y1和y2异号
                return false;
            }
            intersection.x = A.x + t * (B.x - A.x);
            intersection.y = 0.0;
            return true;
        }
        else{
            // 线段不跨越x轴
            return false;
        }
    }
}

__device__ void bezier_clipping(int num,CNPB_5 non_parametric_Bezier,double &d_xmin,double &d_xmax)
{
    int id=blockDim.x*blockIdx.x+threadIdx.x;
    if(id>=num)
        return;
    CNPB_5 temp=non_parametric_Bezier;
    DevicePoint d_points[6];
    d_points[0].x=0.0;
    d_points[0].y=temp.coe[0];
    d_points[1].x=0.2;
    d_points[1].y=temp.coe[1];
    d_points[2].x=0.4;
    d_points[2].y=temp.coe[2];
    d_points[3].x=0.6;
    d_points[3].y=temp.coe[3];
    d_points[4].x=0.8;
    d_points[4].y=temp.coe[4];
    d_points[5].x=1.0;
    d_points[5].y=temp.coe[5];

    // if(id==9)
    // {
    //     for(int i=0;i<6;i++)
    //     {
    //         printf("d_points[%d].x=%f d_points[%d].y=%f\n",i,d_points[i].x,i,d_points[i].y);
    //     }
    // }

    DevicePoint hull[12]; // 最大可能的凸包点数（6点的凸包最多12个点理论上，但实际更少）
    int hull_size = 0;
    
    // 构建凸包
    device_convex_hull(d_points, 6, hull, &hull_size);

    // if(id==9)
    // {
    //     for(int i=0;i<hull_size;i++)
    //     {
    //         printf("hull[%d].x=%f hull[%d].y=%f\n",i,hull[i].x,i,hull[i].y);
    //     }
    // }
    
    // 计算与x轴的交点
    double intersections[12];
    int count = 0;
    device_find_intersections(hull, hull_size, intersections, &count);

    // if(id==9)
    // {
    //     for(int i=0;i<count;i++)
    //     {
    //         printf("intersections[%d]=%f\n",i,intersections[i]);
    //     }
    // }
    
    // 计算最小和最大x值
    double xmin = 1e20;
    double xmax = -1e20;
    for(int i = 0; i < count; i++) {
        if(intersections[i] < xmin) xmin = intersections[i];
        if(intersections[i] > xmax) xmax = intersections[i];
    }
    
    // 处理没有交点的情况
    if(count == 0){
        xmin = xmax = 0.0; // 或者其他标识值
    }

    d_xmin=xmin;
    d_xmax=xmax;
}

__device__ void bezier_clipping_v2(int num,CNPB_5 non_parametric_Bezier,double &d_xmin,double &d_xmax)
{
    int id=blockDim.x*blockIdx.x+threadIdx.x;
    if(id>=num)
        return;
    CNPB_5 temp=non_parametric_Bezier;
    DevicePoint d_points[6];
    // d_points[0].x=0.0;
    // d_points[0].y=evaluate_5degree_non_parameterc_bezier(temp.coe,0.0);
    // d_points[1].x=0.2;
    // d_points[1].y=evaluate_5degree_non_parameterc_bezier(temp.coe,0.2);
    // d_points[2].x=0.4;
    // d_points[2].y=evaluate_5degree_non_parameterc_bezier(temp.coe,0.4);
    // d_points[3].x=0.6;
    // d_points[3].y=evaluate_5degree_non_parameterc_bezier(temp.coe,0.6);
    // d_points[4].x=0.8;
    // d_points[4].y=evaluate_5degree_non_parameterc_bezier(temp.coe,0.8);
    // d_points[5].x=1.0;
    // d_points[5].y=evaluate_5degree_non_parameterc_bezier(temp.coe,1.0);
    d_points[0].x=0.0;
    d_points[0].y=temp.coe[0];
    d_points[1].x=0.2;
    d_points[1].y=temp.coe[1];
    d_points[2].x=0.4;
    d_points[2].y=temp.coe[2];
    d_points[3].x=0.6;
    d_points[3].y=temp.coe[3];
    d_points[4].x=0.8;
    d_points[4].y=temp.coe[4];
    d_points[5].x=1.0;
    d_points[5].y=temp.coe[5];

    double xmin=1.0;
    double xmax=0.0;
    for(int i=0;i<6;i++)
    {
        for(int j=i+1;j<6;j++)
        {
            DevicePoint intersection;
            if(findIntersection_segment_line(d_points[i],d_points[j],intersection))
            {
                if(intersection.x<xmin)
                {
                    xmin=intersection.x;
                }
                if(intersection.x>xmax)
                {
                    xmax=intersection.x;
                }
            }
        }
    }
    d_xmin=xmin;
    d_xmax=xmax;

}

__device__ void bezier_clipping_aftersub(int num,CNPB_5 x_subdivided_bezier,CNPB_5 y_subdivided_bezier,double &d_xmin,double &d_xmax)
{
    int id=blockDim.x*blockIdx.x+threadIdx.x;
    if(id>=num)
        return;
    DevicePoint d_points[6];
    d_points[0].x=x_subdivided_bezier.coe[0];
    d_points[0].y=y_subdivided_bezier.coe[0];
    d_points[1].x=x_subdivided_bezier.coe[1];
    d_points[1].y=y_subdivided_bezier.coe[1];
    d_points[2].x=x_subdivided_bezier.coe[2];   
    d_points[2].y=y_subdivided_bezier.coe[2];
    d_points[3].x=x_subdivided_bezier.coe[3];
    d_points[3].y=y_subdivided_bezier.coe[3];
    d_points[4].x=x_subdivided_bezier.coe[4];
    d_points[4].y=y_subdivided_bezier.coe[4];
    d_points[5].x=x_subdivided_bezier.coe[5];
    d_points[5].y=y_subdivided_bezier.coe[5];

    // if(id==9)
    // {
    //     for(int i=0;i<6;i++)
    //     {
    //         printf("d_points[%d].x=%f d_points[%d].y=%f\n",i,d_points[i].x,i,d_points[i].y);
    //     }
    // }

    DevicePoint hull[12]; // 最大可能的凸包点数（6点的凸包最多12个点理论上，但实际更少）
    int hull_size = 0;
    
    // 构建凸包
    device_convex_hull(d_points, 6, hull, &hull_size);

    // if(id==9)
    // {
    //     for(int i=0;i<hull_size;i++)
    //     {
    //         printf("hull[%d].x=%f hull[%d].y=%f\n",i,hull[i].x,i,hull[i].y);
    //     }
    // }
    
    // 计算与x轴的交点
    double intersections[12];
    int count = 0;
    device_find_intersections(hull, hull_size, intersections, &count);

    // if(id==9)
    // {
    //     for(int i=0;i<count;i++)
    //     {
    //         printf("intersections[%d]=%f\n",i,intersections[i]);
    //     }
    // }
    
    // 计算最小和最大x值
    double xmin = 1e20;
    double xmax = -1e20;
    for(int i = 0; i < count; i++) {
        if(intersections[i] < xmin) xmin = intersections[i];
        if(intersections[i] > xmax) xmax = intersections[i];
    }
    
    // 处理没有交点的情况
    if(count == 0){
        xmin = xmax = 0.0; // 或者其他标识值
    }

    d_xmin=xmin;
    d_xmax=xmax;
}

__device__ void subdivide_5th_bezier(CNPB_5 non_parametric_Bezier,CNPB_5 &subdivided_bezier,double t0,double t2)
{
    CNPB_5 temp;
    temp.coe[0]=pow((1-t0),5)*non_parametric_Bezier.coe[0]+5*t0*pow((1-t0),4)*non_parametric_Bezier.coe[1]+10*pow(t0,2)*pow((1-t0),3)*non_parametric_Bezier.coe[2]+10*pow(t0,3)*pow((1-t0),2)*non_parametric_Bezier.coe[3]+5*pow(t0,4)*(1-t0)*non_parametric_Bezier.coe[4]+pow(t0,5)*non_parametric_Bezier.coe[5];
    temp.coe[1]=                                                pow((1-t0),4)*non_parametric_Bezier.coe[1]+4*t0*pow((1-t0),3)*non_parametric_Bezier.coe[2]+6*pow(t0,2)*pow((1-t0),2)*non_parametric_Bezier.coe[3]+4*pow(t0,3)*(1-t0)*non_parametric_Bezier.coe[4]+pow(t0,4)*non_parametric_Bezier.coe[5];
    temp.coe[2]=                                                                                                pow((1-t0),3)*non_parametric_Bezier.coe[2]+3*t0*pow((1-t0),2)*non_parametric_Bezier.coe[3]+3*pow(t0,2)*(1-t0)*non_parametric_Bezier.coe[4]+pow(t0,3)*non_parametric_Bezier.coe[5];
    temp.coe[3]=                                                                                                                                                 pow((1-t0),2)*non_parametric_Bezier.coe[3]+2*t0*(1-t0)*non_parametric_Bezier.coe[4]+pow(t0,2)*non_parametric_Bezier.coe[5];
    temp.coe[4]=                                                                                                                                                                                                (1-t0)*non_parametric_Bezier.coe[4]+t0*non_parametric_Bezier.coe[5];
    temp.coe[5]=                                                                                                                                                                                                                                        non_parametric_Bezier.coe[5];

    double t1=(t2-t0)/(1-t0);
    subdivided_bezier.coe[0]=temp.coe[0];
    subdivided_bezier.coe[1]=(1-t1)*temp.coe[0]+t1*temp.coe[1];
    subdivided_bezier.coe[2]=pow((1-t1),2)*temp.coe[0]+2*t1*(1-t1)*temp.coe[1]+pow(t1,2)*temp.coe[2];
    subdivided_bezier.coe[3]=pow((1-t1),3)*temp.coe[0]+3*t1*pow((1-t1),2)*temp.coe[1]+3*pow(t1,2)*(1-t1)*temp.coe[2]+pow(t1,3)*temp.coe[3];
    subdivided_bezier.coe[4]=pow((1-t1),4)*temp.coe[0]+4*t1*pow((1-t1),3)*temp.coe[1]+6*pow(t1,2)*pow((1-t1),2)*temp.coe[2]+4*pow(t1,3)*(1-t1)*temp.coe[3]+pow(t1,4)*temp.coe[4];
    subdivided_bezier.coe[5]=pow((1-t1),5)*temp.coe[0]+5*t1*pow((1-t1),4)*temp.coe[1]+10*pow(t1,2)*pow((1-t1),3)*temp.coe[2]+10*pow(t1,3)*pow((1-t1),2)*temp.coe[3]+5*pow(t1,4)*(1-t1)*temp.coe[4]+pow(t1,5)*temp.coe[5];


}

__global__ void _3times_bezier_clipping(int num,CNPB_5 *non_parametric_Bezier,double *distance,int *map_seleted_monotone,CBC *d_monotone_segments,Point test_point)
{
    int id=blockDim.x*blockIdx.x+threadIdx.x;
    if(id>=num)
        return;
    CBC cubic_bezier=d_monotone_segments[map_seleted_monotone[id]];
    CNPB_5 temp=non_parametric_Bezier[id];
    double d_xmin,d_xmax;
    bezier_clipping(num,temp,d_xmin,d_xmax);
    // bezier_clipping_v2(num,temp,d_xmin,d_xmax);
    double ans1=evaluate_5degree_non_parameterc_bezier(temp.coe,d_xmin);
    double ans2=evaluate_5degree_non_parameterc_bezier(temp.coe,d_xmax);
    // printf("id=%d xmin=%f xmax=%f ans1=%f ans2=%f\n",id,d_xmin,d_xmax,ans1,ans2);

    CNPB_5 x_temp;
    x_temp.coe[0]=0;
    x_temp.coe[1]=0.2;
    x_temp.coe[2]=0.4;
    x_temp.coe[3]=0.6;
    x_temp.coe[4]=0.8;
    x_temp.coe[5]=1.0;

    CNPB_5 x_subdivided_bezier;
    CNPB_5 y_subdivided_bezier;

    // printf("id=%d x_temp: %f %f %f %f %f %f\n",id,x_temp.coe[0],x_temp.coe[1],x_temp.coe[2],x_temp.coe[3],x_temp.coe[4],x_temp.coe[5]);
    // printf("id=%d temp: %f %f %f %f %f %f\n",id,temp.coe[0],temp.coe[1],temp.coe[2],temp.coe[3],temp.coe[4],temp.coe[5]);


    subdivide_5th_bezier(temp,y_subdivided_bezier,d_xmin,d_xmax);
    subdivide_5th_bezier(x_temp,x_subdivided_bezier,d_xmin,d_xmax);

    // printf("id=%d x_subdivided_bezier: %f %f %f %f %f %f\n",id,x_subdivided_bezier.coe[0],x_subdivided_bezier.coe[1],x_subdivided_bezier.coe[2],x_subdivided_bezier.coe[3],x_subdivided_bezier.coe[4],x_subdivided_bezier.coe[5]);
    // printf("id=%d y_subdivided_bezier: %f %f %f %f %f %f\n",id,y_subdivided_bezier.coe[0],y_subdivided_bezier.coe[1],y_subdivided_bezier.coe[2],y_subdivided_bezier.coe[3],y_subdivided_bezier.coe[4],y_subdivided_bezier.coe[5]);

    bezier_clipping_aftersub(num,x_subdivided_bezier,y_subdivided_bezier,d_xmin,d_xmax);
    double x_left=(d_xmin-x_subdivided_bezier.coe[0])/(x_subdivided_bezier.coe[5]-x_subdivided_bezier.coe[0]);
    double x_right=(d_xmax-x_subdivided_bezier.coe[0])/(x_subdivided_bezier.coe[5]-x_subdivided_bezier.coe[0]);
    ans1=evaluate_5degree_non_parameterc_bezier(y_subdivided_bezier.coe,x_left);
    ans2=evaluate_5degree_non_parameterc_bezier(y_subdivided_bezier.coe,x_right);

    // printf("id=%d xmin=%f xmax=%f ans1=%f ans2=%f\n",id,d_xmin,d_xmax,ans1,ans2);

    CNPB_5 x_twice_subdivided_bezier;
    CNPB_5 y_twice_subdivided_bezier;

    subdivide_5th_bezier(x_subdivided_bezier,x_twice_subdivided_bezier,x_left,x_right);
    subdivide_5th_bezier(y_subdivided_bezier,y_twice_subdivided_bezier,x_left,x_right);

    // printf("id=%d x_twice_subdivided_bezier: %f %f %f %f %f %f\n",id,x_twice_subdivided_bezier.coe[0],x_twice_subdivided_bezier.coe[1],x_twice_subdivided_bezier.coe[2],x_twice_subdivided_bezier.coe[3],x_twice_subdivided_bezier.coe[4],x_twice_subdivided_bezier.coe[5]);
    // printf("id=%d y_twice_subdivided_bezier: %f %f %f %f %f %f\n",id,y_twice_subdivided_bezier.coe[0],y_twice_subdivided_bezier.coe[1],y_twice_subdivided_bezier.coe[2],y_twice_subdivided_bezier.coe[3],y_twice_subdivided_bezier.coe[4],y_twice_subdivided_bezier.coe[5]);

    bezier_clipping_aftersub(num,x_twice_subdivided_bezier,y_twice_subdivided_bezier,d_xmin,d_xmax);
    x_left=(d_xmin-x_twice_subdivided_bezier.coe[0])/(x_twice_subdivided_bezier.coe[5]-x_twice_subdivided_bezier.coe[0]);
    x_right=(d_xmax-x_twice_subdivided_bezier.coe[0])/(x_twice_subdivided_bezier.coe[5]-x_twice_subdivided_bezier.coe[0]);
    ans1=evaluate_5degree_non_parameterc_bezier(y_twice_subdivided_bezier.coe,x_left);
    ans2=evaluate_5degree_non_parameterc_bezier(y_twice_subdivided_bezier.coe,x_right);

    // Point projection_point=evaluate_bezier_curve(d_xmax,cubic_bezier.control_points,3);
    // // printf("id=%d d_xmax=%f projection_point.x=%f projection_point.y=%f projection_point.z=%f projection_point.w=%f\n",id,d_xmax,projection_point.x,projection_point.y,projection_point.z,projection_point.w);

    // distance[id]=calculate_point_distance(projection_point,test_point);
    // double distance_at_0=calculate_point_distance(evaluate_bezier_curve(0.0,cubic_bezier.control_points,3),test_point);
    // double distance_at_1=calculate_point_distance(evaluate_bezier_curve(1.0,cubic_bezier.control_points,3),test_point);


    // printf("id=%d xmin=%f xmax=%f ans1=%f ans2=%f\n",id,d_xmin,d_xmax,ans1,ans2);
    // printf("id=%d distance=%f distance_at_0=%f distance_at_1=%f\n",id,distance[id],distance_at_0,distance_at_1);

    // if(id==0)
    // {
    //     printf("cubic bezier curve:\n");
    //     for(int i=0;i<4;i++)
    //     {
    //         printf("control_points[%d].x=%f control_points[%d].y=%f control_points[%d].z=%f control_points[%d].w=%f\n",i,cubic_bezier.control_points[i].x,i,cubic_bezier.control_points[i].y,i,cubic_bezier.control_points[i].z,i,cubic_bezier.control_points[i].w);
    //     }
    // }

}

__global__ void _3times_bezier_clipping_multipoint(int num_all_segs,int num_point,CNPB_5 *non_parametric_Bezier,double *potential_t_value,int *map_seleted_monotone,CBC *d_monotone_segments,Point *test_point,int *sum_of_cnt_after_elimation)
{
    int id=blockDim.x*blockIdx.x+threadIdx.x;
    if(id>=num_all_segs)
        return;
    
    int idy=binary_search_int(id,sum_of_cnt_after_elimation,num_point);
    int idx;
    if(id==sum_of_cnt_after_elimation[idy])
        idy++;
    if(idy==0)
        idx=id;
    else
        idx=id-sum_of_cnt_after_elimation[idy-1];
    // if(id==10||id==0||id==1||id==2||id==3||id==99)
    // {
    //     CBC cubic_bezier=d_monotone_segments[map_seleted_monotone[idy*100+idx]];
    //     printf("id=%d idy=%d idx=%d\n",id,idy,idx);
    //     printf("original cubic bezier curve id=%d original.left=%f original.right=%f\n",cubic_bezier.original_bezier_curve_id,cubic_bezier.original_left,cubic_bezier.original_right);

    // }

    CBC cubic_bezier=d_monotone_segments[map_seleted_monotone[idy*100+idx]];
    CNPB_5 temp=non_parametric_Bezier[id];
    double d_xmin,d_xmax;
    bezier_clipping(num_all_segs,temp,d_xmin,d_xmax);
    // bezier_clipping_v2(num,temp,d_xmin,d_xmax);
    double ans1=evaluate_5degree_non_parameterc_bezier(temp.coe,d_xmin);
    double ans2=evaluate_5degree_non_parameterc_bezier(temp.coe,d_xmax);
    // printf("id=%d xmin=%f xmax=%f ans1=%f ans2=%f\n",id,d_xmin,d_xmax,ans1,ans2);

    CNPB_5 x_temp;
    x_temp.coe[0]=0;
    x_temp.coe[1]=0.2;
    x_temp.coe[2]=0.4;
    x_temp.coe[3]=0.6;
    x_temp.coe[4]=0.8;
    x_temp.coe[5]=1.0;

    CNPB_5 x_subdivided_bezier;
    CNPB_5 y_subdivided_bezier;

    // printf("id=%d x_temp: %f %f %f %f %f %f\n",id,x_temp.coe[0],x_temp.coe[1],x_temp.coe[2],x_temp.coe[3],x_temp.coe[4],x_temp.coe[5]);
    // printf("id=%d temp: %f %f %f %f %f %f\n",id,temp.coe[0],temp.coe[1],temp.coe[2],temp.coe[3],temp.coe[4],temp.coe[5]);


    subdivide_5th_bezier(temp,y_subdivided_bezier,d_xmin,d_xmax);
    subdivide_5th_bezier(x_temp,x_subdivided_bezier,d_xmin,d_xmax);

    // printf("id=%d x_subdivided_bezier: %f %f %f %f %f %f\n",id,x_subdivided_bezier.coe[0],x_subdivided_bezier.coe[1],x_subdivided_bezier.coe[2],x_subdivided_bezier.coe[3],x_subdivided_bezier.coe[4],x_subdivided_bezier.coe[5]);
    // printf("id=%d y_subdivided_bezier: %f %f %f %f %f %f\n",id,y_subdivided_bezier.coe[0],y_subdivided_bezier.coe[1],y_subdivided_bezier.coe[2],y_subdivided_bezier.coe[3],y_subdivided_bezier.coe[4],y_subdivided_bezier.coe[5]);

    bezier_clipping_aftersub(num_all_segs,x_subdivided_bezier,y_subdivided_bezier,d_xmin,d_xmax);
    double x_left=(d_xmin-x_subdivided_bezier.coe[0])/(x_subdivided_bezier.coe[5]-x_subdivided_bezier.coe[0]);
    double x_right=(d_xmax-x_subdivided_bezier.coe[0])/(x_subdivided_bezier.coe[5]-x_subdivided_bezier.coe[0]);
    ans1=evaluate_5degree_non_parameterc_bezier(y_subdivided_bezier.coe,x_left);
    ans2=evaluate_5degree_non_parameterc_bezier(y_subdivided_bezier.coe,x_right);

    // printf("id=%d xmin=%f xmax=%f ans1=%f ans2=%f\n",id,d_xmin,d_xmax,ans1,ans2);

    CNPB_5 x_twice_subdivided_bezier;
    CNPB_5 y_twice_subdivided_bezier;

    subdivide_5th_bezier(x_subdivided_bezier,x_twice_subdivided_bezier,x_left,x_right);
    subdivide_5th_bezier(y_subdivided_bezier,y_twice_subdivided_bezier,x_left,x_right);

    // printf("id=%d x_twice_subdivided_bezier: %f %f %f %f %f %f\n",id,x_twice_subdivided_bezier.coe[0],x_twice_subdivided_bezier.coe[1],x_twice_subdivided_bezier.coe[2],x_twice_subdivided_bezier.coe[3],x_twice_subdivided_bezier.coe[4],x_twice_subdivided_bezier.coe[5]);
    // printf("id=%d y_twice_subdivided_bezier: %f %f %f %f %f %f\n",id,y_twice_subdivided_bezier.coe[0],y_twice_subdivided_bezier.coe[1],y_twice_subdivided_bezier.coe[2],y_twice_subdivided_bezier.coe[3],y_twice_subdivided_bezier.coe[4],y_twice_subdivided_bezier.coe[5]);

    bezier_clipping_aftersub(num_all_segs,x_twice_subdivided_bezier,y_twice_subdivided_bezier,d_xmin,d_xmax);
    x_left=(d_xmin-x_twice_subdivided_bezier.coe[0])/(x_twice_subdivided_bezier.coe[5]-x_twice_subdivided_bezier.coe[0]);
    x_right=(d_xmax-x_twice_subdivided_bezier.coe[0])/(x_twice_subdivided_bezier.coe[5]-x_twice_subdivided_bezier.coe[0]);
    ans1=evaluate_5degree_non_parameterc_bezier(y_twice_subdivided_bezier.coe,x_left);
    ans2=evaluate_5degree_non_parameterc_bezier(y_twice_subdivided_bezier.coe,x_right);

    int original_bezier_curve_id=cubic_bezier.original_bezier_curve_id;
    double original_left=cubic_bezier.original_left;
    double original_right=cubic_bezier.original_right;

    potential_t_value[id]=d_xmax;
    // if(0<id&&id<100)
    //     printf("id=%d original_bezier_curve_id=%d original_left=%f original_right=%f potential_t_value=%f\n",id,original_bezier_curve_id,original_left,original_right,potential_t_value[id]);

    // Point projection_point=evaluate_bezier_curve(d_xmax,cubic_bezier.control_points,3);
    // if(id==0||id==1||id==2||id==3||id==99)
    // {
    //     printf("id=%d d_xmax=%f projection_point.x=%f projection_point.y=%f projection_point.z=%f projection_point.w=%f\n",id,d_xmax,projection_point.x,projection_point.y,projection_point.z,projection_point.w);
    // }
    // // printf("id=%d d_xmax=%f projection_point.x=%f projection_point.y=%f projection_point.z=%f projection_point.w=%f\n",id,d_xmax,projection_point.x,projection_point.y,projection_point.z,projection_point.w);

    // distance[id]=calculate_point_distance(projection_point,test_point);
    // double distance_at_0=calculate_point_distance(evaluate_bezier_curve(0.0,cubic_bezier.control_points,3),test_point);
    // double distance_at_1=calculate_point_distance(evaluate_bezier_curve(1.0,cubic_bezier.control_points,3),test_point);


    // printf("id=%d xmin=%f xmax=%f ans1=%f ans2=%f\n",id,d_xmin,d_xmax,ans1,ans2);
    // printf("id=%d distance=%f distance_at_0=%f distance_at_1=%f\n",id,distance[id],distance_at_0,distance_at_1);
    // if(id==0||id==1||id==2||id==3||id==99)
    // {
    //     printf("id=%d ans1=%.9f ans2=%.9f\n",id,ans1,ans2);
    // }

    // if(id==0)
    // {
    //     printf("cubic bezier curve:\n");
    //     for(int i=0;i<4;i++)
    //     {
    //         printf("control_points[%d].x=%f control_points[%d].y=%f control_points[%d].z=%f control_points[%d].w=%f\n",i,cubic_bezier.control_points[i].x,i,cubic_bezier.control_points[i].y,i,cubic_bezier.control_points[i].z,i,cubic_bezier.control_points[i].w);
    //     }
    // }

}

// __device__ double atomicMinDouble(double* address, double val) //原版
// {
//     unsigned long long int* address_as_ull = (unsigned long long int*)address;
//     unsigned long long int old = *address_as_ull, assumed;
//     do {
//         assumed = old;
//         old = atomicCAS(address_as_ull, assumed,
//                         __double_as_longlong(fmin(val, __longlong_as_double(assumed))));
//     } while (assumed != old);
//     // 返回更新后的最小值
//     return __longlong_as_double(old);
// }

__device__ double atomicMinDouble(double* address, double val) //安全版
{
    unsigned long long int* address_as_ull = (unsigned long long int*)address;

    // 检查 address 是否 8 字节对齐
    assert(((uintptr_t)address % sizeof(double)) == 0);

    unsigned long long int old = *address_as_ull, assumed;
    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
                        __double_as_longlong(fmin(val, __longlong_as_double(assumed))));
    } while (assumed != old);
    return __longlong_as_double(old);
}


__global__ void get_minumum_distance(int total_cnt_after_elimation,int num_point,int *sum_of_cnt_after_elimation,double *potential_t_value,double *min_distance,double *ans,int *map_seleted_monotone,CBC *d_monotone_segments,Point *test_point,CBC *d_ans_bezier,double *d_parameters,int *d_spline_id,double *knots)
{
    int id=blockDim.x*blockIdx.x+threadIdx.x;
    if(id>=total_cnt_after_elimation)
        return;
    int idy=binary_search_int(id,sum_of_cnt_after_elimation,num_point);
    int idx;
    if(id==sum_of_cnt_after_elimation[idy])
        idy++;
    if(idy==0)
        idx=id;
    else
        idx=id-sum_of_cnt_after_elimation[idy-1];
    // if(id==10||id==0||id==1||id==2||id==3||id==99)
    // {
    //     printf("id=%d idy=%d idx=%d\n",id,idy,idx);
    // }

    CBC cubic_bezier=d_monotone_segments[map_seleted_monotone[idy*100+idx]];
    // d_ans_bezier[id]=cubic_bezier;
    
    Point projection_point=evaluate_cubic_bezier(cubic_bezier.control_points[0],cubic_bezier.control_points[1],cubic_bezier.control_points[2],cubic_bezier.control_points[3],potential_t_value[id]);
    // if(id==0||id==1||id==2||id==3||id==99)
    // {
    //     printf("id=%d projection_point.x=%f projection_point.y=%f projection_point.z=%f projection_point.w=%f\n",id,projection_point.x,projection_point.y,projection_point.z,projection_point.w);
    // }
    double distance=calculate_point_distance(projection_point,test_point[idy]);
    // if(id==0||id==1||id==2||id==3||id==99)
    // {
    //     printf("id=%d distance=%f\n",id,distance);
    // }

    auto global_min_distance=atomicMinDouble(&min_distance[idy],distance);
    // if(id==0||id==1||id==2||id==3||id==4||id==5)
    // {
    //     printf("id=%d global_min_distance=%f min_distance[%d]=%f\n",id,global_min_distance,idy,min_distance[idy]);
    // }
    if(min_distance[idy]+1e-6>=distance)
    {
        ans[idy]=potential_t_value[id];
        d_ans_bezier[idy]=cubic_bezier;
        d_parameters[idy]=(potential_t_value[id])*cubic_bezier.original_right+(1-potential_t_value[id])*cubic_bezier.original_left;
        d_spline_id[idy]=cubic_bezier.original_bezier_curve_id;
        d_parameters[idy]=knots[d_spline_id[idy]]+(knots[d_spline_id[idy]+1]-knots[d_spline_id[idy]])*d_parameters[idy];
    }
    // if(id==0||id==1||id==2||id==3||id==99)
    // {
    //     printf("id=%d min_distance=%f ans=%f\n",id,min_distance[idy],ans[idy]);
    // }

    // int right=sum_of_cnt_after_elimation[id];
    // int left;
    // if(id==0)
    //     left=0;
    // else
    //     left=sum_of_cnt_after_elimation[id-1];

    // for(int i=left;i<right;i++)
    // {
    //     CBC cubic_bezier=d_monotone_segments[map_seleted_monotone[i]];
    //     if(distance<min_distance[id])
    //     {
    //         min_distance[id]=distance;
    //         ans[id]=potential_t_value[i];
    //     }
    // }
}

__global__ void initial_min_distance(int num_point,double *min_distance,Point *test_point,int m, Point *control_points,double *d_parameters)
{
    int id=blockDim.x*blockIdx.x+threadIdx.x;
    if(id>=num_point)
        return;
    double dis1=calculate_point_distance(control_points[0],test_point[id]);
    double dis2=calculate_point_distance(control_points[m-1],test_point[id]);
    if(dis1<dis2)
    {
        min_distance[id]=dis1;
        d_parameters[id]=0.0;
    }
    else
    {
        min_distance[id]=dis2;
        d_parameters[id]=1.0;
    }
}


__global__ void get_minumum_distance_without_atomic(int total_cnt_after_elimation,int num_point,int *sum_of_cnt_after_elimation,double *potential_t_value,double *min_distance,double *ans,int *map_seleted_monotone,CBC *d_monotone_segments,Point *test_point,CBC *d_ans_bezier)
{
    int id=blockDim.x*blockIdx.x+threadIdx.x;
    if(id>=num_point)
        return;
    int left,right;
    if(id==0)
        left=0;
    else
        left=sum_of_cnt_after_elimation[id-1];
    right=sum_of_cnt_after_elimation[id];


    for(int i=0;i<right-left;i++)
    {
        CBC cubic_bezier=d_monotone_segments[map_seleted_monotone[id*100+i]];
        Point projection_point=evaluate_cubic_bezier(cubic_bezier.control_points[0],cubic_bezier.control_points[1],cubic_bezier.control_points[2],cubic_bezier.control_points[3],potential_t_value[left+i]);
        double distance=calculate_point_distance(projection_point,test_point[id]);
        if(distance<min_distance[id])
        {
            ans[id]=potential_t_value[left+i];
            min_distance[id]=distance;
            d_ans_bezier[id]=cubic_bezier;
        }

    }
    // CBC cubic_bezier=d_monotone_segments[map_seleted_monotone[idy*100+idx]];
    // d_ans_bezier[id]=cubic_bezier;
    
    // Point projection_point=evaluate_cubic_bezier(cubic_bezier.control_points[0],cubic_bezier.control_points[1],cubic_bezier.control_points[2],cubic_bezier.control_points[3],potential_t_value[id]);

    // double distance=calculate_point_distance(projection_point,test_point[idy]);


    // auto global_min_distance=atomicMinDouble(&min_distance[idy],distance);

    // if(min_distance[idy]+1e-6>=distance)
    // {
    //     ans[idy]=potential_t_value[id];
    // }

}

// __global__ void bezier_clipping_v1(int num,CNPB_5 *non_parametric_Bezier)
// {
//     int id=blockDim.x*blockIdx.x+threadIdx.x;
//     if(id>=num)
//         return;
//     CNPB_5 temp=non_parametric_Bezier[id];
//     DevicePoint d_points[6];
//     d_points[0].x=0.0;
//     d_points[0].y=evaluate_5degree_non_parameterc_bezier(temp.coe,0.0);
//     d_points[1].x=0.2;
//     d_points[1].y=evaluate_5degree_non_parameterc_bezier(temp.coe,0.2);
//     d_points[2].x=0.4;
//     d_points[2].y=evaluate_5degree_non_parameterc_bezier(temp.coe,0.4);
//     d_points[3].x=0.6;
//     d_points[3].y=evaluate_5degree_non_parameterc_bezier(temp.coe,0.6);
//     d_points[4].x=0.8;
//     d_points[4].y=evaluate_5degree_non_parameterc_bezier(temp.coe,0.8);
//     d_points[5].x=1.0;
//     d_points[5].y=evaluate_5degree_non_parameterc_bezier(temp.coe,1.0);

//     DevicePoint hull[12]; // 最大可能的凸包点数（6点的凸包最多12个点理论上，但实际更少）
//     int hull_size = 0;
    
//     // 构建凸包
//     device_convex_hull(d_points, 6, hull, &hull_size);
    
//     // 计算与x轴的交点
//     double intersections[12];
//     int count = 0;
//     device_find_intersections(hull, hull_size, intersections, &count);
    
//     // 计算最小和最大x值
//     double xmin = 1e20;
//     double xmax = -1e20;
//     for(int i = 0; i < count; i++) {
//         if(intersections[i] < xmin) xmin = intersections[i];
//         if(intersections[i] > xmax) xmax = intersections[i];
//     }
    
//     // 处理没有交点的情况
//     if(count == 0){
//         xmin = xmax = 0.0; // 或者其他标识值
//     }

//     printf("id=%d xmin=%f xmax=%f\n",id,xmin,xmax);
// }