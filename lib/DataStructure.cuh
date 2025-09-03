#pragma once
#include<vector>

#define NUM 100 //max knot num?
#define SUB_NUM 10000 //max cubic bezier num?
#define SIZE 6
#define IDX2C(i,j,ld) (((j)*(ld))+(i))

#define TOTAL_NUM 100000*1

#define eps 1e-6
#define epsl 1e-4
#define numerial_eps 1e-9

#define Max_thread 256
#define warp_size 32
#define active_warp 8

#define UINT32 unsigned int
#define ERR_NO_ERROR 0x00000000
#define ERR_NAN 0x00000001
#define ERR_INF 0x00000002
#define MAX(a, b) ((a) > (b)) ? (a) : (b)
#define MIN(a, b) ((a) < (b)) ? (a) : (b)

struct Point
{
    double x,y,z,w;
    __host__ __device__ Point(): x(0.0), y(0.0), z(0.0), w(0.0) {};
    __host__ __device__ Point(double x_val, double y_val, double z_val, double w_val):x(x_val), y(y_val), z(z_val), w(w_val) {};

    __host__ __device__ Point operator+(const Point& other) const
    {
        Point result;
        result.x = this->x + other.x;
        result.y = this->y + other.y;
        result.z = this->z + other.z;
        result.w = this->w + other.w;
        return result;
    }

    __host__ __device__ Point operator*(const double &other) const
    {
        Point result;
        result.x = this->x * other;
        result.y = this->y * other;
        result.z = this->z * other;
        result.w = this->w * other;
        return result;
    }

    __host__ __device__ Point operator-(const Point& other) const
    {
        Point result;
        result.x = this->x - other.x;
        result.y = this->y - other.y;
        result.z = this->z - other.z;
        result.w = this->w - other.w;
        return result;
    }

    __host__ __device__ double dot(const Point& other) const
    {
        return this->x * other.x + this->y * other.y + this->z * other.z + this->w * other.w;
    }

    __host__ __device__ Point cross(const Point& other) const
    {
        Point result;
        result.x = this->y * other.z - this->z * other.y;
        result.y = this->z * other.x - this->x * other.z;
        result.z = this->x * other.y - this->y * other.x;
        result.w = 1.0;
        return result;
    }
};
typedef Point P;

struct _3_order_bezier_curve
{
    Point control_points[4];
};
typedef _3_order_bezier_curve BC;
std::vector<double> knots;
std::vector<Point> control_points;
// std::vector<double> weights;
std::vector<double> new_knots;
std::vector<Point> new_control_points;

std::vector<BC> bezier_curves;

std::vector<double> X;
std::vector<std::vector<Point>> Qw;

std::vector<double> querys;
std::vector<double> us;

double N[NUM][SIZE][NUM][SIZE];
double B[NUM][SIZE][SIZE][SIZE];//第k段的Bi,j的第x个系数
double transform_N[NUM][SIZE][SIZE][SIZE];//第k个区间的经过系数变换后的Ni，j的第x个系数
double D[NUM][SIZE][SIZE][SIZE];
double Q[NUM][SIZE][4];//第k个区间上的第i个经过节点插入后的控制点
double bigtransN[NUM*SIZE*SIZE];
double bigD[NUM*SIZE*SIZE];
double big_Q[NUM*SIZE*4];
double cubic_Q[NUM][4][4];
double delta[NUM][2][4];
double pre_G[2][4];
// double inverse_B[3][3];

__device__ double d_N[NUM][SIZE][NUM][SIZE];
__device__ double d_B[NUM][SIZE][SIZE][SIZE];
__device__ double d_transform_N[NUM][SIZE][SIZE][SIZE];
__device__ double d_D[NUM][SIZE][SIZE][SIZE];
__device__ double d_Q[NUM][SIZE][4];////第k个区间上的第i个经过节点插入后的控制点的x,y,z,w坐标
__device__ double d_bigtransN[NUM*SIZE*SIZE];
__device__ double d_bigD[NUM*SIZE*SIZE];
__device__ double d_cubic_Q[NUM][4][4];
// __device__ double d_big_Q[NUM*SIZE*4];
__device__ double d_delta[NUM][2][4];
__device__ double d_pre_G[2][4];

__device__ int seg[SUB_NUM];
__device__ double prefix_error_sum[SUB_NUM];

// __device__ double d_inverse_B[SIZE*SIZE]={1,0,0,0,1,0.3333333,0,0,1,0.6666666,0.3333333,0,1,1,1,1};
// __device__ double d_inverse_B[SIZE*SIZE]={1,0,0,0,0,1,0.25,0,0,0,1,0.5,0.1666667,0,0,1,0.75,0.5,0.25,0,1,1,1,1,1};//p=4
__device__ double d_inverse_B[SIZE*SIZE]={1,0,0,0,0,0,   1,0.2,0,0,0,0,  1,0.4,0.1,0,0,0,  1,0.6,0.3,0.1,0,0,  1,0.8,0.6,0.4,0.2,0,  1,1,1,1,1,1};//p=5
// __device__ double d_inverse_B[SIZE*SIZE]={1,0,0,0,0,0,0,   1, (double)1.0/6.0 ,0,0,0,0,0,  1, (double)1.0/3.0, (double)1.0/15.0 ,0,0,0,0,  1, 0.5, 0.2, 0.05, 0,0,0,  1, (double)2.0/3.0, 0.4, 0.2,(double)1.0/15.0, 0,0,   1, (double)5.0/6.0, (double)2.0/3.0, 0.5, (double)1.0/3.0, (double)1.0/6.0, 0,  1,1,1,1,1,1,1};//p=6
double inverse_B[SIZE*SIZE];

__device__ int C[SIZE*2][SIZE*2];

__device__ double G_m_n[4][SIZE];
__device__ double G_m_m[4][4];
// __device__ Point d_monotone_segments_tr[NUM][NUM][4];


struct cubic_bezier_curve
{
    Point control_points[4];
    int dep;
    int seg;
    int id;
    double error;
    double left;
    double right;
    double original_left;
    double original_right;
    int parent_id;
    int original_bezier_curve_id;
    cubic_bezier_curve **child;
    double sub_value;
};
typedef cubic_bezier_curve CBC;

struct Coefficients_of_5_degree_non_parametric_Bezier
{
    double coe[6];
};

typedef Coefficients_of_5_degree_non_parametric_Bezier CNPB_5;

int p = 5;
const int cnt_point=1000*1;
double ground_truth[cnt_point];