#pragma once
#include"DataStructure.cuh"

int is_number(double x) 
{ 
	return (x == x);
}

/*************************************************
Function: is_finite_number
Description: 判断浮点数是否为inf
Input: 浮点数x
Output: 无
Return: 若浮点数x为inf返回0，否则返回1
Author: Marc Pony(marc_pony@163.com)
*************************************************/
int is_finite_number(double x) 
{ 
	return (x >= -FLT_MAX && x <= FLT_MAX);
}

/*************************************************
Function: solve_quadratic_equation
Description: 求一元二次方程(a*x^2 + b*x + c = 0)的所有实数根
Input: 方程的系数 p = {c, b, a}
Output: 方程的所有实数根x, 实数根的个数rootCount
Return: 错误号
Author: Marc Pony(marc_pony@163.com)
*************************************************/
UINT32 solve_quadratic_equation(double p[], double x[], int *rootCount)
{
	int i;
	double a, b, c, delta, sqrtDelta;
	const double ZERO = FLT_MIN;  // min normalized positive value（1.175494351e-38F）
	const double EPS = FLT_MIN;
	UINT32 errNo = ERR_NO_ERROR;

	*rootCount = 0;
	for (i = 0; i < 3; i++)
	{
		if (!is_number(p[i]))
		{
			errNo = ERR_NAN;
			return errNo;
		}
		if (!is_finite_number(p[i]))
		{
			errNo = ERR_INF;
			return errNo;
		}
	}

	a = p[2];
	b = p[1];
	c = p[0];
	
	if (fabs(a - 0.0) < EPS)
	{
		if (fabs(b - 0.0) > EPS)
		{
			x[0] = -c / b;
			*rootCount = 1;
		}
	}
	else
	{
		b /= a;
		c /= a;
		a = 1.0;

		delta = b * b - 4.0 * a * c;
		if (delta > ZERO)
		{
			if (fabs(c - 0.0) < EPS)	//若c = 0,由于计算误差,sqrt(b*b - 4*a*c）不等于|b|
			{
				x[0] = 0.0;
				x[1] = -b / a;
			}
			else
			{
				sqrtDelta = sqrt(delta);
				if (b > 0.0)
				{
					x[0] = (-2.0 * c) / (b + sqrtDelta);	//避免两个很接近的数相减,导致精度丢失
					x[1] = (-b - sqrtDelta) / (2.0 * a);
				}
				else
				{
					x[0] = (-b + sqrtDelta) / (2.0 * a);
					x[1] = (-2.0 * c) / (b - sqrtDelta);	//避免两个很接近的数相减,导致精度丢失
				}
			}
			*rootCount = 2;
		}
		else if (fabs(delta - 0.0) < EPS)
		{
			x[0] = x[1] = -b / (2.0 * a);
			*rootCount = 2;
		}
		else
		{
			*rootCount = 0;
		}
	}
	return errNo;
}


/*************************************************
Function: solve_cubic_equation
Description: 盛金公式求一元三次方程(a*x^3 + b*x^2 + c*x + d = 0)的所有实数根
			 A = b * b - 3.0 * a * c;
			 B = b * c - 9.0 * a * d;
			 C = c * c - 3.0 * b * d;
			 (1)当A = B = 0时，方程有一个三重实根
			 (2)当Δ = B^2－4 * A * C > 0时，方程有一个实根和一对共轭虚根
			 (3)当Δ = B^2－4 * A * C = 0时，方程有三个实根，其中有一个两重根
			 (4)当Δ = B^2－4 * A * C < 0时，方程有三个不相等的实根
Input: 方程的系数 p = {d, c, b, a}
Output: 方程的所有实数根x，实数根的个数rootCount
Return: 错误号
Author: Marc Pony(marc_pony@163.com)
*************************************************/
UINT32 solve_cubic_equation(double p[], double x[], int *rootCount)
{
	int i;
	double a, b, c, d, A, B, C, delta;
	double Y1, Y2, Z1, Z2, K, parm[3], roots[2], theta, T;
	const double ZERO = FLT_MIN;  // min normalized positive value（1.175494351e-38F）
	const double EPS = FLT_MIN;
	const double CALCULATE_ERROR = 1.0e-7;
	UINT32 errNo = ERR_NO_ERROR;

	*rootCount = 0;
	for (i = 0; i < 4; i++)
	{
		if (!is_number(p[i]))
		{
			errNo = ERR_NAN;
			return errNo;
		}
		if (!is_finite_number(p[i]))
		{
			errNo = ERR_INF;
			return errNo;
		}
	}

	a = p[3];
	b = p[2];
	c = p[1];
	d = p[0];

	if (fabs(a - 0.0) < EPS)
	{
		parm[2] = b;
		parm[1] = c;
		parm[0] = d;

		errNo = solve_quadratic_equation(parm, x, rootCount);
	}
	else
	{
		b /= a;
		c /= a;
		d /= a;
		a = 1.0;

		A = b * b - 3.0 * a * c;
		B = b * c - 9.0 * a * d;
		C = c * c - 3.0 * b * d;

		delta = B * B - 4.0 * A * C;

		if (fabs(A - 0.0) < EPS && fabs(B - 0.0) < EPS)
		{
			x[0] = x[1] = x[2] = -b / (3.0 * a);
			*rootCount = 3;
			return errNo;
		}

		if (delta > ZERO)
		{
			parm[2] = 1.0;
			parm[1] = B;
			parm[0] = A * C;

			errNo = solve_quadratic_equation(parm, roots, rootCount);
			if (errNo != ERR_NO_ERROR)
			{
				return errNo;
			}
			Z1 = roots[0];
			Z2 = roots[1];

			Y1 = A * b + 3.0 * a * Z1;
			Y2 = A * b + 3.0 * a * Z2;

			if (Y1 < 0.0 && Y2 < 0.0)	//pow函数的底数必须为非负数,必须分类讨论
			{
				x[0] = (-b + pow(-Y1, 1.0 / 3.0) + pow(-Y2, 1.0 / 3.0)) / (3.0*a);
			}
			else if (Y1 < 0.0 && Y2 > 0.0)
			{
				x[0] = (-b + pow(-Y1, 1.0 / 3.0) - pow(Y2, 1.0 / 3.0)) / (3.0*a);
			}
			else if (Y1 > 0.0 && Y2 < 0.0)
			{
				x[0] = (-b - pow(Y1, 1.0 / 3.0) + pow(-Y2, 1.0 / 3.0)) / (3.0*a);
			}
			else
			{
				x[0] = (-b - pow(Y1, 1.0 / 3.0) - pow(Y2, 1.0 / 3.0)) / (3.0*a);
			}
			*rootCount = 1;
		}
		else if (fabs(delta - 0.0) < EPS)
		{
			if (fabs(A - 0.0) > EPS)
			{
				K = B / A;
				x[0] = -b / a + K;
				x[1] = x[2] = -0.5 * K;
				*rootCount = 3;
			}
		}
		else
		{
			if (A > 0.0)
			{
				T = (2.0 * A * b - 3.0 * a * B) / (2.0 * pow(A, 3.0 / 2.0));
				if (T > 1.0)	//由于计算误差,T的值可能略大于1(如1.0000001)
				{
					if (T < 1.0 + CALCULATE_ERROR)
					{
						T = 1.0;
					}
					else
					{
						return errNo;
					}
				}
				if (T < -1.0)
				{
					if (T > -1.0 - CALCULATE_ERROR)
					{
						T = -1.0;
					}
					else
					{
						return errNo;
					}
				}
				theta = acos(T);
				x[0] = (-b - 2.0 * sqrt(A) * cos(theta / 3.0)) / (3.0 * a);
				x[1] = (-b + sqrt(A) * (cos(theta / 3.0) + sqrt(3.0) * sin(theta / 3.0))) / (3.0 * a);
				x[2] = (-b + sqrt(A) * (cos(theta / 3.0) - sqrt(3.0) * sin(theta / 3.0))) / (3.0 * a);
				*rootCount = 3;
			}
		}
	}
	return errNo;
}


/*************************************************
Function: solve_quartic_equation
Description: 费拉里法求一元四次方程(a*x^4 + b*x^3 + c*x^2 + d*x + e = 0)的所有实数根
Input: 方程的系数 p = {e, d, c, b, a}
Output: 方程的所有实数根x,实数根的个数rootCount
Return: 错误号
Author: Marc Pony(marc_pony@163.com)
*************************************************/
UINT32 solve_quartic_equation(double p[], double x[], int *rootCount)
{
	double a, b, c, d, e;
	double parm[4], roots[3];
	double y, M, N;
	double x1[2], x2[2];
	int rootCount1, rootCount2, i;
	double MSquare, MSquareTemp, temp, yTemp;
	const double EPS = FLT_MIN;  //min normalized positive value（1.175494351e-38F）
	UINT32 errNo = ERR_NO_ERROR;

	*rootCount = 0;
	for (i = 0; i < 5; i++)
	{
		if (!is_number(p[i]))
		{
			errNo = ERR_NAN;
			return errNo;
		}
		if (!is_finite_number(p[i]))
		{
			errNo = ERR_INF;
			return errNo;
		}
	}

	a = p[4];
	b = p[3];
	c = p[2];
	d = p[1];
	e = p[0];

	if (fabs(a - 0.0) < EPS)
	{
		if (fabs(b - 0.0) < EPS)
		{
			parm[2] = c;
			parm[1] = d;
			parm[0] = e;
			errNo = solve_quadratic_equation(parm, x, rootCount);
		}
		else
		{
			parm[3] = b;
			parm[2] = c;
			parm[1] = d;
			parm[0] = e;
			errNo = solve_cubic_equation(parm, x, rootCount);
		}
	}
	else
	{
		b /= a;
		c /= a;
		d /= a;
		e /= a;

		parm[3] = 1.0;
		parm[2] = -c;
		parm[1] = b * d - 4.0 * e;
		parm[0] = (4 * c - b * b) * e - d * d;

		errNo = solve_cubic_equation(parm, roots, rootCount);
		if (*rootCount != 0)
		{
			y = roots[0];
			MSquare = b * b + 4.0 * (y - c);
			for (i = 1; i < *rootCount; i++)
			{
				yTemp = roots[i];
				MSquareTemp = b * b + 4.0 * (yTemp - c);
				if (MSquareTemp > MSquare)
				{
					MSquare = MSquareTemp;
					y = yTemp;
				}
			}

			if (MSquare > 0.0)
			{
				if (MSquare > 1.0e-8)
				{
					M = sqrt(MSquare);
					N = b * y - 2.0 * d;
					parm[2] = 2.0;
					parm[1] = b + M;
					parm[0] = y + N / M;
					errNo = solve_quadratic_equation(parm, x1, &rootCount1);

					parm[2] = 2.0;
					parm[1] = b - M;
					parm[0] = y - N / M;
					errNo = solve_quadratic_equation(parm, x2, &rootCount2);
				}
				else
				{
					temp = y * y - 4.0 * e;
					if (temp >= 0.0)
					{
						parm[2] = 2.0;
						parm[1] = b;
						parm[0] = y + sqrt(temp);
						errNo = solve_quadratic_equation(parm, x1, &rootCount1);

						parm[2] = 2.0;
						parm[1] = b;
						parm[0] = y - sqrt(temp);
						errNo = solve_quadratic_equation(parm, x2, &rootCount2);
					}
					else
					{
						*rootCount = 0;
						return errNo;
					}
				}

				if (rootCount1 == 2)
				{
					x[0] = x1[0];
					x[1] = x1[1];
					x[2] = x2[0];
					x[3] = x2[1];
				}
				else
				{
					x[0] = x2[0];
					x[1] = x2[1];
					x[2] = x1[0];
					x[3] = x1[1];
				}
				*rootCount = rootCount1 + rootCount2;
			}
			else
			{
				*rootCount = 0;
				return errNo;
			}
		}
	}
	return errNo;
}

std::pair<int, std::vector<double>> cal_quartic_ik(const std::vector<double>& args) {
    double a = args[0], b = args[1], c = args[2], d = args[3], e = args[4];

    double D = 3 * pow(b, 2) - 8 * a * c;
    double E = -pow(b, 3) + 4 * a * b * c - 8 * pow(a, 2) * d;
    double F = 3 * pow(b, 4) + 16 * pow(a, 2) * pow(c, 2) - 16 * a * pow(b, 2) * c + 16 * pow(a, 2) * b * d - 64 * pow(a, 3) * e;

    double A = D * D - 3 * F;
    double B = D * F - 9 * pow(E, 2);
    double C = F * F - 3 * D * pow(E, 2);

    double delta = B * B - 4 * A * C;  // 总判别式

    std::vector<double> roots;

    if (D == 0 && E == 0 && F == 0) {
        // 四重实根
        double x = -b / (4 * a);
        roots.push_back(x);
        return {1, roots};
    }
    if (A == 0 && B == 0 && C == 0 && D * E * F != 0) {
        // 两个实根，其中一个三重实根
        double x1 = (-b * D + 9 * E) / (4 * a * D);
        double x234 = (-b * D - 3 * E) / (4 * a * D);
        roots.push_back(x1);
        roots.push_back(x234);
        return {2, roots};
    }
    if (E == 0 && F == 0 && D != 0) {
        // 一对二重根
        if (D > 0) {  // 根为实数
            double x13 = (-b + sqrt(D)) / (4 * a);
            double x24 = (-b - sqrt(D)) / (4 * a);
            roots.push_back(x13);
            roots.push_back(x24);
            return {2, roots};
        } else {  // 根为虚数
            return {0, {}};
        }
    }
    if (A * B * C != 0 && delta == 0) {
        // 一对二重实根
        double x3 = (-b - std::copysign(1.0, A * B * E) * sqrt(D - B / A)) / (4 * a);
        double x4 = x3; // 重根
        if (A * B > 0) {  // 其余两根为不等实根
            double x1 = (-b + std::copysign(1.0, A * B * E) * sqrt(D - B / A) + sqrt(2 * B / A)) / (4 * a);
            double x2 = (-b + std::copysign(1.0, A * B * E) * sqrt(D - B / A) - sqrt(2 * B / A)) / (4 * a);
            roots.push_back(x1);
            roots.push_back(x2);
            roots.push_back(x3);
            roots.push_back(x4);
            return {4, roots};
        } else {  // 其余两根为共轭虚根
            return {2, {x3, x4}};
        }
    }
    if (delta > 0) {
        // 两个不等实根和一对共轭虚根
        double z1 = A * D + 3 * ((-B + sqrt(delta)) / 2.0);
        double z2 = A * D + 3 * ((-B - sqrt(delta)) / 2.0);

        double z = D * D - D * (std::copysign(1.0, z1) * pow(fabs(z1), 1.0 / 3.0) + std::copysign(1.0, z2) * pow(fabs(z2), 1.0 / 3.0)) +
            (std::copysign(1.0, z1) * pow(fabs(z1), 1.0 / 3.0) + std::copysign(1.0, z2) * pow(fabs(z2), 1.0 / 3.0)) * 
            (std::copysign(1.0, z1) * pow(fabs(z1), 1.0 / 3.0) + std::copysign(1.0, z2) * pow(fabs(z2), 1.0 / 3.0)) - 3 * A;

        double x1 = (-b + std::copysign(1.0, E) * sqrt((D + std::copysign(1.0, z1) * pow(fabs(z1), 1.0 / 3.0) + std::copysign(1.0, z2) * pow(fabs(z2), 1.0 / 3.0)) / 3.0) +
            sqrt((2 * D - std::copysign(1.0, z1) * pow(fabs(z1), 1.0 / 3.0) - std::copysign(1.0, z2) * pow(fabs(z2), 1.0 / 3.0) + 2 * sqrt(z)) / 3.0)) / (4 * a);

        double x2 = (-b + std::copysign(1.0, E) * sqrt((D + std::copysign(1.0, z1) * pow(fabs(z1), 1.0 / 3.0) + std::copysign(1.0, z2) * pow(fabs(z2), 1.0 / 3.0)) / 3.0) -
            sqrt((2 * D - std::copysign(1.0, z1) * pow(fabs(z1), 1.0 / 3.0) - std::copysign(1.0, z2) * pow(fabs(z2), 1.0 / 3.0) + 2 * sqrt(z)) / 3.0)) / (4 * a);

        // 虚根忽略
        return {2, {x1, x2}};
    }
    if (delta < 0) {
        if (E == 0) {
            if (D > 0 && F > 0) {
                // 四个不等实根
                double x1 = (-b + sqrt(D + 2 * sqrt(F))) / (4 * a);
                double x2 = (-b - sqrt(D + 2 * sqrt(F))) / (4 * a);
                double x3 = (-b + sqrt(D - 2 * sqrt(F))) / (4 * a);
                double x4 = (-b - sqrt(D - 2 * sqrt(F))) / (4 * a);
                return {4, {x1, x2, x3, x4}};
            } else {
                // 两对不等共轭虚根
                std::cout << "两对不等共轭虚根" << std::endl;
                return {0, {}};
            }
        } else {
            if (D > 0 && F > 0) {
                // 四个不等实根
                double theta = acos((3 * B - 2 * A * D) / (2 * A * sqrt(A)));
                double y1 = (D - 2 * sqrt(A) * cos(theta / 3.0)) / 3.0;
                double y2 = (D + sqrt(A) * (cos(theta / 3.0) + sqrt(3) * sin(theta / 3.0))) / 3.0;
                double y3 = (D + sqrt(A) * (cos(theta / 3.0) - sqrt(3) * sin(theta / 3.0))) / 3.0;

                double x1 = (-b + std::copysign(1.0, E) * sqrt(y1) + (sqrt(y2) + sqrt(y3))) / (4 * a);
                double x2 = (-b + std::copysign(1.0, E) * sqrt(y1) - (sqrt(y2) + sqrt(y3))) / (4 * a);
                double x3 = (-b - std::copysign(1.0, E) * sqrt(y1) + (sqrt(y2) - sqrt(y3))) / (4 * a);
                double x4 = (-b - std::copysign(1.0, E) * sqrt(y1) - (sqrt(y2) - sqrt(y3))) / (4 * a);

                return {4, {x1, x2, x3, x4}};
            } else {
                // 两对不等共轭虚根
                std::cout << "两对不等共轭虚根" << std::endl;
                return {0, {}};
            }
        }
    }
    return {0, {}};
}


void solvePolynomial(float a, float b, float c, float d, float e) 
{
    // 设置多项式系数
    Eigen::VectorXf coeffs(5);
    coeffs(4) = a;
    coeffs(3) = b;
    coeffs(2) = c;
    coeffs(1) = d;
    coeffs(0) = e;

    // 构造矩阵
    Eigen::Matrix4f A;
    A << -coeffs(3) / coeffs(4), -coeffs(2) / coeffs(4), -coeffs(1) / coeffs(4), -coeffs(0) / coeffs(4),
         1, 0, 0, 0,
         0, 1, 0, 0,
         0, 0, 1, 0;

    // 求解特征值
    Eigen::EigenSolver<Eigen::MatrixXf> ES(A);
    std::cout << "Roots (Eigen method):\n" << ES.eigenvalues().transpose() << std::endl;

}


// void solveQuartic22(double a, double b, double c, double d, double e) {
//     // 归一化方程
//     double A = b / a;
//     double B = c / a;
//     double C = d / a;
//     double D = e / a;

//     // 求解辅助变量
//     double p = B - (A * A) / 4;
//     double q = A * A * A / 64 - A * B / 16 + C;
//     double r = -D + A * A * A * A / 256 - A * A * B / 64 + A * C / 16;

//     // 求解三次方程的根
//     double delta0 = p * p * p / 27 + q * q / 4 + r;
//     double delta1 = 3 * p;
//     double C_root = std::cbrt(delta1 / 2 + std::sqrt(delta1 * delta1 / 4 + delta0));

//     // 求解
//     double y1 = -1 * (C_root + delta0 / C_root) / 3 - A / 4;
//     double y2 = -1 * (C_root - delta0 / C_root) / 3 - A / 4;

//     std::vector<double> roots;

//     // 通过求解二次方程来找出四次方程的根
//     double discriminant1 = y1 * y1 - 4 * B;
//     double discriminant2 = y2 * y2 - 4 * B;

//     if (discriminant1 >= 0) {
//         roots.push_back((-y1 + std::sqrt(discriminant1)) / 2);
//         roots.push_back((-y1 - std::sqrt(discriminant1)) / 2);
//     }
    
//     if (discriminant2 >= 0) {
//         roots.push_back((-y2 + std::sqrt(discriminant2)) / 2);
//         roots.push_back((-y2 - std::sqrt(discriminant2)) / 2);
//     }
//     printf("root.size()=%d\n",roots.size());
//     // 输出结果
//     printf("实数根为：\n");
//     for (double root : roots) {
//         printf("%f\n", root);
//     }
// }

// void findRealRoots(double a, double b, double c, double d, double e) {
//     // 规范化系数
//     double A = b / a;
//     double B = c / a;
//     double C = d / a;
//     double D = e / a;

//     // 计算辅方程的系数
//     double p = B - A * A / 3;
//     double q = (2 * A * A * A) / 27 - (A * B) / 3 + C;
//     double r = D - (A * C) / 3 + (A * A * B) / 9 - (A * A * A * A) / 81;

//     // 计算判别式
//     double Q = (p * p * p) / 27 + (q * q) / 4;
//     if (Q > 0) {
//         // 有一个实根
//         double u = cbrt(-q / 2 + sqrt(Q));
//         double v = cbrt(-q / 2 - sqrt(Q));
//         double root1 = u + v;

//         // 计算双重根
//         double A2 = A / 2;
//         double delta = A2 * A2 - p / 3;
//         double root2, root3;
//         if (delta >= 0) {
//             root2 = A2 + sqrt(delta);
//             root3 = A2 - sqrt(delta);
//             std::cout << "Roots: " << root1 << ", " << root2 << ", " << root3 << std::endl;
//         } else {
//             std::cout << "Roots: " << root1 << std::endl;
//         }
//     } else {
//         // 有三个实根
//         double theta = acos(-q / (2 * sqrt(-p * p * p / 27)));
//         double root1 = 2 * sqrt(-p / 3) * cos(theta / 3);
//         double root2 = 2 * sqrt(-p / 3) * cos((theta + 2 * M_PI) / 3);
//         double root3 = 2 * sqrt(-p / 3) * cos((theta + 4 * M_PI) / 3);
//         std::cout << "Roots: " << root1 << ", " << root2 << ", " << root3 << std::endl;
//     }
// }


//一元4次方程求解
void solveQuartic(double a, double b, double c, double d, double e) 
{
    // 创建多项式系数的矩阵
    Eigen::MatrixXd M(5, 5);
    M << a, b, c, d, e,
         0, 1, 0, 0, 0,
         0, 0, 1, 0, 0,
         0, 0, 0, 1, 0,
         0, 0, 0, 0, 1; // 线性项的单位矩阵

    // 特征值分解
    Eigen::EigenSolver<Eigen::MatrixXd> es(M);
    auto roots = es.eigenvalues();

    // 打印实数解
    printf("Real roots:\n");
    for (const auto& root : roots) 
    {
        if (root.imag() == 0) 
        {
            printf("%lf\n", root.real());
        }
    }
}

__device__ Eigen::VectorXd definePolynomial(int n,double *a) 
{
    // 多项式 f(x) = 2x^3 - 3x^2 - 12x + 1
    Eigen::VectorXd coeffs(n);
    for(int i=0;i<n;i++)
        coeffs[i]=a[n-i-1];
    return coeffs;
}

// 求导数
__device__ Eigen::VectorXd derivative(const Eigen::VectorXd& coeffs) {
    Eigen::VectorXd derivativeCoeffs = Eigen::VectorXd::Zero(coeffs.size() - 1);
    for (int i = 0; i < coeffs.size() - 1; ++i) {
        derivativeCoeffs[i] = coeffs[i] * (coeffs.size() - 1 - i);
    }
    return derivativeCoeffs;
}

// 计算多项式值
__device__ double evaluatePolynomial(const Eigen::VectorXd& coeffs, double x) {
    double result = 0.0;
    for (int i = 0; i < coeffs.size(); ++i) {
        result = result * x + coeffs[i];
    }
    return result;
}

// 定义多项式函数和导数函数
__device__ double polynomial(const Eigen::VectorXd& coeffs, double x) {
    double result = 0.0;
    for (int i = 0; i < coeffs.size(); ++i) {
        result = result * x + coeffs[i];
    }
    return result;
}

__device__ double polynomialDerivative(const Eigen::VectorXd& coeffs, double x) {
    double result = 0.0;
    for (int i = 0; i < coeffs.size() - 1; ++i) {
        result = result * x + coeffs[i] * (coeffs.size() - 1 - i);
    }
    return result;
}

// 牛顿迭代法求解多项式根
__device__ double newtonIteration(const Eigen::VectorXd& coeffs, double initial_guess, int max_iter = 100, double tol = 1e-8) 
{
    double x = initial_guess;
    for (int i = 0; i < max_iter; ++i) 
    {
        double fx = polynomial(coeffs, x);
        double fpx = polynomialDerivative(coeffs, x);

        if (std::abs(fpx) < 1e-8) 
        {
            printf("Derivative too small, stopping iteration.\n");
            break;
        }

        double x_new = x - fx / fpx;

        if (std::abs(x_new - x) < tol) 
        {
            if(x_new<0)
                return 0;
            else if(x_new>1)
                return 1;
            else
                return x_new;
        }

        x = x_new;
    }
    printf("Maximum iterations reached, returning last computed value.\n");

    // std::cerr << "Maximum iterations reached, returning last computed value." << std::endl;
    if(x<0)
        return 0;
    else if(x>1)
        return 1;
    else
        return x;
}

// // 构建伴随矩阵
// __device__ Eigen::MatrixXd buildCompanionMatrix(const Eigen::VectorXd& coeffs) {
//     int degree = coeffs.size() - 1;
//     Eigen::MatrixXd companionMatrix = Eigen::MatrixXd::Zero(degree, degree);
    
//     for (int i = 1; i < degree; ++i) {
//         companionMatrix(i, i - 1) = 1.0;
//     }

//     companionMatrix.col(degree - 1) = -coeffs.head(degree) / coeffs(degree);
    
//     return companionMatrix;
// }

// // 求解多项式的根
// __device__ Eigen::VectorXcd findPolynomialRoots(const Eigen::VectorXd& coeffs) {
//     Eigen::MatrixXd companionMatrix = buildCompanionMatrix(coeffs);
//     Eigen::EigenSolver<Eigen::MatrixXd> eigenSolver(companionMatrix);
//     return eigenSolver.eigenvalues();
// }