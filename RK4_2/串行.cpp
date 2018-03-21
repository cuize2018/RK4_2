//#include <iostream>
//using namespace std;
//
//#include <vector>
//#include <math.h>
//#include <fstream>  
//#include<omp.h>
//using namespace std;
//
//#define N 80
//
//class RK
//{
//public:
//	class DiffFunc
//	{
//	public:
//		double operator()(double x, double y)
//		{
//			// y'=cos(x)
//			return cos(x);
//			//return 2 * x;
//			// y'=x*y-1
//			// return x*y-1;
//
//			//y' = x*cos(x)
//			//return x * cos(x);
//		}
//	} m_df;//	  微分方程
//
//		   //步长h
//	double m_h;
//	//初始点x0,y0
//	double m_x0, m_y0;
//	//x范围
//	double m_max_x;
//
//	RK(double xend, double x0, double y0, double h = 1e-3)
//	{
//		m_max_x = xend;
//		m_x0 = x0;
//		m_y0 = y0;
//		m_h = h;
//		m_half_h = h / 2.0;
//	}//构造函数
//
//	void Solve()
//	{
//		ofstream outfile, outfile2;
//		outfile.open("y.txt");
//		outfile2.open("x.txt");
//
//		int i;
//		int x_m = m_max_x / m_h;
//		cout << x_m << endl;
//
//		double yn = m_y0, xstart = m_x0;
//		double* y = new double[x_m];
//		double* x = new double[x_m];
//
//		x[0] = xstart;
//		y[0] = yn;
//
//		double start = omp_get_wtime();
//		//串行计算y轴
//		//#pragma omp parallel for
//		for (i = 0; i < x_m; i++)
//		{
//
//			y[i + 1] = y[i] + K(x[i], y[i])*m_h;//y(n+1)=y(n)+h*y'
//			x[i + 1] = x[i] + m_h;
//		}
//		double end = omp_get_wtime();
//		cout << "time:  " << end - start << endl;
//		////写入文件
//		for (int j = 0; j < x_m; j++)
//		{
//			outfile << y[j] << "," << endl;
//			outfile2 << x[j] << "," << endl;
//		}
//		outfile.close();
//		outfile2.close();
//	}//求解
//
//
//private:
//	double m_half_h;
//	double m_ptx, m_pty;
//	double  K1(double x, double y)
//	{
//		double v = m_df(x, y);
//		m_ptx = x;
//		m_pty = y;
//		return v;
//	}
//	double  K2(double _k1)
//	{
//		double v = m_df(m_ptx + m_half_h, m_pty + m_half_h * _k1);
//		return v;
//	}
//	double  K3(double _k2)
//	{
//
//		double v = m_df(m_ptx + m_half_h, m_pty + m_half_h * _k2);
//		return v;
//	}
//	double  K4(double _k3)
//	{
//		double v = m_df(m_ptx + m_h, m_pty + m_h * _k3);
//		return v;
//	}
//	double  K(double x, double y)
//	{
//		double _k1 = K1(x, y), _k2 = K2(_k1), _k3 = K3(_k2), _k4 = K4(_k3);
//		return  (_k1 + 2.0*_k2 + 2.0*_k3 + _k4) / 6.0;
//	}
//
//};
//
//int main()
//{
//	
//	RK s1(N, 0, 0);//给出起始点和x轴范围
//	s1.Solve();//求解
//	
//
//	
//	return 0;
//}
//
