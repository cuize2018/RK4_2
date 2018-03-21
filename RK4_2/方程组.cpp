////二元一阶微分方程组
//#include <iostream>
//#include <vector>
//#include <math.h>
//#include <fstream>  
//#include<omp.h>
//#include<ctime>
//using namespace std;
//
//#define N 10000
//
//class DiffFunc
//{
//public:
//	double xr, y1r, y2r;
//	DiffFunc(double x, double y1, double y2)
//	{
//		xr = x;
//		y1r = y1;
//		y2r = y2;
//	}
//	double fun1(double x, double y1, double y2) //y1' = y2
//	{
//		return y2;
//	}
//	double fun2(double x, double y1, double y2)//y2' = -y1 + 3
//	{
//		return -1 * y1 + 3;
//	}
//} m_df;//	  微分方程
//
//class RK
//{
//public:	  
//	double m_h; //步长h
//	//初始点x0,y0
//	double m_x0, m_y0,m_y1;
//	//x范围
//	double m_max_x;
//
//	RK(double xend, double x0, double y0,double y1, double h = 1e-3)
//	{
//		m_max_x = xend;
//		m_x0 = x0;
//		m_y0 = y0;
//		m_y1 = y1;
//		m_h = h;
//		m_half_h = h / 2.0;
//	}//构造函数
//
//	void Solve()
//	{
//		/*ofstream outfile, outfile2;
//		outfile.open("y.txt");
//		outfile2.open("x.txt");*/
//
//		int i, t, threa_id;
//		int x_m = m_max_x / m_h;
//		cout << x_m << endl;
//
//		t = omp_get_num_procs();
//		omp_set_num_threads(t);
//
//		double yn1 = m_y0,yn2 = m_y1, xstart = m_x0;
//		
//		double* x = new double[x_m];
//		double* y1 = new double[x_m];
//		double* y2 = new double[x_m];
//
//		x[0] = xstart;
//		y1[0] = yn1;
//		y2[0] = yn2;
//		double start = omp_get_wtime();
//		/*time_t start, end;*/
//		/*start = time(NULL);*/
//		for (int j = 1; j < t; j++)
//		{
//			x[j] = xstart + j * m_h;
//		}
//
//		//并行计算出x轴
//		#pragma omp parallel for private(threa_id)
//		for (int i = 0; i < x_m; i++)
//		{
//			threa_id = omp_get_thread_num();
//			x[threa_id + (i % (x_m / t))*t] = x[threa_id] + m_h * t*(i % (x_m / t));
//		}
//		//串行计算y轴
//		#pragma omp parallel for private(threa_id)
//		for (i = 0; i < x_m; i++)
//		{
//			threa_id = omp_get_thread_num();
//			if (threa_id == 0)
//			{
//				y1[i + 1] = y1[i] + cal_avgk(x[i], y1[i], y2[i], true) * m_h;//y(n+1)=y(n)+h*y'
//			}
//			if (threa_id == 1)
//			{
//				y2[i + 1] = y2[i] + cal_avgk(x[i], y1[i], y2[i], false) * m_h;//y(n+1)=y(n)+h*y'
//			}
//		}
//		double end = omp_get_wtime();
//		/*end = time(NULL);*/
//		cout << "time:  " << end - start << endl;
//
//		////写入文件
//		//for (int j = 0; j < x_m; j++)
//		//{
//		//	outfile << y[j] << "," << endl;
//		//	outfile2 << x[j] << "," << endl;
//		//}
//		//outfile.close();
//		//outfile2.close();
//	}//求解
//
//
//private:
//	double m_half_h;
//	double m_ptx, m_pty1,m_pty2;
//	double m_ptxg, m_ptyg1, m_ptyg2;
//	double  K1(double x ,double y1, double y2)
//	{
//		double v = m_df.fun1(x, y1,y2);
//		m_ptx = x;
//		m_pty1 = y1;
//		m_pty2 = y2;
//
//		return v;
//	}
//	double  K2(double _k1,double _g1)
//	{
//		double v = m_df.fun1(m_ptx + m_half_h, m_pty1 + m_half_h * _k1,m_pty2 + m_half_h * _g1);
//		return v;
//	}
//	double  K3(double _k2,double _g2)
//	{
//
//		double v = m_df.fun1(m_ptx + m_half_h, m_pty1 + m_half_h * _k2, m_pty2 + m_half_h * _g2);
//		return v;
//	}
//	double  K4(double _k3,double _g3)
//	{
//		double v = m_df.fun1(m_ptx + m_h, m_pty1 + m_h * _k3, m_pty2 + m_half_h * _g3);
//		return v;
//	}
//	//-----------------------------------------------
//	double  G1(double x, double y1, double y2)
//	{
//		double v = m_df.fun2(x, y1, y2);
//		m_ptx = x;
//		m_pty1 = y1;
//		m_pty2 = y2;
//
//		return v;
//	}
//	double  G2(double _k1, double _g1)
//	{
//		double v = m_df.fun2(m_ptx + m_half_h, m_pty1 + m_half_h * _k1, m_pty2 + m_half_h * _g1);
//		return v;
//	}
//	double  G3(double _k2, double _g2)
//	{
//
//		double v = m_df.fun2(m_ptx + m_half_h, m_pty1 + m_half_h * _k2, m_pty2 + m_half_h * _g2);
//		return v;
//	}
//	double  G4(double _k3, double _g3)
//	{
//		double v = m_df.fun2(m_ptx + m_h, m_pty1 + m_h * _k3, m_pty2 + m_half_h * _g3);
//		return v;
//	}
//	double  cal_avgk(double x, double y1, double y2,bool flag)
//	{
//		double _g1 = G1(x, y1, y2);
//		double _k1 = K1(x, y1, y2);
//		
//		double _g2 = G2(_k1, _g1);
//		double _k2 = K2(_k1, _g1);
//
//		double _g3 = G3(_k2, _g2);
//		double _k3 = K3(_k2, _g2);
//
//		double _g4 = G4(_k3, _g3);
//		double _k4 = K4(_k3, _g3);
//
//		if(flag == true)
//		{
//			return (_k1 + 2.0*_k2 + 2.0*_k3 + _k4) / 6.0;
//		}
//		else
//		{ 
//			return (_g1 + 2.0*_g2 + 2.0*_g3 + _g4) / 6.0;
//		}
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
//	return 0;
//}
//
