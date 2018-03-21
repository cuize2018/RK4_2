//#include <iostream>
//#include <math.h>
//#include <fstream>  
//#include<omp.h>
//using namespace std;
//
//#define N 400
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
//			// return x*y;
//
//			//y' = x*cos(x)
//			//return x * cos(x);
//		}
//	} m_df;//	  微分方程
//
//	//步长h
//	double  m_h;
//	//初始点x0,y0
//	double m_x0, m_y0;
//	//x范围
//	double m_max_x;
//
//
//	RK(double xend, double x0, double y0, double h = 0.001)
//	{
//		m_max_x = xend;
//		m_x0 = x0;
//		m_y0 = y0;
//		m_h = h;
//	}//构造函数
//
//	void Solve()
//	{
//		ofstream outfile, outfile2,outfile3;
//		outfile.open("y.txt");
//		outfile2.open("x.txt");
//		outfile3.open("wucha.txt");
//
//		int i,threa_id;
//		double k1, k2, k3, k4;
//		int static t;
//		int static x_m = m_max_x / m_h;
//		cout << x_m << endl;
//
//		/*t = omp_get_num_procs();*/
//		t = 120;
//		omp_set_num_threads(t);
//
//		double yn = m_y0, xstart = m_x0;
//		double* y = new double[x_m];
//		double* x = new double[x_m];
//		double* y1 = new double[x_m];
//
//		x[0] = xstart;
//		y[0] = yn;
//
//		double start = omp_get_wtime();
//		for (int j = 1; j < t; j++)
//		{
//			int index = j * (x_m / t);
//			x[index] = xstart + index * m_h;
//			y[index] = sin(x[index]); 
//			//y[index] = exp((0.5) * pow(x[index],2));
//			y1[index] = 0;
//		}
//		
//		//并行计算
//		#pragma offload target(mic:1) inout(x,y:length(x_m) alloc_if(1) free_if(1)) 
//		#pragma omp parallel for private(threa_id,k1,k2,k3,k4)
//		for (i = 0; i < x_m; i++)
//		{
//			threa_id = omp_get_thread_num();
//
//			if ((i + 1) != (threa_id + 1) * (x_m / t))
//			{
//				x[i + 1] = x[i] + m_h;
//					//y[i + 1] = y[i] + K(x[i], y[i]) * m_h;//y(n+1)=y(n)+h*y'
//
//				k1 = m_df(x[i], y[i]);
//				k2 = m_df(x[i] + m_h / 2, y[i] + (m_h / 2) * k1);
//				k3 = m_df(x[i] + m_h / 2, y[i] + (m_h / 2) * k2);
//				k4 = m_df(x[i] + m_h, y[i] + m_h * k3);
//				y[i + 1] = y[i] + ((k1 + 2 * k2 + 2 * k3 + k4) / 6) * m_h;
//					
//			}
//		}
//
//		double end = omp_get_wtime();
//
//		for(int i = 0;i < x_m; i++)
//		{
//			y1[i] = fabs(sin(x[i]) - y[i]);
//		}
//
//		/*end = time(NULL);*/
//		cout << "time:  " << end - start << endl;
//
//		//写入文件
//		for (int j = 0; j < x_m; j++)
//		{
//			outfile << y[j] << "," << endl;
//			outfile2 << x[j] << "," << endl;
//			outfile3 << y1[j] << "," << endl;
//		}
//		outfile.close();
//		outfile2.close();
//		outfile3.close();
//	}//求解
//	 
//};
//
//int main()
//{
//	RK s1(N, 0, 0);//给出起始点和x轴范围
//	s1.Solve();//求解	
//	return 0;
//}
//

//#include<iostream>
//#include<windows.h>
//#include<time.h>
//using namespace std;
//
//int main()
//{
//	double t1 = time(NULL);
//	Sleep(1542);
//	double t2 = time(NULL);
//
//	cout << t1 << " " << t2 << " " << t2 - t1 << endl;
//	return 0;
//}
