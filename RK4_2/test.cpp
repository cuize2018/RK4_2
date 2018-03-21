//#include<iostream>
//#include<omp.h>
//using  namespace std;
//
//int main()
//{
//	int t,threa_id;
//	t = omp_get_num_procs();
//	omp_set_num_threads(t);
//	int sum[8] = {0};
//	
//	#pragma omp parallel for private(threa_id)
//	for (int i = 0; i < 100; i++) 
//	{
//		threa_id = omp_get_thread_num();
//		sum[threa_id] = sum[threa_id] + (i+1);
//	}
//	int sum1 = 0;
//	for (int i = 0; i < 8; i++)
//	{
//		sum1 = sum1 + sum[i];
//	}
//	for (int i = 0; i < 8; i++)
//	{
//		cout << sum[i] << endl;
//	}
//	cout << sum1 << endl;
//
//	return 0;
//}