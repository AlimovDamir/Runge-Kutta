
#include "stdio.h"
#include "stdlib.h"
#include "mpi.h"
#include <iostream>
#include <fstream>
#include <ctime>
#include <cmath>

using namespace std;

double test(double* a, double* b, int n);

int main(int argc, char* argv[])
{
	double** y;
	int n = 400;
	double** b;
	double** buf2; 
	int t_0 = 0;
	int t_1 = 2;
	int nc = 500;
	int p = 4;
	int len;
	ofstream fout;
	MPI_Status status;
	char size[100];
	double h,h1,h2;
	double* k1,*k2,*k3,*k4, *buf;

    p = atoi(argv[1]);
    n = atoi(argv[2]);
	h = t_1 - t_0;
	h = h/nc;
	k1 = new double[n];
	k2 = new double[n];
	k3 = new double[n];
	k4 = new double[n];
	b = new double*[n];
	buf2 = new double*[p];
	buf = new double[n/p*2];
	y = new double*[nc];
	for(int i = 0; i <= n - 1; i++)
	{
		b[i] = new double[n];
		for(int j = 0; j <= n - 1; j++)
	    {
			b[i][j] =  rand()%5 - 5;
		}
	}
	for(int i = 0; i <= nc-1; i++)
	{
		y[i] = new double[n];
	}
	for(int i = 0; i <= n - 1; i++)
	{
		y[0][i] = rand()%10 + 1.5;
	}
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &len);
    char lent[50];
    if (len == 0)
        {sprintf(lent,"test_%d_%d.txt",p,n);
        fout.open(lent); 
        fout<<p<<endl;
        fout<<n<<endl;
    }
	for(int i = 0; i <= p - 1; i++)
	{
		if (i != len)
		{
			buf2[i] = new double[n/p*2];
		}
	}
	for(int i = 0; i <= nc - 2; i++)
	{
			for(int j = n/p*len; j <= n/p*(len+1) - 1; j++)
	        {
				k1[j] = test(y[i],b[j],n);
				y[i+1][j] = y[i][j] + k1[j]*h/6;
				k1[j] = y[i][j] + h*k1[j]/2;
				buf[j - n/p*len] = y[i+1][j];
				buf[j-n/p*(len - 1)] = k1[j];
			}
			for(int kj = 0; kj <= p-1; kj++)
			{
				if(kj != len)
				{
					MPI_Send(buf,n/p*2,MPI_DOUBLE,kj,kj+(len+1)*10 + (p+1)*i,MPI_COMM_WORLD);
					MPI_Recv(buf2[kj],n/p*2,MPI_DOUBLE,kj,(kj+1)*10 + len + (p+1)*i,MPI_COMM_WORLD,&status);
				}
			}
			for(int pi = 0; pi <= p - 1; pi++)
			{
				if(pi != len)
				{
					for(int j = n/p*(pi); j <= n/p*(pi+1) - 1; j++)
					{
						y[i+1][j] = buf2[pi][j-n/p*pi];
						k1[j] = buf2[pi][j-n/p*(pi-1)];
					}
				}
			}
			for(int j = n/p*len; j <= n/p*(len+1) - 1; j++)
	        {
				k2[j] = test(k1,b[j],n);
				y[i+1][j] = y[i+1][j] + 2*k2[j]*h/6;
				k2[j] = y[i][j] + h*k2[j]/2;
				buf[j - n/p*len] = y[i+1][j];
				buf[j-n/p*(len - 1)] = k2[j];
			}
			for(int kj = 0; kj <= p-1; kj++)
			{
				if(kj != len)
				{
					MPI_Send(buf,n/p*2,MPI_DOUBLE,kj,kj+(len+1)*100+(p+1)*i,MPI_COMM_WORLD);
					MPI_Recv(buf2[kj],n/p*2,MPI_DOUBLE,kj,(kj+1)*100 + len+(p+1)*i,MPI_COMM_WORLD,&status);
				}
			}
			for(int pi = 0; pi <= p - 1; pi++)
			{
				if(pi != len)
				{
					for(int j = n/p*(pi); j <= n/p*(pi+1) - 1; j++)
					{
						y[i+1][j] = buf2[pi][j-n/p*pi];
						k2[j] = buf2[pi][j-n/p*(pi-1)];
					}
				}
			}

			for(int j = n/p*len; j <= n/p*(len+1) - 1; j++)
	        {
				k3[j] = test(k2,b[j],n);
				y[i+1][j] = y[i+1][j] + 2*k3[j]*h/6;
				k3[j] = y[i][j] + h*k3[j];
				buf[j - n/p*len] = y[i+1][j];
				buf[j-n/p*(len - 1)] = k3[j];
			}
			for(int kj = 0; kj <= p-1; kj++)
			{
				if(kj != len)
				{
					MPI_Send(buf,n/p*2,MPI_DOUBLE,kj,kj+(len+1)*1000+(p+1)*i,MPI_COMM_WORLD);
					MPI_Recv(buf2[kj],n/p*2,MPI_DOUBLE,kj,(kj+1)*1000 + len+(p+1)*i,MPI_COMM_WORLD,&status);
				}
			}
			for(int pi = 0; pi <= p - 1; pi++)
			{
				if(pi != len)
				{
					for(int j = n/p*(pi); j <= n/p*(pi+1) - 1; j++)
					{
						y[i+1][j] = buf2[pi][j-n/p*pi];
						k3[j] = buf2[pi][j-n/p*(pi-1)];
					}
				}
			}

			for(int j = n/p*len; j <= n/p*(len+1) - 1; j++)
	        {
				k4[j] = test(k3,b[j],n);
			    y[i+1][j] = y[i+1][j] + k4[j]*h/6;
				buf[j - n/p*len] = y[i+1][j];
				buf[j-n/p*(len - 1)] = k4[j];
			}
			for(int kj = 0; kj <= p-1; kj++)
			{
				if(kj != len)
				{
					MPI_Send(buf,n/p*2,MPI_DOUBLE,kj,kj+(len+1)*10000+(p+1)*i,MPI_COMM_WORLD);
					MPI_Recv(buf2[kj],n/p*2,MPI_DOUBLE,kj,(kj+1)*10000 + len+(p+1)*i,MPI_COMM_WORLD,&status);
				}
			}
			for(int pi = 0; pi <= p - 1; pi++)
			{
				if(pi != len)
				{
					for(int j = n/p*(pi); j <= n/p*(pi+1) - 1; j++)
					{
						y[i+1][j] = buf2[pi][j-n/p*pi];
						k4[j] = buf2[pi][j-n/p*(pi-1)];
					}
				}
			}
	}

	char bufout[50];
    if(len == 0)
	{sprintf(bufout,"%g %g %g %g\n",y[0][0],y[0][1],y[0][2],y[0][3]);
	fout<<bufout;
	sprintf(bufout,"%g %g %g %g\n",y[nc-1][0],y[nc-1][1],y[nc -1][2],y[nc -1][3]);
	fout<<bufout;
	sprintf(bufout,"%d\n", len);
	fout<<bufout;}
	for(int i = 0; i <= nc - 1; i++)
	{
		delete(y[i]);
	}
	for(int i = 0; i <= n - 1 ; i++)
	{
		delete(b[i]);
	}
	for(int i = 0; i <= n - 1 ; i++)
	{
		if (i != len)
		{
			delete(buf2[i]);
		}
	}
	delete(b);
	delete(y);
	delete(k1);
	delete(k2);
	delete(k3);
	delete(k4);
	delete(buf);
	delete(buf2);
	MPI_Finalize();
    if(len == 0)
	{fout << "runtime = " << clock()/1000.0 << endl;
	fout.close();}
	return 0;
}

double test(double* a, double* b, int n)
{
	int i;
	double s = 0;
    double (*fun[3])(double) = {exp, sin, cos};
    int fi;
	for(i = 0; i <= n-1; i++)
	{
           fi = rand()%5;
           if (fi == 4)
	   {s += a[i] * b[i];}
           else
           {if (fi == 3)
	   {s += a[i] * b[i];}
           else
           {
             s += (*fun[fi])(-abs(a[i]))*b[i];
           }
           }
	}
	return s;
}
