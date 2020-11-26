#pragma once
#include "Complex.h"
#include <math.h>
#include <iostream>
#include <vector>
#include <thread>

using namespace std;
const double PI = 3.141592653589793238463;
class Fourier {
	static void dftThreadU(unsigned char* img,Complex* complexRow, int w, int h, int begin, int end)
	{
		for (int u = 0;u < w;u++)
			for (int v = begin;v < end;v++)
			{

				for (int x = 0;x < w;x++)
					for (int y = 0;y < h;y++) {
						complexRow[v*w + u].re += img[y*w + x] * cos(2 * PI*(((1.0*u*x) / w) + ((1.0*v*y) / h)));
						complexRow[v*w + u].im -= img[y*w + x] * sin(2 * PI*(((1.0*u*x) / w) + ((1.0*v*y) / h)));
					}
			}
	}
	//////////////////////////////
	static void dftThreadR(double* img, Complex* complexRow, int w, int h, int begin, int end)
	{
		for (int u = 0;u < w;u++)
			for (int v = begin;v < end;v++)
			{

				for (int x = 0;x < w;x++)
					for (int y = 0;y < h;y++) {
						complexRow[v*w + u].re += img[y*w + x] * cos(2 * PI*(((1.0*u*x) / w) + ((1.0*v*y) / h)));
						complexRow[v*w + u].im -= img[y*w + x] * sin(2 * PI*(((1.0*u*x) / w) + ((1.0*v*y) / h)));
					}
			}
	}
	//////////////////////////////

	static void rftThread(Complex* waves, double* values, int w, int h, int begin, int end)
	{

		for (int u = 0;u < w;u++)
			for (int v = begin;v < end;v++)
			{
				for (int x = 0;x < w;x++)
					for (int y = 0;y < h;y++)
					{
						values[v*w + u] += (waves[y*w + x].re*cos(2 * PI*((1.0*u*x) / w + (1.0*v*y) / h))
							- waves[y*w + x].im*sin(2 * PI*((1.0*u*x) / w + (1.0*v*y) / h))) / (w*h);
					}
				values[v*w + u] = abs(values[v*w + u]);
			}
	}
public:
	static Complex* DFT(unsigned char* img, int w, int h)
	{
		Complex *complexRow = new Complex[w*h];

		int numOfThreads = 8;
		int step = h / numOfThreads;
		vector<thread> threads;
		int yWave;

		for (yWave = 0; yWave + step <= h; yWave += step) {
			threads.push_back(thread(dftThreadU, img,complexRow, w, h, yWave, yWave + step));
		}
		if (yWave != h)
			threads.push_back(thread(dftThreadU, img, complexRow, w, h, yWave, h));

		for (int i = 0;i < threads.size();i++)
			threads[i].join();

		//
		Complex* result = new Complex[w*h];
		for (int y = 0;y < h / 2;y++)
			for (int x = 0;x < w / 2;x++)
				result[(y + h / 2)*w + (x + w / 2)] = (complexRow[y*w + x]);

		for (int y = 0;y < h / 2;y++)
			for (int x = w / 2;x < w;x++)
				result[(y + h / 2)*w + (x - w / 2)] = (complexRow[y*w + x]);

		for (int y = h / 2;y < h;y++)
			for (int x = 0;x < w / 2;x++)
				result[(y - h / 2)*w + (x + w / 2)] = (complexRow[y*w + x]);

		for (int y = h / 2;y < h;y++)
			for (int x = w / 2;x < w;x++)
				result[(y - h / 2)*w + (x - w / 2)] = (complexRow[y*w + x]);

		//
		delete[] complexRow;

		return result;
	}
	//////////////////////
	static Complex* DFT(double* img, int w, int h)
	{
		Complex *complexRow = new Complex[w*h];

		int numOfThreads = 8;
		int step = h / numOfThreads;
		vector<thread> threads;
		int yWave;

		for (yWave = 0; yWave + step <= h; yWave += step) {
			threads.push_back(thread(dftThreadR, img, complexRow, w, h, yWave, yWave + step));
		}
		if (yWave != h)
			threads.push_back(thread(dftThreadR, img, complexRow, w, h, yWave, h));

		for (int i = 0;i < threads.size();i++)
			threads[i].join();

		//
		Complex* result = new Complex[w*h];
		for (int y = 0;y < h / 2;y++)
			for (int x = 0;x < w / 2;x++)
				result[(y + h / 2)*w + (x + w / 2)] = (complexRow[y*w + x]);

		for (int y = 0;y < h / 2;y++)
			for (int x = w / 2;x < w;x++)
				result[(y + h / 2)*w + (x - w / 2)] = (complexRow[y*w + x]);

		for (int y = h / 2;y < h;y++)
			for (int x = 0;x < w / 2;x++)
				result[(y - h / 2)*w + (x + w / 2)] = (complexRow[y*w + x]);

		for (int y = h / 2;y < h;y++)
			for (int x = w / 2;x < w;x++)
				result[(y - h / 2)*w + (x - w / 2)] = (complexRow[y*w + x]);

		//
		delete[] complexRow;

		return result;
	}
	///////////////////

	static unsigned char* dftSpectrum(Complex* dftImg, int w, int h)
	{
		unsigned char* result = new unsigned char[w*h];
		double* modules = new double[w*h];
		


		double maxim = -1;
		for (int i = 0;i < w*h;i++)
		{
			modules[i] = log(1+dftImg[i].abs());
			if (modules[i] > maxim) maxim = modules[i];
		}

		for (int i = 0;i < w*h;i++)
		{
			result[i] = (unsigned char)((modules[i] / maxim) * 255);
		}
		delete[] modules;

		return result;
	}
	static unsigned char* RFT(Complex* waves, int w, int h)
	{
		unsigned char* result = new unsigned char[w*h];
		double* values = new double[w*h];
		for (int i = 0;i < w*h;i++)
		{
			values[i] = 0;
		}
		double maxim = -1;

		//
		int numOfThreads = 8;
		int step = h / numOfThreads;
		vector<thread> threads;
		int yWave;

		for (yWave = 0; yWave + step <= h; yWave += step) {
			threads.push_back(thread(rftThread, waves, values, w, h, yWave, yWave + step));
		}
		if (yWave != h)
			threads.push_back(thread(rftThread, waves, values, w, h, yWave, h));

		for (int i = 0;i < threads.size();i++)
			threads[i].join();
		//

		for (int i = 0;i < w*h;i++)
			if (maxim < values[i]) maxim = values[i];
		

		for (int i = 0;i < w*h;i++)
		{
			result[i] = (unsigned char)((values[i] / maxim) * 255);
		}
		delete[] values;

		return result;
	}
};

double* lowPassFilter(int w, int h, double d0)
{
	double* result = new double[w*h];

	for (int x = 0;x < w;x++)
		for (int y = 0;y<h;y++)
			if (sqrt((x-w/2.0)*(x-w/2.0)+(y-h/2.0)*(y-h/2.0)) > d0) result[y*w+x] = 0;
			else result[y*w+x] = 1;

	return result;
}
double* lowPassButterworth(int w, int h, double d0, int n)
{
	double* result = new double[w*h];
	double dist;
	for (int y = 0;y < h;y++)
		for (int x = 0;x < w;x++) {
			dist = sqrt((x - w / 2.0)*(x - w / 2.0) + (y - h / 2.0)*(y - h / 2.0));
			result[y*w + x] = 1.0 / (1+pow(dist/d0,2*n));
		}
	return result;
}

double* lowPassGausian(int w, int h, double d0)
{

	double* result = new double[w*h];
	double dist;
	for (int y = 0;y < h;y++)
		for (int x = 0;x < w;x++) {
			dist = sqrt((x - w / 2.0)*(x - w / 2.0) + (y - h / 2.0)*(y - h / 2.0));
			result[y*w + x] = exp((-1 * dist*dist) / (2 * d0*d0));
		}
	return result;
}

Complex* kernelFT(double* kernel, int kw, int kh, int newW, int newH)
{
	double* resizedKernel = new double[newW*newH];

	for (int x = 0;x < newW;x++)
		for (int y = 0;y < newH;y++)
			if (x < kw && y < kh) resizedKernel[y*newW + x] = kernel[y*newW + x];
			else resizedKernel[y*newW + x] = 0;
	
	Complex* result = Fourier::DFT(resizedKernel,newW,newH);

	return result;
}

Complex* multiply1D(Complex* ftImg, double* filter, int w, int h)
{
	Complex* result = new Complex[w*h];
	for (int i = 0;i < w*h;i++)
	{
		result[i] = ftImg[i]*filter[i];
	}

	return result;
}

Complex* multiply1D(Complex* ftImg, Complex* filter, int w, int h)
{
	Complex* result = new Complex[w*h];
	for (int i = 0;i < w*h;i++)
	{
		result[i] = ftImg[i] * filter[i];
	}

	return result;
}