#pragma once
#include <complex>
#include <cmath>
#include <thread>
#include <vector>
#include <iostream>

using namespace std;

const double PI = 3.141592653589793238463;

class FourierTransform {
	unsigned char* img;
	double* realOut;
	double* imgOut;
	int w;
	int h;

	static void dftThread(unsigned char *img, int w, int h, int begin, int end, double* realOut, double* imgOut)	// lateral effect
	{

		for (int yWave = begin; yWave < end; yWave++) 
		{
			cout << "DFT merge\n";
			for (int xWave = 0; xWave < w; xWave++) 
			{


				//Complexity O(n^2)
				for (int ySpace = 0; ySpace < h; ySpace++)
					for (int xSpace = 0; xSpace < w; xSpace++) {

						realOut[yWave*w + xWave] += (img[ySpace*w + xSpace] * cos(2 * PI * (((1.0 * xWave * xSpace) / w) + ((1.0 * yWave * ySpace )/ h))));// / sqrt(w * h);
						imgOut[yWave*w + xWave] -= (img[ySpace*w + xSpace] * sin(2 * PI * (((1.0 * xWave * xSpace) / w) + ((1.0 * yWave * ySpace) / h))));// / sqrt(w * h);
					}

			}
		}
	}

	static void rftThread(double* realIn, double* imgIn, int w, int h, int begin, int end,double* realOut, double *imgOut)	// lateral effect
	{

		for (int yWave = begin; yWave < end; yWave++)
		{
			cout << "RFT merge\n";
			for (int xWave = 0; xWave < w; xWave++)
			{


				//Complexity O(n^2)
				for (int ySpace = 0; ySpace < h; ySpace++)
					for (int xSpace = 0; xSpace < w; xSpace++) 
					{

						double lambda = 2 * PI * ((1.0 * xWave * xSpace / w) + (1.0 * yWave * ySpace / h));
						
						realOut[yWave*w + xWave] += (realIn[ySpace*w + xSpace] * cos(lambda) - imgIn[ySpace*w + xSpace] * sin(lambda)) / (w * h);// sqrt(w * h);
						imgOut[yWave*w + xWave] += (realIn[ySpace*w + xSpace] * sin(lambda) + imgIn[ySpace*w + xSpace] * cos(lambda)) / (w * h);// sqrt(w * h);
					}
			}
		}
	}

public:
	FourierTransform(unsigned char* img,int w, int h)
	{
		this->w = w;
		this->h = h;
		this->img = new unsigned char[w*h];
		this->realOut = new double[w*h];
		this->imgOut = new double[w*h];

		for (int i = 0;i < w*h;i++)
		{
			this->img[i] = img[i];
			realOut[i] = 0;
			imgOut[i] = 0;
		}
	}
	~FourierTransform()
	{
		delete[] img;
		delete[] imgOut;
		delete[] realOut;
	}

	void makeDFT()//CALL THIS FUNCTION FIRST
	{
	
		int numOfThreads = 8;
		int step = h / numOfThreads;
		vector<thread> threads;
		int yWave;

		for (yWave = 0; yWave+step <= h; yWave+=step) {
			threads.push_back(thread(dftThread, img, w, h, yWave, yWave + step, realOut, imgOut));
		}
		if(yWave!=h)
			threads.push_back(thread(dftThread, img, w, h, yWave, h, realOut, imgOut));

		for (int i = 0;i < threads.size();i++)
			threads[i].join();

		

	}

	unsigned char* getDFT()//CALL THIS FUNCTION AFTER makeDFT()
	{
		double* ampOut = new double[w*h];


		//Normalize
		double maxVal = 0;
		for (int i = 0;i < w*h;i++)
		{
			ampOut[i] = log(1 + sqrt(realOut[i] * realOut[i] + imgOut[i] * imgOut[i]));
			if (maxVal < ampOut[i])
				maxVal = ampOut[i];
		}

		for (int i = 0;i < h*w;i++) {
			ampOut[i] = (ampOut[i] / maxVal);
		}
		for (int i = 0;i < h*w;i++) {
			ampOut[i] = ampOut[i] * 255.0;
		}

		//Save the result
		unsigned char* result = new unsigned char[w*h];
		for (int y = 0;y < h / 2;y++)
			for (int x = 0;x < w / 2;x++)
				result[(y + h / 2)*w + (x + w / 2)] = (unsigned char)(ampOut[y*w + x]);

		for (int y = 0;y < h / 2;y++)
			for (int x = w / 2;x < w;x++)
				result[(y + h / 2)*w + (x - w / 2)] = (unsigned char)(ampOut[y*w + x]);

		for (int y = h / 2;y < h;y++)
			for (int x = 0;x < w / 2;x++)
				result[(y - h / 2)*w + (x + w / 2)] = (unsigned char)(ampOut[y*w + x]);

		for (int y = h / 2;y < h;y++)
			for (int x = w / 2;x < w;x++)
				result[(y - h / 2)*w + (x - w / 2)] = (unsigned char)(ampOut[y*w + x]);

		delete[] ampOut;

		return result;
	}

	unsigned char* reverseDFT() 
	{
		double* result = new double[w*h];
		double* resultReal = new double[w*h];
		double* resultImg = new double[w*h];

		unsigned char* image = new unsigned char[w*h];

		for (int i = 0;i < w*h;i++)
		{
			resultReal[i] = 0;
			resultImg[i] = 0;
		}

		int numOfThreads = 8;
		int step = h / numOfThreads;
		vector<thread> threads;
		int yWave;

		for (yWave = 0; yWave + step <= h; yWave += step) {
			threads.push_back(thread(rftThread, realOut, imgOut, w, h, yWave, yWave+step, resultReal, resultImg));
		}
		if (yWave != h)
			threads.push_back(thread(rftThread, realOut, imgOut, w, h, yWave, h, resultReal, resultImg));

		for (int i = 0;i < threads.size();i++)
			threads[i].join();




		double maxVal = 0;
		for (int i = 0;i < w*h;i++)
		{
			result[i] = sqrt(resultReal[i] * resultReal[i]);//Simplifica partea img
			if (maxVal < result[i])
				maxVal = result[i];
		}

		for (int i = 0;i < h*w;i++) {
			result[i] = (result[i] / maxVal);
		}
		for (int i = 0;i < h*w;i++) {
			image[i] = (unsigned char)(result[i] * 255.0);
		}
		delete[] result;
		delete[] resultReal;
		delete[] resultImg;


		return image;
	}
};
