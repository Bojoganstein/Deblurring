#pragma once
#include <complex>
#include <cmath>
#include <thread>
#include <vector>


using namespace std;

const double PI = 3.141592653589793238463;

class FourierTransform {
	static void forThread(unsigned char *img, int w, int h, int begin, int end, unsigned char* ampOut)	// lateral effect
	{
		double* realOut = new double[w*h];
		double* imgOut = new double[w*h];

		for (int i = 0;i < w*h;i++)
		{
			realOut[i] = 0;
			imgOut[i] = 0;
		}

		for (int yWave = begin; yWave < end; yWave++) {
			cout << "merge\n";
			for (int xWave = 0; xWave < w; xWave++) {

				for (int ySpace = 0; ySpace < h; ySpace++)
					for (int xSpace = 0; xSpace < w; xSpace++) {

						realOut[yWave*w + xWave] += (img[ySpace*w + xSpace] * cos(2 * PI * ((1.0 * xWave * xSpace / w) + (1.0 * yWave * ySpace / h)))) / sqrt(w * h);
						imgOut[yWave*w + xWave] -= (img[ySpace*w + xSpace] * sin(2 * PI * ((1.0 * xWave * xSpace / w) + (1.0 * yWave * ySpace / h)))) / sqrt(w * h);
					}

				double val = sqrt(realOut[yWave*w + xWave] * realOut[yWave*w + xWave] + imgOut[yWave*w + xWave] * imgOut[yWave*w + xWave]);
				if (val > 255.0)
					ampOut[yWave*w + xWave] = 255;
				else
					ampOut[yWave*w + xWave] = (unsigned char)(val);
			}
		}

		delete[] realOut;
		delete[] imgOut;
	}

public:
	static unsigned char* directDFT(unsigned char *img, int w, int h)
	{

		double* realOut = new double[w*h];
		double* imgOut = new double[w*h];
		unsigned char* ampOut = new unsigned char[w*h];


		int numOfThreads = 32;
		int step = h / numOfThreads;
		vector<thread> threads;
		int yWave;

		for (yWave = 0; yWave+step <= h; yWave+=step) {
			threads.push_back(thread(forThread, img, w, h, yWave, yWave + step, ampOut));
		}
		if(yWave!=h)
			threads.push_back(thread(forThread, img, w, h, yWave, h, ampOut));

		for (int i = 0;i < threads.size();i++)
			threads[i].join();

		return ampOut;
	}

	static unsigned char* reverseDFT(complex<double> * img, int w, int h) 
	{
		unsigned char* result = new unsigned char[w*h];
		complex<double> I(0, 1);

		for (int y = 0;y < h;y++)
			for (int x = 0;x < w;x++)
			{
				complex<double> sum = 0;

				for (int m = 0; m < h;m++)
					for (int n = 0;n < w;n++)
						sum += img[m*h + n] * exp(I*2.0*PI*((double)y / (double)(h)*m + (double)x / (double)(w)*n));

				result[y*w + x] = (unsigned char)(sum.real());
			}
		return result;
	}
};
