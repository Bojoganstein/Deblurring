#pragma once
#include <string>
#include <iostream>
#include <fstream>
using namespace std;
class Filter {
	int width;
	int heigth;
	double **kernel;

	unsigned char convolution(unsigned char* img, int w, int h, int center_x, int center_y)
	{
		double sum = 0;
		for (int y = -heigth / 2;y <= heigth / 2;y++)
			for (int x = -width / 2;x <= width / 2;x++)
				if (center_x + x >= 0 && center_y + y >= 0 && center_x + x < w && center_y + y < h)
				{
					sum += (double)(img[(center_x + x) + (center_y + y)*w]) * kernel[(width / 2) + x][(heigth / 2) + y];
				}
		return (unsigned char)(sum);
	}
	
public:
	Filter(string fileName)
	{
		fstream my_file;
		my_file.open(fileName, ios::in);
		if (!my_file) {
			cout << "No such file";
		}
		else 
		{
			my_file >> heigth >> width;
			kernel = new double*[heigth];
			for (int i = 0;i < heigth;i++)
				kernel[i] = new double[width];

			int index = 0;
			while (index<width*heigth)
			{
				my_file >> kernel[index/heigth][index%width];
				index++;
			}
		}
		my_file.close();
	}

	void printKernel()
	{
		for (int i = 0;i < heigth;i++)
		{
			for (int j = 0;j < width;j++)
				cout << kernel[i][j] << " ";
			cout << "\n";
		}
	}

	unsigned char* filterImage(unsigned char *img, int w, int h)
	{
		unsigned char *result = new unsigned char[w*h];

		for (int y = 0;y < h;y++)
			for (int x = 0;x < w;x++)
				result[w*y + x] = convolution(img, w, h,x,y);
		return result;
	}

	void multiply(double scalar)
	{
		for (int i = 0;i < heigth;i++)
			for (int j = 0;j < width;j++)
			{
				kernel[i][j] *= scalar;
			}
	}
	double* getKernel()
	{
		double* result = new double[width*heigth];

		for (int i = 0;i < heigth;i++)
			for (int j = 0;j < width;j++)
			{
				result[i*width + j] = kernel[i][j];
			}
		return result;
	}
	int getWidth() 
	{
		return width;
	}

	int getHeight()
	{
		return heigth;
	}

	~Filter()
	{
		for (int i = 0;i < heigth;i++)
			delete[] kernel[i];
		delete[] kernel;
	}
};