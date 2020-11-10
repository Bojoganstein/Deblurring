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

	void sort(unsigned char *pixels, int size)//WARNING SIDE EFFECTS (SORTED VECTOR)
	{
		unsigned char aux;
		for (int i = 0;i < size - 1;i++)
			for (int j = i;j < size;j++)
			{
				if (pixels[i] > pixels[j])
				{
					aux = pixels[i];
					pixels[i] = pixels[j];
					pixels[j] = aux;
				}
			}
	}

	unsigned char medFilter(unsigned char* img, int w, int h, int center_x, int center_y)
	{
		unsigned char * pixels = new unsigned char[width*heigth];
		int size = 0;
		unsigned char median = 0;

		for (int y = -heigth / 2;y <= heigth / 2;y++)
			for (int x = -width / 2;x <= width / 2;x++)
				if (center_x + x >= 0 && center_y + y >= 0 && center_x + x < w && center_y + y < h)
				{
					pixels[size] = (img[(center_x + x) + (center_y + y)*w]);
					size++;
				}

		//sorting
		sort(pixels, size);
		median = pixels[size / 2];
		delete[] pixels;

		return median;
	}

	void createHistogram(
		unsigned char* img,
		int w, int h,
		int center_x, int center_y,
		int *histogram)
	{
		for (int i = 0;i < 256;i++)
			histogram[i] = 0;

		for (int y = -heigth / 2;y <= heigth / 2;y++)
			for (int x = -width / 2;x <= width / 2;x++)
				if (center_x + x >= 0 && center_y + y >= 0 && center_x + x < w && center_y + y < h)
				{
					histogram[img[(center_x + x) + (center_y + y)*w]] ++;
				}
	}

	void updateToRight(				//apelez functia cand ma misc pe axa OX in sens pozitiv
		unsigned char* img,			//imi asum ca functia x>0
		int w, int h,
		int center_x, int center_y,
		int *histogram)
	{
		if(center_x- width/2 -1>=0)													//daca am careul intreg la stanga
			for (int y = -heigth / 2;y <= heigth / 2;y++)	//ma deplasez pe fosta coloana din stanga
				if(center_y+y>=0 && center_y+y<h) 
					histogram[img[center_x - width / 2 -1 + (center_y + y)*w]]--;	//decrementez pixelii
		
		if (center_x + width / 2 <= w)												//daca am careul intreg la dreapta
			for (int y = -heigth / 2;y <= heigth / 2;y++)//ma deplasez noua coloana din dreapta
				if (center_y + y >= 0 && center_y + y<h)
					histogram[img[center_x + width / 2 + (center_y + y)*w]]++;		//incrementez pixelii
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

	unsigned char* medianFilter(unsigned char *img, int w, int h)
	{
		unsigned char *result = new unsigned char[w*h];

		for (int y = 0;y < h;y++)
			for (int x = 0;x < w;x++)
				result[w*y + x] = medFilter(img, w, h, x, y);
		return result;
	}

	unsigned char* medianFilterHuang(unsigned char *img, int w, int h)
	{
		unsigned char *result = new unsigned char[w*h];
		int *frequency=new int[256];
	
		for (int y = 0;y < h;y++)
		{
			//creez histograma
			createHistogram(img, w, h, 0, y, frequency);

			//shiftez coloanele
			for (int x = 0;x < w;x++) 
			{
				if (x > 0) updateToRight(img, w, h, x, y, frequency);
				
				//fac suma valorilor din histograma si ma opresc la width*heigth/2
				int sum = 0;
				int i =256;
				for (i = 0;i < 256;i++)
				{
					sum += frequency[i];
					if (sum > (width*heigth) / 2)	//ma opresc cand depasesc width*heigth/2
						break;
				}
				result[x + y*w] = (unsigned char)i;
			}
		}
		delete[] frequency;

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

	~Filter()
	{
		for (int i = 0;i < heigth;i++)
			delete[] kernel[i];
		delete[] kernel;
	}
};