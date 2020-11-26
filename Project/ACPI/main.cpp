#include <QApplication>
#include "ImageGrid.h"
#include <iostream>
#include "Fourier.h"
#include "Filter.h"

using namespace std;

int main(int argc, char *argv[])
{
	
	QApplication a(argc, argv);

	ImageGrid *grid = new ImageGrid("Prelucrarea imaginilor");

	QString imageDir = "Images/";
	QString imageFile(imageDir + "football_64.jpg");

	grid->addImage(imageFile, 0, 0);

	int w, h;
	unsigned char* img = Tools::readImageGray8(imageFile, w, h);
	
	Filter motionBlur("motion_blur_kernel9x9.txt");
	motionBlur.multiply(1.0 / 9);
	unsigned char* blurredImg = motionBlur.filterImage(img, w, h);
	
	//Original Image FT
	Complex* waves = Fourier::DFT(img, w, h);
	unsigned char* spectrum = Fourier::dftSpectrum(waves, w, h);
	unsigned char* rftImage = Fourier::RFT(waves, w, h);

	//Blurred Image FT
	Complex* blurredFT = Fourier::DFT(blurredImg, w, h);
	unsigned char* blurrySpec = Fourier::dftSpectrum(blurredFT, w, h);
	
	//Blurred in Frequency Domain
	/*
	Complex* blurrKernel = kernelFT(motionBlur.getKernel(), motionBlur.getWidth(), motionBlur.getHeight(), w, h);
	Complex* blurredWave = multiply1D(waves, blurrKernel, w, h);
	unsigned char* blurryImg = Fourier::RFT(blurredWave, w, h);
	*/

	grid->addImage(spectrum, w, h, 0, 1, "Spectrum");
	grid->addImage(rftImage, w, h, 0, 2, "Reversed");
	grid->addImage(blurredImg, w, h, 1, 0, "Blurred Image");
	grid->addImage(blurrySpec, w, h, 1, 1, "Blurred Image Spectrum");
	//grid->addImage(blurryImg, w, h, 1, 2, "Blurred in Freq. Domain");
	
	grid->show();
	
	return a.exec();
	
}

