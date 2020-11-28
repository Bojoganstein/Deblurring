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
	
	Filter motionBlur("blur_filter.txt");
	motionBlur.multiply(1.0 / 15.0);
	unsigned char* blurred_in_space = motionBlur.filterImage(img, w, h);
	
	//Original Image FT
	Complex* waves = Fourier::DFT(img, w, h);
	//unsigned char* spectrum = Fourier::dftSpectrum(waves, w, h);
	//unsigned char* rftImage = Fourier::RFT(waves, w, h);

	//Blurred Image FT
	Complex* bis_FT = Fourier::DFT(blurred_in_space, w, h);
	//unsigned char* bis_spec = Fourier::dftSpectrum(bis_FT, w, h);
	
	//Blurred in Frequency Domain
	Complex* blurrKernel = Fourier::filter_DFT(motionBlur.getKernel(), motionBlur.getWidth(), motionBlur.getHeight(), w, h);
	Complex* blurredWave = multiply1D(waves, blurrKernel, w, h);
	unsigned char* blurred_in_freq= Fourier::RFT(blurredWave, w, h);
	unsigned char* bif_spec = Fourier::dftSpectrum(blurredWave, w, h);
	
	
	//Aplicarea filtrului Wiener peste o imagine blurata in domeniul frecventa
	Complex* wiener_filter = wienerFilterApprox(blurrKernel, 0.01, w, h);
	Complex* filtered_image1 = multiply1D(blurredWave, wiener_filter, w, h);
	unsigned char* end_result1 = Fourier::RFT(filtered_image1, w, h);
	
	//Aplicarea filtrului Wiener peste o imagine blurata in domeniul spatiu
	Complex* filtered_image2 = multiply1D(bis_FT, wiener_filter, w, h);
	unsigned char* end_result2 = Fourier::RFT(filtered_image2, w, h);

	//cout << sum(blurred_in_space, w, h) /(double)sum(blurred_in_freq, w, h)<<"\n";

	//Questions to ask for another time
	/*
	grid->addImage(blurred_in_space, w, h, 0, 1, "blurred in space");
	grid->addImage(bis_spec, w, h, 0, 2, "Blurred in space Spectrum");
	grid->addImage(blurred_in_freq, w, h, 1, 0, "Blurred in Freq. Domain");
	grid->addImage(dark_img, w, h, 1, 1, "Blurred in Freq. Domain and darkened");
	grid->addImage(bif_spec, w, h, 1, 2, "Blurred in frequency Spectrum");
	grid->addImage(end_result, w, h, 2, 0, "Result");
	*/
	grid->addImage(blurred_in_space, w, h, 1, 0, "blurred in space");
	grid->addImage(blurred_in_freq, w, h, 2, 0, "Blurred in Freq. Domain");
	grid->addImage(end_result1, w, h, 2, 1, "blurred in frequency result");
	grid->addImage(end_result2, w, h, 1, 1, "blurred in space result");
	grid->show();
	

	//Hard memory leak
	/*
	delete[] grid;
	delete[] img;
	delete[] blurred_in_space;
	delete[] blurred_in_freq;
	delete[] waves;
	delete[] bis_FT;
	delete[] blurrKernel;
	delete[] blurredWave;
	delete[] bif_spec;
	delete[] wiener_filter;
	delete[] filtered_image1;
	delete[] filtered_image2;
	delete[] end_result1;
	delete[] end_result2;
	*/
	return a.exec();
	
}

