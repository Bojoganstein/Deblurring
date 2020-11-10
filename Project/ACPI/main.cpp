#include <QApplication>
#include "ImageGrid.h"
#include "operatii.h"
#include <iostream>
#include "FourierTransform.h"
using namespace std;

int main(int argc, char *argv[])
{
	
	QApplication a(argc, argv);

	ImageGrid *grid = new ImageGrid("Prelucrarea imaginilor");

	QString imageDir = "Images/";
	QString imageFile(imageDir + "lena512.bmp");

	grid->addImage(imageFile, 0, 0);

	int w, h;
	unsigned char* img = Tools::readImageGray8(imageFile, w, h);
	
	unsigned char * fftImg = FourierTransform::directDFT(img, w, h);
	cout << "I'm here\n";
	//unsigned char* rftImg = FourierTransform::reverseDFT(fftImg, w, h);

	cout << "now i'm here\n";

	grid->addImage(fftImg, w, h, 0, 1, "fft");
	grid->show();
	
	return a.exec();
	
}

