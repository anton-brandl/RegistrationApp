#pragma once


#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"


class ImageShort2D
{
public:
	ImageShort2D();
	~ImageShort2D();

private:
	const int dimension = 2;
	typedef short PixelType;

};

