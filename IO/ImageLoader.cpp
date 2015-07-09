/*#include "ImageLoader.h"

template <class PixelType>
itk::Image<PixelType, 2>* ImageLoader::loadDicomImage(std::string src)
{
	const    unsigned int    ImageDimension = 2;
	typedef itk::Image< PixelType, ImageDimension >  ImageType;
	typedef itk::ImageFileReader< ImageType  > ImageReaderType;

	ImageReaderType::Pointer  imageReader = ImageReaderType::New();

	imageReader->SetFileName(src);

	return imageReader->GetOutput();
}

template <class PixelType>
itk::Image<PixelType, 3>* ImageLoader::loadDicomVolume(std::string src)
{

}*/