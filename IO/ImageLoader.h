#include <memory>
#include <string>
#include "itkImage.h"
#include "itkImageSeriesReader.h"
#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkSmartPointer.h"
#include <stdio.h>
#include <vector>

//typedef itk::Image<float, 2> Image2dType;
//typedef itk::Image<float, 3> Volume3dType;
class ImageLoader
{
public:
	
	template <typename PixelType>
	typename itk::SmartPointer<itk::Image<PixelType, 2>> loadDicomImage(std::string src)
	{
		const    unsigned int    ImageDimension = 2;
		typedef itk::Image< float, ImageDimension > ImageType;
		typedef itk::ImageFileReader< ImageType  > ImageReaderType;

		ImageReaderType::Pointer  imageReader = ImageReaderType::New();

		imageReader->SetFileName(src);
		imageReader->Update();
		return imageReader->GetOutput();
	}
	

	template <typename PixelType>
	typename itk::SmartPointer<itk::Image<PixelType, 3>> loadDicomVolume(std::string src) 
	{
		const unsigned int      Dimension = 3;
		typedef itk::Image< PixelType, Dimension > ImageType;
		typedef itk::ImageSeriesReader< ImageType > ReaderType;
		
		ReaderType::Pointer reader = ReaderType::New();
		
		typedef itk::GDCMImageIO       ImageIOType;
		ImageIOType::Pointer dicomIO = ImageIOType::New();

		reader->SetImageIO(dicomIO);
		
		typedef itk::GDCMSeriesFileNames NamesGeneratorType;
		NamesGeneratorType::Pointer nameGenerator = NamesGeneratorType::New();

		nameGenerator->SetUseSeriesDetails(true);
		//nameGenerator->AddSeriesRestriction("0008|0021");

		nameGenerator->SetDirectory(src);
		// Software Guide : EndCodeSnippet


		try
		{
			std::cout << std::endl << "The directory: " << std::endl;
			std::cout << std::endl << src << std::endl << std::endl;
			std::cout << "Contains the following DICOM Series: ";
			std::cout << std::endl << std::endl;

			// Software Guide : BeginCodeSnippet
			typedef std::vector< std::string >    SeriesIdContainer;

			const SeriesIdContainer & seriesUID = nameGenerator->GetSeriesUIDs();

			SeriesIdContainer::const_iterator seriesItr = seriesUID.begin();
			SeriesIdContainer::const_iterator seriesEnd = seriesUID.end();
		/*	while (seriesItr != seriesEnd)
			{
				std::cout << seriesItr->c_str() << std::endl;
				++seriesItr;
			}*/
		
			std::string seriesIdentifier;

			seriesIdentifier = seriesUID.begin()->c_str();
			
			// Software Guide : EndCodeSnippet


			std::cout << std::endl << std::endl;
			std::cout << "Now reading series: " << std::endl << std::endl;
			std::cout << seriesIdentifier << std::endl;
			std::cout << std::endl << std::endl;

		
			typedef std::vector< std::string >   FileNamesContainer;
			FileNamesContainer fileNames;

			fileNames = nameGenerator->GetFileNames(seriesIdentifier);
			
			reader->SetFileNames(fileNames);
			
			try
			{
				reader->Update();
			}
			catch (itk::ExceptionObject &ex)
			{
				std::cout << ex << std::endl;
				return NULL;
			}
			itk::Image<PixelType, 3>::Pointer result = reader->GetOutput();
			return result;
		}
		catch (itk::ExceptionObject &ex)
		{
			std::cout << ex << std::endl;
			return NULL;
		}
	}
		
};