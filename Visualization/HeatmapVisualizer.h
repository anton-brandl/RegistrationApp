#pragma once


#include "Visualizer.h"
#include <assert.h>

#include "itkImage.h"
#include "itkRGBPixel.h"

#include "itkRescaleIntensityImageFilter.h"
#include "itkImageToVTKImageFilter.h"
#include "itkVectorMagnitudeImageFilter.h"
#include "itkScalarToRGBColormapImageFilter.h"

#include <vtkSmartPointer.h>
#include <vtkImageViewer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>

#include <vtkImageBlend.h>

#include <vtkImageFlip.h>
#include <vtkArrowSource.h>
#include <vtkGlyph2D.h>
#include <vtkGlyph3D.h>
#include <vtkPolyDataMapper.h>
#include <vtkImageSlice.h>
#include <vtkImageSliceMapper.h>
#include <vtkRenderer.h>
#include <vtkActor.h>
#include <vtkInteractorStyleImage.h>

#include "itkThresholdImageFilter.h"

template <class ImageType, class FieldType, unsigned int Dimension>
class HeatmapVisualizer : public Visualizer<ImageType, FieldType, Dimension> 
{
	public:
		HeatmapVisualizer(){};
		~HeatmapVisualizer(){};
		void visualize() {
			assert(_image != nullptr && _field != nullptr);

			typedef itk::VectorMagnitudeImageFilter<
				FieldType, ImageType >  VectorMagnitudeFilterType;
			VectorMagnitudeFilterType::Pointer magnitudeFilter =
				VectorMagnitudeFilterType::New();
			magnitudeFilter->SetInput(_field);

			// To write the magnitude image file, we should rescale the gradient values
			// to a reasonable range

			typedef itk::RescaleIntensityImageFilter<
				ImageType, ImageType > rescaleFilterType;

			rescaleFilterType::Pointer rescaler = rescaleFilterType::New();
			rescaler->SetOutputMinimum(0);
			rescaler->SetOutputMaximum(255);
			rescaler->SetInput(magnitudeFilter->GetOutput());
			rescaler->Update();
			
			typedef itk::ThresholdImageFilter<ImageType> ThresholdFilterType;
			ThresholdFilterType::Pointer thresholdFilter = ThresholdFilterType::New();
			thresholdFilter->SetInput(rescaler->GetOutput());
			thresholdFilter->ThresholdOutside(100, 255);
			thresholdFilter->SetOutsideValue(0);
			
			typedef itk::RGBPixel<unsigned char>    RGBPixelType;
			typedef itk::Image<RGBPixelType, Dimension>  RGBImageType;
  
			typedef itk::ScalarToRGBColormapImageFilter<ImageType, RGBImageType> RGBFilterType;
			RGBFilterType::Pointer rgbfilter = RGBFilterType::New();
			rgbfilter->SetInput(thresholdFilter->GetOutput());
			rgbfilter->SetColormap(RGBFilterType::Hot);

			//viewer.AddRGBImage(rgbfilter->GetOutput(), true, "Displacement field magnitudes heatmap");

			typedef itk::ImageToVTKImageFilter<RGBImageType> RGBConnectorType;
			RGBConnectorType::Pointer rgbconnector = RGBConnectorType::New();
			rgbconnector->SetInput(rgbfilter->GetOutput());
			rgbconnector->Update();
			/*
			typedef itk::ImageToVTKImageFilter<DisplacementFieldImageType> VectorfieldConnectorType;
			VectorfieldConnectorType::Pointer vFieldConnector = VectorfieldConnectorType::New();
			vFieldConnector->SetInput(dispfieldGenerator->GetOutput());
			vFieldConnector->Update();
			*/

			/*typedef itk::ImageToVTKImageFilter<FixedImageType> BWConnectorType;
			BWConnectorType::Pointer bwConnector = BWConnectorType::New();
			bwConnector->SetInput(movingImageReader->GetOutput());
			bwConnector->Update();
			*/

			vtkSmartPointer<vtkImageData> heatmap = vtkSmartPointer<vtkImageData>::New();
			heatmap->DeepCopy(rgbconnector->GetOutput());


			rgbfilter->SetInput(_image);
			rgbfilter->SetColormap(RGBFilterType::Grey);
			rgbfilter->Update();
			rgbconnector->SetInput(rgbfilter->GetOutput());
			rgbconnector->Update();

			vtkSmartPointer<vtkImageData> anatomy = rgbconnector->GetOutput();


			vtkSmartPointer<vtkImageBlend> blend =
				vtkSmartPointer<vtkImageBlend>::New();
			blend->AddInputData(anatomy);
			blend->AddInputData(heatmap);
			blend->SetOpacity(0, 1);
			blend->SetOpacity(1, .5);

			vtkSmartPointer<vtkImageFlip> flipYFilter =
				vtkSmartPointer<vtkImageFlip>::New();
			flipYFilter->SetFilteredAxis(1); // flip y axis
			flipYFilter->SetInputConnection(blend->GetOutputPort());


			vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
				vtkSmartPointer<vtkRenderWindowInteractor>::New();

			vtkSmartPointer<vtkImageViewer> imageViewer =
				vtkSmartPointer<vtkImageViewer>::New();
			imageViewer->SetInputConnection(flipYFilter->GetOutputPort());
			imageViewer->SetupInteractor(renderWindowInteractor);
			imageViewer->Render();
			imageViewer->SetColorWindow(255);
			imageViewer->SetColorLevel(128);
			renderWindowInteractor->Start();

		};
};

