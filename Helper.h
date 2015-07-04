#include "itkImage.h"


#include <vtkSmartPointer.h>
#include <vtkImageViewer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include "itkImageToVTKImageFilter.h"

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



class Helper {

public:

	template<class ImageType, unsigned int Dimension>
	int PaintImageWithHeatmap(itk::Image<ImageType, Dimension>::Pointer image, itk::Image<ImageType, Dimension>::Pointer map);
	#include "Helper.cpp"
};