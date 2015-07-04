// Software Guide : BeginCodeSnippet
#include "itkImageRegistrationMethodv4.h"
#include "itkTranslationTransform.h"
#include "itkMeanSquaresImageToImageMetricv4.h"
#include "itkRegularStepGradientDescentOptimizerv4.h"
// Software Guide : EndCodeSnippet


#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkResampleImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkSubtractImageFilter.h"

#include <itkDisplacementFieldTransformParametersAdaptor.h>



//DeformableRegistrationTest:


#include "itkCheckerBoardImageFilter.h"
#include "QuickView.h"

//#include "itkImageRegistrationMethodv4.h"
//#include "itkMeanSquaresImageToImageMetricv4.h"


#include "itkBSplineTransform.h"
#include "itkLBFGSOptimizerv4.h"
#include "itkBSplineTransformParametersAdaptor.h"


//#include "itkImageFileReader.h"
//#include "itkImageFileWriter.h"

//#include "itkResampleImageFilter.h"
//#include "itkCastImageFilter.h"
#include "itkSquaredDifferenceImageFilter.h"

#include "itkIdentityTransform.h"

#include "itkBSplineTransformInitializer.h"
#include "itkTransformToDisplacementFieldFilter.h"


#include "itkVectorImage.h"

#include "itkThresholdImageFilter.h"



#include "itkLabelOverlayImageFilter.h"

#include "itkCastImageFilter.h"

#include "Visualization/HeatmapVisualizer.h"


int registerTest2D(std::string image1, std::string image2, std::string output, std::string differenceAfter = "", std::string differenceBefore = ""); 
int deformableRegistrationTest(std::string image1, std::string image2, std::string outputImage, std::string outputDisplacement = "", std::string outputDeformationField = "");