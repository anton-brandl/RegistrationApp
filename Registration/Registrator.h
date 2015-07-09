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
//#include "QuickView.h"

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
#include "IO/ImageLoader.h"
#include <itkSmartPointer.h>

class Registrator {

public: 
int translationRegistration2DTest(std::string image1, std::string image2, std::string output, std::string differenceAfter = "", std::string differenceBefore = "");

template<class ImageType, class FieldType, unsigned int Dimension>
typename itk::SmartPointer<FieldType> deformableRegistrationTest(itk::SmartPointer<ImageType> fixed, itk::SmartPointer<ImageType> moving){
	/*
	const    unsigned int    ImageDimension = 2;
	typedef  float   PixelType;

	ImageLoader<PixelType> imageLoader;
	typedef itk::Image< PixelType, ImageDimension >  FixedImageType;
	typedef itk::Image< PixelType, ImageDimension >  MovingImageType;

	FixedImageType::Pointer fixedImage = imageLoader.loadDicomImage(image1);
	MovingImageType::Pointer movingImage = imageLoader.loadDicomImage(image2);

	//  Software Guide : BeginLatex
	//
	//  We instantiate the type of the \code{BSplineTransform} using
	//  as template parameters the type for coordinates representation, the
	//  dimension of the space, and the order of the BSpline.
	//
	//  \index{BSplineTransform!New}
	//  \index{BSplineTransform!Instantiation}
	//
	//  Software Guide : EndLatex

	// Software Guide : BeginCodeSnippet

	const unsigned int SpaceDimension = ImageDimension;
	*/

	const unsigned int SplineOrder = 3;
	typedef double CoordinateRepType;

	typedef itk::BSplineTransform<
		CoordinateRepType,
		Dimension,
		SplineOrder >     TransformType;
	// Software Guide : EndCodeSnippet

	typedef itk::LBFGSOptimizerv4       OptimizerType;


	typedef itk::MeanSquaresImageToImageMetricv4<
		ImageType,
		ImageType >    MetricType;

	typedef itk::ImageRegistrationMethodv4<
		ImageType,
		ImageType >    RegistrationType;

	MetricType::Pointer         metric = MetricType::New();
	OptimizerType::Pointer      optimizer = OptimizerType::New();
	RegistrationType::Pointer   registration = RegistrationType::New();


	registration->SetMetric(metric);
	registration->SetOptimizer(optimizer);

	//typedef itk::ImageFileReader< ImageType  > FixedImageReaderType;
	//typedef itk::ImageFileReader< ImageType > MovingImageReaderType;

	//	FixedImageReaderType::Pointer  fixedImageReader = FixedImageReaderType::New();
	//	MovingImageReaderType::Pointer movingImageReader = MovingImageReaderType::New();

	//	fixedImageReader->SetFileName(image1);
	//	movingImageReader->SetFileName(image2);

	//	FixedImageType::ConstPointer fixedImage = fixedImageReader->GetOutput();
	//	FixedImageType::ConstPointer movingImage = movingImageReader->GetOutput();

	registration->SetFixedImage(fixed);
	registration->SetMovingImage(moving);


	//  Software Guide : BeginLatex
	//
	//  We construct the transform object, initialize its parameters and
	//  connect that to the registration object.
	//
	//  \index{itk::RegistrationMethod!SetTransform()}
	//
	//  Software Guide : EndLatex

	// Software Guide : BeginCodeSnippet
	TransformType::Pointer  outputBSplineTransform = TransformType::New();

	// Initialize the fixed parameters of transform (grid size, etc).
	//
	typedef itk::BSplineTransformInitializer<
		TransformType,
		ImageType> InitializerType;

	InitializerType::Pointer transformInitializer = InitializerType::New();

	unsigned int numberOfGridNodesInOneDimension = 25;

	TransformType::MeshSizeType             meshSize;
	meshSize.Fill(numberOfGridNodesInOneDimension - SplineOrder);

	transformInitializer->SetTransform(outputBSplineTransform);
	transformInitializer->SetImage(fixed);
	transformInitializer->SetTransformDomainMeshSize(meshSize);
	transformInitializer->InitializeTransform();

	// Set transform to identity
	//
	typedef TransformType::ParametersType     ParametersType;
	const unsigned int numberOfParameters =
		outputBSplineTransform->GetNumberOfParameters();
	ParametersType parameters(numberOfParameters);
	parameters.Fill(0.0);
	outputBSplineTransform->SetParameters(parameters);

	registration->SetInitialTransform(outputBSplineTransform);
	registration->InPlaceOn();
	//  Software Guide : EndCodeSnippet

	//  Software Guide : BeginLatex
	//
	//  The registration process is run in three levels. The shrink factors
	//  and smoothing sigmas are set for each level.
	//
	//  Software Guide : EndLatex

	// Software Guide : BeginCodeSnippet
	const unsigned int numberOfLevels = 1;

	RegistrationType::ShrinkFactorsArrayType shrinkFactorsPerLevel;
	shrinkFactorsPerLevel.SetSize(numberOfLevels);
	shrinkFactorsPerLevel[0] = 1;
	//	shrinkFactorsPerLevel[1] = 2;
	//	shrinkFactorsPerLevel[2] = 1;

	RegistrationType::SmoothingSigmasArrayType smoothingSigmasPerLevel;
	smoothingSigmasPerLevel.SetSize(numberOfLevels);
	smoothingSigmasPerLevel[0] = 0;
	//	smoothingSigmasPerLevel[1] = 1;
	//	smoothingSigmasPerLevel[2] = 0;

	registration->SetNumberOfLevels(numberOfLevels);
	registration->SetSmoothingSigmasPerLevel(smoothingSigmasPerLevel);
	registration->SetShrinkFactorsPerLevel(shrinkFactorsPerLevel);
	//  Software Guide : EndCodeSnippet

	//  Software Guide : BeginLatex
	//
	//  Create the transform adaptors to modify the flexibility
	//  of the deformable transform for each level of this
	//  multi-resolution scheme.
	//
	//  Software Guide : EndLatex

	// Software Guide : BeginCodeSnippet
	RegistrationType::TransformParametersAdaptorsContainerType adaptors;

	// First, get fixed image physical dimensions
	TransformType::PhysicalDimensionsType             fixedPhysicalDimensions;
	for (unsigned int i = 0; i< Dimension; i++)
	{
		fixedPhysicalDimensions[i] = fixed->GetSpacing()[i] *
			static_cast<double>(
			fixed->GetLargestPossibleRegion().GetSize()[i] - 1);
	}

	// Create the transform adaptors specific to B-splines
	for (unsigned int level = 0; level < numberOfLevels; level++)
	{
		typedef itk::ShrinkImageFilter<
			ImageType,
			ImageType> ShrinkFilterType;
		ShrinkFilterType::Pointer shrinkFilter = ShrinkFilterType::New();
		shrinkFilter->SetShrinkFactors(shrinkFactorsPerLevel[level]);
		shrinkFilter->SetInput(fixed);
		shrinkFilter->Update();

		// A good heuristic is to double the b-spline mesh resolution at each level
		//
		TransformType::MeshSizeType requiredMeshSize;
		for (unsigned int d = 0; d < Dimension; d++)
		{
			requiredMeshSize[d] = meshSize[d] << level;
		}

		typedef itk::BSplineTransformParametersAdaptor<TransformType>
			BSplineAdaptorType;
		BSplineAdaptorType::Pointer bsplineAdaptor = BSplineAdaptorType::New();
		bsplineAdaptor->SetTransform(outputBSplineTransform);
		bsplineAdaptor->SetRequiredTransformDomainMeshSize(requiredMeshSize);
		bsplineAdaptor->SetRequiredTransformDomainOrigin(
			shrinkFilter->GetOutput()->GetOrigin());
		bsplineAdaptor->SetRequiredTransformDomainDirection(
			shrinkFilter->GetOutput()->GetDirection());
		bsplineAdaptor->SetRequiredTransformDomainPhysicalDimensions(
			fixedPhysicalDimensions);

		adaptors.push_back(bsplineAdaptor.GetPointer());
	}

	registration->SetTransformParametersAdaptorsPerLevel(adaptors);
	//  Software Guide : EndCodeSnippet

	// Scale estimator
	typedef itk::RegistrationParameterScalesFromPhysicalShift<MetricType> ScalesEstimatorType;
	ScalesEstimatorType::Pointer scalesEstimator = ScalesEstimatorType::New();
	scalesEstimator->SetMetric(metric);
	scalesEstimator->SetTransformForward(true);
	scalesEstimator->SetSmallParameterVariation(1.0);

	// Set Optimizer
	optimizer->SetScalesEstimator(scalesEstimator);
	optimizer->SetGradientConvergenceTolerance(0.05);
	optimizer->SetLineSearchAccuracy(0.9);
	optimizer->SetDefaultStepLength(1.5);
	optimizer->TraceOn();
	optimizer->SetMaximumNumberOfFunctionEvaluations(20);

	std::cout << "Starting Registration "
		<< std::endl;

	try
	{
		registration->Update();
		std::cout << "Optimizer stop condition = "
			<< registration->GetOptimizer()->GetStopConditionDescription()
			<< std::endl;
	}
	catch (itk::ExceptionObject & err)
	{
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
		return NULL;
	}

	// Finally we use the last transform in order to resample the image.
	//
	typedef itk::ResampleImageFilter<ImageType, ImageType >    ResampleFilterType;

	ResampleFilterType::Pointer resample = ResampleFilterType::New();

	resample->SetTransform(outputBSplineTransform);
	resample->SetInput(moving);

	resample->SetSize(fixed->GetLargestPossibleRegion().GetSize());
	resample->SetOutputOrigin(fixed->GetOrigin());
	resample->SetOutputSpacing(fixed->GetSpacing());
	resample->SetOutputDirection(fixed->GetDirection());
	resample->SetDefaultPixelValue(100);

	typedef  unsigned char  OutputPixelType;

	typedef itk::Image< OutputPixelType, Dimension > OutputImageType;

	typedef itk::CastImageFilter<
		ImageType,
		OutputImageType > CastFilterType;

	typedef itk::ImageFileWriter< OutputImageType >  WriterType;


	WriterType::Pointer      writer = WriterType::New();
	CastFilterType::Pointer  caster = CastFilterType::New();


	writer->SetFileName("output.mha");


	caster->SetInput(resample->GetOutput());
	writer->SetInput(caster->GetOutput());
	//QuickView viewer;
	//viewer.AddImage(resample->GetOutput(), true, "Resampled Image");


	try
	{
		writer->Update();
	}
	catch (itk::ExceptionObject & err)
	{
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
		return NULL;
	}

	typedef itk::SquaredDifferenceImageFilter<
		ImageType,
		ImageType,
		OutputImageType > DifferenceFilterType;

	DifferenceFilterType::Pointer difference = DifferenceFilterType::New();

	WriterType::Pointer writer2 = WriterType::New();
	writer2->SetInput(difference->GetOutput());


	// Compute the difference image between the
	// fixed and resampled moving image.

	difference->SetInput1(fixed);
	difference->SetInput2(resample->GetOutput());
	writer2->SetFileName("differenceFixedResMov.mha");
	try
	{
		writer2->Update();
	}
	catch (itk::ExceptionObject & err)
	{
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
		return NULL;
	}

	typedef itk::CheckerBoardImageFilter< ImageType > CheckerBoardFilterType;
	typedef CheckerBoardFilterType::PatternArrayType PatternArrayType;
	CheckerBoardFilterType::Pointer checkerBoardFilter = CheckerBoardFilterType::New();
	checkerBoardFilter->SetInput1(fixed);
	checkerBoardFilter->SetInput2(resample->GetOutput());

	PatternArrayType pattern;
	pattern[0] = 10; // number of checkers along X
	pattern[1] = 10; // number of checkers along Y
	checkerBoardFilter->SetCheckerPattern(pattern);
	checkerBoardFilter->Update();
	//viewer.AddImage(checkerBoardFilter->GetOutput(), true, "Difference image between fixed and resampled moving");

	// Compute the difference image between the
	// fixed and moving image before registration.
	writer2->SetFileName("differenceBeforeReg.mha");
	difference->SetInput1(fixed);
	difference->SetInput2(moving);
	difference->Update();
	checkerBoardFilter->SetInput2(moving);
	checkerBoardFilter->Update();
	//viewer.AddImage(checkerBoardFilter->GetOutput(), true, "Difference image between fixed and resampled moving");
	try
	{
		writer2->Update();
	}
	catch (itk::ExceptionObject & err)
	{
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
		return NULL;
	}

	// Generate the explicit deformation field resulting from
	// the registration.
	typedef itk::Vector< float, Dimension >          VectorPixelType;
	typedef itk::Image< VectorPixelType, Dimension > DisplacementFieldImageType;

	typedef itk::TransformToDisplacementFieldFilter<
		DisplacementFieldImageType,
		CoordinateRepType >             DisplacementFieldGeneratorType;

	// Create an setup displacement field generator. 
	DisplacementFieldGeneratorType::Pointer dispfieldGenerator =
		DisplacementFieldGeneratorType::New();
	dispfieldGenerator->UseReferenceImageOn();
	dispfieldGenerator->SetReferenceImage(moving);
	dispfieldGenerator->SetTransform(outputBSplineTransform);

	try
	{
		dispfieldGenerator->Update();
	}
	catch (itk::ExceptionObject & err)
	{
		std::cerr << "Exception detected while generating deformation field";
		std::cerr << " : " << err << std::endl;
		return NULL;
	}

	typedef itk::ImageFileWriter< DisplacementFieldImageType >  FieldWriterType;
	FieldWriterType::Pointer fieldWriter = FieldWriterType::New();

	fieldWriter->SetInput(dispfieldGenerator->GetOutput());

	fieldWriter->SetFileName("outputDeformationField.mha");
	try
	{
		fieldWriter->Update();
	}
	catch (itk::ExceptionObject & excp)
	{
		std::cerr << "Exception thrown " << std::endl;
		std::cerr << excp << std::endl;
		return NULL;
	}


	typedef itk::VectorMagnitudeImageFilter<
		DisplacementFieldImageType, ImageType >  VectorMagnitudeFilterType;
	VectorMagnitudeFilterType::Pointer magnitudeFilter =
		VectorMagnitudeFilterType::New();
	magnitudeFilter->SetInput(dispfieldGenerator->GetOutput());

	// To write the magnitude image file, we should rescale the gradient values
	// to a reasonable range

	typedef itk::RescaleIntensityImageFilter<
		ImageType, ImageType > rescaleFilterType;

	rescaleFilterType::Pointer rescaler =
		rescaleFilterType::New();
	rescaler->SetOutputMinimum(0);
	rescaler->SetOutputMaximum(255);
	rescaler->SetInput(magnitudeFilter->GetOutput());
	rescaler->Update();

	//viewer.AddImage(rescaler->GetOutput(), true, "Displacement field magnitudes");


	typedef itk::RGBPixel<unsigned char>    RGBPixelType;
	typedef itk::Image<RGBPixelType, Dimension>  RGBImageType;


	return dispfieldGenerator->GetOutput();
}
};