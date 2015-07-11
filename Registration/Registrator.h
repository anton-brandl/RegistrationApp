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
#include "itkLBFGSBOptimizerv4.h"
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

#include "../Visualization/HeatmapVisualizer.h"
#include "../IO/ImageLoader.h"
#include <itkSmartPointer.h>
#include "itkCommand.h"


#include "itkImageRegionIterator.h"

#include "itkDemonsRegistrationFilter.h"
#include "itkHistogramMatchingImageFilter.h"
#include "itkWarpImageFilter.h"

class CommandIterationUpdateBSpline : public itk::Command
{
public:
	typedef  CommandIterationUpdateBSpline   Self;
	typedef  itk::Command             Superclass;
	typedef itk::SmartPointer<Self>   Pointer;
	itkNewMacro(Self);

protected:
	CommandIterationUpdateBSpline() {};

public:
	typedef itk::LBFGSBOptimizerv4     OptimizerType;
	typedef   const OptimizerType *    OptimizerPointer;

	void Execute(itk::Object *caller, const itk::EventObject & event) ITK_OVERRIDE
	{
		Execute((const itk::Object *)caller, event);
	}

		void Execute(const itk::Object * object, const itk::EventObject & event) ITK_OVERRIDE
	{
		OptimizerPointer optimizer = static_cast< OptimizerPointer >(object);
		if (!(itk::IterationEvent().CheckEvent(&event)))
		{
			return;
		}
		std::cout << optimizer->GetCurrentIteration() << "   ";
		std::cout << optimizer->GetCurrentMetricValue() << "   ";
		std::cout << optimizer->GetInfinityNormOfProjectedGradient() << std::endl;
	}
};

class CommandIterationUpdateDemons : public itk::Command
{
public:
	typedef  CommandIterationUpdateDemons                     Self;
	typedef  itk::Command                               Superclass;
	typedef  itk::SmartPointer<CommandIterationUpdateDemons>  Pointer;
	itkNewMacro(CommandIterationUpdateDemons);
protected:
	CommandIterationUpdateDemons() {};

	typedef itk::Image< float, 2 >            InternalImageType;
	typedef itk::Vector< float, 2 >           VectorPixelType;
	typedef itk::Image<  VectorPixelType, 2 > DisplacementFieldType;

	typedef itk::DemonsRegistrationFilter<
		InternalImageType,
		InternalImageType,
		DisplacementFieldType>   RegistrationFilterType;

public:

	void Execute(itk::Object *caller, const itk::EventObject & event) ITK_OVERRIDE
	{
		Execute((const itk::Object *)caller, event);
	}

		void Execute(const itk::Object * object, const itk::EventObject & event) ITK_OVERRIDE
	{
		const RegistrationFilterType * filter = static_cast< const RegistrationFilterType * >(object);
		if (!(itk::IterationEvent().CheckEvent(&event)))
		{
			return;
		}
		std::cout << filter->GetMetric() << std::endl;
	}
};
class Registrator {

public: 
int translationRegistration2DTest(std::string image1, std::string image2, std::string output, std::string differenceAfter = "", std::string differenceBefore = "");

template<class ImageType, class FieldType, unsigned int Dimension>
typename itk::SmartPointer<FieldType> bsplinesRegistration(itk::SmartPointer<ImageType> fixed, itk::SmartPointer<ImageType> moving){


	const unsigned int SplineOrder = 3;
	typedef double CoordinateRepType;

	typedef itk::BSplineTransform<
		CoordinateRepType,
		Dimension,
		SplineOrder >     TransformType;
	// Software Guide : EndCodeSnippet

	typedef itk::LBFGSBOptimizerv4       OptimizerType;


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

	unsigned int numberOfGridNodesInOneDimension = 20;

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

	const unsigned int numParameters =
		outputBSplineTransform->GetNumberOfParameters();
	OptimizerType::BoundSelectionType boundSelect(numParameters);
	OptimizerType::BoundValueType upperBound(numParameters);
	OptimizerType::BoundValueType lowerBound(numParameters);

	boundSelect.Fill(OptimizerType::UNBOUNDED);
	upperBound.Fill(0.0);
	lowerBound.Fill(0.0);

	optimizer->SetBoundSelection(boundSelect);
	optimizer->SetUpperBound(upperBound);
	optimizer->SetLowerBound(lowerBound);

	optimizer->SetCostFunctionConvergenceFactor(1e+12);
	optimizer->SetGradientConvergenceTolerance(1.0e-35);
	optimizer->SetNumberOfIterations(100);
	optimizer->SetMaximumNumberOfFunctionEvaluations(100);
	optimizer->SetMaximumNumberOfCorrections(5);
	// Software Guide : EndCodeSnippet

	// Create the Command observer and register it with the optimizer.
	//
	CommandIterationUpdateBSpline::Pointer observer = CommandIterationUpdateBSpline::New();
	optimizer->AddObserver(itk::IterationEvent(), observer);

	std::cout << "Starting Registration "
		<< std::endl;


//VERSION 2 END



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


template<class ImageType, class FieldType, unsigned int Dimension>
typename itk::SmartPointer<FieldType> demonsAsymetricRegistration(itk::SmartPointer<ImageType> fixed, itk::SmartPointer<ImageType> moving){

		// Software Guide : BeginCodeSnippet
		typedef float                                      InternalPixelType;
		typedef itk::Image< InternalPixelType, Dimension > InternalImageType;
		typedef itk::CastImageFilter< ImageType,
			InternalImageType >  ImageCasterType;

		ImageCasterType::Pointer fixedImageCaster = ImageCasterType::New();
		ImageCasterType::Pointer movingImageCaster = ImageCasterType::New();

		fixedImageCaster->SetInput(fixed);
		movingImageCaster->SetInput(moving);
		// Software Guide : EndCodeSnippet

		// Software Guide : BeginLatex
		//
		// The demons algorithm relies on the assumption that pixels representing the
		// same homologous point on an object have the same intensity on both the
		// fixed and moving images to be registered. In this example, we will
		// preprocess the moving image to match the intensity between the images
		// using the \doxygen{HistogramMatchingImageFilter}.
		//
		// \index{itk::HistogramMatchingImageFilter}
		//
		// The basic idea is to match the histograms of the two images at a
		// user-specified number of quantile values. For robustness, the histograms
		// are matched so that the background pixels are excluded from both
		// histograms.  For MR images, a simple procedure is to exclude all gray
		// values that are smaller than the mean gray value of the image.
		//
		// Software Guide : EndLatex

		// Software Guide : BeginCodeSnippet
		typedef itk::HistogramMatchingImageFilter<
			InternalImageType,
			InternalImageType >   MatchingFilterType;
		MatchingFilterType::Pointer matcher = MatchingFilterType::New();
		// Software Guide : EndCodeSnippet


		// Software Guide : BeginLatex
		//
		// For this example, we set the moving image as the source or input image and
		// the fixed image as the reference image.
		//
		// \index{itk::HistogramMatchingImageFilter!SetInput()}
		// \index{itk::HistogramMatchingImageFilter!SetSourceImage()}
		// \index{itk::HistogramMatchingImageFilter!SetReferenceImage()}
		//
		// Software Guide : EndLatex

		// Software Guide : BeginCodeSnippet
		matcher->SetInput(movingImageCaster->GetOutput());
		matcher->SetReferenceImage(fixedImageCaster->GetOutput());
		// Software Guide : EndCodeSnippet


		// Software Guide : BeginLatex
		//
		// We then select the number of bins to represent the histograms and the
		// number of points or quantile values where the histogram is to be
		// matched.
		//
		// \index{itk::HistogramMatchingImageFilter!SetNumberOfHistogramLevels()}
		// \index{itk::HistogramMatchingImageFilter!SetNumberOfMatchPoints()}
		//
		// Software Guide : EndLatex

		// Software Guide : BeginCodeSnippet
		matcher->SetNumberOfHistogramLevels(1024);
		matcher->SetNumberOfMatchPoints(7);
		// Software Guide : EndCodeSnippet


		// Software Guide : BeginLatex
		//
		// Simple background extraction is done by thresholding at the mean
		// intensity.
		//
		// \index{itk::HistogramMatchingImageFilter!ThresholdAtMeanIntensityOn()}
		//
		// Software Guide : EndLatex

		// Software Guide : BeginCodeSnippet
		matcher->ThresholdAtMeanIntensityOn();
		// Software Guide : EndCodeSnippet


		// Software Guide : BeginLatex
		//
		// In the \doxygen{DemonsRegistrationFilter}, the deformation field is
		// represented as an image whose pixels are floating point vectors.
		//
		// \index{itk::DemonsRegistrationFilter}
		//
		// Software Guide : EndLatex

		// Software Guide : BeginCodeSnippet
//		typedef itk::Vector< float, Dimension >           VectorPixelType;
	//	typedef itk::Image<  VectorPixelType, Dimension > DisplacementFieldType;
		typedef itk::DemonsRegistrationFilter<
			InternalImageType,
			InternalImageType,
			FieldType> RegistrationFilterType;
		RegistrationFilterType::Pointer filter = RegistrationFilterType::New();
		// Software Guide : EndCodeSnippet


		// Create the Command observer and register it with the registration filter.
		//
		CommandIterationUpdateDemons::Pointer observer = CommandIterationUpdateDemons::New();
		filter->AddObserver(itk::IterationEvent(), observer);


		// Software Guide : BeginLatex
		//
		// The input fixed image is simply the output of the fixed image casting
		// filter.  The input moving image is the output of the histogram matching
		// filter.
		//
		// \index{itk::DemonsRegistrationFilter!SetFixedImage()}
		// \index{itk::DemonsRegistrationFilter!SetMovingImage()}
		//
		// Software Guide : EndLatex

		// Software Guide : BeginCodeSnippet
		filter->SetFixedImage(fixedImageCaster->GetOutput());
		filter->SetMovingImage(matcher->GetOutput());
		// Software Guide : EndCodeSnippet


		// Software Guide : BeginLatex
		//
		// The demons registration filter has two parameters: the number of
		// iterations to be performed and the standard deviation of the Gaussian
		// smoothing kernel to be applied to the deformation field after each
		// iteration.
		// \index{itk::DemonsRegistrationFilter!SetNumberOfIterations()}
		// \index{itk::DemonsRegistrationFilter!SetStandardDeviations()}
		//
		// Software Guide : EndLatex

		// Software Guide : BeginCodeSnippet
		filter->SetNumberOfIterations(50);
		filter->SetStandardDeviations(1.0);
		// Software Guide : EndCodeSnippet


		// Software Guide : BeginLatex
		//
		// The registration algorithm is triggered by updating the filter. The
		// filter output is the computed deformation field.
		//
		// Software Guide : EndLatex

		// Software Guide : BeginCodeSnippet
		filter->Update();
		// Software Guide : EndCodeSnippet


		// Software Guide : BeginLatex
		//
		// The \doxygen{WarpImageFilter} can be used to warp the moving image with
		// the output deformation field. Like the \doxygen{ResampleImageFilter},
		// the \code{WarpImageFilter} requires the specification of the input image to be
		// resampled, an input image interpolator, and the output image spacing and
		// origin.
		//
		// \index{itk::WarpImageFilter}
		// \index{itk::WarpImageFilter!SetInput()}
		// \index{itk::WarpImageFilter!SetInterpolator()}
		// \index{itk::WarpImageFilter!SetOutputSpacing()}
		// \index{itk::WarpImageFilter!SetOutputOrigin()}
		//
		// Software Guide : EndLatex

		// Software Guide : BeginCodeSnippet
		typedef itk::WarpImageFilter<
			ImageType,
			ImageType,
			FieldType  >     WarperType;
		typedef itk::LinearInterpolateImageFunction<
			ImageType,
			double          >  InterpolatorType;
		WarperType::Pointer warper = WarperType::New();
		InterpolatorType::Pointer interpolator = InterpolatorType::New();

		warper->SetInput(moving);
		warper->SetInterpolator(interpolator);
		warper->SetOutputSpacing(fixed->GetSpacing());
		warper->SetOutputOrigin(fixed->GetOrigin());
		warper->SetOutputDirection(fixed->GetDirection());
		// Software Guide : EndCodeSnippet


		// Software Guide : BeginLatex
		//
		// Unlike \code{ResampleImageFilter}, \code{WarpImageFilter}
		// warps or transforms the input image with respect to the deformation field
		// represented by an image of vectors.  The resulting warped or resampled
		// image is written to file as per previous examples.
		//
		// \index{itk::WarpImageFilter!SetDisplacementField()}
		//
		// Software Guide : EndLatex

		// Software Guide : BeginCodeSnippet
		warper->SetDisplacementField(filter->GetOutput());
		// Software Guide : EndCodeSnippet


		// Write warped image out to file
		typedef  unsigned char                           OutputPixelType;
		typedef itk::Image< OutputPixelType, Dimension > OutputImageType;
		typedef itk::CastImageFilter<
			ImageType,
			OutputImageType >          CastFilterType;
		typedef itk::ImageFileWriter< OutputImageType >  WriterType;

		WriterType::Pointer      writer = WriterType::New();
		CastFilterType::Pointer  caster = CastFilterType::New();

		writer->SetFileName("field.mha");

		caster->SetInput(warper->GetOutput());
		writer->SetInput(caster->GetOutput());
		writer->Update();

		
		if (Dimension == 2) {
			typedef FieldType            VectorImage2DType;
			typedef FieldType::PixelType Vector2DType;

			VectorImage2DType::ConstPointer vectorImage2D = filter->GetOutput();

			VectorImage2DType::RegionType  region2D = vectorImage2D->GetBufferedRegion();
			VectorImage2DType::IndexType   index2D = region2D.GetIndex();
			VectorImage2DType::SizeType    size2D = region2D.GetSize();


			typedef itk::Vector< float, 3 >  Vector3DType;
			typedef itk::Image< Vector3DType, 3 >  VectorImage3DType;

			typedef itk::ImageFileWriter< VectorImage3DType > VectorImage3DWriterType;

			VectorImage3DWriterType::Pointer writer3D = VectorImage3DWriterType::New();

			VectorImage3DType::Pointer vectorImage3D = VectorImage3DType::New();

			VectorImage3DType::RegionType  region3D;
			VectorImage3DType::IndexType   index3D;
			VectorImage3DType::SizeType    size3D;

			index3D[0] = index2D[0];
			index3D[1] = index2D[1];
			index3D[2] = 0;

			size3D[0] = size2D[0];
			size3D[1] = size2D[1];
			size3D[2] = 1;

			region3D.SetSize(size3D);
			region3D.SetIndex(index3D);

			vectorImage3D->SetRegions(region3D);
			vectorImage3D->Allocate();

			typedef itk::ImageRegionConstIterator< VectorImage2DType > Iterator2DType;

			typedef itk::ImageRegionIterator< VectorImage3DType > Iterator3DType;

			Iterator2DType  it2(vectorImage2D, region2D);
			Iterator3DType  it3(vectorImage3D, region3D);

			it2.GoToBegin();
			it3.GoToBegin();

			Vector2DType vector2D;
			Vector3DType vector3D;

			vector3D[2] = 0; // set Z component to zero.

			while (!it2.IsAtEnd())
			{
				vector2D = it2.Get();
				vector3D[0] = vector2D[0];
				vector3D[1] = vector2D[1];
				it3.Set(vector3D);
				++it2;
				++it3;
			}


			writer3D->SetInput(vectorImage3D);

			writer3D->SetFileName("field3D.mha");

			try
			{
				writer3D->Update();
			}
			catch (itk::ExceptionObject & excp)
			{
				std::cerr << excp << std::endl;
				return NULL;
			}
		}
			
		return filter->GetOutput();
	



		}
};