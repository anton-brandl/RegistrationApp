#include "configwin.h"

itk::SmartPointer<itk::Image<float, 2>> loadDicomImage(std::string src)
{
	const    unsigned int    ImageDimension = 2;
	typedef itk::Image< float, ImageDimension > ImageType;
	typedef itk::ImageFileReader< ImageType  > ImageReaderType;

	ImageReaderType::Pointer  imageReader = ImageReaderType::New();

	imageReader->SetFileName(src);
	imageReader->Update();
	return imageReader->GetOutput();
}

ConfigWin::ConfigWin(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::ConfigWin)
{
	ui->setupUi(this);
	image_view_fixed = vtkSmartPointer<vtkResliceImageViewer>::New();
	image_view_moving = vtkSmartPointer<vtkResliceImageViewer>::New();

	ui->txtFixedImagePath->setText("endExhale.IMA");
	ui->txtMovingImagePath->setText("endInhale.IMA");
	ui->txtOutputPath->setText("output.mha");
	ui->comboTransformType->addItem("Translation");
	ui->comboMetricType->addItem("MeanSquares"); 
	ui->comboOptimizerType->addItem("RegularStepGradientDescent");
	

	connect(ui->btnTest, SIGNAL(released()), this, SLOT(testOutput()));
	connect(ui->btnRegister, SIGNAL(released()), this, SLOT(registerImages()));
	connect(ui->btnSelectMovingImage, SIGNAL(released()), this, SLOT(browseMovingImage()));
	connect(ui->btnSelectFixedImage, SIGNAL(released()), this, SLOT(browseFixedImage()));
	connect(ui->btnSelectOutputPath, SIGNAL(released()), this, SLOT(browseOutput()));
	connect(ui->txtFixedImagePath, SIGNAL(editingFinished()), this, SLOT(updateShowFixedImage()));
	connect(ui->txtMovingImagePath, SIGNAL(editingFinished()), this, SLOT(updateShowMovingImage()));
}
 
ConfigWin::~ConfigWin()
{
    delete ui;
}


void ConfigWin::testOutput()
{
	std::cout << "Test";
}
 
void ConfigWin::registerImages()
{
	
	QString output(ui->txtOutputPath->text());
	std::string after = "after.png";
	std::string before = "before.png";

	
//	ImageLoader loader;
	typedef itk::Image<PixelType, 2> ImageType2D;
	typedef itk::Vector< PixelType, 2 >          VectorPixelType2D;
	typedef itk::Image< VectorPixelType2D, 2 > DisplacementFieldImageType2D;

	typedef itk::Image<PixelType, 3> VolumeType3D;
	typedef itk::Vector< PixelType, 3 > VectorPixelType3D;
	typedef itk::Image< VectorPixelType3D, 3 > DisplacementFieldImageType3D;

	if (ui->rbImage->isChecked()) {
		Registrator reg;
		itk::SmartPointer<DisplacementFieldImageType2D> field = reg.demonsAsymetricRegistration<ImageType2D, DisplacementFieldImageType2D, 2>(fixedImage, movingImage);

		HeatmapVisualizer<ImageType2D, DisplacementFieldImageType2D, 2> heatmapVis;

		heatmapVis.setImage(fixedImage);
		heatmapVis.setField(field);
		heatmapVis.visualize();
	}
	else if (ui->rbVolume->isChecked()) {
		Registrator reg;
		DisplacementFieldImageType3D::Pointer field = reg.demonsAsymetricRegistration<VolumeType3D, DisplacementFieldImageType3D, 3>(fixedVolume, movingVolume);

		HeatmapVisualizer<VolumeType3D, DisplacementFieldImageType3D, 3> heatmapVis;

		heatmapVis.setImage(fixedVolume);
		heatmapVis.setField(field);
		heatmapVis.visualize();
	}
	//deformableRegistrationTest(image1.toStdString(), image2.toStdString(), output.toStdString(), "outputDisplacement.mha", "outputField.mha");

}

void ConfigWin::browseMovingImage()
{
	QString filename;
	if (ui->rbVolume->isChecked()) {
		filename = QFileDialog::getExistingDirectory(this, tr("Open Moving Volume"));
	}
	else if (ui->rbImage->isChecked()) {
		filename = QFileDialog::getOpenFileName(this, tr("Open Moving Image"), "", "All Files (*.*);;Siemens Dicom (*.IMA);;Dicom (*.dcm)");
	}
	
	ui->txtMovingImagePath->setText(filename);
	updateShowMovingImage();

}

void ConfigWin::browseFixedImage()
{
	QString filename;
	if (ui->rbVolume->isChecked()) {
		filename = QFileDialog::getExistingDirectory(this, tr("Open Moving Volume"));
	}
	else if (ui->rbImage->isChecked()) {
		filename = QFileDialog::getOpenFileName(this, tr("Open Moving Image"), "", "All Files (*.*);;Siemens Dicom (*.IMA);;Dicom (*.dcm)");
	}
	ui->txtFixedImagePath->setText(filename);
	updateShowFixedImage();
}

void ConfigWin::browseOutput()
{
	QString filename = QFileDialog::getSaveFileName(this, tr("Open Fixed Image"), "", "All Files (*.*);;Siemens Dicom (*.IMA);;Dicom (*.dcm)");
	ui->txtOutputPath->setText(filename);
	updateShowMovingImage();
}

void ConfigWin::updateShowFixedImage() {
	updateShowImage(Imagetype::Fixed);
}

void ConfigWin::updateShowMovingImage() {
	updateShowImage(Imagetype::Moving);
}

void ConfigWin::updateShowImage(Imagetype type)
{
	loadImage(type);
	typedef  float           PixelType;
	typedef itk::Image< PixelType, 2 >  ImageType2D;
	typedef itk::Image< PixelType, 3 >  ImageType3D;


	//typedef itk::ImageFileReader< ImageType  >   ImageReaderType;
	//ImageReaderType::Pointer   imageReader = ImageReaderType::New();
	typedef itk::ImageToVTKImageFilter<ImageType2D> ConnectorType2D;
	typedef itk::ImageToVTKImageFilter<ImageType3D> ConnectorType3D;
	
	std::string filename = "";
	QVTKWidget* widget;
	vtkSmartPointer<vtkResliceImageViewer> image_view;
	ConnectorType2D::Pointer connector2d = ConnectorType2D::New();
	ConnectorType3D::Pointer connector3d = ConnectorType3D::New();
	if (type == Imagetype::Fixed)
	{
		
		widget = ui->qvtkFixed;
		image_view = image_view_fixed;
		if (ui->rbVolume->isChecked()) {
			connector3d->SetInput(fixedVolume);
		} else if (ui->rbImage->isChecked()) {
			connector2d->SetInput(fixedImage);
		}
	}
	else if (type == Imagetype::Moving)
	{
		
		widget = ui->qvtkMoving;		
		image_view = image_view_moving; 
		if (ui->rbVolume->isChecked()) {
			connector3d->SetInput(movingVolume);
		}
		else if (ui->rbImage->isChecked()) {
			connector2d->SetInput(movingImage);
		}
	}

	vtkSmartPointer<vtkImageActor> actor =
		vtkSmartPointer<vtkImageActor>::New();

	vtkImageData * image = vtkImageData::New();
	if (ui->rbVolume->isChecked()) {
		connector3d->Update();
		image->DeepCopy(connector3d->GetOutput());
	}
	else if (ui->rbImage->isChecked()){
		connector2d->Update();
		image->DeepCopy(connector2d->GetOutput());
	}
	
	vtkSmartPointer<vtkImageFlip> flipYFilter =
		vtkSmartPointer<vtkImageFlip>::New();
	flipYFilter->SetFilteredAxis(1); // flip y axis
	flipYFilter->SetInputData(image);
	flipYFilter->Update();
	//set VTK Viewer to QVTKWidget in Qt's UI
	widget->SetRenderWindow(image_view->GetRenderWindow());
	image_view->SetupInteractor(widget->GetRenderWindow()->GetInteractor());
	//Set input image to VTK viewer
	image_view->SetInputData(flipYFilter->GetOutput());
	image_view->SetSlice(image_view->GetSliceMax() / 2);
	image_view->GetRenderer()->ResetCamera();
	image_view->Render();


	widget->update();
}

void ConfigWin::loadImage(Imagetype type) {
	QString fixed(ui->txtFixedImagePath->text());
	QString moving(ui->txtMovingImagePath->text());

	ImageLoader loader;
	typedef itk::Image<PixelType, 2> ImageType2D;
	typedef itk::Vector< PixelType, 2 >          VectorPixelType2D;
	typedef itk::Image< VectorPixelType2D, 2 > DisplacementFieldImageType2D;

	typedef itk::Image<PixelType, 3> VolumeType3D;
	typedef itk::Vector< PixelType, 3 > VectorPixelType3D;
	typedef itk::Image< VectorPixelType3D, 3 > DisplacementFieldImageType3D;

	if (ui->rbImage->isChecked()) {
		if (type==Imagetype::Fixed)
			fixedImage = loadDicomImage(fixed.toStdString());
		else if (type==Imagetype::Moving)
			movingImage = loadDicomImage(moving.toStdString());
	}
	else if (ui->rbVolume->isChecked()) {
		if (type == Imagetype::Fixed)
			fixedVolume = loader.loadDicomVolume<PixelType>(fixed.toStdString());
		else if (type == Imagetype::Moving)
			movingVolume = loader.loadDicomVolume<PixelType>(moving.toStdString());
	}
}