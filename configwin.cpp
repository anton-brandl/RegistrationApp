#include "configwin.h"

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
	QString image1(ui->txtFixedImagePath->text());
	QString image2(ui->txtMovingImagePath->text());
	QString output(ui->txtOutputPath->text());
	std::string after = "after.png";
	std::string before = "before.png";

	/*const    unsigned int    Dimension = 2;
	typedef  float           PixelType;

	typedef itk::Image< PixelType, Dimension >  FixedImageType;
	typedef itk::Image< PixelType, Dimension >  MovingImageType;

	typedef itk::ImageFileReader< FixedImageType  >   FixedImageReaderType;
	typedef itk::ImageFileReader< MovingImageType >   MovingImageReaderType;
	FixedImageReaderType::Pointer   fixedImageReader = FixedImageReaderType::New();
	MovingImageReaderType::Pointer  movingImageReader = MovingImageReaderType::New();

	fixedImageReader->SetFileName(image1.toStdString());
	movingImageReader->SetFileName(image2.toStdString());
*/
	
	
	//registerTest2D(image1.toStdString(), image2.toStdString(), output.toStdString(), after, before);
	deformableRegistrationTest(image1.toStdString(), image2.toStdString(), output.toStdString(), "outputDisplacement.mha", "outputField.mha");
}

void ConfigWin::browseMovingImage()
{
	QString filename = QFileDialog::getOpenFileName(this, tr("Open Moving Image"), "", "All Files (*.*);;Siemens Dicom (*.IMA);;Dicom (*.dcm)");
	ui->txtMovingImagePath->setText(filename);


}

void ConfigWin::browseFixedImage()
{
	QString filename = QFileDialog::getOpenFileName(this, tr("Open Fixed Image"), "", "All Files (*.*);;Siemens Dicom (*.IMA);;Dicom (*.dcm)");
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
	const    unsigned int    Dimension = 2;
	typedef  float           PixelType;

	typedef itk::Image< PixelType, Dimension >  ImageType;

	typedef itk::ImageFileReader< ImageType  >   ImageReaderType;
	ImageReaderType::Pointer   imageReader = ImageReaderType::New();
	typedef itk::ImageToVTKImageFilter<ImageType> ConnectorType;
	
	std::string filename = "";
	QVTKWidget* widget;
	vtkSmartPointer<vtkResliceImageViewer> image_view;
	if (type == Imagetype::Fixed)
	{
		filename = ui->txtFixedImagePath->text().toStdString();
		widget = ui->qvtkFixed;
		image_view = image_view_fixed;

	}
	else if (type == Imagetype::Moving)
	{
		filename = ui->txtMovingImagePath->text().toStdString();
		widget = ui->qvtkMoving;		
		image_view = image_view_moving;
	}

	imageReader->SetFileName(filename);
	ConnectorType::Pointer connector = ConnectorType::New();

	connector->SetInput(imageReader->GetOutput());

	vtkSmartPointer<vtkImageActor> actor =
		vtkSmartPointer<vtkImageActor>::New();

	connector->Update();
	vtkImageData * image = vtkImageData::New();
	image->DeepCopy(connector->GetOutput());

	//set VTK Viewer to QVTKWidget in Qt's UI
	widget->SetRenderWindow(image_view->GetRenderWindow());
	image_view->SetupInteractor(widget->GetRenderWindow()->GetInteractor());
	//Set input image to VTK viewer
	image_view->SetInputData(image);
	image_view->SetSlice(image_view->GetSliceMax() / 2);
	image_view->GetRenderer()->ResetCamera();
	image_view->Render();


	widget->update();

	
}