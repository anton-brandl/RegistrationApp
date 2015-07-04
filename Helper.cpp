
template<class ImageType , unsigned int Dimension>
int Helper::PaintImageWithHeatmap(itk::Image<ImageType, Dimension>::Pointer image, itk::Image<ImageType, Dimension>::Pointer map)
{

	typedef itk::ScalarToRGBColormapImageFilter<FixedImageType, RGBImageType> RGBFilterType;
	RGBFilterType::Pointer rgbfilter = RGBFilterType::New();
	rgbfilter->SetInput(map);
	rgbfilter->SetColormap(RGBFilterType::Hot);

	//viewer.AddRGBImage(rgbfilter->GetOutput(), true, "Displacement field magnitudes heatmap");

	typedef itk::ImageToVTKImageFilter<HeatmapType> RGBConnectorType;
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


	rgbfilter->SetInput(image);
	rgbfilter->SetColormap(RGBFilterType::Grey);
	rgbfilter->Update();
	rgbconnector->SetInput(rgbfilter->GetOutput());
	rgbconnector->Update();

	vtkSmartPointer<vtkImageData> anatomy = rgbconnector->GetOutput();



	//GLYPHS

	// Setup the arrows

	/*
	// Create actors
	vtkSmartPointer<vtkImageSliceMapper> imageMapper = vtkSmartPointer<vtkImageSliceMapper>::New();

	imageMapper->SetInputData(anatomy);


	vtkSmartPointer<vtkImageSlice> imageSlice = vtkSmartPointer<vtkImageSlice>::New();
	imageSlice->SetMapper(imageMapper);

	vtkSmartPointer<vtkArrowSource> arrowSource = vtkSmartPointer<vtkArrowSource>::New();
	arrowSource->Update();
	vtkSmartPointer<vtkGlyph3D> glyphFilter = vtkSmartPointer<vtkGlyph3D>::New();
	glyphFilter->SetSourceConnection(arrowSource->GetOutputPort());
	glyphFilter->OrientOn();
	glyphFilter->SetVectorModeToUseVector();
	glyphFilter->SetInputData(vFieldConnector->GetOutput());
	glyphFilter->SetScaleModeToDataScalingOff();
	glyphFilter->SetScaleFactor(3);
	glyphFilter->Update();

	vtkSmartPointer<vtkPolyDataMapper> vectorMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	vectorMapper->SetInputConnection(glyphFilter->GetOutputPort());
	vectorMapper->ScalarVisibilityOn();

	vtkSmartPointer<vtkActor> vectorActor = vtkSmartPointer<vtkActor>::New();
	vectorActor->SetMapper(vectorMapper);

	// Setup renderer
	vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();


	// Setup render window
	vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->AddRenderer(renderer);

	// Setup render window interactor
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor2 = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	// Render and start interaction
	renderWindowInteractor2->SetRenderWindow(renderWindow);

	renderer->AddViewProp(imageSlice);
	renderer->AddActor(vectorActor);

	renderWindow->Render();
	//renderWindowInteractor2->Initialize();

	renderWindowInteractor2->Start();

	*/






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

	return EXIT_SUCCESS;
}