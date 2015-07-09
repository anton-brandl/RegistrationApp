#ifndef CONFIGWIN_H
#define CONFIGWIN_H
 
#include "ui_configwin.h"
#include <iostream> 
#include <QFileDialog>

#include "itkImageToVTKImageFilter.h"


#include "vtkSmartPointer.h"
#include "vtkImageActor.h"
#include "vtkImageMapper3D.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkResliceImageViewer.h"

#include "IO/LoadSeriesSafeAsVolume.h"
#include "Registration/Registrator.h" 

namespace Ui {
    class ConfigWin;
}
 
class ConfigWin : public QDialog
{
    Q_OBJECT
 
public:
    explicit ConfigWin(QWidget *parent = 0);
	~ConfigWin();
	enum class Imagetype { Fixed, Moving };

private slots:
	void testOutput();
	void registerImages();
	void browseMovingImage();
	void browseFixedImage();
	void browseOutput();
	void updateShowFixedImage();
	void updateShowMovingImage();

private:
	void updateShowImage(ConfigWin::Imagetype type);
	Ui::ConfigWin *ui;
	vtkSmartPointer<vtkResliceImageViewer> image_view_moving;
	vtkSmartPointer<vtkResliceImageViewer> image_view_fixed;
};
 
#endif // CONFIGWIN_H
