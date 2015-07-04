#pragma once

template <class ImageType, class FieldType, unsigned int Dimension>
class Visualizer
{
public:
	virtual void visualize() = 0;

	ImageType* getImage() const { return _image; }
	FieldType* getField() const { return _field; }

	void setImage(ImageType* image) { _image = image; }
	void setField(FieldType* field) { _field = field; }

protected:
ImageType* _image;
FieldType* _field;
Visualizer(){};
~Visualizer(){};


};
