#include "image.h"
#include <math.h>
#include <stdlib.h>
#include <algorithm>
#include <string.h>
#include <float.h>
#include <iostream>

/**
 * Image
 **/
Image::Image (int width_, int height_){

    assert(width_ > 0);
    assert(height_ > 0);

    width           = width_;
    height          = height_;
    num_pixels      = width * height;
    sampling_method = IMAGE_SAMPLING_POINT;
    
    data.raw = new uint8_t[num_pixels*4];
	int b = 0; //which byte to write to
	for (int j = 0; j < height; j++){
		for (int i = 0; i < width; i++){
			data.raw[b++] = 0;
			data.raw[b++] = 0;
			data.raw[b++] = 0;
			data.raw[b++] = 0;
		}
	}

    assert(data.raw != NULL);
}

Image::Image (const Image& src){
	
	width           = src.width;
    height          = src.height;
    num_pixels      = width * height;
    sampling_method = IMAGE_SAMPLING_POINT;
    
    data.raw = new uint8_t[num_pixels*4];
    
    //memcpy(data.raw, src.data.raw, num_pixels);
    *data.raw = *src.data.raw;
}

Image::Image (char* fname){

	int numComponents; //(e.g., Y, YA, RGB, or RGBA)
	data.raw = stbi_load(fname, &width, &height, &numComponents, 4);
	
	if (data.raw == NULL){
		printf("Error loading image: %s", fname);
		exit(-1);
	}
	

	num_pixels = width * height;
	sampling_method = IMAGE_SAMPLING_POINT;
	
}

Image::~Image (){
    delete data.raw;
    data.raw = NULL;
}

void Image::Write(char* fname){
	
	int lastc = strlen(fname);

	switch (fname[lastc-1]){
	   case 'g': //jpeg (or jpg) or png
	     if (fname[lastc-2] == 'p' || fname[lastc-2] == 'e') //jpeg or jpg
	        stbi_write_jpg(fname, width, height, 4, data.raw, 95);  //95% jpeg quality
	     else //png
	        stbi_write_png(fname, width, height, 4, data.raw, width*4);
	     break;
	   case 'a': //tga (targa)
	     stbi_write_tga(fname, width, height, 4, data.raw);
	     break;
	   case 'p': //bmp
	   default:
	     stbi_write_bmp(fname, width, height, 4, data.raw);
	}
}

void Image::AddNoise (double factor)
{
	int x, y;
	for (x = 0; x < Width(); x++) {
		for (y = 0; y < Height(); y++) {
			Pixel p = GetPixel(x, y);
			GetPixel(x, y) = PixelLerp(p, PixelRandom(), factor);
		}
	}
}

void Image::Brighten (double factor)
{
	int x,y;
	for (x = 0 ; x < Width() ; x++)
	{
		for (y = 0 ; y < Height() ; y++)
		{
			Pixel p = GetPixel(x, y);
			Pixel scaled_p = p*factor;
			GetPixel(x,y) = scaled_p;
		}
	}
}


void Image::ChangeContrast (double factor)
{
	int sum = 0;
	for (int x = 0; x < Width(); ++x) {
		for (int y = 0; y < Height(); ++y) {
			sum += GetPixel(x, y).Luminance();
		}
	}
	sum /= Width() * Height();
	Pixel temp;
	temp.Set(sum, sum, sum);

	for (int x = 0; x < Width(); ++x) {
		for (int y = 0; y < Height(); ++y) {
			GetPixel(x, y) = PixelLerp(temp, GetPixel(x, y), factor);
		}
	}
}


void Image::ChangeSaturation(double factor)
{
	for (int x = 0; x < Width(); ++x) {
		for (int y = 0; y < Height(); ++y) {
			Pixel p = GetPixel(x, y);
			int temp = p.Luminance();
			Pixel localGray(temp, temp, temp,255);
		    GetPixel(x,y) = PixelLerp(localGray, p, factor);
		}
	}
}


Image* Image::Crop(int x, int y, int w, int h)
{
	if (x + w > Width()) {
		//reduce until an image is cropable
		if (x < Width()) {
			w = Width() - x;
		}
		else { //make nothing
			w = 0;
		}
	}
	if (y + h > Height()) {
		if (y < Height()) {
			y = Height() - y;
		}
		else {
			h = 0;
		}
	}
	Image* img = new Image(w, h);
	for (int i = 0; i < w; ++i) {
		for (int j = 0; j < h; ++j) {
			img->GetPixel(i, j) = GetPixel(x + i, y + j);
		}
	}
	return img;
}


void Image::ExtractChannel(int channel)
{
	for (int i = 0; i < Width(); ++i) {
		for (int j = 0; j < Height(); ++j) {
			Pixel p = GetPixel(i, j);
			switch (channel) {
				case(IMAGE_CHANNEL_RED) :
					p.Set(p.r, 0, 0, 0);
					break;
				case(IMAGE_CHANNEL_BLUE) :
					p.Set(0, 0, p.b, 0);
					break;
				case(IMAGE_CHANNEL_GREEN) :
					p.Set(0, p.g, 0, 0);
					break;
				case(IMAGE_CHANNEL_ALPHA) :
					p.Set(0, 0, 0, p.a);
					break;
				default :
					return;
			}
			GetPixel(i, j) = p;
		}
	}
}


void Image::Quantize (int nbits)
{
	for (int i = 0; i < Height(); ++i) {
		for (int j = 0; j < Width(); ++j) {
			SetPixel(j, i, PixelQuant(GetPixel(j, i), nbits));
		}
	}
}

void Image::RandomDither (int nbits)
{
	this->AddNoise(0.3);
	this->Quantize(nbits);
}


static int Bayer4[4][4] =
{
    {15,  7, 13,  5},
    { 3, 11,  1,  9},
    {12,  4, 14,  6},
    { 0,  8,  2, 10}
};


void Image::OrderedDither(int nbits)
{
	// 15  7 13  5 15...
	//  3 11  1  9  3
	// 12  4 14  6 12
	//  0  8  2 10  0
	// 15  7 13  5 15

	double error = 0;
	int scale = int(pow(2, 8 - nbits));
	for (int i = 0; i < Height(); ++i) {
		for (int j = 0; j < Width(); ++j) {
			Pixel p = PixelQuant(GetPixel(j, i), nbits);
			error = GetPixel(j, i).r - p.r;
			/*if (j == 1 && i == 1) {
				std::cout << "error is " << error << std::endl;
				std::cout << "Get pixel" << int(GetPixel(j, i).r) << std::endl;
				std::cout << "Quantized " << int(p.r) << std::endl;
			}*/

			if (error > Bayer4[j % 4][i % 4] * scale / 16) {
				p.r += scale;
			}
			// else do nothing i.e. floor(Img)
			error = GetPixel(j, i).g - p.g;
			if (error > Bayer4[j % 4][i % 4] * scale / 16) {
				p.g += scale;
			}
			error = GetPixel(j, i).b - p.b;
			if (error > Bayer4[j % 4][i % 4] * scale / 16) {
				p.b += scale;
			}
			SetPixel(j, i, p);
		}
	}

}

/* Error-diffusion parameters */
const double
    ALPHA = 7.0 / 16.0,
    BETA  = 3.0 / 16.0,
    GAMMA = 5.0 / 16.0,
    DELTA = 1.0 / 16.0;

void Image::FloydSteinbergDither(int nbits)
{
	double error = 0;
	double* cumError = new double[Height() * Width() * 3];
	//std::cout << Height()*Width() * 3 << std::endl;
	for (int i = 0; i < Height(); ++i) {
		for (int j = 0; j < Width(); ++j) {
			for (int k = 0; k < 3; ++k) {
				//std::cout << "Which index did it go wrong " << k + j * 3 + i*Width() * 3 << std::endl;
				cumError[ k + j*3 + i*Width()*3] = 0;
			}
		}
	}
	for (int i = 0; i < Height(); ++i) {
		//std::cout << "Which i did it go wrong " << i << std::endl;
		for (int j = 0; j < Width(); ++j) {
			GetPixel(j, i).r += int(cumError[0 + j * 3 + i*Width() * 3]);
			GetPixel(j, i).g += int(cumError[1 + j * 3 + i*Width() * 3]);
			GetPixel(j, i).b += int(cumError[2 + j * 3 + i*Width() * 3]);
			Pixel p = GetPixel(j, i);
			Pixel newPixel = PixelQuant(p, nbits);
			SetPixel(j, i, newPixel);

			error = p.r - newPixel.r;
			/*if (i == 0 && j < 10){
				std::cout << "error in red " << error << std::endl;
				std::cout << "p.r and newPixel.r" << double(p.r) <<' '<< double(newPixel.r)<<std::endl;
			}*/
			if (j + 1 < Width()) { //ignore out of bounds
				cumError[0 + (j + 1) * 3 + i*Width() * 3] += error *ALPHA;
				if (i + 1 < Height()) {
					cumError[0 + (j + 1) * 3 + (i + 1)*Width() * 3] += error *DELTA;
				}
			}
			if (i + 1 < Height()) {
				cumError[0 + (j) * 3 + (i + 1)*Width() * 3] += error *GAMMA;
				if (j - 1 > 0) {
					cumError[0 + (j - 1) * 3 + (i + 1)*Width() * 3] += error *BETA;
				}
			}

			error = p.g - newPixel.g;

			if (j + 1 < Width()) { //ignore out of bounds
				cumError[1 + (j + 1) * 3 + i*Width() * 3] += error *ALPHA;
				if (i + 1 < Height()) {
					cumError[1 + (j + 1) * 3 + (i + 1)*Width() * 3] += error *DELTA;
				}
			}
			if (i + 1 < Height()) {
				cumError[1 + (j) * 3 + (i + 1)*Width() * 3] += error *GAMMA;
				if (j - 1 > 0) {
					cumError[1 + (j - 1) * 3 + (i + 1)*Width() * 3] += error *BETA;
				}
			}

			error = p.b - newPixel.b;

			if (j + 1 < Width()) { //ignore out of bounds
				cumError[2 + (j + 1) * 3 + i*Width() * 3] += error *ALPHA;
				if (i + 1 < Height()) {
					cumError[2 + (j + 1) * 3 + (i + 1)*Width() * 3] += error *DELTA;
				}
			}
			if (i + 1 < Height()) {
				cumError[2 + (j) * 3 + (i + 1)*Width() * 3] += error *GAMMA;
				if (j - 1 > 0) {
					cumError[2 + (j - 1) * 3 + (i + 1)*Width() * 3] += error *BETA;
				}
			}

		}

	}

	delete[] cumError;

}

void Image::Blur(int n)
{
	//If mask is odd, pixel will center on the middle mask
	//If mask is even, pixel will center on the top left off center.
	// |  |  |  |  |
	// |  | x|  |  |
	// |  |  |  |  |
	// |  |  |  |  |
	//create a mask of size n by n
	double* mask = new double[n * n];
	double midPoint;
	if (n % 2) {//if odd
		midPoint = n / 2;
	}
	else {
		midPoint = (n - 1) / 2.0;
	}
	double sigma = (double)n / 4;
	double sum = 0;
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			double temp = exp(-(i - midPoint)*(i - midPoint) / sigma / sigma / 2) *
				exp(-(j - midPoint)*(j - midPoint) / sigma / sigma / 2);
			temp /= (2 * 3.141592654 * sigma * sigma);
			mask[i + j * n] = temp;
			sum += temp;
		}
	}
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			mask[i + j*n] /= sum;
		}
	}

	//debug output
	/*for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			printf("  %f  ", mask[i + j * n]);
		}
		printf("\n");
	}*/

	//one day I will write a helper function
	// Convolute(Height(), Width(), n, a, b);
	//convolute
	Image newImage(Width(), Height());
	double r = 0;
	double g = 0;
	double b = 0;
	int lagY = 0;
	int lagX = 0;
	//std::cout << "Width and Height " << Width() << ' ' << Height() << std::endl;
	for (int i = 0; i < Height(); ++i) {
		for (int j = 0; j < Width(); ++j) {
			/*if (i == 0 && j == 0) {
				std::cout << "R G B" << double(GetPixel(j, i).r) 
					<< ' ' << double(GetPixel(j, i).g) << ' ' << double(GetPixel(j, i).b)<<std::endl;
			}*/

			if (n % 2) { //if n is odd
				lagY = i - (n / 2);
			}
			else { //if n is even
				lagY = i - ((n - 1) / 2);
			}
			int pixelX;
			int pixelY;
			r = 0; g = 0; b = 0;
			for (int k = 0; k < n; ++k) { //height of convolution matrix
				pixelY = lagY;
				//Account for edge cases in Y direction
				if (lagY < 0 || lagY >= Height()) {
					if (lagY < 0) {
						//lagY = lagY % -Height(); //returns (-Height(),0)
						//lagY = -lagY; //returns (0,n) mirror operation
						pixelY = -(lagY % -Height());
					}
					else { //lag > Height()
						//lagY = lagY % Height();
						//lagY = Height() - lagY;
						pixelY = (Height()-2) - (lagY % Height() );
					}
				}

				if (n % 2) {
					lagX = j - (n / 2);
				}
				else {
					lagX = j - ((n - 1) / 2);
				}
				for (int l = 0; l < n; ++l) { //width of convolution matrix
					pixelX = lagX;
					//Account for edge cases in X direction
					if (lagX < 0 || lagX >= Width()) {
						if (lagX < 0) {
							pixelX = -(lagX % -Width());
						}
						else {
							pixelX = (Width()-2) - (lagX % Width());
						}
					}
					
					/*if (i == 0 && j == 0) {
						std::cout << "lag values" << lagX << std::endl;
						std::cout << "pixel x and y" << pixelX << ' ' << pixelY << std::endl;
						std::cout << "R G B " << double(GetPixel(j, i).r) << ' ' 
							<< double(GetPixel(j, i).g) << ' ' << double(GetPixel(j, i).b) << std::endl;
						std::cout << mask[k*n + l] << std::endl;
						std::cout << double(GetPixel(pixelX, pixelY).r) << ' ' << double(GetPixel(pixelX, pixelY).g)<<' '
							<< double(GetPixel(pixelX, pixelY).b) << std::endl;
						std::cout << std::endl;
					}*/
					//converting double to component possible loss of data
					r += mask[k * n + l] * double(GetPixel(pixelX, pixelY).r);
					//if (i == 0 && j == 0)
					//	std::cout << "value of g " << g << std::endl;
					g += mask[k * n + l] * double(GetPixel(pixelX, pixelY).g);
					b += mask[k * n + l] * double(GetPixel(pixelX, pixelY).b);
					lagX++;
				}
				lagY++;
			}
			/*if (i == 0 && j == 0) {
				std::cout << "r g b " << r << ' '
					<< g << ' ' << b << std::endl;
			}*/
			Pixel p(int(r), int(g), int(b),255);
			newImage.SetPixel(j, i, p);
		}

	}

	for (int i = 0; i < Height(); ++i) {
		for (int j = 0; j < Width(); ++j) {
			SetPixel(j, i, newImage.GetPixel(j,i));
			/*if (i == 0 && j == 0) {
				std::cout << "R G B" << double(GetPixel(j, i).r) << ' ' << double(GetPixel(j, i).g) << ' ' << double(GetPixel(j, i).b);
			}*/
		}
	}

	delete[] mask;
}

void Image::Sharpen(int n)
{
	Image temp = *this;
	temp.Blur(n); //get blurred version
	//At this point I decided to follow the row major method to access the data
	//Also I decided to use SetPixel for ease of readibility. Stylistic wise GetPixel
	//should adhere to const correctness and not be use as in Brighten(). But I am too lazy
	//to change things
	for (int i = 0; i < Height(); ++i) {
		for (int j = 0; j < Width(); ++j) {
			SetPixel(j, i, PixelLerp(temp.GetPixel(j, i), GetPixel(j, i), 2.0));
		}
	}

}

void Image::EdgeDetect()
{
	this->ChangeContrast(0.5);
	this->Brighten(0.1);
	this->ChangeContrast(2);
	//this->ChangeSaturation(0);
	this->Brighten(0.1);
	//this->ChangeContrast(0.5);
	//edge detect kernel obtained from Wikipedia
	//https://en.wikipedia.org/wiki/Kernel_(image_processing)
	double* kernel = new double[3 * 3];
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			kernel[i * 3 + j] = -1;
		}
	}
	kernel[4] = 8;

	Image newImage(Width(), Height());
	int lagY = 0;
	int lagX = 0;
	//double r;
	//double g;
	//double b;
	//std::cout << "Width and Height " << Width() << ' ' << Height() << std::endl;
	for (int i = 0; i < Height(); ++i) {
		//std::cout << "Height changes" << i << std::endl;
		for (int j = 0; j < Width(); ++j) {
			//std::cout << "Width changes " << j << std::endl;
			lagY = i - 1;
			int pixelX;
			int pixelY;
			//r = 0; g = 0; b = 0;
			for (int k = 0; k < 3; ++k) { //height of convolution matrix
				pixelY = lagY;
				//Account for edge cases in Y direction
				if (lagY < 0 || lagY >= Height()) {
					if (lagY < 0) {
						//lagY = lagY % -Height(); //returns (-Height(),0)
						//lagY = -lagY; //returns (0,n) mirror operation
						pixelY = -(lagY % -Height());
					}
					else { //lag > Height()
						   //lagY = lagY % Height();
						   //lagY = Height() - lagY;
						pixelY = (Height() - 2) - (lagY % Height());
					}
				}

				lagX = j - 1;
				for (int l = 0; l < 3; ++l) { //width of convolution matrix
					pixelX = lagX;
					//Account for edge cases in X direction
					/*if (i == 0 && j == 898) {
						std::cout << "lagX is " << lagX << std::endl;
					}*/
					if (lagX < 0 || lagX >= Width()) {
						if (lagX < 0) {
							pixelX = -(lagX % -Width());
						}
						else {
							pixelX = (Width() - 2) - (lagX % Width());
						}
					}
					/*if (i == 0 && j == 898) {
						std::cout << "pixel x and y" << pixelX << ' ' << pixelY << std::endl;
					}*/
					newImage.GetPixel(j,i).r += kernel[k * 3 + l] * GetPixel(pixelX, pixelY).r;
					newImage.GetPixel(j,i).g += kernel[k * 3 + l] * GetPixel(pixelX, pixelY).g;
					newImage.GetPixel(j,i).b += kernel[k * 3 + l] * GetPixel(pixelX, pixelY).b;
					lagX++;
				}
				lagY++;
			}
			//Pixel p(int(r), int(g), int(b),255);
			//newImage.SetPixel(j, i, p);
		}
	}
	//std::cout << "Maybe I am broken here" << std::endl;
	for (int i = 0; i < Height(); ++i) {
		for (int j = 0; j < Width(); ++j) {
			SetPixel(j, i, newImage.GetPixel(j,i));
		}
	}

	delete[] kernel;

	//this->ChangeContrast(2);
	this->ChangeSaturation(0);
	//this->ChangeContrast(2);
	//this->Quantize(1);
}

Image* Image::Scale(double sx, double sy)
{
	//scale to the nearest int
	int nWidth = int(Width() * sx);
	int nHeight = int(Height() * sy);
	//std::cout << "ori width and height " << Width() << ' ' << Height() << std::endl;
	//std::cout << "scaled width and height " << nWidth << ' ' << nHeight << std::endl;
	Image* newImage = new Image(nWidth, nHeight);
	for (int i = 0; i < nHeight; ++i) {
		//std::cout << "I am wrong when " << i << std::endl;
		for (int j = 0; j < nWidth; ++j) {
			//if(j >= 1797)
			//	std::cout << "J am wrong when " << j << std::endl;
			double iWidth = double(j) / sx;
			double iHeight = double(i) / sy;
			if (i == 0 && j == 0) {
				//std::cout << "iWidth and iHeight " << iWidth << ' ' << iHeight << std::endl;
				Pixel p = Sample(iWidth, iHeight);
				//std::cout << "Sampled pixel " << int(p.r) << ' ' << int(p.g) << ' ' << int(p.b) << std::endl;
			}
			newImage->SetPixel(j, i, Sample(iWidth, iHeight));
		}
	}
	return newImage;
}

Image* Image::Rotate(double angle)
{
	//assume angle is in radians
	//remember original image
	double vertices[] = {
		0, double(Height()),//top left
		double(Width()), double(Height()),//top right
		0, 0,//bottom left
		double(Width()), 0//bottom right
	};

	int minX = 0;
	int maxX = 0;
	int minY = 0;
	int maxY = 0;

	//rotate the image
	for (int i = 0; i < 8; i += 2) {
		double oldx = vertices[i];
		double oldy = vertices[i + 1];
		vertices[i] = oldx * cos(angle) - oldy * sin(angle);
		vertices[i + 1] = oldx * sin(angle) + oldy * cos(angle);
		//std::cout << "some vertice value x " << vertices[i] << std::endl;
		//std::cout << "some vertice value y " << vertices[i + 1] << std::endl;
		//std::cout << "the min values for y " << minY << std::endl;
		if (vertices[i] > maxX || vertices[i] < minX) {
			if (vertices[i] > maxX) {
				maxX = int(vertices[i] + 0.5);
			}
			else {
				minX = int(vertices[i] + 0.5);
			}
		}
		if (vertices[i + 1] > maxY || vertices[i + 1] < minY) {
			if (vertices[i + 1] > maxY) {
				maxY = int(vertices[i + 1] + 0.5);
			}
			else {
				minY = int(vertices[i + 1] + 0.5);
			}
		}
	}

	//std::cout << "Width and Height " << Width() << ' ' << Height() << std::endl;
	int nWidth = int(maxX - minX);
	int nHeight = int(maxY - minY);
	//std::cout << "ori width and height " << Width() << ' ' << Height() << std::endl;
	//std::cout << "nWidth and nHeight " << nWidth << ' ' << nHeight << std::endl;
	//std::cout << "minx and min y " << minX << ' ' << minY << std::endl;

	Image* newImage = new Image(nWidth, nHeight); 
	for (int i = 0; i < nHeight; ++i) {
		for (int j = 0; j < nWidth; ++j) {
			double newx = minX + j;
			double newy = maxY - i;
			/*if (i == 334 && j == 0) {
				std::cout << "new x and y " << newx << ' ' << newy << std::endl;
			}*/
			double x = newx * cos(-angle) - newy * sin(-angle);
			double y = newx * sin(-angle) + newy * cos(-angle);
			/*if (i == 334 && j == 0) {
				std::cout << "This is double " << x << ' ' << y << std::endl;
			}*/
			if (x < 0 || x >= Width() || y < 0 || y >= Height()){
				newImage->GetPixel(j, i).Set(0, 0, 0, 0);
			}
			else {
				/*Pixel p;
				if (i < 10) {
					p.Set(255, 0, 255, 255);
				}
				else {
					p.Set(255, 255, 255, 255);
				}
				newImage->SetPixel(j, i, p);*/
				newImage->SetPixel(j, i, Sample(x, Height() - y));
				/*if ( j == 0) {
					std::cout << "We have color " << int(newImage->GetPixel(j, i).r) 
						<<" at newx newy" <<newx<<' '<<newy<< std::endl;
				}*/
			}
		}
	}
	return newImage;
}

void Image::Fun()
{
	Image newImage(Width(), Height());
	double a = std::min(Height(), Width())/4;
	double r;
	int newX = Width() / 2;
	int newY = Height() / 2;
	//std::cout << "ori width and height " << Width() << ' ' << Height() << std::endl;
	for (int i = 0; i < Height(); ++i) {
		//std::cout << "THe mistake in I " << i << std::endl;
		for (int j = 0; j < Width(); ++j) {
			//if(i == 141)
			//	std::cout << "THe mistake in j " << j << std::endl;
			int tempX = j - newX;
			int tempY = i - newY;
			r = tempX * tempX + tempY * tempY;
			r = sqrt(r);
			
			if (r >= std::min(Height(), Width())) {
				Pixel p(0, 0, 0, 0);
				newImage.SetPixel(j, i, p);
			}
			else {
				double angle = r / a;
				/*if (i == 141 && j < 200) {
					std::cout << angle << std::endl;
				}*/
				double x = tempX * cos(-angle) - tempY * sin(-angle);
				double y = tempX * sin(-angle) + tempY * cos(-angle);
				x += newX;
				y += newY;

				/*if (i == 141 && j == 635) {
					std::cout << "x and y" << x << ' ' << y << std::endl;
				}*/
				if (int(x) < 0 || int(x) >= Width() || int(y) < 0 || int(y) >= Height()) {
					if (int(x) >= Width() && int(y) >= Height()) {
						newImage.SetPixel(j, i, Sample(Width() - 1, Height() - 1));
					}
					else if (int(x) >= Width()) {
						if (int(y) < 0) {
							newImage.SetPixel(j, i, Sample(Width() - 1, 0));
						}
						else {
							newImage.SetPixel(j, i, Sample(Width() - 1, y));
						}
					}
					else if (int(y) >= Height()) {
						if (int(x) < 0) {
							newImage.SetPixel(j, i, Sample(0, Height() - 1));
						} else {
							newImage.SetPixel(j, i, Sample(x, Height() - 1));
						}
					}
					else if (int(x) < 0) {
						if (int(y) > Height()) {
							newImage.SetPixel(j, i, Sample(0, Height() - 1));
						}
						else {
							newImage.SetPixel(j, i, Sample(0, y));
						}
					}
					else {
						if (int(x) > Width()) {
							newImage.SetPixel(j, i, Sample(Width() - 1, 0));
						}
						else {
							newImage.SetPixel(j, i, Sample(x, 0));
						}
					}
						
				}
				else {
					//std::cout << "something worked" << std::endl;
					newImage.SetPixel(j, i, Sample(x, y));
				}
				
			}
		}
	}

	for (int i = 0; i < Height(); ++i) {
		for (int j = 0; j < Width(); ++j) {
			SetPixel(j, i, newImage.GetPixel(j, i));
		}
	}
}

/**
 * Image Sample
 **/
void Image::SetSamplingMethod(int method)
{
    assert((method >= 0) && (method < IMAGE_N_SAMPLING_METHODS));
    sampling_method = method;
}


Pixel Image::Sample (double u, double v){
	if (sampling_method == IMAGE_SAMPLING_POINT) {
		return GetPixel(int(u), int(v));
	}
	if (sampling_method == IMAGE_SAMPLING_BILINEAR) {
		double x1 = floor(u);
		double x2 = ceil(u);
		if (int(x2) >= Width()) {
			x2 = x1;
		}
		double y1 = floor(v);
		double y2 = ceil(v);
		if (int(y2) >= Height()) {
			y2 = y1;
		}

		double midx1_r, midx1_g, midx1_b, midx2_r, midx2_g, midx2_b;
		if (y2 == y1) {
			midx1_r = double(GetPixel(x1, y1).r);
			midx1_g = double(GetPixel(x1, y1).g);
			midx1_b = double(GetPixel(x1, y1).b);

			midx2_r = double(GetPixel(x2, y1).r);
			midx2_g = double(GetPixel(x2, y1).g);
			midx2_b = double(GetPixel(x2, y1).b);
		}
		else {
			midx1_r = (v - y1) / (y2 - y1) * double(GetPixel(x1, y1).r)
				+ (y2 - v) / (y2 - y1) * double(GetPixel(x1, y2).r);
			midx1_g = (v - y1) / (y2 - y1) * double(GetPixel(x1, y1).g)
				+ (y2 - v) / (y2 - y1) * double(GetPixel(x1, y2).g);
			midx1_b = (v - y1) / (y2 - y1) * double(GetPixel(x1, y1).b)
				+ (y2 - v) / (y2 - y1) * double(GetPixel(x1, y2).b);
			midx2_r = (v - y1) / (y2 - y1) * double(GetPixel(x2, y1).r)
				+ (y2 - v) / (y2 - y1) * double(GetPixel(x2, y2).r);
			midx2_g = (v - y1) / (y2 - y1) * double(GetPixel(x2, y1).g)
				+ (y2 - v) / (y2 - y1) * double(GetPixel(x2, y2).g);
			midx2_b = (v - y1) / (y2 - y1) * double(GetPixel(x2, y1).b)
				+ (y2 - v) / (y2 - y1) * double(GetPixel(x2, y2).b);
		}
		
		double final_r, final_g, final_b;
		if (x1 == x2) {
			final_r = midx1_r;
			final_g = midx1_g;
			final_b = midx1_b;
		}
		else {
			final_r = (u - x1) / (x2 - x1) * midx1_r + (x2 - u) / (x2 - x1) *midx2_r;
			final_g = (u - x1) / (x2 - x1) * midx1_g + (x2 - u) / (x2 - x1) *midx2_g;
			final_b = (u - x1) / (x2 - x1) * midx1_b + (x2 - u) / (x2 - x1) *midx2_b;
		}
		

		Pixel p(int(final_r), int(final_g), int(final_b), 255);
		return p;
	}
	if (sampling_method == IMAGE_SAMPLING_GAUSSIAN) {
		//This method is stupid slow because I wanted to have a scale/rotate 
		//function that is agnostic of the sampling method.
		//so it must be done with an n^{2} algorithm rather than an 2n algorithm
		double radius = 2; //width of gaussian to take
		double sigma = 2 * radius / 4;
		double sum = 0;
		double r = 0;
		double g = 0;
		double b = 0;

		for (int x = floor(u - radius); x < ceil(u + radius); ++x) {
			int pixelx = x;
			if (x < 0 || x >= Width()) {
				if (x < 0) {
					pixelx = -(x % -Width());
				}
				else {
					pixelx = (Width() - 2) - (x % Width());
				}
			}
			if (pixelx == -1) {
				std::cout << "pixelx " << pixelx << std::endl;
				std::cout << "x " << x << std::endl;
			}
			for (int y = floor(v - radius); y < ceil(v + radius); ++y) {
				//check oob
				int pixely = y;
				if (y < 0 || y >= Height()) {
					if (y < 0) {
						pixely = -(y % -Height());
					}
					else {
						pixely = (Height() - 2) - (y % Height());
					}
				}
				
				if ( pixely == -1) {
					std::cout << "pixely " << pixely << std::endl;
					std::cout << "y is " << y << std::endl;
					std::cout << "v is " << v << std::endl;
				}
				
				double weight = exp(-(pixelx - u)*(pixelx - u) / sigma / sigma / 2) *
					exp(-(pixely - v)*(pixely - v) / sigma / sigma / 2);
				sum += weight;
				
				r += double(GetPixel(pixelx, pixely).r) * weight;
				g += double(GetPixel(pixelx, pixely).g) * weight;
				b += double(GetPixel(pixelx, pixely).b) * weight;
			}
		}

		Pixel p(int(r / sum),int (g / sum),int( b / sum), 255);
		return p;
	}
	return Pixel();
}