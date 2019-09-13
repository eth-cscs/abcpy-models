/*
TEphra TRAsport Simulator (TETRAS)
Copyright (C) 2015  University of Geneva, Switzerland

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#include "Display2D.hpp"

/*
 * Implementation of constructor and destructor
 */

Display2D::Display2D(Domain *dom, GridTerrain *terrain, int ySlice, bool saveVideo, int maxParticle):domain_(dom),
terrain_(terrain),
width_(dom->getXSize()),
height_(dom->getZSize()),
ySlice_(ySlice),
img_(Mat(Size (width_ , height_ ), CV_8UC3 , Scalar(0 ,0 ,0) )),
saveVideo_(saveVideo),
maxParticle_(maxParticle)
{
	namedWindow ("Window", CV_WINDOW_KEEPRATIO );
	if(saveVideo_){
		outputVideo_.open("outputvideo.avi",CV_FOURCC('F','F','V','1'),200.0,cv::Size (width_, height_ ),true);
		if (!outputVideo_.isOpened())
		{
			std::cout  << "Could not open the output video for write" << std::endl;
		}
	}
}

Display2D::~Display2D(){
	destroyWindow("Window");
}

/*
 * Implementation of public methods
 */

void Display2D::show(){
	// compute for each site a number between 0 and 250 (corresponding to the particle density)
	// and convert it into RGB
	//int maxParticle = domain_->getMaxNParticles();
	double gradient;
	shared_array<double> color(new double[3]);
	for(int j=0; j<height_; j++){
		for(int i=0; i<width_; i++){
			if(domain_->isInGroundContactAt({int(i),0,int(j)})){
				// color brown
				color[2] = 139;
				color[1] = 69;
				color[0] = 19;
			}
			else{
				gradient = (double(domain_->getNParticleAt({int(i),int(ySlice_),int(j)}))/double(maxParticle_))*250.0;
				if(gradient>250.0)gradient=250.0;
				groundColorMix(color, gradient, 0.0, 255.0);
			}
			img_.data[(i*3)+(height_-1-j)*width_*3] = char(color[0]);
			img_.data[(i*3)+(height_-1-j)*width_*3+1] = char(color[1]);
			img_.data[(i*3)+(height_-1-j)*width_*3+2] = char(color[2]);
		}
	}

	if(terrain_ != NULL){
		// ...
	}

	imshow ("Window", img_ );
	if(saveVideo_)outputVideo_<<img_;
	waitKey(1);
}

void Display2D::showArray(){
	std::cout<<"[";
	for(int j=0; j<height_; j++){
		std::cout<<"[";
		for(int i=0; i<width_; i++){
			std::cout<<domain_->getNParticleAt({int(i),int(ySlice_),int(j)})<<",";
		}
		std::cout<<"],"<<std::endl;
	}
	std::cout<<"]";
}


void Display2D::groundColorMix(shared_array<double> color, double x, double min, double max)
{
	/*
	 * Red = 0
	 * Green = 1
	 * Blue = 2
	 */
	double posSlope = (max-min)/60;
	double negSlope = (min-max)/60;

	if( x < 60 )
	{
		color[0] = max;
		color[1] = posSlope*x+min;
		color[2] = min;
		return;
	}
	else if ( x < 120 )
	{
		color[0] = negSlope*x+2*max+min;
		color[1] = max;
		color[2] = min;
		return;
	}
	else if ( x < 180  )
	{
		color[0] = min;
		color[1] = max;
		color[2] = posSlope*x-2*max+min;
		return;
	}
	else if ( x < 240  )
	{
		color[0] = min;
		color[1] = negSlope*x+4*max+min;
		color[2] = max;
		return;
	}
	else if ( x < 300  )
	{
		color[0] = posSlope*x-4*max+min;
		color[1] = min;
		color[2] = max;
		return;
	}
	else
	{
		color[0] = max;
		color[1] = min;
		color[2] = negSlope*x+6*max;
		return;
	}
}
