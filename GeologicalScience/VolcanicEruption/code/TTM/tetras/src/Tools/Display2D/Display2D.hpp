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


#ifndef Display2D_H_
#define Display2D_H_

#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <iostream>
#include <boost/shared_array.hpp>

//#include "../../Simulator/GridTerrain.hpp"
//#include "../../Simulator/Domain.hpp"

#include "piaf.hpp"

using namespace cv;
using namespace std;
using namespace boost;

using namespace piaf;

/*!
 * \class Display2D
 * \brief Class used to Display2D the domain
 */
class Display2D {
public:

	/*!
	 * \fn Display2D(Domain &dom)
	 * \brief Constructor of the classe Display2D, the window appear when this object is constructed
	 * \param dom : The Domain to Display2D
	 */
	Display2D(Domain *dom, GridTerrain *terrain, int ySlice, bool saveVideo, int maxParticle);

	/*!
	 * \fn ~Display2D()
	 * \brief Destructor of the classe Display2D
	 */
	~Display2D();

	/*!
	 * \fn show()
	 * \brief show the domain into the window, colors represents particle densities into the domain
	 */
	void show();

	void showArray();


	/*!
	* \fn groundColorMix(shared_array<double> color, double x, double min, double max)
	* \brief Computes the color gradiant (source : http://stackoverflow.com/a/7139909)
	* \param color : The output vector
	* \param x : the gradiant (beetween 0 and 360)
	* \param min : Minimum of the RGB channel
	* \param max : Maximum of the RGB channel
	*/
	static void groundColorMix(shared_array<double> color, double x, double min, double max);

private:

	/*! \brief Reference to the domain to show */
	Domain *domain_;
	/*! \brief Reference to the terrain */
	GridTerrain *terrain_;
	/* \brief Width of the domain */
	int width_;
	/* \brief Height of the domain */
	int height_;
	/* \brief the slice to display */
	int ySlice_;
	/* \brief The Mat OpenCV object used to write RGB color values */
	Mat img_;

	cv::VideoWriter outputVideo_;
	bool saveVideo_;
	int maxParticle_;
};

#endif /* Display2D_H_ */
