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

#include "VolcanoFileManager.hpp"
#include "../../version.h"

VolcanoFileManager::VolcanoFileManager()
{
	// TODO Auto-generated constructor stub
}

VolcanoFileManager::~VolcanoFileManager()
{
	// TODO Auto-generated destructor stub
}

void VolcanoFileManager::write(string filename, GridTerrain *terrain, std::vector<double> factors,
							   EruptionParameters parameters, std::vector<int> injectedParticles, double dx, double dt, ESimulatorType simType,
							   double eruptionTotalDuration, std::vector<std::string> vNames)
{

	//std::cout << "avant recuperation du depot" << std::endl;
	std::vector<std::vector<int>> deposition = terrain->getParticleDeposition();
	//std::cout << "apres recuperation du depot" << std::endl;

	double dxterrain = terrain->getDx();
	Int2 size = terrain->getDiscreteSize();
	std::vector<std::shared_ptr<GenericParticleFamily>> families = terrain->getParticleFamilies();

	hid_t file_id;
	hsize_t dims[2] = {hsize_t(size.x_), hsize_t(size.y_)};
	//int data[6] = {1,2,3,4,5,6};
	double *data = new double[size.x_ * size.y_];
	//herr_t status;
	// create the file
	file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	// create one dataset per particle family, each dataset take his name from the id of the particle family
	BOOST_FOREACH (std::shared_ptr<GenericParticleFamily> fam, families)
	{
		//for( std::shared_ptr< GenericParticleFamily > fam: families ){
		if (typeid(*fam) == typeid(TephraParticleFamily))
		{
			// convert data from deposition, which are absolute number of particles, to kg/m^2
			/*
			* \todo this syntax to transform the shared_ptr<GenericParticleFamily> to TephraParticleFamily* is weird...
			* \try to find something else
			*/
			TephraParticleFamily *family = (TephraParticleFamily *)(&(*fam));
			double mass = family->mass_;
			double totalMass;
			double massPerSurfaceUnit;
			for (int y = 0; y < size.y_; y++)
			{
				for (int x = 0; x < size.x_; x++)
				{
					// the total mass of the current particle class for this specific site (times the multiplication factor)
					totalMass = (deposition[terrain->index2d({x, y})])[family->familyId_] * mass * factors[family->familyId_];
					// we want the mass for 1 m^2, so we divide the total mass by the surface of the site
					massPerSurfaceUnit = totalMass / (dxterrain * dxterrain);
					//data[terrain->index2d({x,y})] = massPerSurfaceUnit;
					// stored as y,x in order to read x,y from data
					data[y + x * terrain->getYSize()] = massPerSurfaceUnit;
				}
			}
		}
		// else, the particle family type is not known... if you add some kind of particle family, its up to you to add the
		// support here
		else
		{
			throw - 1;
		}
		//status = H5LTmake_dataset(file_id,("f"+boost::lexical_cast<std::string>(fam->familyId_)).c_str(),2,dims,H5T_NATIVE_DOUBLE,data);
		H5LTmake_dataset(file_id, ("f" + boost::lexical_cast<std::string>(fam->familyId_)).c_str(), 2, dims, H5T_NATIVE_DOUBLE, data);
		//H5LTset_attribute_double( file_id, ("f"+boost::lexical_cast<std::string>(fam->familyId_)).c_str(), "dx", dx, 1 );
		//H5LTset_attribute_double( file_id, ("f"+boost::lexical_cast<std::string>(fam->familyId_)).c_str(), "position", position, 2 );
		if (typeid(*fam) == typeid(TephraParticleFamily))
		{
			TephraParticleFamily *family = (TephraParticleFamily *)(&(*fam));
			double phi[1] = {family->grainSize_};
			double density[1] = {family->density_};
			double weightPct[1] = {family->originalWeightPercent_};
			int injectedPart[1] = {injectedParticles[fam->familyId_]};
			H5LTset_attribute_double(file_id, ("f" + boost::lexical_cast<std::string>(fam->familyId_)).c_str(), "phi", phi, 1);
			H5LTset_attribute_double(file_id, ("f" + boost::lexical_cast<std::string>(fam->familyId_)).c_str(), "density", density, 1);
			H5LTset_attribute_double(file_id, ("f" + boost::lexical_cast<std::string>(fam->familyId_)).c_str(), "weighPct", weightPct, 1);
			H5LTset_attribute_int(file_id, ("f" + boost::lexical_cast<std::string>(fam->familyId_)).c_str(), "injectedParticles", injectedPart, 1);
		}
	}
	// set attributes
	double ht[1];
	double hb[1];
	double duration[1];
	double U0[1];
	double L0[1];
	double n0[1];
	double eccentricity[1];
	double focusSource[2];
	double topSteadyHeight[1];
	double eruptedMass[1];
	double theta0[1];
	double maxWindSpeed[1];
	double windDirection[1];
	if (parameters.useWoodsModel())
	{
		ht[0] = parameters.getHt();
		hb[0] = parameters.getHb();
		duration[0] = parameters.getDuration();
		U0[0] = parameters.getColumnProfile()->getU(0);
		L0[0] = parameters.getColumnProfile()->getL(0);
		n0[0] = parameters.getColumnProfile()->getN0();
		theta0[0] = parameters.getColumnProfile()->getTheta0();
		eccentricity[0] = parameters.getEccentricity();
		focusSource[0] = parameters.getFocusSource().x_;
		focusSource[1] = parameters.getFocusSource().y_;
		topSteadyHeight[0] = parameters.getAtmosphere()->getTopSteadyHeight();
		eruptedMass[0] = parameters.getEruptedMass();
		maxWindSpeed[0] = parameters.getMaxWindSpeed();
		windDirection[0] = parameters.getWindDirection();
	}
	double ventPosition[3] = {parameters.getVentPosition().x_, parameters.getVentPosition().y_, parameters.getVentPosition().z_};
	double columnVerticalDiffusion[1] = {parameters.getColumnVerticalDiffusion()};
	double columnHorizontalDiffusion[1] = {parameters.getColumnHorizontalDiffusion()};
	double tropopause[1] = {parameters.getAtmosphere()->getTropopauseHeight()};
	double T0[1] = {parameters.getAtmosphere()->getT0()};
	double P0[1] = {parameters.getAtmosphere()->getP0()};
	double atmosphereVerticalDiffusion[1] = {parameters.getAtmosphereVerticalDiffusion()};
	double atmosphereHorizontalDiffusion[1] = {parameters.getAtmosphereHorizontalDiffusion()};
	double dxtab[1] = {dx};
	double dttab[1] = {dt};
	std::string simulatorType;
	if (simType == CA)
		simulatorType = "CA";
	if (simType == EXACT)
		simulatorType = "EXACT";
	std::string version = VERSION;
	double dxt[1] = {terrain->getDx()};
	double terrainPosition[2] = {terrain->getPosition().x_, terrain->getPosition().y_};
	double totalDuration[1] = {eruptionTotalDuration};
	std::string speedsStr;
	for (uint i = 0; i < vNames.size(); i++)
	{
		if (i != vNames.size() - 1)
			speedsStr.append(vNames[i].append(", "));
		else
			speedsStr.append(vNames[i]);
	}

	if (parameters.useWoodsModel())
	{
		H5LTset_attribute_double(file_id, "/", "general_ht", ht, 1);
		H5LTset_attribute_double(file_id, "/", "general_hb", hb, 1);
		H5LTset_attribute_double(file_id, "/", "general_duration", duration, 1);
		H5LTset_attribute_double(file_id, "/", "column_U0", U0, 1);
		H5LTset_attribute_double(file_id, "/", "column_L0", L0, 1);
		H5LTset_attribute_double(file_id, "/", "column_n0", n0, 1);
		H5LTset_attribute_double(file_id, "/", "umbrella_eccentricity", eccentricity, 1);
		H5LTset_attribute_double(file_id, "/", "umbrella_focusSource", focusSource, 2);
		H5LTset_attribute_double(file_id, "/", "general_eruptedMass", eruptedMass, 1);
		H5LTset_attribute_double(file_id, "/", "general_totalDuration", totalDuration, 1);
		H5LTset_attribute_double(file_id, "/", "general_ventPosition", ventPosition, 3);
		H5LTset_attribute_double(file_id, "/", "atmosphere_topSteadyHeight", topSteadyHeight, 1);
		H5LTset_attribute_double(file_id, "/", "column_theta0", theta0, 1);
		H5LTset_attribute_double(file_id, "/", "atmosphere_maxWindSpeed", maxWindSpeed, 1);
		H5LTset_attribute_double(file_id, "/", "atmosphere_windDirection", windDirection, 1);
	}

	H5LTset_attribute_double(file_id, "/", "column_vertical_diffusion", columnVerticalDiffusion, 1);
	H5LTset_attribute_double(file_id, "/", "column_horizontal_diffusion", columnHorizontalDiffusion, 1);
	H5LTset_attribute_double(file_id, "/", "atmosphere_tropopause", tropopause, 1);
	H5LTset_attribute_double(file_id, "/", "atmosphere_T0", T0, 1);
	H5LTset_attribute_double(file_id, "/", "atmosphere_P0", P0, 1);
	H5LTset_attribute_double(file_id, "/", "atmosphere_vertical_diffusion", atmosphereVerticalDiffusion, 1);
	H5LTset_attribute_double(file_id, "/", "atmosphere_horizontal_diffusion", atmosphereHorizontalDiffusion, 1);
	H5LTset_attribute_double(file_id, "/", "simulation_dx", dxtab, 1);
	H5LTset_attribute_double(file_id, "/", "simulation_dt", dttab, 1);
	H5LTset_attribute_string(file_id, "/", "simulation_type", simulatorType.c_str());
	H5LTset_attribute_string(file_id, "/", "tetras_version", version.c_str());
	H5LTset_attribute_double(file_id, "/", "terrain_dx", dxt, 1);
	H5LTset_attribute_double(file_id, "/", "terrain_position", terrainPosition, 2);
	H5LTset_attribute_string(file_id, "/", "velocities_name", speedsStr.c_str());
	// close file
	//status = H5Fclose(file_id);
	H5Fclose(file_id);
	delete[] data;
}

EruptionParameters VolcanoFileManager::loadEruptionParameters(string filename, double U0, double L0)
{

	std::string line;
	std::ifstream file(filename.c_str());

	if (file.is_open())
	{
		while (true)
		{
			getline(file, line);
			if (file.eof())
			{
				file.close();
				break;
			}
			line.erase(std::remove(line.begin(), line.end(), '\t'), line.end());
			line.erase(std::remove(line.begin(), line.end(), ' '), line.end());
			if (line.size() > 0)
			{
				if (line.at(0) != '#')
				{
					std::vector<std::string> tokens;
					boost::split(tokens, line, boost::is_any_of(","));
					if (tokens[0].compare("plumeModel") == 0)
					{
						if (tokens[1].compare("woods") == 0)
						{
							file.close();
							return loadWoodsEruptionParameters(filename, U0, L0);
						}
						if (tokens[1].compare("wim") == 0)
						{
							file.close();
							EruptionParameters res = loadWimEruptionParameters(filename);
							return res;
						}
					}
				}
			}
		}
	}
	else
		throw runtime_error("unable to open file " + filename);

	throw runtime_error("eruption file incomplete, plume model type not specified");
}

std::vector<std::shared_ptr<TephraParticleFamily>> VolcanoFileManager::interpolateParticleFamilies(std::vector<std::shared_ptr<TephraParticleFamily>> families, uint n)
{
	// create the vector of original densities, grainsize and weight percent
	std::vector<double> densities, grainSizes, wtpcts;
	for (auto partFam : families)
	{
		densities.push_back(partFam->density_);
		grainSizes.push_back(partFam->grainSize_);
		wtpcts.push_back(partFam->originalWeightPercent_);
	}

	auto grainSizesInterp = createVec(grainSizes.front(), grainSizes.back(), n);
	auto densitiesInterp = interpolate(grainSizes, densities, n).back();
	auto wtpctsInterp = interpolate(grainSizes, wtpcts, n).back();
	wtpctsInterp = preciseNormalize(wtpctsInterp, 100.0);

	// create interpolated particle families
	std::vector<std::shared_ptr<TephraParticleFamily>> familiesInterp;
	for (int id = 0; id < int(grainSizesInterp.size()); id++)
	{
		familiesInterp.push_back(std::make_shared<TephraParticleFamily>(
			TephraParticleFamily(id, grainSizesInterp[id], densitiesInterp[id], wtpctsInterp[id])));
	}

	return familiesInterp;
}

void VolcanoFileManager::loadCommonEruptionParameters(string filename,
													  std::vector<std::shared_ptr<TephraParticleFamily>> &families, Double3 &ventPosition,
													  double &columnVerticalDiffusion, double &columnHorizontalDiffusion, double &tropopause,
													  double &stratosphere, double &T0, double &P0, double &atmosphereVerticalDiffusion,
													  double &atmosphereHorizontalDiffusion, std::vector<std::string> &velocitiesNames, std::string &windModel, Int3 eruptionDate)
{

	std::ifstream file(filename.c_str());
	std::string line;
	int id = 0;
	double values[3];
	bool readSuccess = true;
	uint numParticleFamiliesMult = 1;
	std::vector<std::shared_ptr<TephraParticleFamily>> tmpFamilies;

	if (file.is_open())
	{
		while (true)
		{
			getline(file, line);
			if (file.eof())
				break;
			// remove spaces and tab from line
			line.erase(std::remove(line.begin(), line.end(), '\t'), line.end());
			line.erase(std::remove(line.begin(), line.end(), ' '), line.end());
			try
			{
				if (line.at(0) != '#')
				{
					std::vector<std::string> tokens;
					boost::split(tokens, line, boost::is_any_of(","));

					if (tokens[0].compare("family") == 0)
					{
						values[0] = stod(tokens[1]);
						values[1] = stod(tokens[2]);
						values[2] = stod(tokens[3]);
						tmpFamilies.push_back(std::make_shared<TephraParticleFamily>(
							TephraParticleFamily(id, values[0], values[2], values[1])));
						id++;
					}

					else if (tokens[0].compare("velocities") == 0)
					{
						for (uint i = 1; i < tokens.size(); i++)
						{
							velocitiesNames.push_back(tokens[i]);
						}
					}

					else if (tokens[0].compare("ventPosition") == 0)
					{
						ventPosition = {stod(tokens[1]), stod(tokens[2]), stod(tokens[3])};
					}

					else if (tokens[0].compare("numParticleFamiliesMult") == 0)
					{
						numParticleFamiliesMult = stod(tokens[1]);
					}

					else if (tokens[0].compare("columnVerticalDiffusion") == 0)
					{
						columnVerticalDiffusion = stod(tokens[1]);
					}

					else if (tokens[0].compare("columnHorizontalDiffusion") == 0)
					{
						columnHorizontalDiffusion = stod(tokens[1]);
					}

					else if (tokens[0].compare("tropopause") == 0)
					{
						tropopause = stod(tokens[1]);
					}

					else if (tokens[0].compare("stratosphere") == 0)
					{
						stratosphere = stod(tokens[1]);
					}

					else if (tokens[0].compare("temperatureASL") == 0)
					{
						T0 = stod(tokens[1]);
					}

					else if (tokens[0].compare("pressureASL") == 0)
					{
						P0 = stod(tokens[1]);
					}

					else if (tokens[0].compare("atmosphereHorizontalDiffusion") == 0)
					{
						atmosphereHorizontalDiffusion = stod(tokens[1]);
					}

					else if (tokens[0].compare("atmosphereVerticalDiffusion") == 0)
					{
						atmosphereVerticalDiffusion = stod(tokens[1]);
					}

					else if (tokens[0].compare("windModel") == 0)
					{
						windModel = tokens[1];
					}

					else if (tokens[0].compare("eruptionDate") == 0)
					{
					}

					// else {
					// 	readSuccess = false;
					// 	std::cout<<"error on line : [ " + line + " ]"<<std::endl;
					// }
				}
			}
			catch (out_of_range &e)
			{
			}
		}
	}
	else
		throw std::runtime_error("unable to open file : " + filename);

	if (!readSuccess)
		throw std::runtime_error("syntax error in file : " + filename);

	// values[ 0 ] = stod(tokens[1]);
	// values[ 1 ] = stod(tokens[2]);
	// values[ 2 ] = stod(tokens[3]);
	// families.push_back(std::make_shared<TephraParticleFamily>(
	// 	TephraParticleFamily(id,values[0],values[2],values[1]))
	// );
	// id++;

	// sort particles by size before finishing
	//std::sort(tmpFamilies.begin(), tmpFamilies.end(), [](std::shared_ptr<TephraParticleFamily> i, std::shared_ptr<TephraParticleFamily> j) { return i->diameter_ < j->diameter_; });

	std::sort(tmpFamilies.begin(), tmpFamilies.end(), [](std::shared_ptr<TephraParticleFamily> i, std::shared_ptr<TephraParticleFamily> j) { return i->grainSize_ < j->grainSize_; });	

	if (numParticleFamiliesMult > 1)
		tmpFamilies = interpolateParticleFamilies(tmpFamilies, tmpFamilies.size()*numParticleFamiliesMult);

	// update ids
	id = 0;
	for (auto family : tmpFamilies)
	{
		families.push_back(std::make_shared<TephraParticleFamily>(
			TephraParticleFamily(id, family->grainSize_, family->density_, family->originalWeightPercent_)));
		id++;
	}

	// for( auto family : families ){
	// 	std::cout << "size : " << family->diameter_ << std::endl;
	// }
}

template <class T>
std::vector<T> VolcanoFileManager::readVector(std::vector<std::string> tokens)
{
	std::vector<T> res;
	for (std::string &item : tokens)
		res.push_back(boost::lexical_cast<T>(item));
	return res;
}

void VolcanoFileManager::readWindData(WindProfile &wind, string filename)
{

	double t = -1.0;
	vector<double> height;
	vector<double> directionTo;
	vector<double> speed;

	std::ifstream file(filename.c_str());
	std::string line;
	bool readSuccess = true;

	if (file.is_open())
	{
		while (true)
		{
			getline(file, line);
			if (file.eof())
				break;

			// remove spaces and tab from line
			line.erase(std::remove(line.begin(), line.end(), '\t'), line.end());
			line.erase(std::remove(line.begin(), line.end(), ' '), line.end());
			try
			{
				if (line.at(0) != '#')
				{
					std::vector<std::string> tokens;
					boost::split(tokens, line, boost::is_any_of(","));
					if (tokens[0].compare("time") == 0)
					{
						t = stod(tokens[1]);
					}
					else if (tokens[0].compare("height") == 0)
					{
						height = readVector<double>(std::vector<std::string>(tokens.begin() + 1, tokens.end()));
					}
					else if (tokens[0].compare("speed") == 0)
					{
						speed = readVector<double>(std::vector<std::string>(tokens.begin() + 1, tokens.end()));
					}
					else if (tokens[0].compare("directionTo") == 0)
					{
						directionTo = readVector<double>(std::vector<std::string>(tokens.begin() + 1, tokens.end()));
					}
				}
			}
			catch (out_of_range &e)
			{
			}
		}
	}
	else
		throw std::runtime_error("unable to open file : " + filename);

	if (!readSuccess)
		throw std::runtime_error("syntax error in file : " + filename);

	if (t == -1.0 || height.size() == 0 || directionTo.size() == 0)
	{
		throw std::runtime_error("missing data in file : " + filename);
	}

	wind.addDirectionTo(directionTo, t);
	wind.addHeight(height, t);
	wind.addSpeed(speed, t);
}

void VolcanoFileManager::readPlumeData(WimPlume &plume, string filename)
{

	double t = -1.0;
	double directionTo = -1.0;
	vector<double> angle;
	vector<double> z;
	vector<double> x;
	double ms = -1.0;
	vector<double> r;
	vector<double> u;
	vector<double> m;

	std::ifstream file(filename.c_str());
	std::string line;
	bool readSuccess = true;

	if (file.is_open())
	{
		while (true)
		{
			getline(file, line);
			if (file.eof())
				break;

			// remove spaces and tab from line
			line.erase(std::remove(line.begin(), line.end(), '\t'), line.end());
			line.erase(std::remove(line.begin(), line.end(), ' '), line.end());
			try
			{
				if (line.at(0) != '#')
				{
					std::vector<std::string> tokens;
					boost::split(tokens, line, boost::is_any_of(","));

					if (tokens[0].compare("time") == 0)
					{
						t = stod(tokens[1]);
					}

					else if (tokens[0].compare("directionTo") == 0)
					{
						directionTo = stod(tokens[1]);
					}

					else if (tokens[0].compare("angle") == 0)
					{
						angle = readVector<double>(std::vector<std::string>(tokens.begin() + 1, tokens.end()));
					}

					else if (tokens[0].compare("z") == 0)
					{
						z = readVector<double>(std::vector<std::string>(tokens.begin() + 1, tokens.end()));
					}

					else if (tokens[0].compare("x") == 0)
					{
						x = readVector<double>(std::vector<std::string>(tokens.begin() + 1, tokens.end()));
					}

					else if (tokens[0].compare("m_s") == 0)
					{
						ms = stod(tokens[1]);
					}

					else if (tokens[0].compare("r") == 0)
					{
						r = readVector<double>(std::vector<std::string>(tokens.begin() + 1, tokens.end()));
					}

					else if (tokens[0].compare("u") == 0)
					{
						u = readVector<double>(std::vector<std::string>(tokens.begin() + 1, tokens.end()));
					}

					else if (tokens[0].compare("m") == 0)
					{
						m = readVector<double>(std::vector<std::string>(tokens.begin() + 1, tokens.end()));
					}
				}
			}
			catch (out_of_range &e)
			{
			}
		}
	}
	else
		throw std::runtime_error("unable to open file : " + filename);

	if (!readSuccess)
		throw std::runtime_error("syntax error in file : " + filename);

	if (t == -1.0 || directionTo == -1.0 || angle.size() == 0 || z.size() == 0 ||
		x.size() == 0 || ms == -1.0 || r.size() == 0 || u.size() == 0 || m.size() == 0)
	{
		throw std::runtime_error("missing data in file : " + filename);
	}

	plume.addDirectionTo(directionTo, t);
	plume.addAngle(angle, t);
	plume.addZ(z, t);
	plume.addX(x, t);
	plume.addMs(ms, t);
	plume.addR(r, t);
	plume.addU(u, t);
	plume.addM(m, t);
}

EruptionParameters VolcanoFileManager::loadWimEruptionParameters(string filename)
{

	// common properties
	std::vector<std::shared_ptr<TephraParticleFamily>> families;
	Double3 ventPosition = {-1.0, -1.0, -1.0};
	double columnVerticalDiffusion = -1.0;
	double columnHorizontalDiffusion = -1.0;
	double tropopause = -1.0;
	double stratosphere = -1.0;
	double T0 = -1.0;
	double P0 = -1.0;
	double atmosphereVerticalDiffusion = -1.0;
	double atmosphereHorizontalDiffusion = -1.0;

	std::vector<std::string> velocitiesNames;
	std::string windModel;
	Int3 eruptionDate;

	double eccentricity = -1.0;
	Double2 focusSource = {-1.0, -1.0};

	std::vector<std::tuple<double, double>> eruptionEvents;
	loadCommonEruptionParameters(filename, families, ventPosition, columnVerticalDiffusion, columnHorizontalDiffusion,
								 tropopause, stratosphere, T0, P0, atmosphereVerticalDiffusion, atmosphereHorizontalDiffusion, velocitiesNames,
								 windModel, eruptionDate);

	WimPlume plume(ventPosition);
	WindProfile wind;

	bool readSuccess = true;

	std::ifstream file(filename.c_str());
	string line;

	if (file.is_open())
	{
		while (true)
		{
			getline(file, line);
			if (file.eof())
				break;

			// remove spaces and tab from line
			line.erase(std::remove(line.begin(), line.end(), '\t'), line.end());
			line.erase(std::remove(line.begin(), line.end(), ' '), line.end());
			try
			{
				if (line.at(0) != '#')
				{
					std::vector<std::string> tokens;
					boost::split(tokens, line, boost::is_any_of(","));

					if (tokens[0].compare("eruptionEvent") == 0)
					{
						eruptionEvents.push_back(make_tuple(stod(tokens[1]),
															stod(tokens[2])));
					}

					else if (tokens[0].compare("plume") == 0)
					{
						std::vector<std::string> tokensPath;
						boost::split(tokensPath, filename, boost::is_any_of("/"));
						tokensPath.pop_back();
						std::string path;
						for (std::string item : tokensPath)
							path += item + "/";
						readPlumeData(plume, path + tokens[2]);
					}

					else if (tokens[0].compare("wind") == 0)
					{
						if (windModel.compare("muscle") == 0)
						{
							throw std::runtime_error("Unable to load wind data from files while using muscle for wind model");
						}
						std::vector<std::string> tokensPath;
						boost::split(tokensPath, filename, boost::is_any_of("/"));
						tokensPath.pop_back();
						std::string path;
						for (std::string item : tokensPath)
							path += item + "/";
						readWindData(wind, path + tokens[2]);
					}

					else if (tokens[0].compare("eccentricity") == 0)
					{
						eccentricity = stod(tokens[1]);
					}
					else if (tokens[0].compare("focusSource") == 0)
					{
						focusSource = {stod(tokens[1]), stod(tokens[2])};
					}
				}
			}
			catch (out_of_range &e)
			{
			}
		}
	}
	else
		throw std::runtime_error("unable to open file : " + filename);

	if (!readSuccess)
		throw std::runtime_error("syntax error in file : " + filename);

	if (ventPosition.x_ < 0.0 || ventPosition.y_ < 0.0 || ventPosition.z_ < 0.0 || columnVerticalDiffusion < 0.0 ||
		columnHorizontalDiffusion < 0.0 || families.size() < 1 || tropopause < 0.0 || stratosphere < 0.0 || P0 < 0.0 ||
		T0 < 0.0 || atmosphereVerticalDiffusion < 0.0 || atmosphereHorizontalDiffusion < 0.0)
	{
		throw std::runtime_error("error - some parameters invalid or undefined in file (e1) : " + filename);
	}

	return EruptionParameters::buildWimEruption(families, ventPosition, columnVerticalDiffusion, columnHorizontalDiffusion,
												tropopause, stratosphere, T0, P0, atmosphereVerticalDiffusion, atmosphereHorizontalDiffusion, eruptionEvents, plume,
												wind, eccentricity, focusSource, velocitiesNames, windModel, eruptionDate);
}

EruptionParameters VolcanoFileManager::loadWoodsEruptionParameters(string filename, double U0p, double L0p)
{
	std::string line;
	std::vector<std::shared_ptr<TephraParticleFamily>> families;
	//double values[3];
	double ht = -1.0;
	double hb = -1.0;
	Double3 ventPosition = {-1.0, -1.0, -1.0};
	double duration = -1.0;
	double eruptedMass = -1.0;
	double U0 = U0p;
	double L0 = L0p;
	double n0 = -1.0;
	double theta0 = -1.0;
	double columnVerticalDiffusion = -1.0;
	double columnHorizontalDiffusion = -1.0;
	double eccentricity = -1.0;
	Double2 focusSource = {-1.0, -1.0};
	double tropopause = -1.0;
	double stratosphere = -1.0;
	double T0 = -1.0;
	double P0 = -1.0;
	double maxWindSpeed = -1.0;
	double windDirection = -1.0;
	double atmosphereVerticalDiffusion = -1.0;
	double atmosphereHorizontalDiffusion = -1.0;
	std::vector<std::string> velocitiesNames;
	std::string windModel;
	Int3 eruptionDate;


	loadCommonEruptionParameters(filename, families, ventPosition, columnVerticalDiffusion, columnHorizontalDiffusion,
								 tropopause, stratosphere, T0, P0, atmosphereVerticalDiffusion, atmosphereHorizontalDiffusion, velocitiesNames,
								 windModel, eruptionDate);


	bool readSuccess = true;
	//int id = 0;

	std::ifstream file(filename.c_str());

	if (file.is_open())
	{
		while (true)
		{
			getline(file, line);
			if (file.eof())
				break;

			// remove spaces and tab from line
			line.erase(std::remove(line.begin(), line.end(), '\t'), line.end());
			line.erase(std::remove(line.begin(), line.end(), ' '), line.end());
			try
			{
				if (line.at(0) != '#')
				{
					std::vector<std::string> tokens;
					boost::split(tokens, line, boost::is_any_of(","));
					if (tokens[0].compare("ht") == 0)
					{
						if (tokens[1].compare("unknown") != 0)
							ht = boost::lexical_cast<double>(tokens[1]);
					}
					else if (tokens[0].compare("hb") == 0)
					{
						if (tokens[1].compare("unknown") != 0)
							hb = boost::lexical_cast<double>(tokens[1]);
					}

					// else if( tokens[ 0 ].compare( "ventPosition" ) == 0 ){
					// 	ventPosition = { boost::lexical_cast< double >( tokens[ 1 ] ), boost::lexical_cast< double >( tokens[ 2 ] ), boost::lexical_cast< double >( tokens[ 3 ] ) };
					// }

					// else if( tokens[ 0 ].compare( "family" ) == 0 ){
					// 	values[ 0 ] = stod(tokens[1]);
					// 	values[ 1 ] = stod(tokens[2]);
					// 	values[ 2 ] = stod(tokens[3]);
					// 	families.push_back(std::make_shared<TephraParticleFamily>(TephraParticleFamily(id,values[0],values[2],values[1])));
					// 	id++;
					// }

					else if (tokens[0].compare("duration") == 0)
					{
						duration = stod(tokens[1]);
					}
					else if (tokens[0].compare("eruptedMass") == 0)
					{
						eruptedMass = stod(tokens[1]);
					}
					else if (tokens[0].compare("U0") == 0)
					{
						if (tokens[1].compare("unknown") != 0 && tokens[1].compare("0") != 0)
							U0 = stod(tokens[1]);
					}
					else if (tokens[0].compare("L0") == 0)
					{
						if (tokens[1].compare("unknown") != 0 && tokens[1].compare("0") != 0)
							L0 = stod(tokens[1]);
					}
					else if (tokens[0].compare("n0") == 0)
					{
						if (tokens[1].compare("unknown") != 0)
							n0 = stod(tokens[1]);
					}
					else if (tokens[0].compare("theta0") == 0)
					{
						if (tokens[1].compare("unknown") != 0)
							theta0 = stod(tokens[1]);
					}

					// else if( tokens[ 0 ].compare( "columnVerticalDiffusion" ) == 0 ){
					// 	columnVerticalDiffusion = stod(tokens[1]);
					// }

					// else if( tokens[ 0 ].compare( "columnHorizontalDiffusion" ) == 0 ){
					// 	columnHorizontalDiffusion = stod(tokens[1]);
					// }

					else if (tokens[0].compare("eccentricity") == 0)
					{
						eccentricity = stod(tokens[1]);
					}
					else if (tokens[0].compare("focusSource") == 0)
					{
						focusSource = {stod(tokens[1]), stod(tokens[2])};
					}

					// else if( tokens[ 0 ].compare( "tropopause" ) == 0 ){
					// 	tropopause = stod(tokens[1]);
					// }

					// else if( tokens[ 0 ].compare( "stratosphere" ) == 0 ){
					// 	stratosphere = stod(tokens[1]);
					// }

					// else if( tokens[ 0 ].compare( "temperatureASL" ) == 0 ){
					// 	T0 = stod(tokens[1]);
					// }

					// else if( tokens[ 0 ].compare( "pressureASL" ) == 0 ){
					// 	P0 = stod(tokens[1]);
					// }

					else if (tokens[0].compare("maxWindSpeed") == 0)
					{
						maxWindSpeed = stod(tokens[1]);
					}
					else if (tokens[0].compare("windDirection") == 0)
					{
						windDirection = stod(tokens[1]);
					}

					// else if( tokens[ 0 ].compare( "atmosphereHorizontalDiffusion" ) == 0 ){
					// 	atmosphereHorizontalDiffusion = stod(tokens[1]);
					// }

					// else if( tokens[ 0 ].compare( "atmosphereVerticalDiffusion" ) == 0 ){
					// 	atmosphereVerticalDiffusion = stod(tokens[1]);
					// }

					// else {
					// 	readSuccess = false;
					// 	std::cout<<"error on line : [ " + line + " ]"<<std::endl;
					// }
				}
			}
			catch (out_of_range &e)
			{
			}
		}
	}
	else
		throw std::runtime_error("unable to open file : " + filename);

	if (!readSuccess)
		throw std::runtime_error("syntax error in file : " + filename);

	if (ventPosition.x_ < 0.0 || ventPosition.y_ < 0.0 || ventPosition.z_ < 0.0 || duration < 0.0 || eruptedMass < 0.0 ||
		columnVerticalDiffusion < 0.0 || columnHorizontalDiffusion < 0.0 || eccentricity < 0.0 || focusSource.x_ < 0.0 ||
		focusSource.y_ < 0.0 || families.size() < 1 || tropopause < 0.0 || stratosphere < 0.0 || P0 < 0.0 || T0 < 0.0 ||
		maxWindSpeed < 0.0 || windDirection < 0.0 || atmosphereVerticalDiffusion < 0.0 || atmosphereHorizontalDiffusion < 0.0)
	{
		throw std::runtime_error("error - some parameters invalid or undefined in file (e2) : " + filename);
	}

	if (ht < 0.0 && (U0 < 0.0 || L0 < 0.0 || n0 < 0.0 || theta0 < 0.0))
	{
		throw std::runtime_error("error - some parameters invalid or undefined in file (e3) : " + filename);
	}

	std::cout << "When calling buildWoodsEruption U0 = " << U0 << ", L0 = " << L0 << std::endl;

	return EruptionParameters::buildWoodsEruption(ht, hb, ventPosition, duration, eruptedMass, U0, L0, n0, theta0,
												  columnVerticalDiffusion, columnHorizontalDiffusion, eccentricity, focusSource, families, tropopause, stratosphere,
												  T0, P0, maxWindSpeed, windDirection, atmosphereVerticalDiffusion, atmosphereHorizontalDiffusion, velocitiesNames, eruptionDate);
}
