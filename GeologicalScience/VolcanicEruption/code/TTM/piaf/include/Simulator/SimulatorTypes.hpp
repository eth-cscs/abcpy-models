/*
Particles in Advection Field (piaf)
Copyright (C) 2015  University of Geneva, Switzerland

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef SIMULATORTYPES_H_
#define SIMULATORTYPES_H_

#define UNIT_M 0
#define UNIT_MM 1

#include <boost/serialization/access.hpp>
#include <cmath>

namespace piaf{

	// taken from http://www.cplusplus.com/forum/articles/3638/
	template <typename FloatType>
	FloatType roundhalfup( const FloatType& value )
	{
		return std::floor( value +0.5 );
	}

	template <typename FloatType>
	FloatType roundhalfdown( const FloatType& value )
	{
		return std::ceil( value -0.5 );
	}

	template <typename FloatType>
	FloatType roundhalftowardzero( const FloatType& value )
	{
		return ( value >= 0.0 ) ? roundhalfdown( value ) : roundhalfup( value );
	}

	/*! \todo : utiliser des templates pour les types xxx2 et xxx3 ? */

	/*!
	* \class Uint
	* \brief unsigned int
	*/
	typedef unsigned int Uint;

	//------ Added by Mohamed------------
	/// Integer type for Tetra (taken from palabos)
	/** On some architectures, this type is larger
	*  than int. Using plint instead of int ensures
	*  64-bit compatibility of the code.
	*/
	#ifdef ARC_BGP
	typedef long long int plint;
	#else
	typedef ptrdiff_t plint;
	#endif

	/// Unsigned integer type Tetra (taken from palabos)
	/** On some architectures, this type is larger
	*  than int. Using fluplint instead of unsigned plint
	*  ensures 64-bit compatibility of the code.
	**/
	#ifdef ARC_BGP
	typedef unsigned long long int pluint;
	#else
	typedef size_t pluint;
	#endif
	//------------------------------------

	/*!
	* \struct Int2
	* \brief Vector of 2 int items
	*/
	class Int2{
	private:

		friend class boost::serialization::access;

		template< class Archive >
		void serialize( Archive & ar, const unsigned int version )
		{
			ar & x_;
			ar & y_;
		}

	public:

		int x_;
		int y_;
	};

	/*!
	* \fn Int2 operator+(Int2 const& a, Int2 const& b)
	* \brief Adds two Int2 vectors
	* \param a : first vector
	* \param b : second vector
	* \return the sum of the two vectors (a + b)
	*/
	inline Int2 operator+(Int2 const& a, Int2 const& b){return {a.x_+b.x_, a.y_+b.y_ };}

	/*!
	* \fn Int2 operator-(Int2 const& a, Int2 const& b)
	* \brief Substract two vectors
	* \param a : first vector
	* \param b : second vector
	* \return a - b
	*/
	inline Int2 operator-(Int2 const& a, Int2 const& b){return {a.x_-b.x_, a.y_-b.y_};}

	/*!
	* \fn int operator*(Int2 const& a, Int2 const& b)
	* \brief Scalar product of two vectors
	* \param a : first vector
	* \param b : second vector
	* \return scalar product (element by element) (a.x * b.x + a.y * b.y)
	*/
	inline int operator*(Int2 const& a, Int2 const& b){return a.x_*b.x_ + a.y_*b.y_;}

	/*!
	* \fn Int2 operator*(int const& a, Int2 const& b)
	* \brief Scalar*vector product
	* \param a : an integer
	* \param b : a vector
	* \return a*b
	*/
	inline Int2 operator*(int const& a, Int2 const& b){return {a*b.x_, a*b.y_};}

	/*!
	* \fn Int2 operator*(Int2 const& b, int const& a)
	* \brief Vector*scalar product
	* \param b : a vector
	* \param a : an integer
	* \return a*b
	*/
	inline Int2 operator*(Int2 const& b, int const& a){return {a*b.x_, a*b.y_};}

	/*!
	* \fn Int2 operator/(Int2 const& b, int const& a)
	* \brief Vector/scalar product
	* \param b : a vector
	* \param a : an integer
	* \return b/a
	*/
	inline Int2 operator/(Int2 const& b, int const& a){return {b.x_/a, b.y_/a};}


	/*!
	* \class Uint2
	* \brief Vector of 2 Uint items
	*/
	class Uint2{
	private:

		friend class boost::serialization::access;

		template<class Archive>
		void serialize(Archive & ar, const unsigned int version)
		{
			ar & x_;
			ar & y_;
		}

	public:

		Uint x_;
		Uint y_;
	};

	/*!
	* \fn Uint2 operator+(Uint2 const& a, Uint2 const& b)
	* \brief Adds two vectors
	*/
	inline Uint2 operator+(Uint2 const& a, Uint2 const& b){return {a.x_+b.x_, a.y_+b.y_};}

	/*!
	* \fn Uint2 operator-(Uint2 const& a, Uint2 const& b)
	* \brief Substract two vectors
	*/
	inline Uint2 operator-(Uint2 const& a, Uint2 const& b){return {a.x_-b.x_, a.y_-b.y_};}

	/*!
	* \fn Uint2 operator*(Uint2 const& a, Uint2 const& b)
	* \brief Vector product
	*/
	//inline Uint2 operator^(Uint2 const& a, Uint2 const& b){return {a.y_*b.z_ - a.z_*b.y_, a.z_*b.x_ - a.x_*b.z_, a.x_*b.y_ - a.y_*b.x_};}

	/*!
	* \fn Uint operator*(Uint2 const& a, Uint2 const& b)
	* \brief Scalar product
	*/
	inline Uint operator*(Uint2 const& a, Uint2 const& b){return a.x_*b.x_ + a.y_*b.y_;}

	/*!
	* \fn Uint2 operator*(Uint const& a, Uint2 const& b)
	* \brief Scalar*vector product
	*/
	inline Uint2 operator*(Uint const& a, Uint2 const& b){return {a*b.x_, a*b.y_};}

	/*!
	* \fn Uint2 operator*(Uint2 const& b, Uint const& a)
	* \brief Vector*scalar product
	*/
	inline Uint2 operator*(Uint2 const& b, Uint const& a){return {a*b.x_, a*b.y_};}

	/*!
	* \fn Uint2 operator/(Uint2 const& b, Uint const& a)
	* \brief Vector/scalar product
	*/
	inline Uint2 operator/(Uint2 const& b, Uint const& a){return {b.x_/a, b.y_/a};}




	/*!
	* \class Double2
	* \brief Vector of 2 double items
	*/
	class Double2{
	private:

		friend class boost::serialization::access;

		template<class Archive>
		void serialize(Archive & ar, const unsigned int version)
		{
			ar & x_;
			ar & y_;
		}

	public:

		double x_;
		double y_;
	};

	/*!
	* \fn Double2 operator+(Double2 const& a, Double2 const& b)
	* \brief Adds two vectors
	*/
	inline Double2 operator+(Double2 const& a, Double2 const& b){return {a.x_+b.x_, a.y_+b.y_};}

	/*!
	* \fn Double2 operator-(Double2 const& a, Double2 const& b)
	* \brief Substract two vectors
	*/
	inline Double2 operator-(Double2 const& a, Double2 const& b){return {a.x_-b.x_, a.y_-b.y_};}

	/*!
	* \fn Double2 operator*(Double2 const& a, Double2 const& b)
	* \brief Vector product
	*/
	//inline Double2 operator^(Double2 const& a, Double2 const& b){return {a.y_*b.z_ - a.z_*b.y_, a.z_*b.x_ - a.x_*b.z_, a.x_*b.y_ - a.y_*b.x_};}

	/*!
	* \fn double operator*(Double2 const& a, Double2 const& b)
	* \brief Scalar product
	*/
	inline double operator*(Double2 const& a, Double2 const& b){return a.x_*b.x_ + a.y_*b.y_;}

	/*!
	* \fn Double2 operator*(double const& a, Double2 const& b)
	* \brief Scalar*vector product
	*/
	inline Double2 operator*(double const& a, Double2 const& b){return {a*b.x_, a*b.y_};}

	/*!
	* \fn Double2 operator*(Double2 const& b, double const& a)
	* \brief Vector*scalar product
	*/
	inline Double2 operator*(Double2 const& b, double const& a){return {a*b.x_, a*b.y_};}

	/*!
	* \fn Double2 operator/(Double2 const& b, double const& a)
	* \brief Vector/scalar product
	*/
	inline Double2 operator/(Double2 const& b, double const& a){return {b.x_/a, b.y_/a};}



	/*!
	* \class Int3
	* \brief Vector of 3 int items
	*/
	class Int3{
	private:

		friend class boost::serialization::access;

		template<class Archive>
		void serialize(Archive & ar, const unsigned int version)
		{
			ar & x_;
			ar & y_;
			ar & z_;
		}

	public:

		int x_;
		int y_;
		int z_;
	};

	/*!
	* \fn Int3 operator+(Int3 const& a, Int3 const& b)
	* \brief Adds two vectors
	*/
	inline Int3 operator+(Int3 const& a, Int3 const& b){return {a.x_+b.x_, a.y_+b.y_, a.z_+b.z_};}

	/*!
	* \fn Int3 operator-(Int3 const& a, Int3 const& b)
	* \brief Substract two vectors
	*/
	inline Int3 operator-(Int3 const& a, Int3 const& b){return {a.x_-b.x_, a.y_-b.y_, a.z_-b.z_};}

	/*!
	* \fn Int3 operator*(Int3 const& a, Int3 const& b)
	* \brief Vector product
	*/
	inline Int3 operator^(Int3 const& a, Int3 const& b){return {a.y_*b.z_ - a.z_*b.y_, a.z_*b.x_ - a.x_*b.z_, a.x_*b.y_ - a.y_*b.x_};}

	/*!
	* \fn int operator*(Int3 const& a, Int3 const& b)
	* \brief Scalar product
	*/
	inline int operator*(Int3 const& a, Int3 const& b){return a.x_*b.x_ + a.y_*b.y_ + a.z_*b.z_;}

	/*!
	* \fn Int3 operator*(int const& a, Int3 const& b)
	* \brief Scalar*vector product
	*/
	inline Int3 operator*(int const& a, Int3 const& b){return {a*b.x_, a*b.y_, a*b.z_};}

	/*!
	* \fn Int3 operator*(Int3 const& b, int const& a)
	* \brief Vector*scalar product
	*/
	inline Int3 operator*(Int3 const& b, int const& a){return {a*b.x_, a*b.y_, a*b.z_};}

	/*!
	* \fn Int3 operator/(Int3 const& b, int const& a)
	* \brief Vector/scalar product
	*/
	inline Int3 operator/(Int3 const& b, int const& a){return {b.x_/a, b.y_/a, b.z_/a};}


	/*!
	* \class Uint3
	* \brief Vector of 3 Uint items
	*/
	class Uint3{
	private:

		friend class boost::serialization::access;

		template<class Archive>
		void serialize(Archive & ar, const unsigned int version)
		{
			ar & x_;
			ar & y_;
			ar & z_;
		}

	public:

		Uint x_;
		Uint y_;
		Uint z_;
	};

	/*!
	* \fn Uint3 operator+(Uint3 const& a, Uint3 const& b)
	* \brief Adds two vectors
	*/
	inline Uint3 operator+(Uint3 const& a, Uint3 const& b){return {a.x_+b.x_, a.y_+b.y_, a.z_+b.z_};}

	/*!
	* \fn Uint3 operator-(Uint3 const& a, Uint3 const& b)
	* \brief Substract two vectors
	*/
	inline Uint3 operator-(Uint3 const& a, Uint3 const& b){return {a.x_-b.x_, a.y_-b.y_, a.z_-b.z_};}

	/*!
	* \fn Uint3 operator*(Uint3 const& a, Uint3 const& b)
	* \brief Vector product
	*/
	inline Uint3 operator^(Uint3 const& a, Uint3 const& b){return {a.y_*b.z_ - a.z_*b.y_, a.z_*b.x_ - a.x_*b.z_, a.x_*b.y_ - a.y_*b.x_};}

	/*!
	* \fn Uint operator*(Uint3 const& a, Uint3 const& b)
	* \brief Scalar product
	*/
	inline Uint operator*(Uint3 const& a, Uint3 const& b){return a.x_*b.x_ + a.y_*b.y_ + a.z_*b.z_;}

	/*!
	* \fn Uint3 operator*(Uint const& a, Uint3 const& b)
	* \brief Scalar*vector product
	*/
	inline Uint3 operator*(Uint const& a, Uint3 const& b){return {a*b.x_, a*b.y_, a*b.z_};}

	/*!
	* \fn Uint3 operator*(Uint3 const& b, Uint const& a)
	* \brief Vector*scalar product
	*/
	inline Uint3 operator*(Uint3 const& b, Uint const& a){return {a*b.x_, a*b.y_, a*b.z_};}

	/*!
	* \fn Uint3 operator/(Uint3 const& b, Uint const& a)
	* \brief Vector/scalar product
	*/
	inline Uint3 operator/(Uint3 const& b, Uint const& a){return {b.x_/a, b.y_/a, b.z_/a};}


	/*!
	* \class Double3
	* \brief Vector of 3 double items
	*/
	class Double3{
	private:

		friend class boost::serialization::access;

		template<class Archive>
		void serialize(Archive & ar, const unsigned int version)
		{
			ar & x_;
			ar & y_;
			ar & z_;
		}

	public:

		double x_;
		double y_;
		double z_;
	};

	/*!
	* \fn Double3 operator+(Double3 const& a, Double3 const& b)
	* \brief Adds two vectors
	*/
	inline Double3 operator+(Double3 const& a, Double3 const& b){return {a.x_+b.x_, a.y_+b.y_, a.z_+b.z_};}

	/*!
	* \fn Double3 operator-(Double3 const& a, Double3 const& b)
	* \brief Substract two vectors
	*/
	inline Double3 operator-(Double3 const& a, Double3 const& b){return {a.x_-b.x_, a.y_-b.y_, a.z_-b.z_};}

	/*!
	* \fn Double3 operator*(Double3 const& a, Double3 const& b)
	* \brief Vector product
	*/
	inline Double3 operator^(Double3 const& a, Double3 const& b){return {a.y_*b.z_ - a.z_*b.y_, a.z_*b.x_ - a.x_*b.z_, a.x_*b.y_ - a.y_*b.x_};}

	/*!
	* \fn double operator*(Double3 const& a, Double3 const& b)
	* \brief Scalar product
	*/
	inline double operator*(Double3 const& a, Double3 const& b){return a.x_*b.x_ + a.y_*b.y_ + a.z_*b.z_;}

	/*!
	* \fn Double3 operator*(double const& a, Double3 const& b)
	* \brief Scalar*vector product
	*/
	inline Double3 operator*(double const& a, Double3 const& b){return {a*b.x_, a*b.y_, a*b.z_};}

	/*!
	* \fn Double3 operator*(Double3 const& b, double const& a)
	* \brief Vector*scalar product
	*/
	inline Double3 operator*(Double3 const& b, double const& a){return {a*b.x_, a*b.y_, a*b.z_};}

	/*!
	* \fn Double3 operator/(Double3 const& b, double const& a)
	* \brief Vector/scalar product
	*/
	inline Double3 operator/(Double3 const& b, double const& a){return {b.x_/a, b.y_/a, b.z_/a};}


	/*!
	* \enum ESpeedRepresentation
	* \brief used to specifie the type of the speed
	*/
	enum ESpeedType {
		/*! \brief eulerian speed (aplies to a whole site) computed once at the begining of the simulation ( a static speed must be time independant ! ) */
		EULERIAN_STATIC,
		/*! \brief eulerian speed (aplies to a whole site) computed at each step for each site comtaining a particle */
		EULERIAN_DYNAMIC,
		/*! \brief lagrangian speed (aplies to a unique particle, has to be re-calculated at each step for each particle) */
		LAGRANGIAN_DYNAMIC,
		/* \brief undefined */
		//UNDEFINED
	};

	/*!
	* \enum ESimulatorType
	* \brief used to specifie the simulator type
	*/
	enum ESimulatorType {
		/*! \brief simulator based on cellular automata */
		CA,
		/* \brief simulator that keep exact position of particles (lagrangian on grid model) */
		EXACT,
		/* \brief undefined simulator type */
		UNDEFINED_SIM
	};


	/*!
	* \class Particle
	* \brief A Particle is represented by the id from which its belongs and its "displacement_" regarding the site to which it is linked.
	* Moreover, we add a "lagrangian" speed to each Particle. The domain holds the advecting field and we can add a particular speed to each Particle
	*/
	class Particle{
	protected:

		/*! \todo : possible d'optimiser la serialization ? voir tuto boost mpi ? */

		friend class boost::serialization::access;

		template<class Archive>
		void serialize(Archive & ar, const unsigned int version)
		{
			ar & particleId_;
			ar & familyId_;
			ar & displacement_;
			ar & lagrangianSpeed_;
			ar & diffusionSpeed_;
            ar & scaling_;
		}

		static long long unsigned int getNextParticleId(){
			static long long unsigned int currentParticleId = 0;
			return ++currentParticleId;
		}

	public:

		/*! \fn Particle( int familyId, Double3 disp, Double3 speed )
		* \brief constructor of the struct Particle
		* \param familyId : the id of the corresponding particle family
		* \param disp : the initial displacement (regarding the current site)
		* \param speed : the initial (lagrangian) speed of the particle
		*/

		Particle( int familyId, Double3 disp, Double3 speed ):
		familyId_			(familyId),
		displacement_		(disp),
		lagrangianSpeed_	(speed),
        diffusionSpeed_  ({0.0, 0.0, 0.0}),
        scaling_         (1.0){ particleId_ = getNextParticleId(); }

		Particle( int familyId, Double3 disp, Double3 speed, Double3 diffusionSpeed ):
		familyId_			(familyId),
		displacement_		(disp),
		lagrangianSpeed_	(speed),
        diffusionSpeed_  (diffusionSpeed),
        scaling_         (1.0){ particleId_ = getNextParticleId(); }

        Particle( int familyId, Double3 disp, Double3 speed, double scaling ):
		familyId_			(familyId),
        displacement_		(disp),
        lagrangianSpeed_	(speed),
        diffusionSpeed_  ({0.0, 0.0, 0.0}),
        scaling_         (scaling){ particleId_ = getNextParticleId(); }

        Particle( int familyId, Double3 disp, Double3 speed, Double3 diffusionSpeed, double scaling ):
		familyId_			(familyId),
        displacement_		(disp),
        lagrangianSpeed_	(speed),
        diffusionSpeed_  (diffusionSpeed),
        scaling_         (scaling){ particleId_ = getNextParticleId(); }


		/*! \fn Particle(  )
		* \brief constructor of the struct Particle
		*/
		Particle( ):
		familyId_		( 0 ),
		displacement_	( { 0.0, 0.0, 0.0 } ),
		lagrangianSpeed_	( { 0.0, 0.0, 0.0 } ),
        diffusionSpeed_  ({0.0, 0.0, 0.0}),
        scaling_         (1.0){ particleId_ = getNextParticleId(); }

		long long unsigned int particleId_;

		/*!
		* \brief the id of the corresponding particle family
		*/
		int familyId_;

		/*!
		* \brief displacement (regarding the current site)
		*/
		Double3 displacement_;

		/*!
		* \brief (lagrangian) speed of the particle
		*/
		Double3 lagrangianSpeed_;

		/*!
		* \brief to store diffusion speed of the particle
		* !!! the diffusion speed is already included in the lagrangian speed !!!
		*/
		Double3 diffusionSpeed_;

        double scaling_;

		long long unsigned int getParticleId(){ return particleId_; }

		int getFamilyId(){ return familyId_; }

		Double3 getDisplacement(){ return displacement_; }

		Double3 getLagrangianSpeed(){ return lagrangianSpeed_; }

		Double3 getDiffusionSpeed(){ return diffusionSpeed_; }

		void setFamilyId( int fid ){ familyId_ = fid; }

		void setDisplacement( Double3 disp ){ displacement_ = disp; }

		void setLagrangianSpeed( Double3 ls ){ lagrangianSpeed_ = ls; }

		void setDiffusionSpeed( Double3 ds ){ diffusionSpeed_ = ds; }

	};

	class BoundaryParticle: public Particle{

	protected:

		friend class boost::serialization::access;

		template<class Archive>
		void serialize(Archive & ar, const unsigned int version)
		{
			ar & particleId_;
			ar & familyId_;
			ar & displacement_;
			ar & lagrangianSpeed_;
			ar & timeStamp_;
		}

		double timeStamp_;

	public:

		BoundaryParticle( Particle part, double timeStamp ):
		Particle( part.familyId_, part.displacement_, part.lagrangianSpeed_, part.diffusionSpeed_ ),
		timeStamp_( timeStamp ){}

		BoundaryParticle( ):
		Particle( 0, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0} ),
		timeStamp_ ( 0.0 ){}

		double getTimeStamp(){ return timeStamp_; }

		void setTimeStamp( double ts ){ timeStamp_ = ts; }

	};

} // namespace piaf

#endif /* SIMULATORTYPES_H_ */
