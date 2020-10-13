#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <string>
#include <cmath>
#include <sstream>
#include <random>
#include <list>
#include <sys/stat.h>
#include <sys/types.h>
#include <assert.h>
#include <math.h>

using namespace std;

std::vector<double> MSD;
std::vector<double> t;

struct Point_
{
    double x , y , z ;
    double x0, y0, z0;
};

inline double logistic(double x, double loc, double scale)
{
    return loc + scale*log(x/(1.-x));
}

inline double pareatoValue(double r, double xmin, double alpha)
{
    double x = log(xmin) - (1./alpha)*log(1.-r);
    return exp(x);
}

inline void periodicity(double& pos, double domainSize)
{
    if (pos > domainSize)
    {
        pos = pos - (double)((int)(std::abs(pos) / domainSize)) * domainSize;
    }
    else if (pos < 0.)
    {
        pos = pos + (double)((int)(std::abs(pos) / domainSize) + 1) * domainSize;
    }
}

inline void upper_bounce_back(double& pos, double domainSize)
{
    if (pos > domainSize)
    {
        double jump = pos / domainSize;
        if (jump > 2.)
        {
            // Add it in the boundary layer.
            pos = -1.;
            return;
        }
            
        pos = domainSize - (jump - 1.) * domainSize;
    }
}

inline bool boundary_layer(double x, double y, double z, double lx, double ly,
                           int** candidate_deposition, int nx, int ny)
{
    if (z > 0.)
        return false;
    
    int i = (int)((x/lx)*nx);
    int j = (int)((y/ly)*ny);

    // account for more particles at the same position (stacking)
    candidate_deposition[i][j] += 1;
    return true;
}

class flight
{
    public:
        flight(double lx_, double ly_, double lz_,
               size_t num_particles_,
               int nx_, int ny_,
               double shear_rate_x_, double shear_rate_y_,
               double dt_,
               double v_xy_, double v_z_,
               double alpha_xy_, double alpha_z_,
               std::default_random_engine* gen_):
            lx(lx_), ly(ly_), lz(lz_),
            num_particles(num_particles_),
            nx(nx_), ny(ny_),
            shear_rate_x(shear_rate_x_), shear_rate_y(shear_rate_y_),
            dt(dt_),
            v_xy(v_xy_), v_z(v_z_),
            alpha_xy(alpha_xy_), alpha_z(alpha_z_),
            gen(gen_)
        {
            particles = std::list<Point_>(num_particles);

            std::uniform_real_distribution<double> distribution_x(0.0, lx);
            std::uniform_real_distribution<double> distribution_y(0.0, ly);
            std::uniform_real_distribution<double> distribution_z(0.0, lz);

            for (it = particles.begin(); it != particles.end(); ++it)
            {
                it->x = distribution_x(*gen);
                it->y = distribution_y(*gen);
                it->z = distribution_z(*gen);

                it->x0 = it->x;
                it->y0 = it->y;
                it->z0 = it->z;
            }

            // Gaussian Random Walk
            std_gaussian = std::normal_distribution<double>(0.0, 1.0);

            // Uniform distribution for the direction on the xy-plane
            distribution_theta = std::uniform_real_distribution<double>(0.0, 2.*std::acos((double)-1));

            // Up - down & Pareto
            rng = std::uniform_real_distribution<double>(0.0, 1.0);

            // Re-shuffling etc.
            candidate_deposition_tmp = new int *[nx];
            for (int i = 0; i < nx; i++)
            {
                candidate_deposition_tmp[i] = new int[ny];
                for (int j = 0; j < ny; j++)
                    candidate_deposition_tmp[i][j] = 0;
            }
        }

        ~flight()
        {
            for (int i = 0; i < nx; i++)
                delete[] candidate_deposition_tmp[i];
            
            delete[] candidate_deposition_tmp;
        }

        void advance(int** candidate_deposition)
        {            
            for (it = particles.begin(); it != particles.end(); ++it)
            {
                // Move up or down
                direction_z = (2 * (int)(rng(*gen) + 0.5) - 1);

                // direction on the xy-plane
                theta = distribution_theta(*gen);

#ifdef PARETO
                velocity = pareatoValue(rng(*gen), v_xy, alpha_xy);
#endif
#ifdef GAUSSIAN
                velocity = v_xy * std::abs(std_gaussian(*gen));
#endif

                it->x += (velocity * cos(theta) + shear_rate_x * it->z) * dt;
                it->y += (velocity * sin(theta) + shear_rate_y * it->z) * dt;
                periodicity(it->x, lx);
                periodicity(it->y, ly);

#ifdef PARETO
                velocity = pareatoValue(rng(*gen), v_z, alpha_z);
#endif
#ifdef GAUSSIAN
                //velocity = v_z * std::abs(std_gaussian(*gen));
                velocity = std::abs(logistic(rng(*gen), 2.1896808363537395e-06, 1.e-4));
#endif
                
                it->z += (double)(direction_z) * velocity * dt;
                upper_bounce_back(it->z, lz);
                if (boundary_layer(it->x, it->y, it->z, lx, ly, candidate_deposition, nx, ny))
                {
                    it = particles.erase(it);
                    --it;
                }
            }
        }
        
        void reinject(int** candidate_deposition)
        {
            // Do nothing
            // If in Boundary Layer then deposit it
            // Interested on min val of D
        }

        double computeNoPart()
        {
            return (double)(particles.size());
        }

    public:
        std::list<Point_> particles;
        std::list<Point_>::iterator it;

        double lx, ly, lz;
        size_t num_particles;
        int nx, ny;
        double shear_rate_x, shear_rate_y;
        double dt;
        
        double v_xy, v_z;
        double alpha_xy, alpha_z;

        std::default_random_engine* gen;

        std::normal_distribution<double> std_gaussian;

        double theta;
        std::uniform_real_distribution<double> distribution_theta;

        int direction_z;
        std::uniform_real_distribution<double> rng;
        
        double velocity;

        int** candidate_deposition_tmp;
};


double slope(const vector<double>& x, const vector<double>& y)
{
    if(x.size() != y.size()){
        exit(EXIT_FAILURE);
    }
    double n = x.size();

    double avgX = accumulate(x.begin(), x.end(), 0.0) / n;
    double avgY = accumulate(y.begin(), y.end(), 0.0) / n;

    double numerator = 0.0;
    double denominator = 0.0;

    for(int i=0; i<n; ++i){
        numerator += (x[i] - avgX) * (y[i] - avgY);
        denominator += (x[i] - avgX) * (x[i] - avgX);
    }

    if(denominator == 0){
        exit(EXIT_FAILURE);
    }

    return numerator / denominator;
}


namespace Maths
{
namespace Regression
{
//! Given a set of points, this class calculates the linear regression parameters and evaluates the regression line at arbitrary abscissas.
class Linear
{
public:
 
    //! Class constructor
    Linear(int n, double *x, double *y)
    {

        // calculate the averages of arrays x and y
        double xa = 0, ya = 0;
        for (int i = 0; i < n; i++) {
            xa += x[i];
            ya += y[i];
        }
        xa /= n;
        ya /= n;

        // calculate auxiliary sums
        double xx = 0, yy = 0, xy = 0;
        for (int i = 0; i < n; i++) {
            double tmpx = x[i] - xa, tmpy = y[i] - ya;
            xx += tmpx * tmpx;
            yy += tmpy * tmpy;
            xy += tmpx * tmpy;
        }

        // calculate regression line parameters

        // make sure slope is not infinite
        assert(fabs(xx) != 0);

            m_b = xy / xx;
            m_a = ya - m_b * xa;
        m_coeff = (fabs(yy) == 0) ? 1 : xy / sqrt(xx * yy);
    }
 
    //! Evaluates the linear regression function at the given abscissa.
    double getValue(double x)
    {
        return m_a + m_b * x;
    }
 
    //! Returns the slope of the regression line
    double getSlope()
    {
    return m_b;
    }
 
    //! Returns the intercept on the Y axis of the regression line
    double getIntercept()
    {
    return m_a;
    }
 
    //! Returns the linear regression coefficient
    double getCoefficient()
    {
    return m_coeff;
    }
 
private:
    double m_a, m_b, m_coeff;
}; 
}
}


void compute_MSD(flight& AP, flight& NAP)
{
    double MSD_ = 0;

    std::list<Point_>::iterator it;
    for (it = AP.particles.begin(); it != AP.particles.end(); ++it)
        MSD_ += (it->z - it->z0)*(it->z - it->z0);
    for (it = NAP.particles.begin(); it != NAP.particles.end(); ++it)
        MSD_ += (it->z - it->z0)*(it->z - it->z0);

    MSD_ /= (AP.computeNoPart() + NAP.computeNoPart());
    MSD.push_back(MSD_);
}


int main()
{
    // initial conditions for Activated PLTs (AP) & nonActivated PLTs (NAP)
    // refer to ul (micro-liter)
    int noAP   = 0;
    int noNAP  = 172200;

    // For the distributions of the particles
    // Same for both AP & NAP
    double v_xy = 0.0;
    double v_z = 0.0;
    double alpha_xy = 0.0;
    double alpha_z = 0.0;

    //-----------

    std::random_device rd;
    std::mt19937 rng(rd());
    std::default_random_engine generator(rng());

    double Lz = 1.e-3; // thickness of the blood layer
    double dt = 1.e-2;  // time step
    int nx = 447;
    int ny = 447;

    flight activatePlatelets(1.e-3, 1.e-3, Lz, // double lx_, double ly_, double lz_
                            noAP, // size_t num_particles_
                            nx, ny, // int nx_, int ny_
                            100., 0., // double shear_rate_x_, double shear_rate_y_
                            dt, // double dt_
                            v_xy, v_z,
                            alpha_xy, alpha_z,
                            &generator);
    
    int** candidate_deposit_AP;
    // allocation memory for the arrays
    candidate_deposit_AP = new int *[nx];
    for (int i = 0; i < nx; i++)
    {
        candidate_deposit_AP[i] = new int[ny];
        for (int j = 0; j < ny; j++)
            candidate_deposit_AP[i][j] = 0;
    }

    //d1q3 nonActivatePlatelets(nz, rhoNonActivated, D, dt, dz, W);
    flight nonActivatePlatelets(1.e-3, 1.e-3, Lz, // double lx_, double ly_, double lz_
                                noNAP, // size_t num_particles_
                                nx, ny, // int nx_, int ny_
                                100., 0., // double shear_rate_x_, double shear_rate_y_
                                dt, // double dt_
                                v_xy, v_z,
                                alpha_xy, alpha_z,
                                &generator);

    int** candidate_deposit_NAP;
    // allocation memory for the arrays
    candidate_deposit_NAP = new int *[nx];
    for (int i = 0; i < nx; i++)
    {
        candidate_deposit_NAP[i] = new int[ny];
        for (int j = 0; j < ny; j++)
            candidate_deposit_NAP[i][j] = 0;
    }
    

    double tMax = 20; // the simulation will run for thant many seconds
    int nbIter = tMax / dt; // the number of iterations

    int iter = 0;
    double t_ = 0.;
    t.push_back(0.);
    MSD.push_back(0.);

    while (iter < nbIter)
    {
        activatePlatelets.advance(candidate_deposit_AP);
        nonActivatePlatelets.advance(candidate_deposit_NAP);

        t_ += dt;
        t.push_back(t_);
        compute_MSD(activatePlatelets, nonActivatePlatelets);
        
        iter++;
    }

    Maths::Regression::Linear A(MSD.size(), t.data(), MSD.data());

    std::cout << "        Diffsion Coefficient (eqn) : " << (v_z*v_z) * dt * 0.5 << std::endl;
    //std::cout << "        Diffsion Coefficient (MSD) : " << slope(t, MSD) * 0.5 << std::endl;
    std::cout << "        Diffsion Coefficient (MSD) : " << A.getSlope() * 0.5 << std::endl;
    std::cout << "max num of PLTs deposited (per ul) : " << ((noAP  - activatePlatelets.computeNoPart()) +
                                                            (noNAP - nonActivatePlatelets.computeNoPart()))
                                                         << std::endl;

    std::string filename = "MSD.csv";
    std::ofstream MSD_file(filename.c_str());
    for (size_t i = 0; i < MSD.size(); ++i)
        MSD_file << t[i] << "," << MSD[i] << std::endl;
    MSD_file.close();

    return 0;
}