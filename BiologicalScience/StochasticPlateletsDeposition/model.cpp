#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <cmath>
#include <sstream>
#include <random>
#include <list>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

using namespace std;

// eqn 4.4 - slightly modified
inline double P(double p, double maxAlbumine, double albumine)
{
    if (albumine >= maxAlbumine)
        return 0.;
    else
        return p * (1.0 - albumine / maxAlbumine);
}

// eqn 4.5 & 4.6 - slightly modified
inline double QR(double p, double attenuation, double maxAlbumine, double albumine)
{
    if (albumine >= maxAlbumine)
        return 0.;
    else
        return p * exp(-attenuation * (albumine / maxAlbumine));
}

// pgm image
void writePgm(int** d, int nx, int ny, int imax, std::string filename)
{
    std::ofstream file;
    file.open (filename);
    file << "P2" << std::endl;
    file << nx << " " << ny << std::endl;
    file << imax << std::endl;
    for(int j=0; j<ny; j++){
        for(int i=0; i<nx; i++){
            file << d[i][j] << " ";
        }
        file << std::endl;
    }
    file.close();
}

// pgm image : inverse version
void writePgm_inv(int** d, int nx, int ny, int imax, std::string filename)
{
    std::ofstream file;
    file.open (filename);
    file << "P2" << std::endl;
    file << nx << " " << ny << std::endl;
    file << imax << std::endl;
    for(int j=0; j<ny; j++){
        for(int i=0; i<nx; i++){
            file << (imax - d[i][j]) << " ";
        }
        file << std::endl;
    }
    file.close();
}

struct Point_
{
    double x , y , z ;
    // For MSD
    //double x0, y0, z0;
};

// Inverse of the pareto CDF
inline double pareato_(double r, double xmin, double alpha)
{
    double x = log(xmin) - (1./alpha)*log(1.-r);
    return exp(x);
}

// Inverse of the logistic CDF
inline double logistic_(double x, double loc, double scale)
{
    return loc + scale*log(x/(1.-x));
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

            // for MSD calculation
            //it->x0 = it->x;
            //it->y0 = it->y;
            //it->z0 = it->z;
        }

        // Gaussian Random Walk
        std_gaussian = std::normal_distribution<double>(0.0, 1.0);

        // Uniform distribution for the direction on the xy-plane
        distribution_theta = std::uniform_real_distribution<double>(0.0, 2.*std::acos((double)-1));

        // Up or down, Pareto, logistic
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

            // same for z & xy (just an approximation for now)
            velocity = v_z * std::abs(std_gaussian(*gen));

            it->x += (velocity * cos(theta) + shear_rate_x * it->z) * dt;
            it->y += (velocity * sin(theta) + shear_rate_y * it->z) * dt;
            periodicity(it->x, lx);
            periodicity(it->y, ly);
            
            it->z += (double)(direction_z) * velocity * dt;
            upper_bounce_back(it->z, lz);
            if (boundary_layer(it->x, it->y, it->z, lx, ly, candidate_deposition, nx, ny))
            {
                it = particles.erase(it);
                --it;
            }
        }
    }
    
    void reshuffleSubstrate(int** candidate_deposition)
    {
        // No need for re-injection since the PLTs are trapped in the CFL
        // At least the majority of them

        int size_;
        double x, y;

        for (int i = 0; i < nx; i++)
            for (int j = 0; j < ny; j++)
                candidate_deposition_tmp[i][j] = candidate_deposition[i][j];

        // Transport/ Reshuffling inside the Boundary Layer
        int i_, j_;
        for (int i = 0; i < nx; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                size_ = candidate_deposition_tmp[i][j];
                for (int k = 0; k < size_; ++k)
                {
                    x = (((double)(i) + 0.5)/(double)(nx))*lx;
                    y = (((double)(j) + 0.5)/(double)(ny))*ly;

                    theta = distribution_theta(*gen);
                    
                    velocity = v_z * std::abs(std_gaussian(*gen));
                
                    x += (velocity * cos(theta)) * dt;
                    y += (velocity * sin(theta)) * dt;
                    periodicity(x, lx);
                    periodicity(y, ly);

                    i_ = (int)((x/lx)*nx);
                    j_ = (int)((y/ly)*ny);

                    candidate_deposition[i ][j ] -= 1;
                    candidate_deposition[i_][j_] += 1;
                }
            }
        }
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


class Substrate
{
public:
    int nx, ny;

    std::default_random_engine* gen;
    std::uniform_real_distribution<double> rng;

    int *sizeOfCluster;
    int *volumeOfCluster;
    int *maxHeightOfCluster;
    int **deposit; // refers to PLTs on the Surface S
    double **avail; // array to hold if site is available for adhesion
                    // refers to albumin
    int **label; // label visited clusters for counting purposes

    int noCluster;
    int noColl;
    int noPart;         // total no of particle created in space
    int aggregateIndex; // to enumerate the clusters
    double averageSize;
    double averageVol;
    double averageHeight;

    Substrate(int nnx, int nny, std::default_random_engine* gen_) :
        nx(nnx),
        ny(nny),
        gen(gen_)
    {
        rng = std::uniform_real_distribution<double>(0, 1);

        sizeOfCluster = new int[nx * ny]; //  max num of clusters
        volumeOfCluster = new int[nx * ny];
        maxHeightOfCluster = new int[nx * ny];
        
        // allocation memory for the arrays
        deposit = new int *[nx];
        for (int i = 0; i < nx; i++)
            deposit[i] = new int[ny];

        avail = new double *[nx];
        for (int i = 0; i < nx; i++)
            avail[i] = new double[ny];

        label = new int *[nx];
        for (int i = 0; i < nx; i++)
            label[i] = new int[ny];

        // initialized variables to zeros:
        for (int i = 0; i < nx; i++)
            for (int j = 0; j < ny; j++)
            {
                deposit[i][j] = 0;
                avail[i][j] = 0.;
                label[i][j] = 0;
            }
        
        noPart = 0;
        noCluster = 0;
    } // end constructor

    ~Substrate()
    {
        for (int i = 0; i < nx; i++)
        {
            delete[] deposit[i];
            delete[] avail[i];
            delete[] label[i];
        }

        delete[] deposit;
        delete[] avail;
        delete[] label;
        delete[] sizeOfCluster;
        delete[] volumeOfCluster;
        delete[] maxHeightOfCluster;
    }

    int countDeposit()
    {
        int b = 0;
        for (int i = 0; i < nx; i++)
            for (int j = 0; j < ny; j++)
                b += deposit[i][j];

        return b;
    }

    double countDepositedAlbumine()
    {
        double b = 0;
        for (int i = 0; i < nx; i++)
            for (int j = 0; j < ny; j++)
                b += avail[i][j];
        
        return b;
    }

    // Counts the no of clusters, label each pixel of each cluster with a
    // a number going from 1 to no of cluster.
    int countAggregate()
    {
        int ii, jj;
        
        int *myStack;
        myStack = new int[2 * nx * ny]; // pointer to hold pixels (i,j) that still need
                                        // to be explored for a continuation of the cluster
        int top = 0;                    // position in stack

        for (int i = 0; i < nx; i++)
            for (int j = 0; j < ny; j++)
                label[i][j] = 0;

        aggregateIndex = 0; // to enumerate the clusters

        for (int i = 0; i < nx; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                if (deposit[i][j] == 0)
                    continue; // this is not a cluster
                if (label[i][j] != 0)
                    continue; // this is an already visited cluster

                myStack[top] = i;
                top++;
                myStack[top] = j;
                top++;
                
                // Get the same label as your neighbors, if any.
                // Else get a new label given by newIndex
                // Set a unique number to all pixels belonging to the same cluster
                aggregateIndex++;
                
                int clusterSize = 0;
                int clusterVolume = 0;
                int maxHeight = 0.;

                while (top > 0)
                {
                    top--;
                    jj = myStack[top];
                    top--;
                    ii = myStack[top];
                    
                    label[ii][jj] = aggregateIndex;
                    
                    clusterSize++;
                    clusterVolume += deposit[ii][jj];
                    if (deposit[ii][jj] > maxHeight)
                        maxHeight = deposit[ii][jj];

                    //  j,y
                    //  ^
                    //  |
                    //  |
                    //  o------->i,x
                    //
                    // visit all 8 neighbors of (ii,jj) to see if the cluster continues
                    //
                    // i-1,j+1    i,j+1    i+1,j+1
                    //             ^
                    //             | 
                    // i-1,j <--- i,j ---> i+1,j
                    //             |
                    //             v
                    // i-1,j-1    i,j-1    i+1,j-1

                    // Use the modulo technique to take proper care of the boundary nodes
                    // For example to get the i+1 neighbors:
                    // Do (i+1+nx)%nx
                    // 0 <= i < nx, nx <= i + nx < 2*nx so i + nx + 1 at max 2*nx

                    // i,j-1
                    if (deposit [ii][(jj + ny - 1) % ny] != 0 &&
                        label   [ii][(jj + ny - 1) % ny] == 0)
                    {
                        // to avoid to push it again on stack when seen by another neighbor
                        label[ii][(jj + ny - 1) % ny] = -1;
                        myStack[top] = ii;
                        top++;
                        myStack[top] = (jj + ny - 1) % ny;
                        top++;
                    }

                    // i-1,j+1
                    if (deposit [(ii + nx - 1) % nx][(jj + 1) % ny] != 0 &&
                        label   [(ii + nx - 1) % nx][(jj + 1) % ny] == 0)
                    {
                        label[(ii + nx - 1) % nx][(jj + 1) % ny] = -1;
                        myStack[top] = (ii + nx - 1) % nx;
                        top++;
                        myStack[top] = (jj + 1) % ny;
                        top++;
                    }

                    // i-1,j
                    if (deposit [(ii + nx - 1) % nx][jj] != 0 &&
                        label   [(ii + nx - 1) % nx][jj] == 0)
                    {
                        label[(ii + nx - 1) % nx][jj] = -1;
                        myStack[top] = (ii + nx - 1) % nx;
                        top++;
                        myStack[top] = jj;
                        top++;
                    }

                    // i-1,j-1
                    if (deposit [(ii + nx - 1) % nx][(jj + ny - 1) % ny] != 0 &&
                        label   [(ii + nx - 1) % nx][(jj + ny - 1) % ny] == 0)
                    {
                        label[(ii + nx - 1) % nx][(jj + ny - 1) % ny] = -1;
                        myStack[top] = (ii + nx - 1) % nx;
                        top++;
                        myStack[top] = (jj + ny - 1) % ny;
                        top++;
                    }

                    // i+1,j-1
                    if (deposit [(ii + nx + 1) % nx][(jj + ny - 1) % ny] != 0 &&
                        label   [(ii + nx + 1) % nx][(jj + ny - 1) % ny] == 0)
                    {
                        label[(ii + nx + 1) % nx][(jj + ny - 1) % ny] = -1;
                        myStack[top] = (ii + nx + 1) % nx;
                        top++;
                        myStack[top] = (jj + ny - 1) % ny;
                        top++;
                    }

                    // i+1,j
                    if (deposit [(ii + nx + 1) % nx][jj] != 0 &&
                        label   [(ii + nx + 1) % nx][jj] == 0)
                    {
                        label[(ii + nx + 1) % nx][jj] = -1;
                        myStack[top] = (ii + nx + 1) % nx;
                        top++;
                        myStack[top] = jj;
                        top++;
                    }

                    // i+1,j+1
                    if (deposit [(ii + nx + 1) % nx][(jj + ny + 1) % ny] != 0 &&
                        label   [(ii + nx + 1) % nx][(jj + ny + 1) % ny] == 0)
                    {
                        label[(ii + nx + 1) % nx][(jj + ny + 1) % ny] = -1;
                        myStack[top] = (ii + nx + 1) % nx;
                        top++;
                        myStack[top] = (jj + ny + 1) % ny;
                        top++;
                    }

                    // i,j+1
                    if (deposit [ii][(jj + ny + 1) % ny] != 0 &&
                        label   [ii][(jj + ny + 1) % ny] == 0)
                    {
                        label[ii][(jj + ny + 1) % ny] = -1;
                        myStack[top] = ii;
                        top++;
                        myStack[top] = (jj + ny + 1) % ny;
                        top++;
                    }
                } // end while

                sizeOfCluster[aggregateIndex] = clusterSize;
                volumeOfCluster[aggregateIndex] = clusterVolume;
                maxHeightOfCluster[aggregateIndex] = maxHeight;
            } // end j loop
        }     // end i loop

        averageSize = 0.;
        averageVol = 0.;
        averageHeight = 0.;

        if (aggregateIndex > 0)
        {
            for (int i = 1; i <= aggregateIndex; i++)
            {
                averageSize += sizeOfCluster[i];
                averageVol += volumeOfCluster[i];
                averageHeight += maxHeightOfCluster[i];
            }

            averageSize = averageSize / (double)aggregateIndex;
            averageVol = averageVol / (double)aggregateIndex;
            averageHeight = averageHeight / (double)aggregateIndex;
        }

        delete[] myStack;
        return aggregateIndex;
    }

    // adhesion process refers to Activated PLTs
    double adhesion(double pAdhesion, double attenuation, double maxAlbumine, double dt, int** candidate_deposition)
    {
        int platelets = 0;
        int newdeposit;
        double p;

        for (int i = 0; i < nx; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                // Check the quantity of deposited albumin
                if (avail[i][j] >= maxAlbumine)
                    continue; // site not available (full of albumin)

                // Deposition takes place with prob pAdhesion on existing platelets
                if (deposit[i][j] != 0)
                {
                    newdeposit = rng(*gen) < pAdhesion * dt * (double)candidate_deposition[i][j] ? 1 : 0;

                    candidate_deposition[i][j] -= newdeposit;
                    
                    deposit[i][j] += newdeposit;
                    platelets += newdeposit; // a platelet has been used
                }
                // pure adhesion (cell PLT empty)
                else
                {
                    // adhesion rate is affected by the amount of albumin locally
                    p = QR(pAdhesion, attenuation, maxAlbumine, avail[i][j]);

                    newdeposit = rng(*gen) < p * dt * (double)candidate_deposition[i][j] ? 1 : 0;

                    candidate_deposition[i][j] -= newdeposit;

                    deposit[i][j] = newdeposit;
                    platelets += newdeposit; // a platelet has been used
                    noCluster += newdeposit;
                }
            } // end of j loop
        }     // end of i loop

        return (double)platelets / (double)(nx * ny);
    }

    // aggregation process refers to NON-Activated PLTs
    double aggregation(double pAggreg, double attenuation, double pTop, double maxAlbumine, double dt, int** candidate_deposition)
    {
        int local_density;
        int platelets = 0;
        int newdeposit;
        double p;

        for (int i = 0; i < nx; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                // deposit on the side
                local_density = deposit [(i + 1 + nx) % nx][j                ] +
                                deposit [(i - 1 + nx) % nx][j                ] +
                                deposit [i                ][(j - 1 + ny) % ny] +
                                deposit [i                ][(j + 1 + ny) % ny];

                // adsorb next to an existing deposit
                if (local_density > 0 && deposit[i][j] == 0)
                {
                    // aggregation rate is affected as well by the amount of albumin locally
                    p = QR(pAggreg, attenuation, maxAlbumine, avail[i][j]);

                    newdeposit = rng(*gen) < p * dt * (double)candidate_deposition[i][j] ? 1 : 0;
                    
                    candidate_deposition[i][j] -= newdeposit;

                    deposit[i][j] = newdeposit;
                    platelets += newdeposit; // the platelets has been consumed
                }

                // deposit on top
                if (deposit[i][j] > 0)
                {
                    newdeposit =  rng(*gen) < pTop * dt * (double)candidate_deposition[i][j] ? 1 : 0;

                    candidate_deposition[i][j] -= newdeposit;

                    deposit[i][j] += newdeposit;
                    platelets += newdeposit; // the platelets has been consumed
                }
            }
        }

        return (double)platelets / (double)(nx * ny);
    }

    double competingAdhesion(double pfilling, double maxAlbumine, double flux)
    {
        double albumine = 0;
        double newPart;
        double p;

        for (int i = 0; i < nx; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                if (avail[i][j] >= maxAlbumine || deposit[i][j] > 0)
                    continue; // site already filled (either al or PLTs)

                // The deposition rate is max for avail[i][j] == 0 &
                // min for avail[i][j] == maxAlbumine
                p = P(pfilling, maxAlbumine, avail[i][j]);
                
                // deposited particles = pd(i,j,t) * N(t) * dt (paper)
                // For albumin this formula has physical correspondence.
                // For AP & NAP this formula is interpreted as a probability 
                // and based on this probability there is/ is not deposition.
                // For ALBUMIN, there is no stochasticity for the phenomenon.
                // This explains why the deposition rates vary a lot as absolute numbers.
                newPart = flux * p;
                
                avail[i][j] += newPart;
                albumine += newPart;
            }
        }

        return (double)albumine / (double)(nx * ny);
    }
}; // -----------------------------end class particles

void Deposition( double *results, unsigned int nrows, unsigned int ncols,
                 int noAP, int noNAP, // number of PLTs per micro-liter
                 double SR_x, // shear rate (100 s^-1 or 400 s^-1)
                 double pAd, double pAg, double pT, double pF, double aT, 
                 double v_z_AP, double v_z_NAP,
                 int seed )
{
    // Everything in SI Units

#ifdef STAND_ALONE
    std::string OutputDir = "./tmp_/";
    mkdir(OutputDir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
#endif

    // Fix or not the output
    std::default_random_engine generator(seed);

    // Substrate discretization
    int nx = 447;
    int ny = 447;
    // Time step
    double dt = 1.e-2;
    // Thickness of the blood layer
    double Lz = 0.82e-3;

    // surface of a dx*dy cell in m^2
    double ds = 5.e-12;
    double mm2 = nx * ny * ds * 1.e6;

    ///////////////////////////////////////////////////////////////////////////
    // Activated PLTs
    ///////////////////////////////////////////////////////////////////////////
    flight activatePlatelets(1.e-3, 1.e-3, Lz, // double lx_, double ly_, double lz_
                            (int)(noAP*0.82), // size_t num_particles_
                            nx, ny, // int nx_, int ny_
                            SR_x, 0., // double shear_rate_x_, double shear_rate_y_
                            dt, // double dt_
                            0., v_z_AP,
                            0., 0.,
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
    
    ///////////////////////////////////////////////////////////////////////////
    // non-Activated PLTs
    ///////////////////////////////////////////////////////////////////////////
    flight nonActivatePlatelets(1.e-3, 1.e-3, Lz, // double lx_, double ly_, double lz_
                                (int)(noNAP*0.82), // size_t num_particles_
                                nx, ny, // int nx_, int ny_
                                SR_x, 0., // double shear_rate_x_, double shear_rate_y_
                                dt, // double dt_
                                0., v_z_NAP,
                                0., 0.,
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

    ///////////////////////////////////////////////////////////////////////////
    // Albumin
    ///////////////////////////////////////////////////////////////////////////
    double rhoAlbumine = 2.69e13 * 1.e9; // albumin / m^3
    double maxAlbumine = 100000; // max number of albumin to saturate a cell
    // Deposition substrate
    // Almost 1 RBC diameter (~10 um)
    double DZ = 1.e-5;
    double N_al = rhoAlbumine * ds * DZ; // initial values of the boundary layer

    ///////////////////////////////////////////////////////////////////////////

    Substrate deposition(nx, ny, &generator);

    ///////////////////////////////////////////////////////////////////////////

    double tMax = 300; // the simulation will run for that many seconds

    int iter = 0;
    int nbIter = tMax / dt; // the number of iterations

    int outputIndex = 0;
    results[nrows * outputIndex + 0] = iter * dt;
    results[nrows * outputIndex + 1] = deposition.countAggregate() / mm2;
    results[nrows * outputIndex + 2] = deposition.averageSize * ds * 1e12;
    results[nrows * outputIndex + 3] = nonActivatePlatelets.computeNoPart() / 0.82; // per micro-liter
    results[nrows * outputIndex + 4] = activatePlatelets.computeNoPart() / 0.82; // per micro-liter

    while (iter < nbIter)
    {
        iter++;

        // Activated PLTs
        activatePlatelets.advance(candidate_deposit_AP);
        deposition.adhesion(pAd, aT, maxAlbumine, dt, candidate_deposit_AP);
        //activatePlatelets.reshuffleSubstrate(candidate_deposit_AP);

        // Albumin
        deposition.competingAdhesion(pF, maxAlbumine, N_al * dt);

        // Non-Activated PLTs
        nonActivatePlatelets.advance(candidate_deposit_NAP);
        deposition.aggregation(pAg, aT, pT, maxAlbumine, dt, candidate_deposit_NAP);
        //nonActivatePlatelets.reshuffleSubstrate(candidate_deposit_NAP);
        
        //------------- measurement branch
        if (iter * dt == 20 || iter * dt == 60 || iter * dt == 120 || iter * dt == 300)
        {
#ifdef STAND_ALONE
            int imax = 0;
            for(int j=0; j<ny; j++)
                for(int i=0; i<nx; i++)
                    if (deposition.deposit[i][j] > imax)
                        imax = deposition.deposit[i][j];

            writePgm    (deposition.deposit, nx, ny, imax, OutputDir + "normal_"+ std::to_string(iter)+ ".pgm");
            writePgm_inv(deposition.deposit, nx, ny, imax, OutputDir + "inv_"   + std::to_string(iter)+ ".pgm");
#endif

            outputIndex++;
            results[nrows * outputIndex + 0] = iter * dt;
            results[nrows * outputIndex + 1] = deposition.countAggregate() / mm2;
            results[nrows * outputIndex + 2] = deposition.averageSize * ds * 1e12;
            results[nrows * outputIndex + 3] = nonActivatePlatelets.computeNoPart() / 0.82; // per micro-liter
            results[nrows * outputIndex + 4] = activatePlatelets.computeNoPart() / 0.82; // per micro-liter
        }
        //---------- end measurmement  branch
    }
    //------------ end time iteration
} // end Deposition


void model( double *results, unsigned int rsize, unsigned int k,
            int noAP, int noNAP,
            double SR_x,
            double pAd, double pAg, double pT, double pF, double aT, 
            double v_z_AP, double v_z_NAP,
            int seed )
{
    std::mt19937 rng(seed);
    
    for (unsigned int i = 0; i < k; ++i)
    {
        int msize = rsize / k;
        int pos = msize * i;
        double *ires = results + pos;

        Deposition( ires, 5, 5,
                    noAP, noNAP,
                    SR_x,
                    pAd, pAg, pT, pF, aT,
                    v_z_AP, v_z_NAP,
                    rng() );
    }
}

#ifdef STAND_ALONE
int main()
{
    //--- declare the matrix to get the results
    int nrows = 5; // time measurements at time 0, 20 , 60, 120 and 300 s
    int ncols = 5; // time, #clusters, average-size, #non-activated p, #activated p
    int k = 1;
    int rsize = nrows * ncols * k;

    double *results = new double[rsize];
    
    //-----------

    // initial conditions for Activated PLTs (AP) & nonActivated PLTs (NAP)
    // ATTENTION: Values refer to ul (micro-liter)
    int noAP   = 4808;
    int noNAP  = 172200;

    // Shear rate
    double SR_x = 100.0;

    // params
    double pAd = 110.0;
    double pAg = 14.6;
    double pT = 0.6;
    double pF = 1.7e-3;
    double aT = 6.0;

    double v_z_AP = 3e-3;
    double v_z_NAP = 3e-4;

    //-----------
    
    model( results, rsize, k,
           noAP, noNAP,
           SR_x,
           pAd, pAg, pT, pF, aT,
           v_z_AP, v_z_NAP,
           1 );

    //-----------

    ofstream output_file;
    output_file.open("output.csv");

    for (int h = 0; h < k; ++h)
    {
        for (int i = 0; i < nrows; i++)
        {
            for (int j = 0; j < ncols; j++)
            {
                int pos = h * nrows * ncols + nrows * i + j;
                cout << results[pos] << " ";
                output_file << results[pos] << ", ";
            }
            cout << endl;
            output_file << endl;
        }
        cout << "------" << endl;
        if (k > 1)
            output_file << "------" << endl;
    }

    delete[] results;
    return 0;
}
#endif