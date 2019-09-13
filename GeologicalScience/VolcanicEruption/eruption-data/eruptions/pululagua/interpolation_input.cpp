#include "gsl/gsl_interp.h"
#include "mpfr.h"
#include <iostream>
#include <vector>

std::vector<double> createVec(double start, double stop, int n)
{
    std::vector<double> res(n);
    int i = 0;
    for (double v = start; v <= stop; v += (stop - start) / (n - 1))
    {
        res[i] = v;
        i++;
    }
    if (res.back() != stop)
        res[res.size() - 1] = stop;
    return res;
}

void printVec(std::vector<double> v)
{
    for (auto i : v)
        std::cout << i << ", ";
    std::cout << std::endl;
}

std::vector<std::vector<double>> interpolate(std::vector<double> xs, std::vector<double> ys, int n)
{
    auto xs_interp = createVec(xs.front(), xs.back(), n);
    std::vector<double> ys_interp(xs_interp.size());

    auto interp = gsl_interp_alloc(gsl_interp_linear, xs.size());
    gsl_interp_init(interp, xs.data(), ys.data(), xs.size());
    auto acc = gsl_interp_accel_alloc();

    for (int i = 0; i < ys_interp.size(); i++)
    {
        ys_interp[i] = gsl_interp_eval(interp, xs.data(), ys.data(), xs_interp[i], acc);
    }

    gsl_interp_free(interp);
    gsl_interp_accel_free(acc);

    return std::vector<std::vector<double>>{xs_interp, ys_interp};
}

double sum(std::vector<double> v)
{
    mpfr_t res_mpfr;
    mpfr_init2 (res_mpfr, 200);
    mpfr_set_d (res_mpfr, 0.0, MPFR_RNDD);

    mpfr_t e_mpfr;
    mpfr_init2 (e_mpfr, 200);

    for (auto e : v){
        mpfr_set_d (e_mpfr, e, MPFR_RNDD);
        mpfr_add (res_mpfr, res_mpfr, e_mpfr, MPFR_RNDD);
    }

    return mpfr_get_d (res_mpfr, MPFR_RNDD);
}

// normalize v to make its sum = sum
std::vector<double> normalize(std::vector<double> v, double targetsum){
    mpfr_t e_mpfr, sum_mpfr, targetsum_mpfr, sumdivtargetsum_mpfr;
    mpfr_init2 (e_mpfr, 200);
    mpfr_init2 (sum_mpfr, 200);
    mpfr_init2 (targetsum_mpfr, 200);
    mpfr_init2 (sumdivtargetsum_mpfr, 200);
    mpfr_set_d (sum_mpfr, sum(v), MPFR_RNDD);
    mpfr_set_d (targetsum_mpfr, targetsum, MPFR_RNDD);
    mpfr_div (sumdivtargetsum_mpfr, sum_mpfr, targetsum_mpfr, MPFR_RNDD);

    for (auto &e : v){
        mpfr_set_d (e_mpfr, e, MPFR_RNDD);
        mpfr_div (e_mpfr, e_mpfr, sumdivtargetsum_mpfr, MPFR_RNDD);
        e = mpfr_get_d (e_mpfr, MPFR_RNDD);
    }

    return v;
}

int main()
{
    auto xs = createVec(-7, 10, 18);
    std::vector<double> densities = {600, 600, 600, 600, 600, 600, 600, 837.5, 1075, 1312.5, 1550, 1787.5, 2025, 2262.5, 2500, 2500, 2500, 2500};
    std::vector<double> wtpct = {0.139, 0.4, 1.65, 3.27, 5.96, 9.346, 12.756, 16.72, 21.037, 18.871, 7.089, 1.082, 0.603, 0.716, 0.26, 0.082, 0.016, 0.003};

    //printVec(xs);
    //printVec(densities);

    auto interpolated = interpolate(xs, densities, 10 * xs.size());
    auto xs_inter = interpolated.front();
    auto densities_inter = interpolated.back();
    interpolated = interpolate(xs, wtpct, 10 * xs.size());
    auto wtpct_inter = interpolated.back();
    wtpct_inter = normalize(wtpct_inter, sum(wtpct));

    //printVec(interpolated.front());
    //printVec(interpolated.back());

    //printVec(wtpct);

    std::cout << "sum of wtpct : " << sum(wtpct) << std::endl;

    //printVec(wtpct);
    
    std::cout << "sum of wtpct_inter : " << sum(wtpct_inter) << std::endl;

    //printVec(wtpct_inter);
}