#include <algorithm>  // std::random_shuffle
#include<random>
#include<iostream>
#include <RcppEigen.h>
using namespace Rcpp;  // call any Rcpp constructs by their given name without the Rcpp:: prefix.
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::SparseMatrix;
using Eigen::SparseVector;
typedef Eigen::MappedSparseMatrix<double> MapSMatXd;
typedef Map<VectorXd> MapVecXd;


// TODO:
//  - currently all betas output are NaN... need to debug...
//  - support random sampling data points.
//  - call inline functions [in current version, calling inline functions are
//  very slow.]

// inline float wProb(SparseVector<double> x, VectorXd betas) {
    // // int nsamples = x.cols();
    // // VectorXd W(nsamples);
    // // for (int i=0; i<nsamples; i++) {
    // float value = -x.dot(betas);
    // float w = 1 / (1 + exp(value));
        // // W(i) = w;
    // // }
    // // return W;
    // return w;
// }

// inline double logLik(SparseVector<double> x, double y, VectorXd betas, double lambda=0) {
    // // int nsamples = x.cols();
    // // VectorXd m(nsamples, 1);
    // int m = 1;
    // // double l = 0;
    // // for (int i=0; i<nsamples; i++) {
    // float value = x.dot(betas);
    // double l = y * value - m * log(1 + exp(value));
    // // }
    // // double penalty = lambda * betas.dot(betas);
    // // double l_penalty = l + penalty;
    // return -l;
// }

// inline VectorXd logLikPrime(SparseVector<double> x, double y, VectorXd betas, double lambda=0) {
    // double w = wProb(x, betas);
    // int nfeatures = x.size();
    // VectorXd l_prime(nfeatures);
    // for (SparseVector<double>::InnerIterator it(x); it; ++it) {
        // int i = it.index();
        // double value = -it * (y - 1 * w) + lambda * 2 * betas(i);
        // l_prime(i) = value;
    // }
    // // VectorXd l_prime = -x * (y - 1 * w) + lambda * 2 * betas;
    // return l_prime;
// }

// inline VectorXd samplingData(int y_size, int maxiter, bool replace=true) {
    // std::vector<int> y_index(y_size);
    // for (int i=0; i<y_size; i++) y_index[i] = i;  // vector including index.
    // std::vector<int> index;  //initialize a vector
    // if (!replace and maxiter > y_size) {
        // int a = maxiter / y_size;
        // int b = maxiter % y_size;
        // for (int i=0; i<a; i++) {
            // std::random_shuffle(y_index.begin(), y_index.end());
            // index.insert(index.end(), y_index.begin(), y_index.end());
        // }
        // std::random_shuffle(y_index.begin(), y_index.end());
        // y_index.resize(b);
        // index.insert(index.end(), y_index.begin(), y_index.end());
    // }
    // else {
        // std::random_device rd;
        // std::mt19937 gen(rd());
        // std::uniform_int_distribution<int> dis(0, y_index.back());
        // for (int n=0; n<maxiter; n++) index.push_back(dis(gen));
    // }
    // VectorXd index_eigen(index.size());
    // for (int i=0; i<index.size(); i++) {
        // index_eigen(i) = index[i];
    // }
    // return (index_eigen);
// }


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
SEXP AdaGradSGD(MapSMatXd x, VectorXd y, VectorXd betas0, double step_size, int rep, double lambda, bool replace=true, double weight=0.5) {


    // // generate sample index
    // VectorXd sample_index = sgdSubFn::samplingData(y_size, maxiter, replace);
    // // sample data
    // SparseMatrix<double> x_sample;
    // x_sample.setZero(nfeatures, maxiter);
    // VectorXd y_sample(maxiter);
    // for (int i=0; i<maxiter; i++) {
        // y_sample(i) = y(sample_index(i));
    // }
    // for (SparseVector<double>::InnerIterator it(x_sample); it; ++it) {
        // int i = it.index();
        // // double value = -it * (y - 1 * w) + lambda * 2 * betas(i);
        // // l_prime(i) = value;
        // it = it + x.innerVector(sample_index(i));
    // }
    // initialization
    int nsamples, nfeatures;
    SparseMatrix<double> x_sample;
    VectorXd y_sample, betas, gd_cum, adjusted_gd, last_update;
    x_sample = x;
    y_sample = y;
    nsamples = y_sample.size();
    nfeatures = betas0.size();
    betas = betas0;
    gd_cum.setZero(nfeatures);
    adjusted_gd.setZero(nfeatures);
    last_update.setZero(nfeatures);

    // int iter = 1;
    int iter = 0;
    double w, l_average, l_weighted;
    l_average = 0;
    l_weighted = 0;
    VectorXd l_avg_tracking(nsamples*rep);
    VectorXd l_weighted_tracking(nsamples*rep);

    // evaluate betas
    for (int n=0; n<rep; n++) {
        for (int i=0; i<nsamples; i++) {
            SparseVector<double> x_single = x_sample.innerVector(i);
            double y_single = y_sample(i);

            // calculate log likelihood for each single sample.
            // double l_single = logLik(x_single, y_single, betas);  // when
            // calling inline function, it becomes very slow...
            // calculate w.
            // w = sgdSubFn::wProb(x_single, betas);
            double w_ = -x_single.dot(betas);
            w = 1 / (1 + exp(w_));
            double l_single = - y_single * (-w_) + 1 * log(1 + exp(-w_));
            // // tracking exponentially weighted average
            // l_weighted = weight * l_single + (1-weight) * l_weighted;
            // l_weighted_tracking(iter-1) = l_weighted;
            // tracking avg log likelihood
            if (iter > 0) {
                l_average = (l_single + (iter-1)*l_avg_tracking(iter-1)) / iter;
                l_avg_tracking(iter) = l_average;
            }
            // l_avg_tracking(iter-1) = l_average;

            for (SparseVector<double>::InnerIterator it(x_single); it; ++it) {
                int j = it.index();

                // int skip = iter - 1 - last_update(j);
                int skip = iter - last_update(j);
                last_update(j) = iter;

                double cum_penalty = betas(j) * ((1 - pow(1 + lambda * step_size * adjusted_gd(j), skip))/(1 - lambda * step_size * adjusted_gd(j)));

                betas(j) -= cum_penalty;

                double gd_single = -it * (y_single - 1 * w) + lambda * 2 * betas(j);

                gd_cum(j) += gd_single * gd_single;
                adjusted_gd(j) = gd_single / (std::sqrt(gd_cum(j)) + 1e-06);
                betas(j) -= step_size * adjusted_gd(j);
            }

            iter += 1;
        }
    }

    for (int i=0; i<nfeatures; i++) {
        int skip = iter - 1 - last_update(i);
        double cum_penalty = betas(i) * ((1 - pow(1 + lambda * step_size * adjusted_gd(i), skip))/(1 - lambda * step_size * adjusted_gd(i)));
        betas(i) -= cum_penalty;

    }
    return List::create(Named("betas") = betas, Named("l_weighted") = l_weighted_tracking, Named("l_avg") = l_avg_tracking);
}
