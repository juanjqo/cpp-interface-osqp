typedef double mjtNum;
typedef long long c_int;

#include <iostream>
#include "osqp.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <memory>

using namespace std;
using namespace Eigen;

int main() {
    Eigen::Matrix<mjtNum,2,2> Hessian;
    Hessian << 4, 1,
                1, 2;
    Eigen::Matrix<mjtNum,Eigen::Dynamic,Eigen::Dynamic> tempH = Hessian.triangularView<Eigen::Upper>();
    Eigen::SparseMatrix<mjtNum, Eigen::ColMajor> HessianS = tempH.sparseView();
    HessianS.makeCompressed();

    Eigen::Matrix<mjtNum,1,2> gradient;
    gradient << 1,1;

    Eigen::Matrix<mjtNum,3,2> A;
    A << 1,1,
                1,0,
                0,1;
    Eigen::SparseMatrix<mjtNum, Eigen::ColMajor> A_S = A.sparseView();
    A_S.makeCompressed();

    Eigen::Matrix<mjtNum,3,1> u;
    u << 1,
                0.7,
                0.7;

    Eigen::Matrix<mjtNum,3,1> l;
    l << 1,
                0.0,
                0.0;

    c_int n = 2;
    c_int k = 3;

    // Workspace structures
    OSQPWorkspace *work;
    OSQPSettings  *settings = (OSQPSettings *)c_malloc(sizeof(OSQPSettings));
    OSQPData      *data     = (OSQPData *)c_malloc(sizeof(OSQPData));

    // Populate data
    if (data) {
                data->n = n;
                data->m = k;
                data->A = csc_matrix(k,n, A_S.nonZeros(), A_S.valuePtr(), (c_int*)A_S.innerIndexPtr(), (c_int*)A_S.outerIndexPtr());
                data->P = csc_matrix(n, n, HessianS.nonZeros(), HessianS.valuePtr(), (c_int*)HessianS.innerIndexPtr(), (c_int*)HessianS.outerIndexPtr());
                data->q = gradient.data();//q;
                data->l = l.data();
                data->u = u.data();
            }
    // Define solver settings as default
    if (settings) {
                osqp_set_default_settings(settings);
                settings->alpha = 1.0; // Change alpha parameter
                // settings->verbose = 0;
                //settings->delta = 0.1;
                //settings->time_limit = 0.5;
            }

    // Exitflag
    c_int exitflag = 0;
    // Setup workspace
    exitflag = osqp_setup(&work, data, settings);

    // Solve Problem
    osqp_solve(work);

    // Cleanup
    if (data) {
                if (data->A) c_free(data->A);
                if (data->P) c_free(data->P);
                c_free(data);
            }
    if (settings) c_free(settings);
    return exitflag;
};
