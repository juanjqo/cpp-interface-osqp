#include <iostream>
#include "osqp.h"
#include <Eigen/Dense>
#include <memory>

using namespace std;
using namespace Eigen;



int main(int argc, char **argv)
{
    // Load problem data
        c_float P_x[3] = {4.0, 1.0, 2.0, };
        c_int P_nnz = 3;
        c_int P_i[3] = {0, 0, 1, };
        c_int P_p[3] = {0, 1, 3, };
        c_float q[2] = {1.0, 1.0, };
        c_float A_x[4] = {1.0, 1.0, 1.0, 1.0, };
        c_int A_nnz = 4;
        c_int A_i[4] = {0, 1, 0, 2, };
        c_int A_p[3] = {0, 2, 4, };
        c_float l[3] = {1.0, 0.0, 0.0, };
        c_float u[3] = {1.0, 0.7, 0.7, };
        c_int n = 2;
        c_int m = 3;

        // Exitflag
        c_int exitflag = 0;

        OSQPWorkspace *work;
        std::unique_ptr<OSQPSettings>  settings_sptr(new OSQPSettings());
        std::unique_ptr<OSQPData>      data_sptr(new OSQPData());

        if (data_sptr) {
            data_sptr->n = n;
            data_sptr->m = m;
            data_sptr->P = csc_matrix(data_sptr->n, data_sptr->n, P_nnz, P_x, P_i, P_p);
            data_sptr->q = q;
            data_sptr->A = csc_matrix(data_sptr->m, data_sptr->n, A_nnz, A_x, A_i, A_p);
            data_sptr->l = l;
            data_sptr->u = u;
        }

        if (settings_sptr) {
            osqp_set_default_settings(settings_sptr.get());
            settings_sptr->alpha = 1.0; // Change alpha parameter
            settings_sptr->verbose = 0.0;
        }

        // Setup workspace
        exitflag = osqp_setup(&work, data_sptr.get(), settings_sptr.get());

        // Solve Problem
        osqp_solve(work);
        auto rt0 = Eigen::Map<VectorXd>(work->solution->x, 2);
        std::cout<<"Solution: "<<rt0<<std::endl;
        std::cout<<"exitflag: "<<exitflag<<std::endl;


        osqp_cleanup(work);
        return exitflag;
}
