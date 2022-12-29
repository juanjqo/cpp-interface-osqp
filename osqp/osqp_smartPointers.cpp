#include <iostream>
#include "osqp.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <memory>

using namespace std;
using namespace Eigen;


void eigen_to_csc(const MatrixXd &H, const bool &upper_triangle_flag)
{
    int m = H.rows();
    int n = H.cols();
    Eigen::SparseMatrix<double> mat;
    if (upper_triangle_flag == true)
    {
        mat = MatrixXd(H.triangularView<Upper>()).sparseView();
    }else{
        mat = MatrixXd(H).sparseView();
    }

    std::cout << mat<< "\n\n";
    int nnz = mat.nonZeros();
    std::cout <<"Number of Non zeros entries: "<<nnz<<std::endl;
    double H_x[nnz];
    double H_i[nnz];
    double H_p[n+1];

    int i = 0;
    int prev_col;
    int counter=0;
    for (int k=0; k<mat.outerSize(); ++k)
      for (SparseMatrix<double>::InnerIterator it(mat,k); it; ++it)
      {
        it.value();
        it.row();   // row index
        it.col();   // col index (here it is equal to k)
        it.index(); // inner index, here it is equal to it.row()
        H_x[i] = it.value();
        H_i[i] = it.row();
        std::cout << it.value()<<" ("<<it.row()<<", "<<it.col()<<") . Index: "<<it.index()<<" . Iter: "<<i<<std::endl;


        if (i==0){
            prev_col = it.col();
            H_p[counter] = i;
            counter++;
        }
        else if (it.col()  != prev_col)
        {
            H_p[counter] = i;
            prev_col = it.col();
            counter++;
        }
        if (i==nnz-1)
        {
            H_p[counter] = i+1;
            counter++;

        }








        i++;
      }

    auto H_x_vec = Eigen::Map<Eigen::VectorXd>(H_x, nnz);
    auto H_i_vec = Eigen::Map<Eigen::VectorXd>(H_i, nnz);
    auto H_p_vec = Eigen::Map<Eigen::VectorXd>(H_p, n+1);

    std::cout <<m<<", "<<n<<", "<<nnz<<", "<<"{" <<H_x_vec.transpose()<<"}"<<
                ", {"<<H_i_vec.transpose()<<"},  {"<<H_p_vec.transpose()<<"}"<<std::endl;


    std::cout <<"--------------"<< "\n\n";
}

void eigen2csc(const MatrixXd &H, const bool &upper_triangle_flag)
{
    int m = H.rows();
    int n = H.cols();
    Eigen::SparseMatrix<double, Eigen::ColMajor> mat;

    if (upper_triangle_flag == true)
    {
        mat = MatrixXd(H.triangularView<Upper>()).sparseView();
    }else{
        mat = MatrixXd(H).sparseView();
    }
    mat.makeCompressed();

    std::cout <<m<<", "<<n<<", "<<mat.nonZeros()<<", "<<"{" <<mat.valuePtr()<<"}"<<
                ", {"<<(c_int*)mat.innerIndexPtr()<<"},  {"<<(c_int*)mat.outerIndexPtr()<<"}"<<std::endl;

}

typedef long long c_int;

int main(int argc, char **argv)
{
    MatrixXd H = MatrixXd(2,2);
    H << 4,1,1,2;
    //eigen2csc(H, true);

    //MatrixXd A = MatrixXd(3,2);
    //A << 1,1,1,0,0,1;
    //eigen_to_csc(A, false);


    Matrix<double,Dynamic,Dynamic> tempH = H.triangularView<Upper>();
    SparseMatrix<double, ColMajor> HessianS = tempH.sparseView();
    HessianS.makeCompressed();




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
            //data_sptr->P = csc_matrix(data_sptr->n, data_sptr->n, P_nnz, P_x, P_i, P_p);
            data_sptr->P = csc_matrix(n, n, HessianS.nonZeros(), HessianS.valuePtr(), (c_int*)HessianS.innerIndexPtr(), (c_int*)HessianS.outerIndexPtr());
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
