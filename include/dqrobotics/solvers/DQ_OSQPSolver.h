/**
(C) Copyright 2023 DQ Robotics Developers

This file is part of DQ Robotics.

    DQ Robotics is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    DQ Robotics is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with DQ Robotics.  If not, see <http://www.gnu.org/licenses/>.

Contributors:
- Juan Jose Quiroz Omana       (juanjqo@g.ecc.u-tokyo.ac.jp)
*/

#pragma once

#include <vector>
#include <dqrobotics/solvers/DQ_QuadraticProgrammingSolver.h>
#include <iostream>
#include <osqp.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace Eigen;


namespace DQ_robotics
{
class DQ_OSQPSolver: public DQ_QuadraticProgrammingSolver
{
protected:
    std::unique_ptr<OSQPData> data_sptr;
    std::unique_ptr<OSQPSettings>  settings_sptr;

    /**
     * @brief _dense2sparse
     * @param H
     * @param upper_triangle_flag
     * @return
     */
    SparseMatrix<double> _dense2sparse(const MatrixXd &H, const bool &upper_triangle_flag)
    {
        SparseMatrix<double> mat;
        if (upper_triangle_flag == true)
        {
            mat = MatrixXd(H.triangularView<Upper>()).sparseView();
        }else{
            mat = MatrixXd(H).sparseView();
        }
        mat.makeCompressed();
        return mat;
    }


    void _show_data(const int &m, const int &n, const int &nnz,
                   const Map<VectorXd>& H_x_vec, const Map<VectorXi>& H_i_vec, const Map<VectorXi>& H_p_vec)
    {
        std::cout <<"--------------"<<std::endl;
        std::cout <<m<<", "<<n<<", "<<nnz<<", "<<"{" <<H_x_vec.transpose()<<"}"<<
                    ", {"<<H_i_vec.transpose()<<"},  {"<<H_p_vec.transpose()<<"}"<<std::endl;
        std::cout <<"--------------"<<std::endl;
    }



public:
    DQ_OSQPSolver();
   ~DQ_OSQPSolver()=default;

    /**
     * @brief
     *   Solves the following quadratic program
     *   min(x)  0.5*x'Px + q'x
     *   s.t.    lb <= Ax <= ub
     *
     */
    VectorXd _solve_quadratic_program(const MatrixXd& P, const MatrixXd& A) //, const VectorXd& q, const VectorXd& lb, const VectorXd& ub
    {
        SparseMatrix<double> P_sparse = _dense2sparse(P, true);
        int m_P = P_sparse.rows();
        int n_P = P_sparse.cols();
        std::cout << P_sparse<< "\n";
        int nnz_P = P_sparse.nonZeros();
        auto P_x_vec = Map<VectorXd>(P_sparse.valuePtr(), nnz_P);
        auto P_i_vec = Map<VectorXi>(P_sparse.innerIndexPtr(), nnz_P);
        auto P_p_vec = Map<VectorXi>(P_sparse.outerIndexPtr(), n_P+1);
        _show_data(m_P,n_P,nnz_P, P_x_vec, P_i_vec, P_p_vec);


        SparseMatrix<double> A_sparse = _dense2sparse(P, false);
        int m_A = P_sparse.rows();
        int n_A = P_sparse.cols();
        std::cout << A_sparse<< "\n";
        int nnz_A = P_sparse.nonZeros();
        auto A_x_vec = Map<VectorXd>(P_sparse.valuePtr(), nnz_A);
        auto A_i_vec = Map<VectorXi>(P_sparse.innerIndexPtr(), nnz_A);
        auto A_p_vec = Map<VectorXi>(P_sparse.outerIndexPtr(), n_P+1);
        _show_data(m_P,n_P,nnz_A, A_x_vec, A_i_vec, A_p_vec);

        return VectorXd::Zero(1);
    }

    /**
     * @brief
     *   Solves the following quadratic program
     *   min(x)  0.5*x'Hx + f'x
     *   s.t.    Ax <= b
     *           Aeqx = beq.
     * Method signature is compatible with MATLAB's 'quadprog'.
     * @param H the n x n matrix of the quadratic coeficitients of the decision variables.
     * @param f the n x 1 vector of the linear coeficients of the decision variables.
     * @param A the m x n matrix of inequality constraints.
     * @param b the m x 1 value for the inequality constraints.
     * @param Aeq the m x n matrix of equality constraints.
     * @param beq the m x 1 value for the inequality constraints.
     * @return the optimal x
     */
    VectorXd solve_quadratic_program(const MatrixXd& H, const VectorXd& f, const MatrixXd& A, const VectorXd& b, const MatrixXd& Aeq, const VectorXd& beq) override
    {
        const int PROBLEM_SIZE = H.rows();
        const int INEQUALITY_CONSTRAINT_SIZE = b.size();
        const int EQUALITY_CONSTRAINT_SIZE = beq.size();
        if(EQUALITY_CONSTRAINT_SIZE != 0)
            throw std::runtime_error("DQ_OSQPSolver::solve_quadratic_program(): Equality constraints are not implemented yet.");

        ///Check sizes
        //Objective function
        if(H.rows()!=H.cols())
            throw std::runtime_error("DQ_OSQPSolver::solve_quadratic_program(): H must be symmetric. H.rows()="+std::to_string(H.rows())+" but H.cols()="+std::to_string(H.cols())+".");
        if(f.size()!=H.rows())
            throw std::runtime_error("DQ_OSQPSolver::solve_quadratic_program(): f must be compatible with H. H.rows()=H.cols()="+std::to_string(H.rows())+" but f.size()="+std::to_string(f.size())+".");

        //Inequality constraints
        if(b.size()!=A.rows())
            throw std::runtime_error("DQ_OSQPSolver::solve_quadratic_program(): size of b="+std::to_string(b.size())+" should be compatible with rows of A="+std::to_string(A.rows())+".");

        //Equality constraints
        if(beq.size()!=Aeq.rows())
            throw std::runtime_error("DQ_OSQPSolver::solve_quadratic_program(): size of beq="+std::to_string(beq.size())+" should be compatible with rows of Aeq="+std::to_string(Aeq.rows())+".");

        return VectorXd::Zero(1);
    }

};
}
