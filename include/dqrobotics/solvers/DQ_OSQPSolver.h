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

- Juan Jose Quiroz Omana (juanjqo@g.ecc.u-tokyo.ac.jp)
        - Adapted the file DQ_QPOASESSolver.h implemented by Murilo M. Marinho (murilo@nml.t.u-tokyo.ac.jp)
          in (https://github.com/dqrobotics/cpp-interface-qpoases) to use the osqp solver
          (https://osqp.org/)
*/

#pragma once
#include <vector>
#include <dqrobotics/solvers/DQ_QuadraticProgrammingSolver.h>
#include <iostream>
#include <osqp.h>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>


using namespace Eigen;

namespace DQ_robotics
{
class DQ_OSQPSolver: public DQ_QuadraticProgrammingSolver
{
protected:
    //
    std::shared_ptr<OSQPSettings>  settings_sptr_;
    std::shared_ptr<OSQPWorkspace> work_sptr_;
    std::shared_ptr<OSQPData> data_sptr_;

    VectorXd x_0_;
    VectorXd y_0_;
    int MAX_ITER_ = 300;
    bool SOLVE_FIRST_TIME_;
    VectorXd u_solution_;

    /**
     * @brief returns a sparse matrix from a dense matrix.
     * @param H the dense matrix.
     * @param use_upper_triangle_matrix
     * @return the sparse eigen matrix
     */
    SparseMatrix<double> _dense2sparse(const MatrixXd &H, const bool &use_upper_triangle_matrix)
    {
        SparseMatrix<double> mat;
        if (use_upper_triangle_matrix == true)
        {
            mat = MatrixXd(H.triangularView<Upper>()).sparseView();
        }else{
            mat = MatrixXd(H).sparseView();
        }
        mat.makeCompressed();
        return mat;
    }


public:
    DQ_OSQPSolver():SOLVE_FIRST_TIME_(true)
    {
        work_sptr_     = std::make_shared<OSQPWorkspace>(OSQPWorkspace());
        settings_sptr_ = std::make_shared<OSQPSettings>(OSQPSettings());
        data_sptr_     = std::make_shared<OSQPData>(OSQPData());
    }
   ~DQ_OSQPSolver()=default;

    /**
     * @brief Solves the following quadratic program
     *   min(x)  0.5*x'Px + f'x
     *   s.t.    Ax <= ub
     *
     * @param H the n x n matrix of the quadratic coeficitients of the decision variables.
     * @param f the n x 1 vector of the linear coeficients of the decision variables.
     * @param A the m x n matrix of inequality constraints.
     * @param b the m x 1 value for the inequality constraints.
     * @param Aeq the m x n matrix of equality constraints.
     * @param beq the m x 1 value for the inequality constraints.
     * @return the optimal x
     */
    VectorXd _solve_osqp_quadratic_program(const MatrixXd& H, const VectorXd& f, const MatrixXd& A,  const VectorXd& b,
                                           const MatrixXd& Aeq, const VectorXd& beq)
    {
        SparseMatrix<double> P_sparse = _dense2sparse(H, true);
        //const int m_P = P_sparse.rows();
        const int n_P = P_sparse.cols();
        const int nnz_P =  P_sparse.nonZeros();

        c_float P_x_array[nnz_P];
        c_int   P_i_array[nnz_P];
        c_int   P_p_array[n_P+1];
        std::copy(P_sparse.valuePtr(),      P_sparse.valuePtr()      +nnz_P, P_x_array);
        std::copy(P_sparse.innerIndexPtr(), P_sparse.innerIndexPtr() +nnz_P, P_i_array);
        std::copy(P_sparse.outerIndexPtr(), P_sparse.outerIndexPtr() +n_P+1, P_p_array);

        SparseMatrix<double> A_sparse;
        VectorXd ub_;

        // If there are no constraints the solver fails. A naive way to circumvent this is to impose
        // constraints in the form
        //    Identity*x <= inf.
        if (A.size() == 0)
        {
            A_sparse = _dense2sparse(MatrixXd::Identity(n_P,n_P), false);
            ub_ = VectorXd::Ones(n_P)*OSQP_INFTY;
        }
        else{
            A_sparse = _dense2sparse(A, false);
            ub_ = b;
        }

        /**
         * The solver expects constraints in the form
         *     lb <= Ax <= ub.
         * Therefore, we set lb = [-inf] to implement the constraints
         * in the form:
         *      Ax <= ub.
         */
        VectorXd lb_;
        lb_ = -1*VectorXd::Ones(n_P)*OSQP_INFTY;


        const int m_A = A_sparse.rows();
        const int n_A = A_sparse.cols();
        int nnz_A =  A_sparse.nonZeros();

        c_float A_x_array[nnz_A];
        c_int   A_i_array[nnz_A];
        c_int   A_p_array[n_A+1];

        std::copy(A_sparse.valuePtr(),      A_sparse.valuePtr()      +nnz_A, A_x_array);
        std::copy(A_sparse.innerIndexPtr(), A_sparse.innerIndexPtr() +nnz_A, A_i_array);
        std::copy(A_sparse.outerIndexPtr(), A_sparse.outerIndexPtr() +n_A+1, A_p_array);

        //--------------
        osqp_set_default_settings(settings_sptr_.get());
        settings_sptr_->alpha = 1.5; // Change alpha parameter 1.0
        settings_sptr_->verbose = 0.0;
        settings_sptr_->max_iter = MAX_ITER_;
        settings_sptr_->warm_start = 1.0;
        VectorXd f_ = f;

        data_sptr_->n = n_A;
        data_sptr_->m = m_A;
        data_sptr_->P = csc_matrix(data_sptr_->n, data_sptr_->n, nnz_P, P_x_array, P_i_array , P_p_array);
        data_sptr_->q = f_.data();
        data_sptr_->A = csc_matrix(data_sptr_->m, data_sptr_->n, nnz_A, A_x_array, A_i_array , A_p_array);
        data_sptr_->l = lb_.data();
        data_sptr_->u = ub_.data();

        // Setup workspace
        OSQPWorkspace* work = work_sptr_.get();

        //OSQPWorkspace* work;
        auto exitflag = osqp_setup(&work, data_sptr_.get(), settings_sptr_.get());
        //std::cout<<"exitflag: "<<exitflag<<std::endl;


        if (SOLVE_FIRST_TIME_ == false)
        {
            c_float x_0_array[data_sptr_->n];
            c_float y_0_array[data_sptr_->n];
            std::copy(x_0_.data(), x_0_.data()+data_sptr_->n, x_0_array);
            std::copy(y_0_.data(), y_0_.data()+data_sptr_->n, y_0_array);
            osqp_warm_start_x(work, x_0_array);
            //MAX_ITER_ = 300;

        }
        // Solve Problem
        osqp_solve(work);

        x_0_ = Eigen::Map<VectorXd>(work->solution->x, data_sptr_->n);
        y_0_ = Eigen::Map<VectorXd>(work->solution->y, data_sptr_->n);
        u_solution_ = Eigen::Map<VectorXd>(work->solution->x, data_sptr_->n);

        osqp_cleanup(work);
        SOLVE_FIRST_TIME_ = false;
        return u_solution_;
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
    VectorXd solve_quadratic_program(const MatrixXd& H, const VectorXd& f,
                                     const MatrixXd& A, const VectorXd& b,
                                     const MatrixXd& Aeq, const VectorXd& beq) override
    {
        ///---------------------------------------------------------------------------------------------------
        /// Copy from https://github.com/dqrobotics/cpp-interface-qpoases
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
        ///---------------------------------------------------------------------------------------------------
        //OSQP_INFTY
        return _solve_osqp_quadratic_program(H, f, A, b, Aeq, beq);
    }

};
}


