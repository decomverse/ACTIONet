#include <ACTIONet.h>
//#include <mini_cblas.h>
#include <cassert>

#include <chrono>
#define duration(a) std::chrono::duration_cast<std::chrono::nanoseconds>(a).count()
#define now() std::chrono::high_resolution_clock::now()


// Re-implemented from: Fast and Robust Archetypal Analysis for Representation
// Learning
namespace ACTIONet {

// min(|| AX - B ||) s.t. simplex constraint
mat run_simplex_regression_FW_base(mat& A, mat& B, int max_iter, double min_diff) {
    if(max_iter == -1)
        max_iter = A.n_cols;

    printf("Initializing ... ");
    //mat tmp = trans(A) * B;
    mat tmp = cor(A, B);
    
    mat X = zeros(size(tmp));

    for(int j = 0; j < X.n_cols; j++) {
        uword i = index_max(tmp.col(j));
        X(i, j) = 1;
        //X(0, j) = 1;
    } 
    printf("done\n");

    mat At = trans(A);
    mat AtA = At * A;
    mat AtB = At * B;

    mat old_X = X;
    for(int it = 0; it < max_iter; it++) {
        mat grad = (AtA * X) - AtB;
        mat obj = A * X - B;

        for(int k = 0; k < X.n_cols; k++) {
            vec g = grad.col(k);
            vec x = X.col(k);

            int i1, i2;
            int i1_val, i2_val;
            i1_val = i2_val = 1000;
            i1 = i2 = 0;
            for(int i = 0; i < g.n_elem; i++) {
                if(g(i) < i1_val) {
                    i1_val = g(i);
                    i1 = i;
                }
                if(x(i) > 0) {
                    if(g(i) < i2_val) {
                        i2_val = g(i);
                        i2 = i;
                    }
                }
            }
            
            vec d_FW = -x;
            d_FW(i1) = 1 + d_FW(i1);

            vec d_A = x;
            d_A(i2) = d_A(i2) - 1;

            double alpha_max = 1;
            vec d;
            if( dot(g, d_FW) < dot(g, d_A)) {
                d = d_FW;
                alpha_max = 1;
            } else {
                d = d_A;
                alpha_max = x(i1) / (1 - x(i1));
            }

            /*
            // Backtracking line-search
            vec Ad = A * d;
            double e1 = dot(Ad, Ad);
            double alpha = 0;
            if(e1 != 0) {
                double e2 = 2 * dot(obj.col(k), Ad);
                double e3 = 0.5* dot(g, d); // multiplier can be in (0, 0.5]
                alpha = (e3 - e2) / e1;
            }
            */
            double alpha = 2.0 / (it+ 2);

            X.col(k) = x + alpha*d;
        }

        printf("%d- ", it);
        double res = sum(sum(abs(old_X - X))) / X.n_cols;
        printf("%e\n", res);

        if(res < min_diff) {
            break;
        }
        old_X = X;
    }
    
    X = clamp(X, 0, 1);
    X = normalise(X, 1);

    return (X);
}



mat run_simplex_regression_FW_test1(mat& A, mat& B, int max_iter, double min_diff) {
    if(max_iter == -1)
        max_iter = A.n_cols;

    
    mat X = zeros(A.n_cols, B.n_cols);
    X.row(0).ones();
    
   
   /*
    mat tmp = cor(A, B);
    tmp(tmp < 0).zeros();
    mat X = normalise(tmp, 1, 0);
    */

    /*
    printf("Initializing ... ");
    //mat tmp = trans(A) * B;
    mat tmp = cor(A, B);
    for(int j = 0; j < X.n_cols; j++) {
        uword i = index_max(tmp.col(j));
        X(i, j) = 1;
        //X(0, j) = 1;
    } 
    printf("done\n");
    */

    mat At = trans(A);
    mat AtA = At * A;


    vec d;
    mat old_X = X;
    for(int it = 0; it < max_iter; it++) {
        mat grad = (AtA * X) - At * B;
        mat obj = A * X - B;

        mat mask = X;
        mask.transform( [](double val) { return(val == 0? datum::inf:1.0); } );

        mat masked_grad = grad % mask;
        urowvec ii1 = index_min(grad);
        urowvec ii2 = index_min(masked_grad);

        mat D_FW = -X;
        mat D_A = X;
        vec alpha_caps(X.n_cols);
        for(int j = 0; j < X.n_cols; j++) {
            double x = X(ii1(j), j);
            alpha_caps(j) = x / (1-x);
            D_FW(ii1(j), j)++;
            D_A (ii2(j), j)--;
        }

        for(int k = 0; k < X.n_cols; k++) {
            vec g = grad.col(k);
            vec x = X.col(k);
            
            vec d_FW = D_FW.col(k);
            vec d_A = D_A.col(k);
        
            double alpha_max = 1;
            if( dot(g, d_FW) < dot(g, d_A) ) {
                d = d_FW;
            } else {
                d = d_A;
                alpha_max = alpha_caps(k);
            }

            // Backtracking line-search
            vec Ad = A * d;
            double e1 = dot(Ad, Ad);
            double alpha = 0;
            if(e1 != 0) {
                double e2 = 2 * dot(obj.col(k), Ad);
                double e3 = 0.5* dot(g, d); // multiplier can be in (0, 0.5]
                alpha = (e3 - e2) / e1;
            }
            
           // double alpha = 2.0 / (it+ 2);
            alpha = min(alpha, alpha_max);
            
            X.col(k) = x + alpha*d;
        }

        printf("%d- ", it);
        //double res = norm(abs(old_X - X), "fro") / norm(X, "fro");
        //printf("%e\n", res);

/*
        if(res < min_diff) {
            break;
        }
*/        
        old_X = X;
    }
    
    /*
    X = clamp(X, 0, 1);
    X = normalise(X, 1);
    */

    return (X);
}


/*
mat run_simplex_regression_FW(mat& A, mat& B, int max_iter, double min_diff) {
    if(max_iter == -1)
        max_iter = A.n_cols;

    printf("Initializing ... ");
    mat tmp = cor(A, B);    
    mat X = zeros(size(tmp));

    for(int j = 0; j < X.n_cols; j++) {
        vec v = tmp.col(j);
        int i = index_max(v);
        X(i, j) = 1;
    } 
    printf("done\n");

    mat At = trans(A);
    mat AtA = At * A;
    mat AtB = At * B;

    mat old_X = X;
    for(int it = 0; it < max_iter; it++) {

        mat grad = (AtA * X) - AtB;
        mat obj = A * X - B;

        for(int k = 0; k < X.n_cols; k++) {            
            vec g = grad.col(k);
            vec x = X.col(k);
            vec b = B.col(k);

            int i1, i2;
            int i1_val, i2_val;
            i1_val = i2_val = 1000;
            i1 = i2 = 0;
            for(int i = 0; i < g.n_elem; i++) {
                if(g(i) < i1_val) {
                    i1_val = g(i);
                    i1 = i;
                }
                if(x(i) > 0) {
                    if(g(i) < i2_val) {
                        i2_val = g(i);
                        i2 = i;
                    }
                }
            }
            vec d_FW = -x;
            d_FW(i1) = 1 + d_FW(i1);

            vec d_A = x;
            d_A(i2) = d_A(i2) - 1;
            
            double alpha_max = 1;
            vec d;
            if( dot(g, d_FW) < dot(g, d_A)) {
                d = d_FW;
                alpha_max = 1;
            } else {
                d = d_A;
                alpha_max = x(i2) / (1 - x(i2));
            }

            double alpha = 0;

            vec q = A * d;
            double q_norm = dot(q, q);
            if(q_norm > 0) {
                alpha = dot(q, b - A*x) / q_norm;
            }
            alpha = min(alpha, alpha_max);
            
            X.col(k) = x + alpha*d;
        }

        printf("%d- ", it);
        double res = mean(mean(abs(old_X - X)));
        printf("%e\n", res);

        if(res < min_diff) {
            break;
        }
        old_X = X;
    }
    
    
    X = clamp(X, 0, 1);
    X = normalise(X, 1);
    

    return (X);
}
*/
mat run_simplex_regression_FW_working(mat& A, mat& B, int max_iter, double min_diff) {

    double t1 = 0, t2 = 0, t3 = 0, t4 = 0, t5 = 0, t6 = 0;
    std::chrono::duration<double> elapsed;
    std::chrono::high_resolution_clock::time_point start, finish;

    if(max_iter == -1)
        max_iter = A.n_cols;

    start = now();
    //printf("Initializing ... ");
    mat tmp = cor(A, B);    
    mat X = zeros(size(tmp));

    for(int j = 0; j < X.n_cols; j++) {
        vec v = tmp.col(j);
        int i = index_max(v);
        X(i, j) = 1;
    } 
    //printf("done\n");

    mat At = trans(A);
    mat AtA = At * A;
    mat AtB = At * B;
    finish = now();
    t1 += duration(finish - start);

    mat old_X = X;
    for(int it = 0; it < max_iter; it++) {

        start = now();
        mat grad = (AtA * X) - AtB;
        mat AX = A*X;
        mat obj = AX - B;        
        finish = now();
        t2 += duration(finish - start);

        for(int k = 0; k < X.n_cols; k++) {

            start = now();
            vec g = grad.col(k);
            vec x = X.col(k);
            vec b = B.col(k);

            int i1, i2;
            int i1_val, i2_val;
            i1_val = i2_val = 1000;
            i1 = i2 = 0;
            for(int i = 0; i < g.n_elem; i++) {
                if(g(i) < i1_val) {
                    i1_val = g(i);
                    i1 = i;
                }
                if(x(i) > 0) {
                    if(g(i) < i2_val) {
                        i2_val = g(i);
                        i2 = i;
                    }
                }
            }
            finish = now();
            t3 += duration(finish - start);

            start = now();
            vec d_FW = -x;
            d_FW(i1) = 1 + d_FW(i1);

            vec d_A = x;
            d_A(i2) = d_A(i2) - 1;
            
            double alpha_max = 1;
            vec d;
            int direction = 0;
            if( dot(g, d_FW) < dot(g, d_A)) {
                direction = +1; // Adding element
                d = d_FW;
                alpha_max = 1;
            } else {
                direction = -1; // Removing element
                d = d_A;
                alpha_max = x(i2) / (1 - x(i2));
            }
            finish = now();
            t4 += duration(finish - start);

            start = now();

            vec q;
            if(direction == +1) {
                q = A.unsafe_col(i1) - AX.unsafe_col(k);
            } else {
                q = AX.unsafe_col(k) - A.unsafe_col(i2);
            }
            //vec q = A * d;

/*
            double alpha = 0;
            double q_norm = norm(q, 1);
            if(q_norm > 0) {
                double q_norm_sq = q_norm*q_norm;
                vec delta = -obj.col(k);
                alpha = dot(q, delta) / q_norm_sq;
            }            
*/          
            double alpha = 2 / (it + 2);  

            alpha = min(alpha, alpha_max);
            //X.col(k) = x + alpha*d;
            finish = now();
            t5 += duration(finish - start);

        }

        start = now();
        //("%d- ", it);
        double res = norm(old_X - X, "fro");
        //printf("%e\n", res);
        finish = now();
        t6 += duration(finish - start);
/*
        if(res < min_diff) {
            break;
        }
*/        
        old_X = X;
    }
    
    double total = t1 + t2 + t3 + t4 + t5 + t6;
    //printf("t1 = %3.f, t2 = %3.f, t3 = %3.f, t4 = %3.f, t5 = %3.f, t6 = %3.f\n", 100*t1/total, 100*t2/total, 100*t3/total, 100*t4/total, 100*t5/total, 100*t6/total);

    X = clamp(X, 0, 1);
    X = normalise(X, 1);
    
    return (X);
}




mat run_simplex_regression_FW(mat& A, mat& B, int max_iter, double min_diff) {

    double t1 = 0, t2 = 0, t3 = 0, t4 = 0, t5 = 0, t6 = 0;
    std::chrono::duration<double> elapsed;
    std::chrono::high_resolution_clock::time_point start, finish;

    int m = A.n_rows, n = A.n_cols, k = B.n_cols;
    bool compute_AtA = n < 1000; //(n*n*(m+k)) < (2*m*n*k);

    if(max_iter == -1)
        max_iter = A.n_cols;

    start = now();
    mat tmp = cor(A, B);    
    mat X = zeros(size(tmp));
    for(int j = 0; j < X.n_cols; j++) {
        vec v = tmp.col(j);
        int i = index_max(v);
        X(i, j) = 1;
    } 
    mat mask = X;

    mat At = trans(A);
    mat AtB = At * B;

    mat AtA;
    if( compute_AtA ) {
        AtA = At * A;
    }    

    finish = now();
    t1 += duration(finish - start);

    mat grad;
    mat old_X = X;    
    double min_grad_norm = 1e-3;
    for(int it = 0; it < max_iter; it++) {

        start = now();
        mat AX = A*X;
        mat obj = AX - B;

        if(compute_AtA) {
            grad = (AtA * X) - AtB;
        } else {
            grad = (At * (A*X)) - AtB;
        }
        rowvec grad_norms = sqrt(sum(square(grad)));
        uvec unsaturated_cols = find(grad_norms > min_grad_norm);
        if(unsaturated_cols.n_elem == 0) {
            break;
        }
        //printf("%d- %d unsaturated (mean grad_norm = %2e)\n", it, unsaturated_cols.n_elem, mean(grad_norms));

        urowvec s = index_min(grad.cols(unsaturated_cols), 0);
        urowvec v = index_min(grad.cols(unsaturated_cols) % mask.cols(unsaturated_cols), 0);

        mat S = mat(sp_mat(join_vert(s, regspace<urowvec>(0,  s.n_elem-1)), ones(s.n_elem), X.n_rows, unsaturated_cols.n_elem));
        mat V = mat(sp_mat(join_vert(v, regspace<urowvec>(0,  v.n_elem-1)), ones(v.n_elem), X.n_rows, unsaturated_cols.n_elem));
        mat D_FW = S - X.cols(unsaturated_cols);
        mat D_A = X.cols(unsaturated_cols) - V;

/*
        mat Q_FW = A * D_FW;
        mat Q_A = A * D_A;
*/
        rowvec delta = sum(D_FW % grad.cols(unsaturated_cols)) - sum(D_A % grad.cols(unsaturated_cols));

/*
        res(0) = X;
        res(1) = grad;
        res(2) = D_FW;
        res(3) = D_A;
        return(res);
*/

        finish = now();
        t2 += duration(finish - start);

        for(int k = 0; k < unsaturated_cols.n_elem; k++) {
            int j = unsaturated_cols(k);

            vec d;
            double alpha_max = 1;            
            if( (delta(k) <= 0) | (X(v(k), j) == 0) ) {
                d = D_FW.col(k);
                alpha_max = 1;
                mask(s(k), j) = 1;
            } else {
                d = D_A.col(k);
                alpha_max = X(v(k), j) / (1 - X(v(k), j));
            }

            start = now();

/*
            // From: https://thatdatatho.com/gradient-descent-line-search-linear-regression/
            vec g = grad.col(j);
            double num = dot(g, g), denom, alpha = 0;            
            if( compute_AtA ) {
                vec AtAg = AtA *  g;
                denom = dot(g, AtAg);
            } else {
                vec Ag = A * g;
                denom = dot(Ag, Ag);
            }
            alpha =  num / denom;
            printf("\t%d- num = %.2e, denom = %.2e, Alpha = %.2e\n", it, num, denom, alpha);
*/

            double alpha = 2.0 / (it + 2.0);  

/*
            vec q = A * d;
            double alpha = 0;
            double q_norm = norm(q, 1);
            if(q_norm > 0) {
                double q_norm_sq = q_norm*q_norm;
                vec delta = -obj.col(j);
                alpha = dot(q, delta) / q_norm_sq;
            }            
*/

            /*
            // Backtracking line-search
            vec g = grad.col(j);
            vec Ad = A * d;
            double e1 = dot(Ad, Ad);
            double alpha = 0;
            if(e1 != 0) {
                double e2 = 2 * dot(obj.col(j), Ad);
                double e3 = 0.5* dot(g, d); // multiplier can be in (0, 0.5]
                alpha = (e3 - e2) / e1;
            }
            
*/
            alpha = min(alpha, alpha_max);
            X.col(j) += alpha*d;
            if(0 < delta(k)) {
                if(X(v(k), j) < 1e-6)
                    mask(v(k), j) = 0;
            }
            finish = now();
            t5 += duration(finish - start);

        }

        start = now();
        //("%d- ", it);
        //double res = norm(old_X - X, "fro");
        //printf("%e\n", res);
        finish = now();
        t6 += duration(finish - start);
/*
        if(res < min_diff) {
            break;
        }
*/        
        //old_X = X;
    }
    
    double total = t1 + t2 + t3 + t4 + t5 + t6;
    //printf("t1 = %3.f, t2 = %3.f, t3 = %3.f, t4 = %3.f, t5 = %3.f, t6 = %3.f\n", 100*t1/total, 100*t2/total, 100*t3/total, 100*t4/total, 100*t5/total, 100*t6/total);

    //X = clamp(X, 0, 1);
    //X = normalise(X, 1);
    
    return (X);
}
}  // namespace ACTIONet
