#include <ACTIONet.h>
//#include <mini_cblas.h>
#include <cassert>

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
        printf("\tPre ... ");
        mat grad = (AtA * X) - AtB;
        mat obj = A * X - B;
        printf("done\n");
        
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
    
    /*
    X = clamp(X, 0, 1);
    X = normalise(X, 1);
    */

    return (X);
}


}  // namespace ACTIONet
