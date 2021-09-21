#include "PCHA.h"

void Supdate(mat KC, 
             mat CtKC, 
             float SST, 
             float &muS, 
             float &SSE, 
             mat &S, 
             mat &SSt,
	     int niter)
{
	int noc = S.n_rows;
	int J = S.n_cols;
	vec e = ones(noc, 1);
	for (int k = 0; k < niter; k++)
	{
		float SSE_old = SSE;
		mat g = ((CtKC*S) - KC.t())/(SST/J);
		g = g - e*sum(g % S);
		
		int stop= 0;
		mat Sold = S;
		
		
		while (!stop)
		{
			S = Sold - g * muS;
			S.transform( [](double val) { return (val < 0? 0:val); } );

			rowvec d = sum(S);
			for(int i = 0; i < d.n_elem; i++) {
				if(d(i) != 0)
					S.col(i) = S.col(i) / d(i); 
			}
			
			SSt = S * S.t();
			SSE = SST - 2*sum(sum(S.t() % KC)) + sum(sum(CtKC % SSt));
			if (SSE <= SSE_old*(1+1e-9))
			{
				muS = muS * 1.2;
				stop = 1;
			}
			else {
				muS = muS/2;
			}
		} 
		
	}	
}

void Cupdate(mat K,
	     uvec I,
             uvec U,
             mat S, 
             mat &KC, 
             mat SSt, 
             mat &C, 
             float delta, 
             float &muC, 
             float &mualpha, 
             float SST, 
             float &SSE, 
             float &muS, 
             mat &CtKC,
	     int niter)
{
		
	int JJ = K.n_cols;
	int J = C.n_rows;
	int noc = C.n_cols;
	rowvec alphaC;	
	if(delta != 0)
	{
		alphaC = sum(C);
		rowvec temp; temp.ones(alphaC.n_rows);
		C = C * diagmat(temp/alphaC);
	}
	vec e = ones(J,1);
	mat KSt;
	if (U.n_elem != JJ)
		KSt = K.cols(U) * S.t();
	else
		KSt = K * S.t();
	float SSE_old;
	
	for(int k = 0; k < niter; k++)
	{
		SSE_old = SSE;
		mat g;
		if(I.n_elem != JJ)
			g = (KC.rows(I) * SSt - KSt.rows(I))/SST;
		else
			g = (KC * SSt - KSt)/SST;
		if(delta != 0)
			g = g * diagmat(alphaC);
		g = g - e * sum(g % C);
		int stop = 0;
		mat Cold = C;
		while (!stop)
		{
            C=Cold-muC*g;
			C.transform( [](double val) { return (val < 0? 0:val); } );
			
			rowvec nC = sum(C);
			for(int i=0; i < nC.n_elem; i++) {
				if(nC(i) != 0)
					C.col(i) = C.col(i) / nC(i);
			}
				
			mat Ct;
			if (delta != 0)
				Ct = C * diagmat(alphaC);
			else
				Ct = C;
			if(I.n_rows != JJ)
			{
				KC = K.cols(I) * Ct;
				CtKC = Ct.t() * KC.rows(I);
			}
			else
			{
				KC = K * Ct;
				CtKC = Ct.t() * KC;
			}	
			if (U.n_rows != JJ)
				SSE = SST - 2 * sum(sum(S.t() * K.rows(U))) + sum(sum(CtKC % SSt));
			else
				SSE = SST - 2 * sum(sum(S.t() % KC)) + sum(sum(CtKC % SSt));
			if (SSE <= SSE_old * (1+1e-9))
			{
				muC = muC * 1.2;
				stop = 1;
			}
			else
				muC = muC/2;
			
		}	
	}
	SSE_old = SSE;
	if(delta != 0)
	{
		vec g;
		if(I.n_rows != JJ)
			g = (sum(CtKC % SSt) / alphaC - sum(C%KSt.rows(I)))/(SST*J);
		else
			g = (sum(CtKC % SSt) / alphaC - sum(C%KSt))/(SST*J);
		int stop = 0;
		vec alphaCold = alphaC;
		while (!stop)
		{
			alphaC = alphaCold - mualpha * g;
			uvec indices = find(alphaC < 1 - delta);
			alphaC.rows(indices).fill(1-delta);
			indices = find(alphaC > 1 + delta);
			alphaC.rows(indices).fill(1 + delta);
			mat KCt = KC * diagmat(alphaC / alphaCold);
			mat CtKCt = diagmat(alphaC/alphaCold) * CtKC * diagmat(alphaC/alphaCold);
			if(SSE <= SSE_old * (1+1e-9))
			{
				mualpha = mualpha * 1.2;
				stop = 1;
				KC = KCt;
				CtKC = CtKCt;
			}
			else
				mualpha = mualpha/2;
		}

	}
	if (delta != 0)
		C = C * diagmat(alphaC);		
}	

struct PCHAkernel_ret PCHAkernel(mat Z, mat C, mat S, float delta) {
	mat K;
	if(Z.n_rows == Z.n_cols)
		K = Z;
	else
		K = Z.t()*Z;

	printf("K sum = %e, S sum = %e, C sum = %e\n", sum(sum(K)), sum(sum(C)), sum(sum(S)));
	int noc = C.n_cols;
	int n = K.n_rows;	
	
	printf("noc = %d, n = %d\n", noc, n);
	
	uvec U(n);
	uvec I(n);
	for(int i = 0; i < n; i++) {
		U(i) = i;
		I(i) = i;
	}
	
	struct PCHAkernel_ret ret;
	ret.C = C;
	ret.S = S;
	
	if(delta == 0)
	{
		ret.alphaC.eye(noc, noc);
	}
	float conv_crit = 1e-6; 
	int maxiter = 500;
	float SST = trace(K.submat(U,U));


	//matlab code has something to initialise C here, but I'm assuming it's being passed as an argument

	
	
	mat KC = K.cols(I) * ret.C;
	
	float muS = 1;
	float muC = 1;
	float mualpha = 1;


	//matlab code checks if S is in the opts
	mat CtKC = ret.C.t() * K.submat(I,I) * ret.C;
	mat SSt = ret.S*ret.S.t();
	ret.SSE = SST - 2*(sum(sum(ret.S.t() % KC.rows(U)))) + sum(sum(CtKC%SSt));
	printf("part 1 = %e, part 2 = %e, part 3 = %e, part 4 = %e\n", sum(sum(ret.S.t())), (sum(sum(KC.rows(U)))), sum(sum(CtKC)), sum(sum(SSt))); 
		
	int iter = 0;
	int dSSE = 1e10; //inf
	int t1 = clock();
	float varexpl = (SST - ret.SSE)/SST;
	
	printf("SST = %e, SSE = %e, varexp = %e, err == %e\n", SST, ret.SSE, varexpl, conv_crit * abs(ret.SSE) );
	
	
	while( abs(dSSE) >= conv_crit * abs(ret.SSE) && iter < maxiter && varexpl < 0.9999)
	{		
		int told = t1;
		iter++;
		float SSE_old = ret.SSE;
		
		Cupdate(K,I,U,ret.S, KC, SSt, ret.C, delta, muC, mualpha, SST, ret.SSE, muS, CtKC, 5);

		Supdate(KC.rows(U), CtKC, SST, muS, ret.SSE, ret.S, SSt, 5);
	
		dSSE = SSE_old - ret.SSE;
		t1 = clock();
		//if rem(iter, 1) == 0???
		varexpl = (SST - ret.SSE)/SST;
	
	}		
	
	printf("Done\n");
	//varexpl = (SST - ret.SSE)/SST;
	/*uvec indices = sort_index(sum(ret.S,1), "descend");
	ret.S = ret.S(indices);
	ret.C = ret.C.cols(indices);*/
	return ret;
}

