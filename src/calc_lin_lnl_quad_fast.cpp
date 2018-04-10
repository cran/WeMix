#include <RcppArmadillo.h>
#include <map>
#include <iterator>

using namespace Rcpp;

const double log_2_pi = log(2.0*M_PI);

//'Calculate multivariate normal using the inverse cholesky of sigma
//'@param x data matrix
//'@param mean vector of means of x
//'@param rooti inverse of the cholesky decomposition of the  variance covariance matrix
//'@return vector of densities
//'
//'@export
// [[Rcpp::depends(RcppArmadillo)]]
arma::vec dmvnrm_arma(arma::mat x,  
                      arma::rowvec mean,  
                      arma::mat rooti) { 
 //figure out nrow and ncol
  int n_obs = x.n_rows;
  int n_c = x.n_cols;

  // calcualte the constant parts outside for loop
  double chol_plus_logpi=-2.0* log(prod(diagvec(rooti)))+n_c*log_2_pi;

  arma::vec raw_results(n_obs);

  for (int i=0; i < n_obs;i=i+1){
    arma::vec diff = rooti *(x.row(i) - mean).t();   
    raw_results(i)  =  sum(diff%diff);     
  }  
  
  return(-.5*(raw_results+chol_plus_logpi));
}


// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
arma::mat sumByGroup(arma::colvec index, arma::colvec value) {
  std::map <int, double> sumResult;
  
  if (index.size() == value.size()) {
    int n = index.size();
    for (int i=0; i < n; i++){
      sumResult[index[i]] += value[i];
    }
    
    arma::mat sumResultMat(sumResult.size(), 2);
    int i = 0;
    for (auto it : sumResult) {
      sumResultMat(i, 0) = it.first;
      sumResultMat(i, 1) = it.second;
      i++;
    }
    return(sumResultMat);
  } else {
    Rcout << "Error: "
          << "Index and value has different length."
          << std::endl;
    arma::mat sumResultMat(sumResult.size(), 2);
    return(sumResultMat);
  }
}

// [[Rcpp::depends(RcppArmadillo)]]
arma::mat deduplicate(arma::colvec index, arma::colvec value) {
  std::map <int, double> deDupResult;
  if (index.size() == value.size()) {
    int n = index.size();
    for (int i=0; i < n; i++){
      deDupResult[index[i]] = value[i];
    }
    
    arma::mat deDupResultMat(deDupResult.size(), 2);
    int i = 0;
    for (auto it : deDupResult) {
      deDupResultMat(i, 0) = it.first;
      deDupResultMat(i, 1) = it.second;
      i++;
    }
    return(deDupResultMat);
  } else {
    Rcout << "Error: "
          << "Index and value has different length."
          << std::endl;
    arma::mat deDupResultMat(deDupResult.size(), 2);
    return(deDupResultMat);
  }
}

// [[Rcpp::depends(RcppArmadillo)]]
arma::mat CartesianJoin(arma::colvec x, arma::colvec y, bool product = false) {
  arma::mat result(x.n_rows * y.n_rows, 2);
  int index = 0;
  for (size_t y_i = 0; y_i < y.n_rows; y_i++) {
    for (size_t x_i = 0; x_i < x.n_rows; x_i++) {
      result(index, 0) = x(x_i);
      result(index, 1) = y(y_i);
      index += 1;
    }
  }
  if (product) {
    return(result.cols(0, 0) % result.cols(1, 1));
  } else {
    return(result);
  }
}


// [[Rcpp::depends(RcppArmadillo)]]
//'This function calcuates the liklihood of the model using integration by adaptive quadrature 
//'@param y a numeric vector, the response.
//'@param yhat the current predicted values of the response
//'@param level an integer that respresents the number of levels in the likelihood
//'         that is desired. In a two level model this function will be called
//'         with l=2 and it will recurse and call itself with l=1
//'@param Z a list of the Z matricies where the index of Z indicates the level.
//'          of the Z matri. Z[[1]] is NULL because there is no individual lavel Zs.
//'@param ZFull Z expended such that each Z[[i]] matrix contains  one row for each observation (ie each element has same number of rows as original data)
//'@param Qi the scaling factor for the adaptive quadratures (per group)
//'@param QiFull the scaling factors for adaptive quarature points duplicated so each element has one row per observation (like Zfull)
//'@param omega a list of the b estimates for ech group
//'@param omegaFull numeric, the b estimates for each group duplicated so each element has one row per observation (like Zfull)
//'@param C a list of Cholesky decompositions of the Sigma matricies.
//'         C[[1]] is simply the residual variance (a scalar) while C[[l]] for 
//'         l > 1 is a matrix with the name number of rows and columns as the 
//'         Z matrix for that level.
//'@param qp Gaussian quadrature result from statmod::gauss.quad.
//'@param W list of weight matricies. must have `w` and `index` columns
//'@param top boolean set to TRUE to return a single scalar, otherwise returns a vector
//'@param verbose boolean set to TRUE to get verbose output
//'@param acc numeric, accuracy of the mpfr
//'@param atPoint boolean, indicates likelihood should be calculated at single point
//'               at the top level and then integrated below that. This is useful
//'               for finding the maximum posterior (or likelihood) extimate for 
//'               the random effects
//'@param integralMultiplierExponent a single integer, the function evaluates the
//'                                  integral times the randome effect to this power
//'                                  when set to 0, it is just the log-likelhood
//'                                  when set to 1, this can be used to estimate the
//'                                  expected value.
//'@param integralZColumn is the column index of Z to use integralMultiplierExponent on
//'                       only one random effect at a time can be integrated over, and the integration happens at the top level.
//'
//'@description
//'calculates the log-likelihood of an l level mixed model using adaptive quadrature.
//'the genral model is y = Xb + ZU + e
//'but that values of beta and U are not included in the call. Instead this information
//'is contained in yhat which incorporates Xb at the top level and all relevant 
//'Zu information at lower levels.
//'@author Huade Huo, Paul Bailey, Claire Kelley
//'@export
// [[Rcpp::export]]
arma::mat calc_lin_lnl_quad_fast (arma::vec y,
                                  arma::vec yhat, 
                                  int level,
                                  List Z,
                                  List Qi,
                                  List omega,
                                  List W,
                                  List C,
                                  List qp,
                                  List omegaFull,
                                  List QiFull,
                                  List ZFull,
                                  bool top = true,
                                  bool atPoint = false,
                                  int integralMultiplierExponent = 0,
                                  int integralZColumn = 1,
                                  bool verbose = true,
                                  int acc=120) {
  
  // R index start at 1, while cpp is zero-indexed
  arma::mat data;
  arma::mat W1 = W[0];
  
  data = arma::join_rows(y, W1);
  arma::mat Cl = C[level - 1];
  arma::mat Wl = W[level - 1];
  arma::mat Wlm1 = W[level - 2];
  arma::mat Zl = Z[level - 1];
  arma::mat ZFulll = ZFull[level - 1];
  arma::mat Qil = Qi[level - 1];
  arma::mat QiFulll = QiFull[level - 1];
  
  arma::colvec qp_v;
  arma::colvec qp_w1;
  arma::colvec qp_nodes = qp["nodes"];
  arma::colvec qp_weights = qp["weights"];
  
  // Create zero (one) matrix for v and w if atPoint
  if (atPoint) {
    qp_v = arma::zeros<arma::colvec>(1);
    qp_w1 = arma::ones<arma::colvec>(1);
  } else {
    qp_v = qp_nodes;
    qp_w1 = qp_weights;
  }
  
  // Catersian product
  arma::colvec grdw = qp_w1;
  arma::mat grdv = qp_v;
  
  //  declare i as size_t if they will be compared to sizes
  for (size_t grdv_i = 2; grdv_i <= Cl.n_rows; grdv_i++) {
    grdw = CartesianJoin(grdw, qp_w1, true);
    grdv = CartesianJoin(grdv, qp_v, false);
  }
  
  
  arma::colvec wlzero(Wl.n_rows, arma::fill::zeros);
  Wl = arma::join_rows(Wl, wlzero);
  
  arma::uvec i_idx;
  arma::mat grdv_i;
  arma::mat grdv_i_t;
  arma::mat Wlm1_pts, pts;
  arma::mat omegal = omega[level - 1];
  arma::mat omegaFulll = omegaFull[level - 1];
  arma::mat yyh;
  arma::colvec yyh_colvec;
  arma::mat Wlm1_ll;
  arma::colvec y_dnorm(y.size());
  arma::mat agg;
  arma::colvec agg_li;
  arma::colvec agg_li_pre_exp;
  arma::mat Wlm1NonD;
  arma::colvec Wlm1NonD_g_weight;
  arma::mat Wl_ll;
  
  arma::mat rooti = arma::trans(arma::inv(trimatu(Cl)));
  
  arma::colvec Wl_l = arma::zeros<arma::colvec>(Wl.n_rows);
  
  for (size_t i =  0; i < grdw.n_rows; i++) {
    // Generate a vector with regularly spaced elements
    i_idx = arma::regspace<arma::uvec>(i, i);
    grdv_i = grdv.rows(i_idx);
    grdv_i_t = arma::trans(grdv_i);
    
    if (level == 2) {
      Wlm1_pts = omegal + std::sqrt(2) * arma::trans(arma::reshape(grdv_i * Qil, size(arma::trans(omegal))));
      yyh = yhat + arma::sum(Zl % Wlm1_pts, 1);
      yyh_colvec = yyh.cols(0, 0);
      for (size_t dnorm_i = 0; dnorm_i < y.size(); dnorm_i++) {
        y_dnorm(dnorm_i) = R::dnorm(y(dnorm_i), 
                yyh_colvec(dnorm_i), C[0], TRUE);
      }
      // TODO: Select first column in wlm1 (it's hardcoded now)
      Wlm1_ll = Wlm1.cols(0, 0) % y_dnorm;
    } else {
      Wlm1_pts = omegal + std::sqrt(2) * arma::trans(grdv_i * Qil);
      pts = omegaFulll + std::sqrt(2) * arma::trans(grdv_i * QiFulll);
      Wlm1_ll = calc_lin_lnl_quad_fast(y, yyh, level-1, Z, Qi,
                                       omega, W, C, qp, omegaFull, QiFull, ZFull,
                                       false, false, 0,
                                       verbose, acc);
    }
    agg = sumByGroup(Wlm1.cols(1,1), Wlm1_ll);
    if (atPoint) {
      if (!top) {
        return(agg.cols(1,1));
      } else {
        return(sum(agg.cols(1, 1)));
      }
    } else {
      Wlm1NonD = deduplicate(Wlm1.cols(1,1), Wlm1_pts.cols(0, 0));
      for (size_t Wlm1_pts_i = 1; Wlm1_pts_i < Wlm1_pts.n_cols; Wlm1_pts_i++) {
        Wlm1NonD = arma::join_rows(Wlm1NonD, 
                                   deduplicate(Wlm1.cols(1,1), 
                                               Wlm1_pts.cols(Wlm1_pts_i, Wlm1_pts_i)).cols(1, 1));
      }
      
      Wlm1NonD_g_weight = arma::zeros(Wlm1NonD.n_rows);
      for (size_t g_weight_j = 0; g_weight_j < Wlm1NonD.n_rows; g_weight_j++) {
        //
        //Rcout << "x" << std::endl;
        //Rcout << Wlm1NonD(arma::span(g_weight_j, g_weight_j), arma::span(1,Wlm1NonD.n_cols - 1)) << std::endl;
       // Rcout << "mean" << std::endl;
        //Rcout << arma::zeros<arma::rowvec>(Wlm1NonD.n_cols - 1) << std::endl;
        //Rcout << "rooti" << std::endl;
        //Rcout << rooti << std::endl;
        Wlm1NonD_g_weight(g_weight_j) = dmvnrm_arma(
          Wlm1NonD(arma::span(g_weight_j, g_weight_j), arma::span(1,Wlm1NonD.n_cols - 1)), 
          arma::zeros<arma::rowvec>(Wlm1NonD.n_cols - 1), 
          rooti)[0];
        //Rcout << "results" << std::endl;
        //Rcout << Wlm1NonD_g_weight(g_weight_j)  << std::endl;
      }
      if ((log(grdw(i)) + agg.cols(1,1) + Wlm1NonD_g_weight + arma::accu(grdv_i_t % grdv_i_t)).min() < -700) {
        warning("Warning: The likelihood is probably inaccurate. Set fast = TRUE to get a more accurate estimation.");
      }
      agg_li = exp(log(grdw(i)) + agg.cols(1,1) + Wlm1NonD_g_weight + arma::accu(grdv_i_t % grdv_i_t));
      if (integralMultiplierExponent != 0) {
        agg_li = agg_li % pow(Wlm1NonD.cols(1, 1), integralMultiplierExponent);
      }
      Wl_l += agg_li * pow(2, Zl.n_cols / 2.0) % Wl.cols(2, 2);
    }
  }
  
  if (integralMultiplierExponent != 0) {
    return(Wl_l);
  }
  
  Wl_ll = Wl.cols(1,1) % arma::trunc_log(Wl_l);
  
  if (top) {
    return(sum(Wl_ll.cols(0, 0)));
  }
  
  return(Wl_ll);
}


