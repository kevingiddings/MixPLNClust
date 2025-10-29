#include <RcppArmadillo.h>
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <numeric>
#include <cmath>
#include <algorithm>
#include <set>
#include <limits>
using namespace Rcpp;
using namespace arma;
using MatList = std::vector<arma::mat>&;
using VecList = std::vector<arma::vec>&;
using VecMatList = std::vector<std::vector<arma::mat>>&;
using VecVecList = std::vector<std::vector<arma::vec>>&;



// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
List update_g_params(Rcpp::List& GX_g, Rcpp::List& dGX_g, Rcpp::List& z_S_g, Rcpp::List& m_g, Rcpp::List& O_list, Rcpp::List& start_g, Rcpp::List& S_g, Rcpp::List& iOsigO_g, arma::vec& lib_mat, arma::mat& Y, arma::vec& mu_g, arma::mat& z, const int g, const int N) {

    if (S_g.size() < N) S_g = Rcpp::List(N);  
    if (dGX_g.size() < N) dGX_g = Rcpp::List(N);
    if (z_S_g.size() < N) z_S_g = Rcpp::List(N);
    if (GX_g.size() < N) GX_g = Rcpp::List(N);
    if (m_g.size() < N) m_g = Rcpp::List(N);

    lib_mat = arma::clamp(lib_mat, 1e-16, arma::datum::inf);


    for (int i = 0; i < N; i++) {
        arma::mat O_i = Rcpp::as<arma::mat>(O_list[i]);
        arma::vec start_g_i = Rcpp::as<vec>(start_g[i]);
        arma::mat S_g_i = Rcpp::as<arma::mat>(S_g[i]);
        arma::mat iOsigO_g_i = Rcpp::as<mat>(iOsigO_g[i]);
        
        arma::mat dGX_g_i = diagmat(arma::exp(O_i * vectorise(arma::log(lib_mat)) + start_g_i) + 0.5 * S_g_i.diag()) + iOsigO_g_i;
        S_g_i = inv(dGX_g_i);

        arma::mat z_S_g_i = z(i, (g-1)) * S_g_i;


        arma::vec GX_g_i = O_i * vectorise(Y.row(i)) - vectorise(arma::exp(start_g_i + O_i * vectorise(arma::log(lib_mat)) + 0.5 * S_g_i.diag())) - iOsigO_g_i * (start_g_i - O_i * vectorise(mu_g));        



        arma::vec m_g_i = start_g_i + S_g_i * vectorise(GX_g_i);

        S_g[i] = S_g_i;
        dGX_g[i] = dGX_g_i;
        z_S_g[i] = z_S_g_i;
        GX_g[i] = GX_g_i;
        m_g[i] = m_g_i;
    }

    return List::create(
        Rcpp::Named("S_g") = Rcpp::wrap(S_g),
        Rcpp::Named("dGX_g") = Rcpp::wrap(dGX_g),
        Rcpp::Named("z_S_g") = Rcpp::wrap(z_S_g),
        Rcpp::Named("GX_g") = Rcpp::wrap(GX_g),
        Rcpp::Named("m_g") = Rcpp::wrap(m_g));
}


arma::mat for_mu(int i, int g, Rcpp::List iOsigO_g, arma::mat z, Rcpp::List& O_list, arma::mat m_g_i){
    arma::mat O_i = Rcpp::as<mat>(O_list[i]);
    arma::mat iOsigO_g_i = Rcpp::as<mat>(iOsigO_g[i]);

    arma::mat temp = z(i, (g-1)) * trans(O_i) * iOsigO_g_i;
    arma::vec num = temp * m_g_i;
    arma::mat den = temp * O_i;

    arma::mat result = join_rows(den, num);

    return result;
}

// [[Rcpp::export]]
List update_mu(Rcpp::List& m_g, Rcpp::List& O_list, Rcpp::List& iOsigO_g, arma::vec& mu_g, arma::mat z, const int d, const int g, const int N) {
    arma::mat mu_temp = zeros(d, d + 1);
    for (int i = 0; i < N; i++)
    {
        arma::mat for_mu_i = for_mu(i, g, iOsigO_g, z, O_list, m_g[i]);
        mu_temp += for_mu_i;
    }
    mu_g = inv(mu_temp.submat(0, 0, d - 1, d - 1)) * mu_temp.col(d);
    return List::create(Rcpp::Named("mu_g") = mu_g);
}

arma::mat for_sig(int i, int g, Rcpp::List iOsigO_g, arma::mat z, Rcpp::List O_list, arma::vec m_g_i, arma::vec mu_g, arma::mat S_g_i){
    arma::mat O_i = Rcpp::as<mat>(O_list[i]);
    arma::mat iOsigO_g_i = Rcpp::as<mat>(iOsigO_g[i]);

    arma::mat omega = (m_g_i - O_i * mu_g) * (m_g_i - O_i * mu_g).t() + S_g_i;
    arma::mat forsig1 = z(i, (g-1)) * trans(O_i) * iOsigO_g_i * omega * iOsigO_g_i * O_i;
    arma::mat forsig2 = z(i, (g-1)) * trans(O_i) * iOsigO_g_i * O_i;
    //Rcpp::Rcout << "FOR SIG: " << (forsig1 - forsig2) << "i and g:" << i << " " << g << std::endl;
    return -(forsig1 - forsig2);
}



// [[Rcpp::export]]
List update_sig(Rcpp::List m_g, Rcpp::List sigma_new, arma::mat sigma_g, Rcpp::List O_list, Rcpp::List iOsigO_g, arma::vec mu_g, Rcpp::List S_g, arma::mat z, const int d, const int g, const int N, const double step) {
    arma::mat gr = zeros(d, d);

    for (int i = 0; i < N; i++){
        arma::vec m_g_i = Rcpp::as<arma::vec>(m_g[i]);
        arma::mat S_g_i = Rcpp::as<arma::mat>(S_g[i]);

        arma::mat for_sig_i = for_sig(i, g, iOsigO_g, z, O_list, m_g_i, mu_g, S_g_i);
        gr += for_sig_i;
    }
    arma::mat sigma_new_g = sigma_g - step * gr;

    return List::create(Rcpp::Named("sigma_new_g") = sigma_new_g);
}




arma::mat for_sig2(int i, int g, Rcpp::List iOsigO_g, arma::mat z, Rcpp::List O_list, arma::vec m_g_i, arma::vec mu_g, arma::mat S_g_i) {
    arma::mat O_i = Rcpp::as<mat>(O_list[i]);
    arma::mat iOsigO_g_i = Rcpp::as<mat>(iOsigO_g[i]);

    arma::mat omega = (m_g_i - O_i * mu_g) * (m_g_i - O_i * mu_g).t() + S_g_i;
    arma::mat forsig1 = z(i, (g - 1)) * trans(O_i) * iOsigO_g_i * omega * iOsigO_g_i * O_i;
    arma::mat forsig2 = z(i, (g - 1)) * trans(O_i) * iOsigO_g_i * O_i;
    
    return -(forsig1 - forsig2);
}

// [[Rcpp::export]]
List update_sig2(Rcpp::List m_g, Rcpp::List sigma_new, arma::mat sigma_g, Rcpp::List O_list, Rcpp::List iOsigO_g, arma::vec mu_g, Rcpp::List S_g, arma::mat z, const int d, const int g, const int N, const double step) {
    arma::mat gr = zeros(d, d);

    for (int i = 0; i < N; i++) {
        arma::vec m_g_i = Rcpp::as<arma::vec>(m_g[i]);
        arma::mat S_g_i = Rcpp::as<arma::mat>(S_g[i]);

        arma::mat for_sig_i = for_sig(i, g, iOsigO_g, z, O_list, m_g_i, mu_g, S_g_i);
        gr += for_sig_i;
    }
    arma::mat sigma_new_g = sigma_g - step * gr;

    return List::create(Rcpp::Named("sigma_new_g") = sigma_new_g);
}






arma::mat symmetrize(const arma::mat& M) {
    arma::mat M_sym = 0.5 * (M + M.t());

    return M_sym;
}


// [[Rcpp::export]]
bool PD_check(Rcpp::List& sigma_new, int G) {
    Function eigen_R("eigen");
    for (int g = 0; g < G; ++g) {
        arma::mat sigma_g = Rcpp::as<arma::mat>(sigma_new[g]);



        try { 
            arma::vec eig_val = arma::eig_sym((sigma_g + sigma_g.t()) / 2.0);

            //List eigen_result = eigen_R(sigma_g);
            //NumericVector eig_val = eigen_result["values"];
            //arma::vec eig_vec = as<arma::vec>(eig_val);
            
            //best eigen function needs to be found

            // spd
            if (arma::any(eig_val <= 0)) {
                //Rcpp::Rcout << "Matrix at index " << g << " is not positive definite" << std::endl;
                return false;
            }
        }
        catch (const std::exception& e) {
            Rcpp::Rcout << "Eigen decomposition failed for matrix at index " << g << ": " << e.what() << std::endl;
            //Rcpp::Rcout << "Matrix sigma_g:\n" << sigma_g << std::endl;
            //Rcpp::Rcout << "Determinant is:\n" << arma::det(sigma_g) << std::endl;
            return false; 
        }


    }
    return true;
}


// [[Rcpp::export]]
void invert_matrices(Rcpp::List& sigma, Rcpp::List& i_sigma, int G) {
    for (int g = 0; g < G; ++g) {
        arma::mat sigma_g = Rcpp::as<arma::mat>(sigma[g]);
        i_sigma[g] = arma::inv(sigma_g);
    }
}

// [[Rcpp::export]]
arma::rowvec compute_pi_g(const arma::mat& z, int N) {
    return arma::sum(z, 0) / static_cast<double>(N);
}

// [[Rcpp::export]]
arma::mat create_lib_mat_full(const arma::rowvec& lib_mat, int N) {
    return arma::repmat(lib_mat, N, 1); 
}

// [[Rcpp::export]]
void compute_iOsigO(const Rcpp::List& O_list, const Rcpp::List& sigma, Rcpp::List& iOsigO) {
    int G = sigma.size();
    int N = O_list.size();

    for (int g = 0; g < G; ++g) {
        arma::mat sigma_g = Rcpp::as<arma::mat>(sigma[g]);
        Rcpp::List iOsigO_g(N);

        for (int i = 0; i < N; ++i) {
            arma::mat O_i = Rcpp::as<arma::mat>(O_list[i]);
            arma::mat intermediate = O_i * sigma_g * O_i.t();

            iOsigO_g[i] = arma::inv(intermediate);
        }

        iOsigO[g] = iOsigO_g;
    }
}


void add_smallest_positive(arma::mat& F) {
    arma::vec flat = arma::vectorise(F); 
    arma::vec positives = flat.elem(arma::find(flat > 0));

    if (!positives.is_empty()) {
        double smallest_positive = positives.min();
        F += smallest_positive;
    }
}

void na_omit_mat(arma::mat& mat) {
    arma::uvec good_rows = arma::find_finite(arma::sum(mat, 1));
    mat = mat.rows(good_rows);
}

// [[Rcpp::export]]
Rcpp::List compute_F_matrices(const Rcpp::List& S, const Rcpp::List& m, const Rcpp::List& O_list,
    const Rcpp::List& mu, const Rcpp::List& iOsigO, arma::mat& O_mat,
    const arma::mat& Y, const arma::vec& pi_g) {
    int G = S.size();
    int N = O_list.size();

    arma::mat F_raw(N, G);
    arma::mat F(N, G);

    F_raw.fill(datum::nan);
    F.fill(datum::nan);

    //na_omit_mat(O_mat);

    arma::vec log_pi_g = arma::log(pi_g);

    for (int g = 0; g < G; ++g) {
        Rcpp::List S_g = Rcpp::as<Rcpp::List>(S[g]);
        Rcpp::List m_g = Rcpp::as<Rcpp::List>(m[g]);
        Rcpp::List iOsigO_g = Rcpp::as<Rcpp::List>(iOsigO[g]);
        arma::vec mu_g = Rcpp::as<arma::vec>(mu[g]);

        for (int i = 0; i < N; ++i) {
            int O_sum = 0;
            for (uword x = 0; x < O_mat.n_cols; ++x) {
                if (O_mat(i, x) == O_mat(i, x)) {
                    O_sum++;
                }
            }
            
            
            arma::mat S_gi = Rcpp::as<arma::mat>(S_g[i]);
            arma::mat iOsigO_gi = Rcpp::as<arma::mat>(iOsigO_g[i]);
            arma::mat O_i = Rcpp::as<arma::mat>(O_list[i]);
            arma::vec m_gi = Rcpp::as<arma::vec>(m_g[i]);

            F_raw(i, g) = 0.5 * (log(arma::det(S_gi)) + log(arma::det(iOsigO_gi))) -
                0.5 * arma::as_scalar((m_gi - O_i * mu_g).t() * iOsigO_gi * (m_gi - O_i * mu_g)) -
                arma::trace(iOsigO_gi * S_gi) + 0.5 * O_sum +
                arma::as_scalar(m_gi.t() * O_i * Y.row(i).t()) -
                arma::sum(arma::exp(m_gi + 0.5 * S_gi.diag()) + arma::lgamma(O_i * Y.row(i).t() + 1));


            //F(i, g) = pi_g(g) * std::exp(F_raw(i, g));
        }
    }

    for (int i = 0; i < N; ++i) {
        arma::rowvec log_prob = F_raw.row(i) + log_pi_g.t();

        double max_log_prob = log_prob.max();
        arma::rowvec stabilized_log_prob = log_prob - max_log_prob;

        arma::rowvec exp_probs = arma::exp(stabilized_log_prob);
        double denom = arma::sum(exp_probs);

        F.row(i) = exp_probs / denom;
    }

    double loglik = 0.0;
    for (int i = 0; i < N; ++i) {
        arma::rowvec log_prob = F_raw.row(i) + log_pi_g.t();
        double max_log_prob = log_prob.max();
        arma::rowvec stabilized_log_prob = log_prob - max_log_prob;
        double row_logsumexp = max_log_prob + std::log(arma::sum(arma::exp(stabilized_log_prob)));
        loglik += row_logsumexp;
    }

    arma::mat z = F;

    return Rcpp::List::create(
        Rcpp::Named("F_raw") = F_raw,
        Rcpp::Named("F") = F,
        Rcpp::Named("loglik") = loglik,
        Rcpp::Named("z") = z
    );
}

// [[Rcpp::export]]
Rcpp::List aitkens_accel(int it, arma::vec loglik, arma::vec a_loglik) {
    int checks = 0;
    std::string exit_code = "";
    it = it - 1;

    if (it >= (int)loglik.n_elem) {
        loglik.resize(it + 1);
    }
    if (it >= (int)a_loglik.n_elem) {
        a_loglik.resize(it + 1);
    }


    if (it > 3) {
        if ((loglik[it - 1] - loglik[it - 2]) == 0) {
            checks = 1;
            exit_code = "Log Likelihood equal for two iterations";
        }
        else {
            double a = (loglik[it] - loglik[it - 1]) / (loglik[it - 1] - loglik[it - 2]);
            
            double add_to = (1 / (1 - a) * (loglik[it] - loglik[it - 1]));

            a_loglik[it] = loglik[it - 1] + add_to;

            if (std::abs(a_loglik[it] - loglik[it - 1]) < 0.01) {
                checks = 1;
                exit_code = "Aitken's acceleration converged";
            }
        }
        return List::create(Named("checks") = checks,
            Named("exit_code") = exit_code,
            Named("a_loglik") = a_loglik);
    } else {
        return List::create(Named("checks") = 0,
            Named("exit_code") = "Iteration index too small",
            Named("a_loglik") = a_loglik);
    }
}

