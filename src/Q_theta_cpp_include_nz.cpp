#include <Rcpp.h>
using namespace Rcpp;


double expit(double x) {
  double res=exp(x)/(1+exp(x));
  return res;
}

NumericVector expit(NumericVector x) {
  NumericVector res=exp(x)/(1+exp(x));
  return res;
}

double beta(double a,double b) {
  NumericVector a_vec(1);
  NumericVector b_vec(1);
  a_vec[0]=a;
  b_vec[0]=b;
  double res=beta(a_vec,b_vec)[0];
  return res;

}

NumericVector mmultcpp(NumericMatrix mm,NumericVector vv){
  NumericVector out(mm.nrow());
  NumericVector rm1;
  for (int i=0; i<mm.nrow();i++) {
    rm1=mm(i,_);
    out[i]=std::inner_product(rm1.begin(),rm1.end(),vv.begin(),0.0);
  }
  return out;
}


double expect_M_x(NumericVector psi, double delta_1, NumericVector mu) {
  double res=sum(psi*(1-expit(delta_1))*expit(mu));
  return(res);
}

NumericVector test_inf(NumericVector x, NumericVector tau) {
  NumericVector res_vec=clone(x);
  LogicalVector ind_vec=is_infinite(x) & abs(tau)<1e-8;
  res_vec[ind_vec]=0;
  return res_vec;
}

NumericVector NaN_to_0_vec(NumericVector x) {
  LogicalVector nan_ind=is_nan(x);
  NumericVector res=clone(x);
  res[nan_ind]=0;
  return res;
}


NumericVector impute_small_num(NumericVector x, double del=1e-100) {
  LogicalVector too_small_ind= x<del;
  NumericVector res=clone(x);
  res[too_small_ind]=del;
  return res;
}

NumericVector li_1_func(NumericVector para_vec, double yi, double m_star,
                        double x_i, NumericVector confound_vec) {

  int num_conf=confound_vec.size();
  int k=(para_vec.size()-9-3*num_conf)/3;
  double beta_0=para_vec[0];
  double beta_1=para_vec[1];
  double beta_2=para_vec[2];
  double beta_3=para_vec[3];
  double beta_4=para_vec[4];
  double beta_5=para_vec[5];

  double gamma_0=para_vec[6];
  double gamma_1=para_vec[7];
  double phi=para_vec[8];
  double delta=para_vec[9];

  NumericVector alpha_0_vec=para_vec[seq(10,9+k)];
  NumericVector alpha_1_vec=para_vec[seq(10+k,9+2*k)];
  NumericVector psi_vec=para_vec[seq(10+2*k,8+3*k)];

  NumericVector beta_conf=para_vec[seq(9+3*k,8+3*k+num_conf)];
  NumericVector alpha_conf=para_vec[seq(9+3*k+num_conf,8+3*k+2*num_conf)];
  NumericVector gamma_conf=para_vec[seq(9+3*k+2*num_conf,8+3*k+3*num_conf)];


  psi_vec.push_back(1-sum(para_vec[seq(10+2*k,8+3*k)]));


  NumericVector mu=expit(alpha_0_vec+alpha_1_vec*x_i+sum(confound_vec*alpha_conf));
  double delta_i=gamma_0+gamma_1*x_i+sum(confound_vec*gamma_conf);

  double norm_part=pow(yi-beta_0-beta_1*m_star-beta_2-(beta_3+beta_4)*x_i-
                       beta_5*x_i*m_star-sum(confound_vec*beta_conf),2)/(2*delta*delta);
  NumericVector li_1_value= -0.5*log(2*M_PI)-log(delta)-norm_part+
    (mu*phi-1)*log(m_star)+((1-mu)*phi-1)*log(1-m_star)-log(beta(mu*phi,(1-mu)*phi));

  NumericVector beta_part=log(beta(mu*phi,(1-mu)*phi));
  NumericVector li_1_res;
  if(delta_i>200) {
    li_1_res=log(psi_vec) - delta_i + li_1_value;
  } else {
    li_1_res=log(psi_vec) - log(1+exp(delta_i)) + li_1_value;
  }


  // Rcout << " li_1_value " << li_1_value << " norm_part " << norm_part << " mu " << mu <<
  //   " beta_part " << beta_part << " psi_vec " << psi_vec <<
  //     "delta_i " << delta_i << "\n";
  return li_1_res;
}

double li0_2_func(NumericVector para_vec,double yi,double x_i,NumericVector confound_vec) {
  int num_conf=confound_vec.size();
  int k=(para_vec.size()-9-3*num_conf)/3;

  double beta_0=para_vec[0];
  // double beta_1=para_vec[1];
  // double beta_2=para_vec[2];
  double beta_3=para_vec[3];
  // double beta_4=para_vec[4];
  // double beta_5=para_vec[5];

  double gamma_0=para_vec[6];
  double gamma_1=para_vec[7];
  // double phi=para_vec[8];
  double delta=para_vec[9];

  NumericVector beta_conf=para_vec[seq(9+3*k,8+3*k+num_conf)];
  // NumericVector alpha_conf=para_vec[seq(9+3*k+num_conf,8+3*k+2*num_conf)];
  NumericVector gamma_conf=para_vec[seq(9+3*k+2*num_conf,8+3*k+3*num_conf)];


  // NumericVector alpha_0_vec=para_vec[seq(10,9+k)];
  // NumericVector alpha_1_vec=para_vec[seq(10+k,9+2*k)];
  // NumericVector psi_vec(k);
  // psi_vec[seq(0,k-2)]=para_vec[seq(10+2*k,8+3*k)];
  // psi_vec[k-1]=1-sum(para_vec[seq(10+2*k,8+3*k)]);

  double delta_i=gamma_0+gamma_1*x_i+sum(confound_vec*gamma_conf);

  double norm_part=pow(yi-beta_0-beta_3*x_i-sum(confound_vec*beta_conf),2)/(2*delta*delta);
  double li0_2_value= -0.5*log(2*M_PI)-log(delta)-norm_part;
  double return_value;
  if(delta_i>200) {
    return_value=0;
  } else {
    return_value=delta_i-log(1+exp(delta_i));
  }


  return return_value + li0_2_value;

}



NumericVector li_2_func(NumericVector para_vec,double yi,double x_i,double l_i,
                        NumericVector confound_vec) {
  int num_conf=confound_vec.size();
  int k=(para_vec.size()-9-3*num_conf)/3;


  double beta_0=para_vec[0];
  double beta_1=para_vec[1];
  double beta_2=para_vec[2];
  double beta_3=para_vec[3];
  double beta_4=para_vec[4];
  double beta_5=para_vec[5];

  double gamma_0=para_vec[6];
  double gamma_1=para_vec[7];
  double phi=para_vec[8];
  double delta=para_vec[9];

  NumericVector alpha_0_vec=para_vec[seq(10,9+k)];
  NumericVector alpha_1_vec=para_vec[seq(10+k,9+2*k)];
  NumericVector psi_vec=para_vec[seq(10+2*k,8+3*k)];

  NumericVector beta_conf=para_vec[seq(9+3*k,8+3*k+num_conf)];
  NumericVector alpha_conf=para_vec[seq(9+3*k+num_conf,8+3*k+2*num_conf)];
  NumericVector gamma_conf=para_vec[seq(9+3*k+2*num_conf,8+3*k+3*num_conf)];

  psi_vec.push_back(1-sum(para_vec[seq(10+2*k,8+3*k)]));

  NumericVector mu=expit(alpha_0_vec+alpha_1_vec*x_i+sum(confound_vec*alpha_conf));
  double delta_i=gamma_0+gamma_1*x_i+sum(confound_vec*gamma_conf);
  double int_part_0= -pow(yi-beta_0-beta_2-(beta_3+beta_4)*x_i-
                          sum(confound_vec*beta_conf),2)/(2*delta*delta);
  double int_part_li=-pow(yi-beta_0-beta_1/l_i-beta_2-(beta_3+beta_4)*x_i-
                          beta_5*x_i/l_i-sum(confound_vec*beta_conf),2)/(2*delta*delta);

  double const_part=(int_part_0+int_part_li)/2;
  NumericVector int_value(k);
  NumericVector li_to_vec(1);
  li_to_vec[0]=1/l_i;

  for (int i=0;i<k;i++) {
    int_value[i]=Rcpp::pbeta(li_to_vec,mu[i]*phi,(1-mu[i])*phi)[0];
  }


  NumericVector li1_2_value= -0.5*log(2*M_PI) - log(delta) + const_part + log(int_value);
  // NumericVector Psi1=(1-delta_i)*psi_vec;
  // NumericVector li_2_res=log(Psi1) + li1_2_value;

  NumericVector li_2_res;
  if(delta_i>200) {
    li_2_res=log(psi_vec) - delta_i + li1_2_value;
  } else {
    li_2_res=log(psi_vec) - log(1+exp(delta_i)) + li1_2_value;
  }


  return li_2_res;
}




NumericVector tau_1_func(NumericVector para_vec, double yi, double m_star, double x_i,
                         NumericVector confound_vec) {
  int num_conf=confound_vec.size();
  int k=(para_vec.size()-9-3*num_conf)/3;


  // double beta_0=para_vec[0];
  // double beta_1=para_vec[1];
  // double beta_2=para_vec[2];
  // double beta_3=para_vec[3];
  // double beta_4=para_vec[4];
  // double beta_5=para_vec[5];

  // double gamma_0=para_vec[6];
  // double gamma_1=para_vec[7];
  double phi=para_vec[8];
  // double delta=para_vec[9];

  NumericVector alpha_0_vec=para_vec[seq(10,9+k)];
  NumericVector alpha_1_vec=para_vec[seq(10+k,9+2*k)];
  NumericVector psi_vec=para_vec[seq(10+2*k,8+3*k)];

  // NumericVector beta_conf=para_vec[seq(9+3*k,8+3*k+num_conf)];
  NumericVector alpha_conf=para_vec[seq(9+3*k+num_conf,8+3*k+2*num_conf)];
  // NumericVector gamma_conf=para_vec[seq(9+3*k+2*num_conf,8+3*k+3*num_conf)];

  psi_vec.push_back(1-sum(para_vec[seq(10+2*k,8+3*k)]));

  NumericVector mu=expit(alpha_0_vec+alpha_1_vec*x_i+sum(confound_vec*alpha_conf));
  // double delta_i=expit(gamma_0+gamma_1*x_i);

  NumericVector dbeta_value(k);
  NumericVector m_star_to_vec(1);
  m_star_to_vec[0]=m_star;

  for (int i=0;i<k;i++) {
    dbeta_value[i]=Rcpp::dbeta(m_star_to_vec,mu[i]*phi,(1-mu[i])*phi)[0];
  }


  NumericVector res=psi_vec*impute_small_num(dbeta_value)/
    (sum(psi_vec*impute_small_num(dbeta_value)));

  return res;
}

NumericVector tau_2_func(NumericVector para_vec, double yi, double x_i, double l_i,
                         NumericVector confound_vec) {
  double li0_2_value=li0_2_func(para_vec,yi,x_i,confound_vec);
  NumericVector li_2_value=li_2_func(para_vec,yi,x_i,l_i,confound_vec);
  li_2_value.push_front(li0_2_value);
  NumericVector li_diff(li_2_value.length());
  for (int i=0;i<li_2_value.length();i++) {
    li_diff[i]=1/sum(exp(NaN_to_0_vec(li_2_value-li_2_value[i])));
  }
  return(li_diff);
}

double li_1_1taxon_func(NumericVector para_vec, double yi, double m_star, double x_i,
                        NumericVector confound_vec) {
  int num_conf=confound_vec.size();

  double beta_0=para_vec[0];
  double beta_1=para_vec[1];
  double beta_2=para_vec[2];
  double beta_3=para_vec[3];
  double beta_4=para_vec[4];
  double beta_5=para_vec[5];

  double gamma_0=para_vec[6];
  double gamma_1=para_vec[7];
  double phi=para_vec[8];
  double delta=para_vec[9];

  double alpha_0=para_vec[10];
  double alpha_1=para_vec[11];

  NumericVector beta_conf=para_vec[seq(12,11+num_conf)];
  NumericVector alpha_conf=para_vec[seq(12+num_conf,11+2*num_conf)];
  NumericVector gamma_conf=para_vec[seq(12+2*num_conf,11+3*num_conf)];


  // NumericVector alpha_0_vec=para_vec[seq(10,9+k)];
  // NumericVector alpha_1_vec=para_vec[seq(10+k,9+2*k)];
  // NumericVector psi_vec=para_vec[seq(10+2*k,8+3*k)];
  //
  // psi_vec.push_back(1-sum(para_vec[seq(10+2*k,8+3*k)]));


  double mu=expit(alpha_0+alpha_1*x_i+sum(confound_vec*alpha_conf));
  double delta_i=gamma_0+gamma_1*x_i+sum(confound_vec*gamma_conf);
  double norm_part=pow(yi-beta_0-beta_1*m_star-beta_2-(beta_3+beta_4)*x_i-beta_5*x_i*m_star-
                       sum(confound_vec*beta_conf),2)/(2*delta*delta);
  double li_1_value;
  if(delta_i>200) {
    li_1_value= -0.5*log(2*M_PI)-log(delta)-norm_part-delta_i-log(beta(mu*phi,(1-mu)*phi))+
      (mu*phi-1)*log(m_star)+((1-mu)*phi-1)*log(1-m_star);
  } else {
    li_1_value= -0.5*log(2*M_PI)-log(delta)-norm_part-log(1+exp(delta_i))-log(beta(mu*phi,(1-mu)*phi))+
      (mu*phi-1)*log(m_star)+((1-mu)*phi-1)*log(1-m_star);
  }
  // double li_1_value= -0.5*log(2*M_PI)-log(delta)-norm_part-log(1+exp(delta_i))-log(beta(mu*phi,(1-mu)*phi))+
  //   (mu*phi-1)*log(m_star)+((1-mu)*phi-1)*log(1-m_star);


  // NumericVector beta_part=log(beta(mu*phi,(1-mu)*phi));
  //
  // NumericVector Psi1=(1-delta_i)*psi_vec;
  // NumericVector li_1_res=log(Psi1) + li_1_value;

  // Rcout << " li_1_value " << li_1_value << " norm_part " << norm_part << " mu " << mu <<
  //   " beta_part " << beta_part << " psi_vec " << psi_vec << " PSI " << Psi1 << "\n";
  return li_1_value;
}

double li_2_1taxon_func(NumericVector para_vec,double yi,double x_i,double l_i,
                        NumericVector confound_vec) {
  int num_conf=confound_vec.size();

  double beta_0=para_vec[0];
  double beta_1=para_vec[1];
  double beta_2=para_vec[2];
  double beta_3=para_vec[3];
  double beta_4=para_vec[4];
  double beta_5=para_vec[5];

  double gamma_0=para_vec[6];
  double gamma_1=para_vec[7];
  double phi=para_vec[8];
  double delta=para_vec[9];

  double alpha_0=para_vec[10];
  double alpha_1=para_vec[11];

  NumericVector beta_conf=para_vec[seq(12,11+num_conf)];
  NumericVector alpha_conf=para_vec[seq(12+num_conf,11+2*num_conf)];
  NumericVector gamma_conf=para_vec[seq(12+2*num_conf,11+3*num_conf)];

  double mu=expit(alpha_0+alpha_1*x_i+sum(confound_vec*alpha_conf));
  double delta_i=expit(gamma_0+gamma_1*x_i+sum(confound_vec*gamma_conf));
  double norm_part= -pow(yi-beta_0-beta_3*x_i-sum(confound_vec*beta_conf),2)/(2*delta*delta);

  double int_part_0= exp(-pow(yi-beta_0-beta_2-(beta_3+beta_4)*x_i-
                         sum(confound_vec*beta_conf),2)/(2*delta*delta));
  double int_part_li= exp(-pow(yi-beta_0-beta_1/l_i-beta_2-(beta_3+beta_4)*x_i-
                          beta_5*x_i/l_i-sum(confound_vec*beta_conf),2)/(2*delta*delta));

  double const_part=(int_part_0+int_part_li)/2;
  double int_value;
  NumericVector li_to_vec(1);
  li_to_vec[0]=1/l_i;

  int_value=Rcpp::pbeta(li_to_vec,mu*phi,(1-mu)*phi)[0];
  double li1_2_value= -0.5*log(2*M_PI) - log(delta) +
    log(delta_i*exp(norm_part)+(1-delta_i)*const_part*int_value);

  return li1_2_value;
}

// [[Rcpp::export]]
double Q_theta_cpp(NumericVector para_vec_ori,NumericVector para_vec_0_ori,NumericVector yi_vec,
                   NumericVector m_star_vec,NumericVector x_i_vec,NumericVector l_i_vec,
                   NumericMatrix confound_mat, bool x4_inter, bool x5_inter) {


  // Rcout << "step1";
  // para_vec[5]=0;
  // para_vec_0[5]=0;
  //
  // para_vec[0]=0;
  // para_vec_0[0]=0;
  NumericVector para_vec=clone(para_vec_ori);
  NumericVector para_vec_0=clone(para_vec_0_ori);
  
  if (!x4_inter) {
    para_vec[4] = 0;
  }
  
  if (!x5_inter) {
    para_vec[5] = 0;
  }

  double yi;
  double x_i;
  double m_star;
  double l_i;
  NumericVector li_1_value;
  NumericVector tau1_value;
  NumericVector li_2_value;
  double li0_2_value;
  NumericVector tau_2_value;
  double res;
  NumericVector res_vec(yi_vec.length());
  double res_sum;

  for (int i=0;i<yi_vec.length();i++) {
    yi=yi_vec[i];
    m_star=m_star_vec[i];
    x_i=x_i_vec[i];
    l_i=l_i_vec[i];
    NumericVector confound_vec=confound_mat(i,_);
    // Rcout << i <<"\n";
    if (m_star>0.99999) {
      m_star=0.99999;
    }

    if(m_star>1e-50) {
      li_1_value=li_1_func(para_vec,yi,m_star,x_i,confound_vec);
      tau1_value=tau_1_func(para_vec_0,yi,m_star,x_i,confound_vec);
      // Rcout << i << " " << "li_1_value " << li_1_value << " tau1_value " << tau1_value << "\n";

      res=sum(tau1_value* test_inf(li_1_value,tau1_value));

    } else if (m_star<1e-50) {
      li0_2_value=li0_2_func(para_vec,yi,x_i,confound_vec);
      li_2_value=li_2_func(para_vec,yi,x_i,l_i,confound_vec);
      li_2_value.push_front(li0_2_value);

      tau_2_value=tau_2_func(para_vec_0,yi,x_i,l_i,confound_vec);
      res=sum(tau_2_value*test_inf(li_2_value,tau_2_value));
    } else {Rcout << "M_star less 0";}

    // Rcout << i << " " << res << "\n";
    res_vec[i]=res;
  }

  // Rcout << "li_1_value "<<li_1_value<< " tau1_value " << tau1_value <<"\n";
  res_sum= -sum(res_vec);

  return res_sum;
}


// [[Rcpp::export]]
double Q_theta_cpp_nz(NumericVector para_vec_ori,NumericVector para_vec_0_ori,NumericVector yi_vec,
                      NumericVector m_star_vec,NumericVector x_i_vec,NumericVector l_i_vec,
                      NumericMatrix confound_mat, bool x4_inter, bool x5_inter) {

  // para_vec[5]=0;
  // para_vec_0[5]=0;
  //
  // para_vec[0]=0;
  // para_vec_0[0]=0;
  NumericVector para_vec=clone(para_vec_ori);
  NumericVector para_vec_0=clone(para_vec_0_ori);

  para_vec[2]=0;
  para_vec[4]=0;
  para_vec[6]= -100;
  para_vec[7]=0;
  
  if (!x4_inter) {
    para_vec[4] = 0;
  }
  
  if (!x5_inter) {
    para_vec[5] = 0;
  }

  double yi;
  double x_i;
  double m_star;
  double l_i;
  NumericVector li_1_value;
  NumericVector tau1_value;
  NumericVector li_2_value;
  double li0_2_value;
  NumericVector tau_2_value;
  double res;
  NumericVector res_vec(yi_vec.length());
  double res_sum;

  for (int i=0;i<yi_vec.length();i++) {
    yi=yi_vec[i];
    m_star=m_star_vec[i];
    x_i=x_i_vec[i];
    l_i=l_i_vec[i];
    NumericVector confound_vec=confound_mat(i,_);
    if (m_star>0.99999) {
      m_star=0.99999;
    }
    if(m_star>1e-50) {
      li_1_value=li_1_func(para_vec,yi,m_star,x_i,confound_vec);
      tau1_value=tau_1_func(para_vec_0,yi,m_star,x_i,confound_vec);
      res=sum(tau1_value* test_inf(li_1_value,tau1_value));

    } else if (m_star<1e-50) {
      li0_2_value=li0_2_func(para_vec,yi,x_i,confound_vec);
      li_2_value=li_2_func(para_vec,yi,x_i,l_i,confound_vec);
      li_2_value.push_front(li0_2_value);

      tau_2_value=tau_2_func(para_vec_0,yi,x_i,l_i,confound_vec);
      res=sum(tau_2_value*test_inf(li_2_value,tau_2_value));
    } else {Rcout << "M_star less 0";}
    res_vec[i]=res;
  }

  // Rcout << "li_1_value "<<li_1_value<< " tau1_value " << tau1_value <<"\n";
  res_sum= -sum(res_vec);

  return res_sum;
}



// [[Rcpp::export]]
double Q_theta_cpp_nomix(NumericVector para_vec_ori ,NumericVector yi_vec,
                         NumericVector m_star_vec,NumericVector x_i_vec,NumericVector l_i_vec,
                         NumericMatrix confound_mat, bool x4_inter, bool x5_inter) {

  NumericVector para_vec=clone(para_vec_ori);
  
  if (!x4_inter) {
    para_vec[4] = 0;
  }
  
  if (!x5_inter) {
    para_vec[5] = 0;
  }
  
  double yi;
  double x_i;
  double m_star;
  double l_i;
  double li_1_value;

  double li_2_value;

  double res;
  NumericVector res_vec(yi_vec.length());
  double res_sum;

  for (int i=0;i<yi_vec.length();i++) {
    yi=yi_vec[i];
    m_star=m_star_vec[i];
    x_i=x_i_vec[i];
    l_i=l_i_vec[i];
    NumericVector confound_vec=confound_mat(i,_);

    if (m_star>0.99999) {
      m_star=0.99999;
    }

    if(m_star>1e-50) {
      li_1_value=li_1_1taxon_func(para_vec,yi,m_star,x_i,confound_vec);
      res=li_1_value;

    } else if (m_star<1e-50) {
      li_2_value=li_2_1taxon_func(para_vec,yi,x_i,l_i,confound_vec);

      res=li_2_value;
    } else {Rcout << "M_star less 0";}
    res_vec[i]=res;
  }

  // Rcout << "res_vec "<<res_vec<<"\n";
  res_sum= -sum(res_vec);

  return res_sum;
}

// [[Rcpp::export]]
double Q_theta_cpp_nz_nomix(NumericVector para_vec_ori,NumericVector yi_vec,
                            NumericVector m_star_vec,NumericVector x_i_vec,NumericVector l_i_vec,
                            NumericMatrix confound_mat, bool x4_inter, bool x5_inter) {

  NumericVector para_vec=clone(para_vec_ori);

  para_vec[2]=0;
  para_vec[4]=0;
  para_vec[6]= -100;
  para_vec[7]=0;

  if (!x4_inter) {
    para_vec[4] = 0;
  }
  
  if (!x5_inter) {
    para_vec[5] = 0;
  }
  
  double yi;
  double x_i;
  double m_star;
  double l_i;
  double li_1_value;

  double li_2_value;

  double res;
  NumericVector res_vec(yi_vec.length());
  double res_sum;

  for (int i=0;i<yi_vec.length();i++) {
    yi=yi_vec[i];
    m_star=m_star_vec[i];
    x_i=x_i_vec[i];
    l_i=l_i_vec[i];
    NumericVector confound_vec=confound_mat(i,_);
    if (m_star>0.99999) {
      m_star=0.99999;
    }

    if(m_star>1e-50) {
      li_1_value=li_1_1taxon_func(para_vec,yi,m_star,x_i,confound_vec);
      res=li_1_value;

    } else if (m_star<1e-50) {
      li_2_value=li_2_1taxon_func(para_vec,yi,x_i,l_i,confound_vec);

      res=li_2_value;
    } else {Rcout << "M_star less 0";}
    res_vec[i]=res;
  }

  // Rcout << "li_1_value "<<li_1_value<< " tau1_value " << tau1_value <<"\n";
  res_sum= -sum(res_vec);

  return res_sum;
}
