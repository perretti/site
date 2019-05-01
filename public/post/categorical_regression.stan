data {
  int n_obs;
  int n_cat;
  int n_cov;
  int outcome[n_obs];
  matrix[n_obs, n_cov] X;
}

transformed data {
  vector[n_cov] zeros = rep_vector(0, n_cov);
}

parameters {
  matrix[n_cov, n_cat - 1] B_sub;
}

transformed parameters {
  matrix[n_cov, n_cat] B = append_col(zeros, B_sub);
}

model {
  matrix[n_obs, n_cat] score = X * B;
  
  target += normal_lpdf(to_vector(B_sub) | 0, 10); // prior on B_sub
  
  for(i in 1:n_obs) {
    vector[n_cat] p = softmax(to_vector(score[i,]));
    target += categorical_lpmf(outcome[i] | p);
  }
}

generated quantities {
  matrix[n_obs, n_cat] p_est;
  for(i in 1:n_obs) {
    p_est[i,] = to_row_vector(softmax(to_vector(X[i,] * B)));
  }
}

