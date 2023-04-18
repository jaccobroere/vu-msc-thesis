dyn.load("lib/Rfuncs.so")

linreg_path <- function(Y,X,val,idx,jdx,lambda_graph, 
						lambda_sparse,p,m,rho=1.0,
						eps_abs=1e-5,eps_rel=1e-4,
						max_num_iter=1e4,
						linearized_ADMM=FALSE, standard_ADMM=FALSE,
						nonadaptive_varying_rho=FALSE, constant_rho=FALSE,
						general_A=FALSE, graph_only=FALSE,
						graph_sparse=TRUE, reporting=FALSE)
{
	adaptive_varying_rho=TRUE
    if (nonadaptive_varying_rho && constant_rho) { stop("nonadaptive_varying_rho and constant_rho can't both be TRUE! \n") }
    if (constant_rho) { nonadaptive_varying_rho = FALSE; adaptive_varying_rho=FALSE }
	if (nonadaptive_varying_rho) { constant_rho = FALSE; adaptive_varying_rho=FALSE }
	special_ADMM = TRUE;
	if (standard_ADMM && linearized_ADMM) { stop("standard_ADMM and linearized_ADMM can't both be TRUE! \n") }
	if (standard_ADMM) { linearized_ADMM = FALSE; special_ADMM = FALSE }
	if (linearized_ADMM) { standard_ADMM = FALSE; special_ADMM = FALSE }	 
	
    out <- .Call('linreg',Y,t(X),val,
    as.integer(idx),as.integer(jdx),
	lambda_graph,lambda_sparse, 
	as.integer(p), as.integer(m),rho,
    eps_abs,eps_rel,as.integer(max_num_iter),
	 as.logical(special_ADMM), as.logical(linearized_ADMM), as.logical(standard_ADMM),
	 as.logical(adaptive_varying_rho),as.logical(nonadaptive_varying_rho),as.logical(constant_rho),
	 as.logical(general_A), as.logical(graph_only),
	 as.logical(graph_sparse),as.logical(reporting))
	 out$beta_path <- array(out$beta_path,c(p,length(lambda_sparse),length(lambda_graph)))
	 if (reporting) { out$fun_path <- array(out$fun_path,max_num_iter,c(length(lambda_sparse),length(lambda_graph))) }
	return (out)
}


############ lambda_sparse = gamma*lambda_graph ##########
linreg_path_v2 <- function(Y,X,val,idx,jdx,lambda_graph,gamma,p,m,rho=1.0,
                        eps_abs=1e-5,eps_rel=1e-4,
                        max_num_iter=1e4,
                        linearized_ADMM=FALSE, standard_ADMM=FALSE,
                        nonadaptive_varying_rho=FALSE, constant_rho=FALSE,
                        general_A=FALSE, graph_only=FALSE,
                        graph_sparse=TRUE, reporting=FALSE,diag_mat=NULL)
{
  adaptive_varying_rho=TRUE
  diag_mat_given = FALSE
  if (nonadaptive_varying_rho && constant_rho) { stop("nonadaptive_varying_rho and constant_rho can't both be TRUE! \n") }
  if (constant_rho) { nonadaptive_varying_rho = FALSE; adaptive_varying_rho=FALSE }
  if (nonadaptive_varying_rho) { constant_rho = FALSE; adaptive_varying_rho=FALSE }
  special_ADMM = TRUE;
  if (standard_ADMM && linearized_ADMM) { stop("standard_ADMM and linearized_ADMM can't both be TRUE! \n") }
  if (standard_ADMM) { linearized_ADMM = FALSE; special_ADMM = FALSE }
  if (linearized_ADMM) { standard_ADMM = FALSE; special_ADMM = FALSE }	 
  if (!is.null(diag_mat)) { diag_mat_given = TRUE } else { diag_mat = .0 }
  lambda_pairs <- rep(0,2*length(lambda_graph))
  lambda_pairs[seq(1,2*length(lambda_graph)-1,length.out=length(lambda_graph))] <- lambda_graph
  lambda_pairs[seq(2,2*length(lambda_graph),length.out=length(lambda_graph))] <- gamma*lambda_graph
  out <- .Call('linreg_v2',Y,t(X),val,
               as.integer(idx),as.integer(jdx),
               lambda_pairs,as.integer(p), as.integer(m),rho,
               eps_abs,eps_rel,as.integer(max_num_iter),as.double(diag_mat), as.logical(diag_mat_given),
               as.logical(special_ADMM), as.logical(linearized_ADMM), as.logical(standard_ADMM),
               as.logical(adaptive_varying_rho),as.logical(nonadaptive_varying_rho),as.logical(constant_rho),
               as.logical(general_A), as.logical(graph_only),
               as.logical(graph_sparse),as.logical(reporting))
  out$beta_path <- matrix(out$beta_path,p,length(lambda_pairs)/2)
  if (reporting) { out$fun_path <- matrix(out$fun_path,max_num_iter,length(lambda_pairs)/2) }
  return (out)
}


filtering <- function(Y, val, idx, jdx, lambda, m, rho=1.0, eps=1e-3,max_iter=1e5,
                        variantADMM=TRUE,varying_rho=TRUE,general_A=FALSE,reporting=FALSE)
{
    out <- .Call('filtering',Y,val,as.integer(idx),as.integer(jdx),lambda,as.integer(m),
                rho,eps,as.integer(max_iter),as.logical(variantADMM),as.logical(varying_rho),
                as.logical(general_A),as.logical(reporting))
    out$beta_path <- matrix(out$beta_path,p,length(lambda))
    if (reporting) { out$fun_path <- matrix(out$fun_path, max_iter,length(lambda)) }
    return (out)
}

# \hat{\beta}(\lambda) =
#       \argmin_\beta \|y - X \beta|_2^2 + \lambda * (\|D \beta\|_1 + gamma * \|\beta\|_1),
# or    \argmin_\beta \|y - X \beta|_2^2 + \lambda * \|D \beta\|_1 + lambda_sparse * \|\beta\|_1,
gen_lasso <- function(y,X,D,val,idx,jdx,lambda,gamma,lambda_sparse) {
    if (missing(y)) stop("y is missing.")
    if (!is.numeric(y)) stop("y must be numeric.")
    if (length(y) == 0) stop("There must be at least one data point [must have length(y) > 1].")
    if (missing(X)) stop("X is missing.")
    if (!is.matrix(X)) stop("X must be a matrix.")
    if (dim(X)[1] != length(y)) stop("The row dimension of X must agree with the length of the response y.")
    n = dim(X)[1]
    p = dim(X)[2]
    if (missing(D)) # D is sparse #
    {
        if (missing(val) || missing(idx) || missing(jdx))
        {
            stop("Row index 'idx', column index 'jdx', and entry values 'val' must all be provided when D is a sparse matrix")
        }
        m = max(idx)
    } else {
        if (!is.matrix(D)) { stop("D must be a matrix! ") }
        m = dim(D)[1]
    }
}



