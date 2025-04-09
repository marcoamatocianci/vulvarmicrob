import numpy as np
import pandas as pd
import argparse
import warnings
#import h5py
import os
from typing import Union
# Import functions from  SPARCC code
from core_methods import to_fractions, compute_correlation_pvalues
from compositional_methods import variation_mat
def basis_var(Var_mat, M, V_min=1e-4):
    """Estimate the basis variances of compositional data."""
    V_vec = Var_mat.sum(axis=1)
    V_base = np.linalg.solve(M, V_vec)
    basis_variance = np.where(V_base <= 0, V_min, V_base)
    return basis_variance

def C_from_V(Var_mat, V_base):
    """Compute the basis correlation and covariance matrices."""
    Vi, Vj = np.meshgrid(V_base, V_base)
    Cov_base = 0.5 * (Vi + Vj - Var_mat)
    C_base = Cov_base / np.sqrt(Vi) / np.sqrt(Vj)
    return C_base, Cov_base

def new_excluded_pair(C, previously_excluded=[], th=0.1):
    """Find the next component pair to exclude based on the highest correlation."""
    C_temp = np.triu(np.abs(C), 1).copy()  # Use upper triangle, exclude diagonal
    if len(previously_excluded) > 0:
        C_temp[tuple(zip(*previously_excluded))] = 0
    a = np.unravel_index(np.argmax(C_temp), C_temp.shape)
    cmax = C_temp[a]
    if cmax > th:
        return a
    else:
        return None

def run_sparcc(frame, th=0.1, x_iter=10):
    Var_mat = variation_mat(frame)
    Var_mat_temp = Var_mat.copy()
    D = frame.shape[1]
    M = np.ones((D, D)) + np.diag([D - 2] * D)
    
    V_base = basis_var(Var_mat_temp, M)
    C_base, Cov_base = C_from_V(Var_mat, V_base)
    
    excluded_pairs = []
    excluded_comp = np.array([])

    for xi in range(x_iter):
        to_exclude = new_excluded_pair(C=C_base, th=th, previously_excluded=excluded_pairs)
        if to_exclude is None:
            break
        excluded_pairs.append(to_exclude)
        i, j = to_exclude
        M[i, j] -= 1
        M[j, i] -= 1
        M[i, i] -= 1
        M[j, j] -= 1
        Var_mat_temp[i, j] = 0
        Var_mat_temp[j, i] = 0

        nexcluded = np.bincount(np.ravel(excluded_pairs))
        excluded_comp_prev = set(excluded_comp.copy())
        excluded_comp = np.where(nexcluded >= D - 3)[0]
        excluded_comp_new = set(excluded_comp) - excluded_comp_prev

        if len(excluded_comp_new) > 0:
            if len(excluded_comp) > D - 4:
                warnings.warn("Too many components excluded. Returning CLR result.")
                return run_clr(frame)

            for xcomp in excluded_comp_new:
                Var_mat_temp[xcomp, :] = 0
                Var_mat_temp[:, xcomp] = 0
                M[xcomp, :] = 0
                M[:, xcomp] = 0
                M[xcomp, xcomp] = 1

        V_base = basis_var(Var_mat_temp, M)
        C_base, Cov_base = C_from_V(Var_mat, V_base)

        for xcomp in excluded_comp:
            V_base[xcomp] = np.nan
            C_base[xcomp, :] = np.nan
            C_base[:, xcomp] = np.nan
            Cov_base[xcomp, :] = np.nan
            Cov_base[:, xcomp] = np.nan
    np.fill_diagonal(C_base, 1.0)
    np.fill_diagonal(Cov_base, V_base)
    return C_base, Cov_base

def basic_corr(frame, method='sparcc', th=0.1, x_iter=10):
    method = method.lower()
    if method == 'clr':
        return run_clr(frame)
    elif method == 'sparcc':
        return run_sparcc(frame, th=th, x_iter=x_iter)
    else:
        raise ValueError(f"Unsupported method: {method}")

def permute_w_replacement(frame: Union[pd.DataFrame, np.ndarray], axis=0):
    '''
    Generates a bootstrap sample by resampling rows or columns without replacement.
    '''
    if isinstance(frame, pd.DataFrame):
        frame = frame.values

    rows, cols = frame.shape

    if axis == 0:  # Bootstrap rows
        bootstrap_indices = np.random.choice(rows, size=rows, replace=False)
        bootstrap_sample = frame[bootstrap_indices, :]
        return bootstrap_sample

    elif axis == 1:  # Bootstrap columns
        bootstrap_indices = np.random.choice(cols, size=cols, replace=False)
        bootstrap_sample = frame[:, bootstrap_indices]
        return bootstrap_sample

    else:
        raise ValueError("Axis must be 0 (rows) or 1 (columns).")


### Function to Compute Pseudo P-values for Correlation and Covariance
def compute_pseudo_pvals(frame, original_cor, original_cov, n_bootstrap=1000, method='sparcc', th=0.1, x_iter=10, test_type='two-sided'):
    """Compute pseudo p-values for correlation and covariance matrices using bootstrapping."""
    
   # # Step 1: Compute the original correlation and covariance matrices
  #  original_cor, original_cov = basic_corr(frame, method=method, th=th, x_iter=x_iter)
    
    n_components = frame.shape[1]
    
    # Initialize arrays to store bootstrap correlation and covariance matrices
    cor_bootstrap_samples = np.zeros((n_bootstrap, n_components, n_components))
    cov_bootstrap_samples = np.zeros((n_bootstrap, n_components, n_components))
    
    # Step 2: Perform bootstrapping
    for i in range(n_bootstrap):
        # Generate a bootstrap sample (resample rows with replacement)
        bootstrap_sample = permute_w_replacement(frame, axis=0)
        #bootstrap_sample = bootstrap_with_replacement(frame, axis=0)
        # Compute the correlation and covariance matrices for the bootstrap sample
        cor_bootstrap, cov_bootstrap = basic_corr(bootstrap_sample, method=method, th=th, x_iter=x_iter)
        cor_bootstrap = np.nan_to_num(cor_bootstrap)
        cov_bootstrap = np.nan_to_num(cov_bootstrap)
        # Store the bootstrap results
        cor_bootstrap_samples[i] = cor_bootstrap
        cov_bootstrap_samples[i] = cov_bootstrap
    
    # Step 3: Compute pseudo p-values
    pseudo_pvals_cor = np.zeros((n_components, n_components))
    pseudo_pvals_cov = np.zeros((n_components, n_components))
    
    for i in range(n_components):
        for j in range(n_components):
            # One-sided or two-sided test for correlation
            if test_type == 'two-sided':
                pseudo_pvals_cor[i, j] = np.nanmean(np.abs(cor_bootstrap_samples[:, i, j]) >= np.abs(original_cor[i, j]))
            elif test_type == 'right-tailed':
                pseudo_pvals_cor[i, j] = np.nanmean(cor_bootstrap_samples[:, i, j] >= original_cor[i, j])
            elif test_type == 'left-tailed':
                pseudo_pvals_cor[i, j] = np.nanmean(cor_bootstrap_samples[:, i, j] <= original_cor[i, j])
            else:
                raise ValueError("Invalid test_type. Choose from 'two-sided', 'right-tailed', or 'left-tailed'.")
            #print(pseudo_pvals_cor[i, j])
            # One-sided or two-sided test for covariance
            if test_type == 'two-sided':
                pseudo_pvals_cov[i, j] = np.nanmean(np.abs(cov_bootstrap_samples[:, i, j]) >= np.abs(original_cov[i, j]))
            elif test_type == 'right-tailed':
                pseudo_pvals_cov[i, j] = np.nanmean(cov_bootstrap_samples[:, i, j] >= original_cov[i, j])
            elif test_type == 'left-tailed':
                pseudo_pvals_cov[i, j] = np.nanmean(cov_bootstrap_samples[:, i, j] <= original_cov[i, j])
            else:
                raise ValueError("Invalid test_type. Choose from 'two-sided', 'right-tailed', or 'left-tailed'.")
    
    return pseudo_pvals_cor, pseudo_pvals_cov

def main_alg(frame, method='sparcc', th=0.1, x_iter=10, n_iter=20, output_dir='./', col_names=None, n_bootstrap=10, test_type="right-tailed",pval_method='pseudo'):
    # Step 1: Filter out zero-variance columns
    non_zero_var_mask = np.var(frame, axis=0) > 0
    frame_filtered = frame[:, non_zero_var_mask]
    # dirlichlet transformation
    fracs = to_fractions(frame_filtered)
    # clr transformation
    #fracs = run_clr(frame_filtered)
    # Step 2: Keep track of the filtered column names
    col_names_filtered = col_names[non_zero_var_mask]

    # Step 3: Run SparCC
    results = []
    for i in range(n_iter):
        print(f'Running iteration {i + 1}/{n_iter}')
        cor_sparse, cov_sparse = basic_corr(fracs, method=method, th=th, x_iter=x_iter)
        results.append((cor_sparse, cov_sparse))

    # Step 4: Average the results
    avg_cor = np.nanmean([res[0] for res in results], axis=0)
    avg_cov = np.nanmean([res[1] for res in results], axis=0)

    # Step 5: Convert the correlation and covariance matrices back to DataFrame
    cor_df = pd.DataFrame(avg_cor, index=col_names_filtered, columns=col_names_filtered)
    cov_df = pd.DataFrame(avg_cov, index=col_names_filtered, columns=col_names_filtered)
    if pval_method == 'pseudo':
        # Step 6: Compute pseudo p-values using bootstrapping
        pseudo_pvals_cor, pseudo_pvals_cov = compute_pseudo_pvals(fracs,avg_cor, avg_cov, n_bootstrap=n_bootstrap, method=method, th=th, x_iter=x_iter, test_type=test_type)
        
        # Convert pseudo p-values to DataFrame
        pseudo_pvals_cor_df = pd.DataFrame(pseudo_pvals_cor, index=col_names_filtered, columns=col_names_filtered)
        pseudo_pvals_cov_df = pd.DataFrame(pseudo_pvals_cov, index=col_names_filtered, columns=col_names_filtered)
    else:
        pearson_pvals = compute_correlation_pvalues(pd.DataFrame(avg_cor, index=col_names_filtered, columns=col_names_filtered))
        
        # Convert Pearson p-values to DataFrame
        pearson_pvals_df = pd.DataFrame(index=col_names_filtered, columns=col_names_filtered)
        for pair, (r, p_value) in pearson_pvals.items():
            pearson_pvals_df.loc[pair[0], pair[1]] = p_value
    # Step 7: Save the results
    os.makedirs(output_dir, exist_ok=True)
    cor_file = os.path.join(output_dir, 'sparcc_correlation_matrix.csv')
    cov_file = os.path.join(output_dir, 'sparcc_covariance_matrix.csv')
    if pval_method == 'pseudo':
        pvals_cor_file = os.path.join(output_dir, 'pseudo_pvals_correlation.csv')
        pvals_cov_file = os.path.join(output_dir, 'pseudo_pvals_covariance.csv')
        pseudo_pvals_cor_df.to_csv(pvals_cor_file, index=True)
        pseudo_pvals_cov_df.to_csv(pvals_cov_file, index=True)
    else:
        pvals_cor_file = os.path.join(output_dir, 'pearson_pvals_correlation.csv')
        pvals_cov_file = os.path.join(output_dir, 'pearson_pvals_covariance.csv')
        pearson_pvals_df.to_csv(pvals_cor_file, index=True)
        pearson_pvals_df.to_csv(pvals_cov_file, index=True)
    cor_df.to_csv(cor_file, index=True)
    cov_df.to_csv(cov_file, index=True)
    print(f"Correlation matrix saved to: {cor_file}")
    print(f"Covariance matrix saved to: {cov_file}")
    if pval_method == 'pseudo':
        print(f"Pseudo p-values for correlation saved to: {pvals_cor_file}")
        print(f"Pseudo p-values for covariance saved to: {pvals_cov_file}")
        #return avg_cor, avg_cov, pseudo_pvals_cor, pvals_cov_file
    else:   
        print(f"Pearson p-values for correlation saved to: {pvals_cor_file}")
        print(f"Pearson p-values for covariance saved to: {pvals_cov_file}")
        #return avg_cor, avg_cov, pearson_pvals_df, pvals_cov_file
    

# Argument parsing
def parse_args():
    parser = argparse.ArgumentParser(description="Run SPARCC on compositional data.")
    
    # Input file
    parser.add_argument("input_file", type=str, help="Path to the input CSV file containing compositional data (rows: samples; columns: components).")
    
    # Output directory
    parser.add_argument("-o", "--output_dir", type=str, default="./output", help="Directory to save the output correlation and covariance matrices (default: ./output).")
    
    # Number of iterations
    parser.add_argument("-n", "--n_iter", type=int, default=20, help="Number of iterations for SparCC (default: 20).")
    
    # Exclusion threshold
    parser.add_argument("-t", "--threshold", type=float, default=0.1, help="Exclusion threshold for SparCC (default: 0.1).")
    
    # Exclusion iterations
    parser.add_argument("-x", "--x_iter", type=int, default=10, help="Number of exclusion iterations for SparCC (default: 10).")

    # Number of bootstrap samples
    parser.add_argument("-b", "--n_bootstrap", type=int, default=1000, help="Number of bootstrap samples for pseudo p-value calculation (default: 1000).")

    parser.add_argument("-m", "--method", type=str, default="right-tailed", help="Test type of the bootstrap statistic (default: right-tailed).")
    
    parser.add_argument("-p", "--pval_method", type=str, choices=['pseudo', 'pearson'], default='pseudo', 
                        help="Method to compute p-values: 'pseudo' for pseudo p-values via bootstrapping, 'pearson' for Pearson correlation p-values (default: pseudo).")
    return parser.parse_args()

# Main function to execute SPARCC
def main():
    # Parse command-line arguments
    args = parse_args()
    # set a random seed
    np.random.seed(42)
    # Read the input data
    print(f"Reading input file: {args.input_file}")
    data = pd.read_csv(args.input_file,sep="\t", index_col=0)
    # Save row and column names
    data = data.T
    col_names = data.columns
    # Convert to NumPy array
    data_np = data.to_numpy()

    # Run SPARCC
    main_alg(data_np,
             method='sparcc',
             th=args.threshold,
             x_iter=args.x_iter,
             n_iter=args.n_iter,
             output_dir=args.output_dir,
             col_names=col_names,
             n_bootstrap=args.n_bootstrap,
             test_type=args.method,
             pval_method=args.pval_method)

if __name__ == "__main__":
    main()
