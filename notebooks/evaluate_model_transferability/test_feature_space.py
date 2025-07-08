import pandas as pd
import os
import glob
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pickle
from scipy.spatial.distance import cdist


def load_tsvs(data_dir, suffix=''):
    """Load all .tsv files from directory and concatenate with sample_id column."""
    files = glob.glob(os.path.join(data_dir, "*.tsv"))
    dfs = []
    
    for file in files:
        df = pd.read_csv(file, sep='\t')
        df['chrm'] = df['chrm'].astype(str)
        df['pos'] = df['pos'].astype(int)
        sample_id = os.path.splitext(os.path.basename(file))[0].replace(suffix, '')
        df['sample_id'] = sample_id
        dfs.append(df)
    
    return pd.concat(dfs, ignore_index=True)


def load_training_data(pickle_path):
    """Load precomputed training data from pickle file."""
    with open(pickle_path, 'rb') as f:
        return pickle.load(f)


def plot_feature_distributions(training_merged, user_data, feature_cols):
    """
    Prepare data and plot violin plots for feature distributions.
    """
    # Split training data by label
    training_passing = training_merged[training_merged['assignment'] == 1][feature_cols + ['sample_id']].copy()
    training_passing['data_type'] = 'Training Passing'

    training_artifact = training_merged[training_merged['assignment'] == 0][feature_cols + ['sample_id']].copy()
    training_artifact['data_type'] = 'Training Artifact'

    training_unlabeled = training_merged[training_merged['assignment'] == -1][feature_cols + ['sample_id']].copy()
    training_unlabeled['data_type'] = 'Training Unlabeled'

    # Prepare user data
    user_features = user_data[feature_cols + ['sample_id']].copy()
    user_features['data_type'] = 'User Data'

    # Combine all data
    plot_data = pd.concat([training_passing, training_artifact, training_unlabeled, user_features], ignore_index=True)

    print(f"Data for plotting:")
    print(f"  Training passing: {len(training_passing)}")
    print(f"  Training artifact: {len(training_artifact)}")
    print(f"  Training unlabeled: {len(training_unlabeled)}")
    print(f"  User data: {len(user_features)}")

    # Create violin plots
    n_features = len(feature_cols)
    n_cols = 3
    n_rows = int(np.ceil(n_features / n_cols))

    plt.figure(figsize=(10, 2.5 * n_rows))

    colors = ['#440054', '#FF6600', '#888888', '#4292C6']  # Blue, Purple, Orange, Gray

    for i, feature in enumerate(feature_cols):
        ax = plt.subplot(n_rows, n_cols, i + 1)
        sns.violinplot(
            data=plot_data,
            y='data_type',
            x=feature,
            palette=colors,
            orient='h',
            cut=0,
            linewidth=0
        )
        if i % n_cols != 0:
            ax.set_yticklabels([])
        ax.set_xlabel(feature)
        ax.set_ylabel('')
        plt.xticks(fontsize=8)
        sns.despine()

    plt.tight_layout()
    plt.show()


def plot_umap_distributions(umap_results):
    """Plot UMAP distributions using precomputed embeddings."""

    def plot_umap_highlight(embedding_dict, highlight_category, colors):
        """Plot UMAP with one category highlighted and others in light gray"""
        
        # Plot all categories in light gray first
        for category, (embedding, label) in embedding_dict.items():
            if category != highlight_category:
                plt.scatter(embedding[:, 0], embedding[:, 1], 
                        c='lightgray', alpha=0.3, s=0.2, edgecolors=None, label='_nolegend_')
        
        # Plot highlighted category on top
        embedding, label = embedding_dict[highlight_category]
        plt.scatter(embedding[:, 0], embedding[:, 1], 
                c=colors[highlight_category], alpha=1, s=0.2, edgecolors=None, label=label)

        sns.despine()
        plt.xticks([])
        plt.yticks([])
        plt.xlabel('UMAP 1')
        plt.ylabel('UMAP 2')
        plt.title(label)

    # Prepare embedding dictionary using precomputed results
    lengths = umap_results['training_lengths']
    training_embedding = umap_results['training_embedding']
    user_embedding = umap_results['user_embedding']
    
    embedding_dict = {
        'pass': (training_embedding[:lengths['pass']], 'Training Pass'),
        'artifact': (training_embedding[lengths['pass']:lengths['pass']+lengths['artifact']], 'Training Artifact'),
        'unlabeled': (training_embedding[lengths['pass']+lengths['artifact']:], 'Training Unlabeled'),
        'user': (user_embedding, 'User Data')
    }

    colors = {
        'pass': '#440054',
        'artifact': '#FF6600', 
        'unlabeled': '#888888',
        'user': '#4292C6'
    }

    # Create 4 plots
    fig, axes = plt.subplots(2, 2, figsize=(10, 10))
    plt.subplots_adjust(hspace=0.3, wspace=0.3)

    for i, category in enumerate(['pass', 'artifact', 'unlabeled', 'user']):
        plt.subplot(2, 2, i+1)
        plot_umap_highlight(embedding_dict, category, colors)

    plt.show()


def compute_coverage_metrics(user_data, feature_cols, training_data_info):
    """
    Compute coverage metrics for user data against training data.
    
    Args:
        user_data: DataFrame with user features
        feature_cols: List of feature column names
        training_data_info: Dictionary with precomputed training data
        
    Returns:
        dict: Coverage analysis results including distances, threshold, and coverage percentage
    """
    # Clean user data and standardize
    user_features_clean = user_data[feature_cols].dropna()
    user_features_std = training_data_info['scaler'].transform(user_features_clean)
    
    # Calculate user-to-training distances
    user_to_training_dists = cdist(user_features_std, training_data_info['training_features_std']).min(axis=1)
    
    # Get precomputed threshold and training distances
    threshold = training_data_info['threshold_99']
    training_dists = training_data_info['min_training_dists']
    
    # Calculate coverage percentage
    coverage_pct = (user_to_training_dists <= threshold).mean() * 100
    
    return {
        'user_to_training_dists': user_to_training_dists,
        'training_dists': training_dists,
        'threshold': threshold,
        'coverage_percentage': coverage_pct,
        'user_features_clean': user_features_clean,
        'uncovered_mask': user_to_training_dists > threshold
    }


def plot_coverage_analysis(coverage_results):
    """
    Plot distance distributions for coverage analysis.
    
    Args:
        coverage_results: Dictionary returned by compute_coverage_metrics
    """
    user_to_training_dists = coverage_results['user_to_training_dists']
    training_dists = coverage_results['training_dists']
    threshold = coverage_results['threshold']
    coverage_pct = coverage_results['coverage_percentage']
    
    plt.figure(figsize=(10, 6))
    
    # Calculate shared bin range and bins
    all_distances = np.concatenate([training_dists, user_to_training_dists])
    min_dist, max_dist = all_distances.min(), all_distances.max()
    bins = np.linspace(min_dist, max_dist, 50)
    
    # Plot training-to-training distances
    plt.hist(training_dists, bins=bins, alpha=0.7, label='Training-to-Training Distances', 
             color='#999999', density=True)
    
    # Plot user-to-training distances  
    plt.hist(user_to_training_dists, bins=bins, alpha=0.7, label='User-to-Training Distances', 
             color='#4292C6', density=True)
    
    # Add threshold line
    plt.axvline(threshold, color='red', linestyle='--', linewidth=2, 
               label=f'99th Percentile Threshold ({threshold:.2f})')
    
    plt.xlabel('Distance in Standardized Feature Space')
    plt.ylabel('Density')
    plt.legend()
    
    # Add text annotation with coverage
    plt.text(0.7, 0.9, f'Coverage: {coverage_pct:.1f}%', 
             transform=plt.gca().transAxes, fontsize=12, 
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    sns.despine()
    
    plt.show()


def analyze_uncovered_points(user_data, coverage_results, training_merged, feature_cols):
    """
    Analyze uncovered points by plotting their feature distributions.
    
    Args:
        user_data: Original user data DataFrame
        coverage_results: Dictionary returned by compute_coverage_metrics
        training_merged: Training data DataFrame
        feature_cols: List of feature column names
    """
    uncovered_mask = coverage_results['uncovered_mask']
    user_features_clean = coverage_results['user_features_clean']
    
    print(f"Uncovered points: {uncovered_mask.sum()} ({uncovered_mask.mean()*100:.1f}%)")
    
    if uncovered_mask.sum() > 0:
        # Get uncovered user data
        uncovered_indices = user_features_clean.index[uncovered_mask]
        uncovered_data = user_data.loc[uncovered_indices]
        
        # Plot feature distributions comparing uncovered points to training data
        print("\nFeature distributions: Uncovered user data vs Training data")
        plot_feature_distributions(training_merged, uncovered_data, feature_cols=feature_cols)
    else:
        print("All user data points are covered - no outliers to analyze!")


def run_coverage_analysis(user_data, feature_cols, training_data_info):
    """
    Run coverage analysis including metrics computation and plotting.
    Does NOT include uncovered points analysis - use analyze_uncovered_points separately.
    
    Args:
        user_data: DataFrame with user features
        feature_cols: List of feature column names  
        training_data_info: Dictionary with precomputed training data
        
    Returns:
        dict: Coverage analysis results
    """
    # Compute coverage metrics
    coverage_results = compute_coverage_metrics(user_data, feature_cols, training_data_info)
    
    # Print coverage summary
    print(f"Coverage: {coverage_results['coverage_percentage']:.1f}% of user data within training feature space")
    
    # Plot coverage analysis
    plot_coverage_analysis(coverage_results)
    
    return coverage_results


def run_full_coverage_analysis(user_data, feature_cols, training_data_info, training_merged):
    """
    Run complete coverage analysis including metrics computation, plotting, and uncovered points analysis.
    
    Args:
        user_data: DataFrame with user features
        feature_cols: List of feature column names  
        training_data_info: Dictionary with precomputed training data
        training_merged: Training data DataFrame for uncovered analysis
        
    Returns:
        dict: Coverage analysis results
    """
    # Run main coverage analysis
    coverage_results = run_coverage_analysis(user_data, feature_cols, training_data_info)
    
    # Analyze uncovered points
    analyze_uncovered_points(user_data, coverage_results, training_merged, feature_cols)
    
    return coverage_results


def compute_umap_results(user_data, feature_cols, training_data_info):
    """
    Compute UMAP embedding for user data using precomputed UMAP reducer.
    
    Args:
        user_data: DataFrame with user features
        feature_cols: List of feature column names
        training_data_info: Dictionary with precomputed training data
        
    Returns:
        dict: UMAP results ready for plotting
    """
    user_features_clean = user_data[feature_cols].dropna()
    user_embedding = training_data_info['umap_reducer'].transform(user_features_clean[feature_cols])
    
    return {
        'training_embedding': training_data_info['training_embedding'],
        'user_embedding': user_embedding,
        'training_lengths': training_data_info['training_lengths']
    }
