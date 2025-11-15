# Author: Joelle Habib
# Date: 15 November 2025
# Script to create the final csv for the clusters
# Bin the clusters and plot each cluster as a section
#
# Input:
#   - indiv_float.csv (individual particle data with object IDs)
#   - volumes_float.csv (sampled volumes per profile and depth)
#   - k-means_pca_results_q25.csv (cluster assignments from morphospace analysis)
#   - list_of_profiles.csv (metadata with lat, lon, date for profiles)
#
# Output:
#   - cluster_concentration.csv (binned cluster data with concentrations and biovolume)

import os
os.environ['PROJ_LIB'] = '/Users/joellehabib/anaconda3/pkgs/proj-8.2.1-hd69def0_0/share/proj'


from pathlib import Path

import numpy as np
import pandas as pd
import warnings  # these two lines to remove the annoying warnings from pandas
warnings.simplefilter(action='ignore', category=FutureWarning)
import itertools


#defining the path

Path_to_data = Path("~/GIT/TRATLEQ/Data/Float_data").expanduser()

#open the files where I have the different particle 
# first I need to import the csv indiv_float that contains object id in order to compute the concentration of each cluster considered as species
Path_to_data = Path("~/GIT/TRATLEQ/Data/Float_data/").expanduser()
indiv = pd.read_csv(Path_to_data/"indiv_float.csv")

# and the sampled volume
volumes = pd.read_csv(Path_to_data/'volumes_float.csv')

# and the object with the associated cluster
cluster = pd.read_csv(Path_to_data/"k-means_pca_results_q25.csv")

cluster.drop(['datetime', 'depth_min', 'Dim.1', 'Dim.2', 'Dim.3', 'Dim.4'], inplace=True, axis=1)
cluster.rename(columns = {'id':'object_id'}, inplace = True)

# delete individuals with column 'keep' = False (result from morphospace R script, it's individuals that are far from their cluster center)
# cluster = cluster[cluster['keep'] == True]


# add the column cluster
indiv = pd.merge(indiv, cluster, on = 'object_id')

# Compute abundance and biovolume per cluster and bin for either the rough or medium taxonomic definiton
indiv_binned = indiv.groupby(['Cluster', 'depth_bin', 'profile'], group_keys=False).agg(n = pd.NamedAgg(column = 'Cluster', aggfunc = 'count'),
                                                                                         vol_sph = pd.NamedAgg(column = "vol_sph", aggfunc = 'sum'),
                                                                                         perim = pd.NamedAgg(column = "perim.", aggfunc = 'mean'),
                                                                                         circ = pd.NamedAgg(column = "circ.", aggfunc = 'mean'),
                                                                                         mean = pd.NamedAgg(column = "mean", aggfunc = 'mean'),
                                                                                         kurt = pd.NamedAgg(column = "kurt", aggfunc = 'mean'),
                                                                                         esd = pd.NamedAgg(column = "esd", aggfunc = 'mean'),
                                                                                         fractal = pd.NamedAgg(column = "fractal", aggfunc = 'mean'))



indiv_binned.reset_index(inplace = True) # to keep a column with exports_groups and depth_bin values
indiv_binned.columns = indiv_binned.columns.str.replace('_bin', '')

volumes=volumes.rename(columns={"process_id": "profile","depth_bin":"depth","volume_L":"vol"})


# add the volume and compute the concentrations
obs = pd.merge(indiv_binned, volumes, how="left", on=['profile', 'depth'])
obs["watervolume"] = obs["vol"] / 1000          # volume in m3
obs["conc"] = obs["n"] / obs["watervolume"]          # concentration in n/m3
obs["vol_sph"] = obs["vol_sph"] / obs["watervolume"] # biovolume concentration in mm3/m3


# keep only the 5m depth bins which have a watervolume <= 1000L
obs = obs[obs["vol"] > 0]
obs = obs.dropna(subset=['watervolume'])


# compute all cluster x bins combinations
depth_bins = np.unique(volumes.depth)

# and all the profiles
profile_list = volumes['profile'].unique()

cluster_list = indiv.Cluster.unique()

# regroup them with the list of clusters
list_ = [profile_list, depth_bins, cluster_list]
# compute the unique combinations of the 3 lists
combination  = [p for p in itertools.product(*list_)]

# Convert it into a dataframe
column_names = ["profile", "depth", "Cluster"]
all1 = pd.DataFrame(combination, columns =["profile", "depth", "Cluster"])

# Add the data from obs
full = pd.merge(all1, obs, how="left", on=['profile', 'depth', "Cluster"])
full.head()

full = pd.DataFrame(full)

# remove the volume_L column. We'll keep only watervolume (in m3)
full.drop('vol', inplace = True, axis = 1)

# consider the observations we don't have as zeroes
cols_to_fill = ['n', 'conc', 'vol_sph', 'watervolume']
full[cols_to_fill]=full[cols_to_fill].fillna(0)


#remove nan value for interpolation 
full = full.dropna()
full = full.reset_index()
full.drop(['index'], inplace=True, axis=1)

#ADD  lon lat and date  
profiles = pd.read_csv(Path_to_data/"list_of_profiles.csv")

full_final=pd.merge(full, profiles, how="left", on=['profile'])
full_final.to_csv(Path_to_data/"cluster_concentration.csv", index=False)
