# Author: Joelle Habib
# Date: 15 November 2025
# Script to run after 1_morphospce.R to extract the data from ecopart and ecotaxa
# that you will need for the next one
#
# Input: 
#   - Float_processed.tsv (ecopart data)
#   - float_ecotaxa.tsv (ecotaxa data)
#
# Output:
#   - volumes_float.csv (sampled volumes per profile and depth)
#   - list_of_profiles.csv (metadata for selected profiles)
#   - ecotaxa_float_df.csv (processed ecotaxa data with morphological features)

import os
os.environ['PROJ_LIB'] = '/Users/joellehabib/anaconda3/pkgs/proj-8.2.1-hd69def0_0/share/proj'


from pathlib import Path
import pandas as pd
import numpy as np
import cmocean
import matplotlib
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from seawater import dpth
import calendar
import time
from datetime import date,datetime
from scipy.signal import savgol_filter
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
import itertools

from tabulate import tabulate

# Function to round to the closest 5 value below the number

def rounditup(x, precision, method):
    import math
    if method == "floor":
        return math.floor(x / precision) * precision
    elif method == "ceiling":
        return math.ceil(x / precision) * precision
    else:
        return "give the parameter floor or ceiling"

#defining the path
Path_to_data = Path("~/GIT/TRATLEQ/Data/Float_data").expanduser()


#########################
#1- ecopart data
#########################

df_ecopart1 = pd.read_csv(Path_to_data / "Float_processed3.tsv", sep='\t')

# rename some of the columns
df_ecopart1.rename(columns={'Vol [L] (sampled for this depth bin)': 'volume_L'}, inplace=True)

# select only ascent profiles :
ecopart = df_ecopart1[df_ecopart1['profile'].str.contains('a')]

# compute max and min depth for the profiles
profiles_selected = ecopart.groupby("profile").agg({"Pressure [dbar]": [np.max, np.min]})
profiles_selected.columns = ['max_depth', "min_depth"]
profiles_selected = profiles_selected.reset_index()
profiles_selected["range"] = profiles_selected["max_depth"] - profiles_selected["min_depth"]  # depth range of the layer



# # look at their distribution
# # max depth
# plt.figure()
# fig_depth_max = sns.displot(profiles_selected, x="max_depth", kde=True, binwidth=50)
# plt.show()
# # fig_depth_max.savefig(str(path_to_figures) + '/Distribution_max_depth.pdf')

# # min depth
# plt.figure()
# fig_depth_min = sns.displot(profiles_selected, x="min_depth", kde=True, binwidth=50)
# plt.show()
# # fig_depth_min.savefig(str(path_to_figures) + '/Distribution_min_depth.pdf')

# add the lon lat and date to it
profiles_selected_list = list(profiles_selected["profile"].unique())
len(profiles_selected_list)

# Keep in the ecopart table only the selected profiles
ecopart_selected = ecopart.merge(profiles_selected, on=['profile'])

ecopart_selected = ecopart_selected.filter(items=['profile', 'Pressure [dbar]', 'volume_L'])
ecopart_selected.rename(columns={'Pressure [dbar]': 'depth'}, inplace=True)
print(ecopart_selected)
print(len(np.unique(ecopart_selected['profile'])))


# first compute a new binned_depth column because ecopart data have 5m bins and I want 50m bins
#ecopart_selected['depth_bin'] = ecopart_selected["depth"].apply(lambda x: rounditup(x, 50, "floor")+ 25)
#ecopart_selected = ecopart_selected.drop('depth', axis = 1)
#ecopart_selected.columns = ecopart_selected.columns.str.replace('_bin', '')

# compute the sum of sampled volume for each bin of each profile
#ecopart_selected = ecopart_selected.groupby(['depth', 'profile']).agg({'volume_L': 'sum'})
#ecopart_selected.reset_index(inplace = True)


# save the table
volumes = ecopart_selected
volumes.to_csv(Path_to_data/"volumes_float.csv", index=False)
volumes = pd.read_csv(Path_to_data/"volumes_float.csv")


# 2.2 Ecotaxa

# Import ecotaxa table
df_ = pd.read_csv(Path_to_data/"float_ecotaxa.tsv",sep='\t',low_memory=False)


RAWfilename=df_.process_id
ct=0
sel_filename = [True for i in range(RAWfilename.size)]
for a in RAWfilename:
    if a[4] == 'a':
        sel_filename[ct]=True
    else:
        sel_filename[ct] = False
    ct+=1
    
df_['sel_filename']=sel_filename
# I remove all sel_filename= false 
df_=df_.loc[df_.sel_filename==True]

df_.head()
df_.shape
df_.object_depth_min.describe()

print("Max lat:", np.max(df_['object_lat']))
print("min lat:", np.min(df_['object_lat']))

print("Max lon:", np.max(df_['object_lon']))
print("min lon:", np.min(df_['object_lon']))

# remove "object_" from the columns names
df_.columns = df_.columns.str.replace('object_', '')


# keep only relevant information :
df_ = df_[['id', 'lat', 'lon', 'date', 'time', 'depth_min', 'depth_max', 'sample_id', 'annotation_category',
     'annotation_hierarchy', 'area','area_exc', 'minor', 'major', 'acq_pixel', 'perim.', 'circ.', 'mean', 'kurt', 'fractal']]

print(df_)

df_["Area_mm2"]=(df_.area)*(df_.acq_pixel)**2
df_["Area_exc_mm2"]=(df_.area_exc)*(df_.acq_pixel)**2


df_["Major_mm"]=(df_.major)*(df_.acq_pixel)
df_["Minor_mm"]=(df_.minor)*(df_.acq_pixel)


df_["esd"] = 2*np.sqrt(df_["Area_mm2"]/np.pi) # equivalent spherical diameter
df_["vol_sph"] = 4/3*np.pi*(df_["esd"]/2)**3 # spherical volume


# shift by 1.2m the ecotaxa depth in order for them to coincide with the ecopart data
df_['depth_min'] = df_['depth_min'].apply(lambda x: x + 1.2)
df_['depth_max'] = df_['depth_max'].apply(lambda x: x + 1.2)


# compute mean depth of ecotaxa bin
df_["depth"] = df_.loc[:, ["depth_min","depth_max"]].mean(axis= 1)



# Function to round to the closest 5 value below the number rounditup
# Bin the data using a new column
df_['depth_bin'] = df_["depth"].apply(lambda x: rounditup(x, 5, "floor")+ 2.5)

# rename columns
df_ = df_.rename(columns={"sample_id": "profile", "id": "object_id", "annotation_category": "taxon", "acq_pixel": "px_size",
             "area": "area_px"})

# compute max and min depth for the profiles
profiles_selected = df_.groupby("profile").agg({"depth": [np.max, np.min]})
profiles_selected.columns = ['max_depth', "min_depth"]
profiles_selected = profiles_selected.reset_index()
profiles_selected["range"] = profiles_selected["max_depth"] - profiles_selected["min_depth"]  # depth range of the layer


# # look at their distribution
# # max depth
# plt.figure()
# fig_depth_max = sns.displot(profiles_selected, x="max_depth", kde=True, binwidth=50)
# plt.show()
# # fig_depth_max.savefig(str(path_to_figures) + '/Distribution_max_depth.pdf')

# # min depth
# plt.figure()
# fig_depth_min = sns.displot(profiles_selected, x="min_depth", kde=True, binwidth=50)
# plt.show()
# # fig_depth_min.savefig(str(path_to_figures) + '/Distribution_min_depth.pdf')

# # add the lon lat and date to it
# profiles_selected_list = list(profiles_selected["profile"].unique())
# len(profiles_selected_list)

# get the metadata of the profiles
coordinates_profiles = df_.filter(items=['profile', 'lat', 'lon', 'date'])

# 2.2.2 Add the ecopart data and keep only the selected profiles

# add the ecotaxa lon and lat to the profiles we selected earlier
profiles = pd.DataFrame(np.unique(ecopart_selected.profile), columns=['profile'])
# add the coordinates and date of each profile
profiles = profiles.merge(coordinates_profiles, on=['profile'])

# keep only one line per profile and per date
list_prof = profiles.groupby(['profile']).agg({'lat': 'mean', 'lon': 'mean', 'date': 'mean'}).reset_index() #j'ai enlever le groupby date

print(len(np.unique(list_prof['profile'])))

# check that it was done properly
print(len(np.unique(ecopart_selected.profile)))
print(len(np.unique(list_prof.profile)))


ecopart_selected_unique = ecopart_selected.groupby(['profile']).size().reset_index(name='counts')
list_prof_unique = list_prof.groupby(['profile']).size().reset_index(name='counts')

# check this selection to see if the same profiles are found here
df_1notin2 = ecopart_selected_unique[~(ecopart_selected_unique['profile'].isin(list_prof_unique['profile']))].reset_index(drop=True)
# this shows profiles which are in ecopart_selected_unique (filtered ecopart profiles)
# but not in list_prof_unique (profiles from ecotaxa which are among the filtered ecopart profiles)
len(df_1notin2) == len(ecopart_selected_unique) - len(list_prof_unique)  # True
print(df_1notin2)  # it's empty so we have the same profiles in both tables. Great !


# save the table for later
list_prof.to_csv(Path_to_data/"list_of_profiles.csv", index=False)

# select in ecotaxa only the profiles which are present in list_prof
list_prof_reduced = list_prof.filter(items=['profile'])  # keep only the profile column
# this will allow us to not have two date, lon and lat columns as they are not exactly the same
# between CTD and UVP data
df_before_taxo = df_.merge(list_prof_reduced, on=['profile'])


# save
df_.reset_index().to_csv(Path_to_data/"ecotaxa_float_df.csv")