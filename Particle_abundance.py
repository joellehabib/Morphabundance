#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Joelle Habib
# Date: 15 November 2025
#
# Script to compute particle abundance and biovolume from ECOTAXA and ECOPART data
#
# Input:
#   - Float_processed.tsv (ecopart data with sampled volumes per depth bin)
#   - float_ecotaxa.tsv (ecotaxa data with particle annotations and morphology)
#
# Output:
#   - Processed dataframe with particle abundance (particles/L) and biovolume (mm³/L)
#     per sample_id and depth_bin

from pathlib import Path
import pandas as pd
import numpy as np
import cmocean
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from seawater import dpth

# Function to round to the closest 5 value below the number

def rounditup(x, precision, method):
    import math
    if method == "floor":
        return math.floor(x / precision) * precision
    elif method == "ceiling":
        return math.ceil(x / precision) * precision
    else:
        return "give the parameter floor or ceiling"
    
cmpa=cmocean.cm.algae
cmoxy=cmocean.cm.balance

def nonlinear_colormap():
    import pylab as pyl
    #import numpy as np
    levels1 = [0, 1, 2]

    ####################################################################
    ###                      non linear colormap                     ###
    ####################################################################

    """
    nlcmap - a nonlinear cmap from specified levels

    Copyright (c) 2006-2007, Robert Hetland <hetland@tamu.edu>
    Release under MIT license.

    Some hacks added 2012 noted in code (@MRR)
    """

    from matplotlib.colors import LinearSegmentedColormap


    class nlcmap(LinearSegmentedColormap):
        """Nonlinear colormap.
           Needs the input of a linear colormap e.g. pylab.cm.jet
           and a list of levels e.g. [0,1,2,3,6]
           in order to change the interpolation, change the levels vector,
           e.g ad 20, 50 to make it more unlinear."""
        import numpy as np
        name = 'nlcmap'

        def __init__(self, cmap, levels):
            import numpy as np
            self.cmap = cmap
            # @MRR: Need to add N for backend
            self.N = cmap.N
            self.monochrome = self.cmap.monochrome
            self.levels = np.asarray(levels, dtype='float64')
            self._x = self.levels / self.levels.max()
            self._y = np.linspace(0.0, 1.0, len(self.levels))

        #@MRR Need to add **kw for 'bytes'
        def __call__(self, xi, alpha=1.0, **kw):
            import numpy as np
            """docstring for fname"""
            # @MRR: Appears broken?
            # It appears something's wrong with the
            # dimensionality of a calculation intermediate
            #yi = stineman_interp(xi, self._x, self._y)
            yi = np.interp(xi, self._x, self._y)
            return self.cmap(yi, alpha)


    cmap_nonlin = nlcmap(pyl.cm.CMRmap, levels1)
    return cmap_nonlin



cm = nonlinear_colormap()


def contour_levels_func(min_contour_level, max_contour_level, levels):
    """Function to define contour levels for contourf"""
    distance_levels = max_contour_level / levels
    contour_levels = np.arange(min_contour_level, max_contour_level, distance_levels)
    return contour_levels


def gridding_func(pos_min_max, depth_min_max, pos_array, depth, param):
    grid_method = "linear"  # choose the gridding method here
    # method to do the regridding, can be nearest (for nearest neighbour) or linear

    xi = np.linspace(min(pos_min_max), max(pos_min_max), 1000)
    yi = np.linspace(min(depth_min_max), max(depth_min_max), 200)
    zi = griddata((pos_array, depth), param, (xi[None, :], yi[:, None]), method=grid_method)
    return xi, yi, zi
    
    
    

#defining the path
Path_to_data = Path("~/GIT/TRATLEQ/Data").expanduser()


df_ = pd.read_csv(Path_to_data/"ecotaxa_export_158_20230427.tsv",sep='\t',low_memory=False)

#variables to utilize from the tables 
#for object details 
#df_.sample_id.unique()

df_=df_[["object_lat","object_lon","object_annotation_status","object_annotation_category",
         "object_area_exc","object_area","object_major","object_minor","process_pixel",
         "acq_id","object_depth_min","object_depth_max","sample_id", "object_esd"]]

#object_area_exc=surface area of the objects excluded the empty spaces (holes) in pixel
# object_area = surface area of the objects included the empty spaces (holes) in pixel
#object_major = primary axis of the best fitting ellipse for the object in pixel
#object_minor = secondary axis of the best fitting ellipse for the object in pixel

#I have to concert the pixel variables in metric form 

# object_area=np.array(df_.object_area)
# process_pixel=np.array(df_.process_pixel)

# Area_mm2=object_area*(process_pixel*10e-3)**2

df_["Area_mm2"]=(df_.object_area)*(df_.process_pixel)**2
df_["Area_exc_mm2"]=(df_.object_area_exc)*(df_.process_pixel)**2


df_["Major_mm"]=(df_.object_major)*(df_.process_pixel)
df_["Minor_mm"]=(df_.object_minor)*(df_.process_pixel)


df_["esd"] = 2*np.sqrt(df_["Area_mm2"]/np.pi) # equivalent spherical diameter
df_["vol_sph"] = 4/3*np.pi*(df_["esd"]/2)**3 # spherical volume


# shift by 1.2m the ecotaxa depth in order for them to coincide with the ecopart data
df_['object_depth_min'] = df_['object_depth_min'].apply(lambda x: x + 1.2)
df_['object_depth_max'] = df_['object_depth_max'].apply(lambda x: x + 1.2)


# compute mean depth of ecotaxa bin
df_["depth"] = df_.loc[:, ["object_depth_min","object_depth_max"]].mean(axis= 1)


# Function to round to the closest 5 value below the number rounditup
# Bin the data using a new column
df_['depth_bin'] = df_["depth"].apply(lambda x: rounditup(x, 5, "floor")+ 2.5)



#Take only particles between 1-5mmm from ecotaxa
min_ESD=1
max_ESD=5

df_=df_[(df_["esd"]>=min_ESD) & (df_["esd"]<=max_ESD)]


list = ["sample_id", "depth_bin"]
#groupping based on sampleid and depth bin for each depth
df_grouped_counts = df_.groupby(list)["sample_id"].agg([("count", np.count_nonzero)])

df_grouped_biov = df_.groupby(["sample_id","depth_bin"]).agg(biov=('vol_sph', np.sum))

df_grouped_counts["biov"]=df_grouped_biov
#grouping of only lat/lon
df_grouped_meta = df_.groupby(["sample_id"]).agg(lat=('object_lat', np.mean),
                                                    lon=('object_lon', np.mean))



# I have to extract the volume of water from the ecopart file *aggregate.tsv
# df_ecopart = pd.read_csv(Path_to_data/"M158_M181_detailed_20230128_PAR_Aggregated.tsv",sep='\t',encoding= 'unicode_escape',low_memory=False)
# df_ecopart=df_ecopart[["Project","Sampled volume [L]","Profile"]]
#here £i decided to take the processede file because I already have the same data

df_ecopart = pd.read_csv(Path_to_data / "M158_M181_processed.tsv", sep='\t')

#print(df_.columns.tolist())
df_ecopart =df_ecopart[["Profile", "Pressure [dbar]", "Latitude", "Longitude",
        "Date_Time", "Project", "Vol [L] (sampled for this depth bin)"]]


project_dic={"M158" : "uvp5_sn210_2019_m158"}
df_ecopart=df_ecopart[df_ecopart["Project"] == "uvp5_sn210_2019_m158"]


list2=["Profile",'Pressure [dbar]']

df_grouped_ecopart = df_ecopart.groupby(list2).agg(vol=('Vol [L] (sampled for this depth bin)', np.mean),
                                                   depth_bin=('Pressure [dbar]', np.mean))



df_grouped_counts = df_grouped_counts.reset_index()

df_grouped_meta = df_grouped_meta.reset_index()
df_grouped_ecopart = df_grouped_ecopart.reset_index()

#merge the data 
#a) based on sasmple id 
df_grouped = df_grouped_meta.merge(df_grouped_counts, how="left", on=["sample_id"])

df_grouped_ecopart=df_grouped_ecopart.sort_values(by=['Profile', 'depth_bin'])
df_grouped_ecopart["sample_id"]=df_grouped_ecopart["Profile"]


#b)
df_grouped= df_grouped_ecopart.merge(df_grouped, how="left", on=["depth_bin","sample_id"])


#L/1000=m3
divisor_volume=1000 #if I want it in m3

df_grouped["vol"] = df_grouped["vol"] #* divisor_volume
df_grouped["abundance_per_L"] = df_grouped["count"] / df_grouped["vol"]
df_grouped["BIOV"]=df_grouped["biov"] / df_grouped["vol"]

#put the correct lon,lat when it is nan and then biov and abundance== nan when ther is no data
# assuming that the nan are in the middle of each profile
df_grouped['lat'].bfill(inplace=True)
df_grouped['lon'].bfill(inplace=True)

df_grouped['BIOV'] = df_grouped['BIOV'].fillna(0)
df_grouped['abundance_per_L'] = df_grouped['abundance_per_L'].fillna(0)

# Print the first 5 rows of results
print(f"\n=== Results ({len(df_grouped)} rows) ===")
if len(df_grouped) > 0:
    cols = ['sample_id', 'depth_bin', 'lat', 'lon', 'abundance_per_L', 'BIOV']
    display_cols = [col for col in cols if col in df_grouped.columns]
    print(df_grouped[display_cols].head(5).to_string(index=False))
else:
    print("No data available")